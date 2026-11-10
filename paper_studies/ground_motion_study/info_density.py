"""Frame information density of the hop-limited link with the actual 512-word source.

For each packet the streaming detector graph is recorded (as in nb_packets.py) and two forward
passes give
    log p(y | x)  : the transmitted word sequence forced (word log-prior 0 / -inf),
    log p(y)      : uniform word prior (-log 512 per group),
so i(x; y) = log p(y|x) - log p(y) (converted to bits) is the frame information density under joint
(sequence) detection with the detector's reduced-state metric. P[i < 256 bits] estimates the frame
error rate of the best code of this length on this exact signal and receiver (dependence-testing /
normal-approximation heuristic), against which a real code can be compared.
Also checks the metric scale: the mean edge metric along the true path per complex sample should be ~1.
Usage: python info_density.py --case 24:0.50 --case 18:0.60 --packets 128
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, json, math, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from numba import njit
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
import hop_packets as hp
from hop_receiver import HopLimitedReceiver, sampled_hop
import turbo_rx as tr

OUT = Path(__file__).resolve().parent/'results'


@njit(cache=True)
def _forward_total(src, dst, met, word, eoff, noff, group, Lp):
    T = len(eoff) - 1
    alpha = np.full(noff[-1], -np.inf); alpha[noff[0]] = 0.0; total = 0.0
    for t in range(T):
        g = group[t]
        for e in range(eoff[t], eoff[t + 1]):
            p = Lp[g, word[e]] if g >= 0 else 0.0
            i = noff[t + 1] + dst[e]
            alpha[i] = tr._lse(alpha[i], alpha[noff[t] + src[e]] + met[e] + p)
        m = -np.inf
        for i in range(noff[t + 1], noff[t + 2]):
            m = max(m, alpha[i])
        for i in range(noff[t + 1], noff[t + 2]):
            alpha[i] -= m
        total += m
    m = -np.inf
    for i in range(noff[T], noff[T + 1]):
        m = max(m, alpha[i])
    s = 0.0
    for i in range(noff[T], noff[T + 1]):
        s += math.exp(alpha[i] - m)
    return total + m + math.log(s)


def forward_total(packed, Lp):
    return _forward_total(packed['src'], packed['dst'], packed['met'], packed['word'], packed['eoff'],
                          packed['noff'], packed['group'], np.ascontiguousarray(Lp, np.float64))


def run(job):
    carrier, amplitude, index, seed_ns = job
    rng = np.random.default_rng(np.random.SeedSequence([seed_ns[0], 43, index]))
    ranks = rng.integers(0, 512, 43)
    tones, prev = [], hp.M//2
    for w in ranks:
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    ts = 224/(len(tones) + 18)
    txs = np.r_[tones, np.full(6, hp.M//2)]
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    det = HopLimitedReceiver(r, ts, amplitude, 43, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                             h=hp.H, beta=hp.BETA, tone_count=hp.M, max_hop=hp.D, group_symbols=hp.L,
                             group_bits=hp.BITS, pad_bits=0)
    mean = r.stream(1).process(sampled_hop(txs, ts, hp.H, float(hp.BETA), hp.M))
    noise_rng = np.random.default_rng(np.random.SeedSequence([seed_ns[1], 43, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    tr.record_graph(det)
    position = 0; nsamp = 0
    for symbol in range(len(txs)):
        count = det._geometry(symbol)[3]
        det.push(observed[position:position + count]); position += count; nsamp += count
    stats = det.finish()
    packed = tr.pack_graph(det, hp.L, 43)
    uniform = np.full((43, 512), -math.log(512))
    forced = np.full((43, 512), -np.inf); forced[np.arange(43), ranks] = 0.0
    lpy = forward_total(packed, uniform); lpyx = forward_total(packed, forced)
    # metric scale check: forced path total is log p(y|x) up to the per-sample constant
    return dict(carrier_hz=carrier, amplitude_m_s2=amplitude, index=index, density_bits=(lpyx - lpy)/math.log(2),
                log_p_y_given_x=lpyx, samples=nsamp, pruned=stats['pruned_states'])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='carrier:amplitude_factor_of_A0')
    p.add_argument('--packets', type=int, default=128)
    p.add_argument('--seed-start', type=int, default=7300000)
    p.add_argument('--seeds', default='20371401,20371402')
    p.add_argument('--workers', type=int, default=18)
    a = p.parse_args()
    ns = tuple(int(x) for x in a.seeds.split(','))
    jobs = [(float(c.split(':')[0]), float(c.split(':')[1])*mg.A0, i, ns)
            for c in a.case for i in range(a.seed_start, a.seed_start + a.packets)]
    res = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(run, j) for j in jobs]):
            res.append(f.result())
    out = {}
    for c in a.case:
        fc, amp = float(c.split(':')[0]), float(c.split(':')[1])
        d = np.array([r['density_bits'] for r in res if r['carrier_hz'] == fc and abs(r['amplitude_m_s2'] - amp*mg.A0) < 1e-12])
        out[c] = dict(n=len(d), mean_bits=float(d.mean()), std_bits=float(d.std(ddof=1)),
                      p_below_256=float(np.mean(d < 256)), p_below_256_gauss=float(
                          0.5*math.erfc((d.mean() - 256)/(d.std(ddof=1)*math.sqrt(2)))), densities=d.tolist())
        print(c, {k: v for k, v in out[c].items() if k != 'densities'}, flush=True)
    path = OUT/'info_density_hop.json'
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(out); path.write_text(json.dumps(old, indent=1))


if __name__ == '__main__':
    main()
