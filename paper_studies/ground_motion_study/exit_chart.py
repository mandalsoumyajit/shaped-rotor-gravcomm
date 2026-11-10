"""EXIT-chart design of the iterative CF-FSK receiver (joint choice of labeling and outer code).

Extrinsic-information-transfer analysis (ten Brink 2001) of the two soft-in/soft-out components:
  * inner: the detector graph (turbo_rx.graph_app) of the tone trellis in the colored receiver noise,
    fed with consistent Gaussian a-priori LLRs of mutual information I_A on the coded bits; its
    transfer curve I_E = T_det(I_A; amplitude, labeling) is measured on simulated packets.
    The labeling (word -> bit pattern) only enters the soft-output stage, so many labelings are
    evaluated on the same detector graphs;
  * outer: candidate tail-biting convolutional codes (memory 2..6, punctured to rate 2/3), whose
    transfer curve I_E = T_code(I_A) is measured with Gaussian priors on their coded bits.
Iterative decoding converges (for long frames) when T_det(I) > T_code^{-1}(I) on 0 <= I < 1,
i.e. the tunnel is open. The design rule is the labeling and code whose tunnel opens at the
lowest amplitude. Outputs: results/exit_detector.json, results/exit_codes.json.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import json, math, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from numba import njit
from scipy.optimize import brentq
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver, sampled_cffsk
from gravcomm.streaming_soft import StreamingSoftReceiver
import turbo_rx as tr

RES = Path(__file__).resolve().parent/'results'
IA = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97])


# ------------------------------------------------------------------ mutual information helpers
def mi_llr(llr, bits):
    """Time-average estimate of I(bit; LLR) for LLR = log P(0)/P(1)."""
    x = llr*(1 - 2*bits.astype(float))
    return float(np.mean(1 - np.logaddexp(0, -x)/np.log(2)))


_SIG = np.linspace(1e-3, 60, 6000)


def _J(s):
    if s < 1e-3:
        return 0.0
    y = np.linspace(-8, 8, 801)
    L = s*s/2 + s*y
    return float(1 - np.trapezoid(np.exp(-y*y/2)/np.sqrt(2*np.pi)*np.logaddexp(0, -L)/np.log(2), y))


_JT = np.array([_J(s) for s in _SIG])


def sigma_for(I):
    return float(np.interp(I, _JT, _SIG)) if I > 0 else 0.0


def gaussian_prior(bits, I, rng):
    s = sigma_for(I)
    return (s*s/2)*(1 - 2*bits.astype(float)) + s*rng.normal(size=bits.shape)


# ------------------------------------------------------------------ generic tail-biting code SISO
def code_table(gens, memory):
    S = 1 << memory
    out = np.array([[[((s << 1 | b) & g).bit_count() % 2 for g in gens] for b in range(2)] for s in range(S)],
                   dtype=np.uint8)
    return out


@njit(cache=True)
def _lse(a, b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    m = max(a, b)
    return m + math.log(math.exp(a - m) + math.exp(b - m))


@njit(cache=True)
def _bcjr_generic(L, outputs, laps=2):
    K = L.shape[0]; S = outputs.shape[0]; mask = S - 1; n = outputs.shape[2]
    g = np.empty((K, S, 2))
    for t in range(K):
        for s in range(S):
            for b in range(2):
                v = 0.0
                for i in range(n):
                    v -= outputs[s, b, i]*L[t, i]
                g[t, s, b] = v
    alpha = np.zeros((K + 1, S)); beta = np.zeros((K + 1, S)); a = np.zeros(S)
    for lap in range(laps):
        alpha[0] = a
        for t in range(K):
            nxt = np.full(S, -np.inf)
            for s in range(S):
                for b in range(2):
                    d = ((s << 1) | b) & mask
                    nxt[d] = _lse(nxt[d], alpha[t, s] + g[t, s, b])
            alpha[t + 1] = nxt - nxt.max()
        a = alpha[K].copy()
    bb = np.zeros(S)
    for lap in range(laps):
        beta[K] = bb
        for t in range(K - 1, -1, -1):
            cur = np.full(S, -np.inf)
            for s in range(S):
                for b in range(2):
                    d = ((s << 1) | b) & mask
                    cur[s] = _lse(cur[s], g[t, s, b] + beta[t + 1, d])
            beta[t] = cur - cur.max()
        bb = beta[0].copy()
    info = np.empty(K); app = np.empty((K, n))
    for t in range(K):
        u0 = -np.inf; u1 = -np.inf
        c0 = np.full(n, -np.inf); c1 = np.full(n, -np.inf)
        for s in range(S):
            for b in range(2):
                d = ((s << 1) | b) & mask
                v = alpha[t, s] + g[t, s, b] + beta[t + 1, d]
                if b == 0:
                    u0 = _lse(u0, v)
                else:
                    u1 = _lse(u1, v)
                for i in range(n):
                    if outputs[s, b, i] == 0:
                        c0[i] = _lse(c0[i], v)
                    else:
                        c1[i] = _lse(c1[i], v)
        info[t] = u0 - u1
        for i in range(n):
            app[t, i] = c0[i] - c1[i]
    return info, app


class TBCode:
    """Tail-biting feedforward convolutional code, rate 1/2 mother code punctured by `pattern`."""
    def __init__(self, gens, memory, pattern=((1, 1), (0, 1)), name=''):
        self.gens, self.memory, self.name = gens, memory, name
        self.out = code_table(gens, memory); self.pattern = np.array(pattern, bool)

    def mask(self, K):
        return self.pattern[np.arange(K) % len(self.pattern)]

    def encode(self, bits):
        K = len(bits); m = self.memory
        state = sum(int(bits[(-j) % K]) << (j - 1) for j in range(1, m + 1))
        out = np.empty((K, 2), np.uint8)
        for t, b in enumerate(bits):
            out[t] = self.out[state, b]; state = ((state << 1) | int(b)) & ((1 << m) - 1)
        return out[self.mask(K)]

    def siso(self, llr, K):
        mask = self.mask(K); L = np.zeros(mask.shape); L[mask] = llr
        info, app = _bcjr_generic(L, self.out)
        return info, (app - L)[mask]


CODES = {   # standard best-distance rate-1/2 feedforward codes (octal), punctured to 2/3 as released
    'K3 (5,7)': TBCode((0o5, 0o7), 2), 'K4 (15,17)': TBCode((0o15, 0o17), 3),
    'K5 (23,35)': TBCode((0o23, 0o35), 4), 'K6 (53,75)': TBCode((0o53, 0o75), 5),
    'K7 (171,133) released': TBCode((0o171, 0o133), 6)}


def code_exit(code, K=2560, reps=4, seed=0):
    rng = np.random.default_rng(seed); out = []
    for I in IA:
        vals = []
        for _ in range(reps):
            info = rng.integers(0, 2, K).astype(np.uint8); c = code.encode(info)
            La = gaussian_prior(c, I, rng)
            _i, ext = code.siso(La, K)
            vals.append(mi_llr(ext, c))
        out.append(float(np.mean(vals)))
    return out


# ------------------------------------------------------------------ detector EXIT
def labelings():
    nat = np.arange(256)
    gray = nat ^ (nat >> 1)
    lab = {'natural (released)': nat, 'gray of index': gray}
    p = RES/'label_opt.json'
    if p.exists():
        lab['BICM-optimized'] = np.array(json.loads(p.read_text())['optimized']['label'])
    # per-digit Gray: digits (d1,d2,d3) in base 7, each Gray-coded (3 bits), word index = rank
    g7 = [0, 1, 3, 2, 6, 7, 5, 4]
    pd = []
    for w in range(256):
        d = (w//49, (w//7) % 7, w % 7)
        pd.append((g7[d[0]] << 6 | g7[d[1]] << 3 | g7[d[2]]) & 0xFF)
    if len(set(pd)) == 256:
        lab['per-digit Gray'] = np.array(pd)
    rng = np.random.default_rng(5)
    lab['random'] = rng.permutation(256)
    return lab


def detector_job(args):
    carrier, amplitude, index = args
    rng = np.random.default_rng(np.random.SeedSequence([20371107, index]))
    words = rng.integers(0, 256, 48)
    symbols = np.column_stack((words//49, words//7 % 7, words % 7)).ravel() - 3
    ts = 224/(3*48 + 18)
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    mean = r.stream(1).process(sampled_cffsk(np.r_[symbols, np.zeros(6, int)], ts, q=.08, transition_fraction=.6))
    noise_rng = np.random.default_rng(np.random.SeedSequence([20371108, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    det = StreamingSoftReceiver(r, ts, amplitude, 48, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                                q=.08, transition_fraction=.6)
    tr.record_graph(det); position = 0
    for s in range(len(symbols) + 6):
        c = det._geometry(s)[3]; det.push(observed[position:position + c]); position += c
    det.finish()
    pk = tr.pack_graph(det, 3, 48)
    res = {}
    for name, lab in labelings().items():
        labels = tr.bit_labels(256, 8)[lab]
        bits = labels[words].astype(np.uint8)
        row = []
        for I in IA:
            La = gaussian_prior(bits, I, np.random.default_rng(index*100 + int(I*100)))
            app = tr.graph_app(pk, La, labels)
            row.append((app - La).ravel().tolist())
        res[name] = (bits.ravel().tolist(), row)
    return res


def detector_exit(carrier, amplitude, packets=48, workers=17):
    with ProcessPoolExecutor(workers) as pool:
        outs = list(pool.map(detector_job, [(carrier, amplitude*mg.A0, i) for i in range(packets)]))
    curves = {}
    for name in outs[0]:
        bits = np.concatenate([np.array(o[name][0]) for o in outs])
        curves[name] = [mi_llr(np.concatenate([np.array(o[name][1][k]) for o in outs]), bits) for k in range(len(IA))]
    return curves


def main():
    import sys
    t0 = time.time()
    if 'codes' in sys.argv:
        out = {n: code_exit(c) for n, c in CODES.items()}
        (RES/'exit_codes.json').write_text(json.dumps(dict(IA=IA.tolist(), curves=out), indent=1) + '\n')
        for n, v in out.items():
            print(n, np.round(v, 3), flush=True)
    else:
        grid = ((24., (0.50, 0.55, 0.60, 0.65)), (18., (0.60, 0.65, 0.70, 0.75)))
        target = RES/'exit_detector.json'
        extra = [x for x in sys.argv[1:] if '=' in x]          # e.g. 24=0.40,0.45 18=0.48,0.54
        if extra:
            grid = tuple((float(x.split('=')[0]), tuple(float(v) for v in x.split('=')[1].split(','))) for x in extra)
            target = RES/'exit_detector_low.json'
        out = dict(IA=IA.tolist(), curves={})
        for fc, amps in grid:
            for a in amps:
                c = detector_exit(fc, a)
                out['curves'][f'{fc:g}_{a:.3f}'] = c
                print(fc, a, {k: np.round(v, 3).tolist() for k, v in c.items()}, f'{time.time() - t0:.0f}s', flush=True)
                target.write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
