"""Packet simulation of the hop-limited link with a non-binary LDPC code over GF(512).

One code symbol = one 4-tone hop-limited group (9 bits, enumerative word), so the detector and the
decoder exchange word (symbol) probabilities, not bit LLRs. Transmitter, channel, frame length and
seeding follow turbo_packets.py ('hop' waveform): 256 information bits (16-bit header, 224-bit
payload, CRC-16), here carried by 29 GF(512) information symbols (the last one with 5 known zero
bits) and 14 parity symbols: 43 groups, 172 data tones, the same 224 s frame.
Usage: python nb_packets.py --case nb_f24_a0.60:24:0.60 --packets 64 --iterations 16
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, json, math, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from numba import njit
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
import hop_packets as hp
from hop_receiver import HopLimitedReceiver, sampled_hop
import turbo_rx as tr
import nb_ldpc as nb
import word_tables as wt

OUT = Path(__file__).resolve().parent/'results'
K = 256
FLOOR = -40.0                                               # log-probability floor of exchanged messages


def crc_ok(bits):
    data = np.packbits(bits)
    return bool(bytes(data[:2]) == b'\x00\x01' and
                binascii.crc_hqx(bytes(data[:-2]), 0xffff) == int.from_bytes(bytes(data[-2:]), 'big'))


@njit(cache=True)
def _graph_word_app(src, dst, met, word, eoff, noff, group, Lp):
    """Forward-backward over the detector graph with word log-priors Lp (G, W); word log-APPs (G, W)."""
    T = len(eoff) - 1; G, W = Lp.shape
    alpha = np.full(noff[-1], -np.inf); beta = np.full(noff[-1], -np.inf)
    alpha[noff[0]] = 0.0
    prior = np.zeros(len(met))
    for t in range(T):
        g = group[t]
        if g >= 0:
            for e in range(eoff[t], eoff[t + 1]):
                prior[e] = Lp[g, word[e]]
    for t in range(T):
        for e in range(eoff[t], eoff[t + 1]):
            i = noff[t + 1] + dst[e]
            alpha[i] = tr._lse(alpha[i], alpha[noff[t] + src[e]] + met[e] + prior[e])
        m = -np.inf
        for i in range(noff[t + 1], noff[t + 2]):
            m = max(m, alpha[i])
        for i in range(noff[t + 1], noff[t + 2]):
            alpha[i] -= m
    for i in range(noff[T], noff[T + 1]):
        beta[i] = 0.0
    for t in range(T - 1, -1, -1):
        for e in range(eoff[t], eoff[t + 1]):
            i = noff[t] + src[e]
            beta[i] = tr._lse(beta[i], met[e] + prior[e] + beta[noff[t + 1] + dst[e]])
        m = -np.inf
        for i in range(noff[t], noff[t + 1]):
            m = max(m, beta[i])
        for i in range(noff[t], noff[t + 1]):
            beta[i] -= m
    post = np.full((G, W), -np.inf)
    for t in range(T):
        g = group[t]
        if g < 0:
            continue
        for e in range(eoff[t], eoff[t + 1]):
            v = alpha[noff[t] + src[e]] + met[e] + prior[e] + beta[noff[t + 1] + dst[e]]
            post[g, word[e]] = tr._lse(post[g, word[e]], v)
    return post


def graph_word_app(packed, Lp):
    return _graph_word_app(packed['src'], packed['dst'], packed['met'], packed['word'], packed['eoff'],
                           packed['noff'], packed['group'], np.ascontiguousarray(Lp, np.float64))


def normalize_log(L):
    L = L - L.max(1, keepdims=True)
    L = L - np.log(np.exp(L).sum(1, keepdims=True))
    return np.maximum(L, FLOOR)


_CODES = {}
_TABLES = {}


def get_table(name):
    if name not in _TABLES:
        _TABLES[name] = wt.TableHopCode(json.loads((OUT/f'words_{name}.json').read_text()))
    return _TABLES[name]


def get_code(seed, dv=2):
    key = (seed, dv)
    if key not in _CODES:
        _CODES[key] = nb.make_code(43, 29, dv, seed)
    return _CODES[key]


def info_symbols(info):
    bits = np.r_[info, np.zeros(29*9 - len(info), int)]
    return bits.reshape(29, 9) @ (1 << np.arange(8, -1, -1))


def symbols_info(u):
    return ((np.asarray(u)[:, None] >> np.arange(8, -1, -1)) & 1).ravel()[:K].astype(np.uint8)


def run(job, return_graph=False):
    carrier, amplitude, index, iterations, seed_ns, code_seed, bp_iters = job[:7]
    dv = job[7] if len(job) > 7 else 2
    words = job[8] if len(job) > 8 else 'lex'
    device = job[9] if len(job) > 9 else 'cpu'
    table = None if words == 'lex' else get_table(words)
    enc = hp.CODE if table is None else table
    code = get_code(code_seed, dv)
    start = time.perf_counter()
    rng = np.random.default_rng(np.random.SeedSequence([seed_ns[0], K, index]))
    payload = rng.integers(0, 256, (K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + binascii.crc_hqx(body, 0xffff).to_bytes(2, 'big'), dtype=np.uint8))
    x = code.encode(info_symbols(info))
    perm = np.random.default_rng(4100 + code_seed).permutation(code.n)     # channel group i carries x[perm[i]]
    ranks = x[perm]
    tones, prev = [], hp.M//2
    for w in ranks:
        s = enc.encode(int(w), prev); tones += s; prev = s[-1]
    ngroups = len(ranks)
    ts = 224/(len(tones) + 18)
    txs = np.r_[tones, np.full(6, hp.M//2)]
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    Det = HopLimitedReceiver
    if device == 'gpu':
        from hop_receiver_gpu import HopLimitedReceiverGPU as Det, pack_graph_gpu, GraphGPU
    det = Det(r, ts, amplitude, ngroups, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                             h=hp.H, beta=hp.BETA, tone_count=hp.M, max_hop=hp.D, group_symbols=hp.L,
                             group_bits=hp.BITS, pad_bits=0, word_table=table)
    mean = r.stream(1).process(sampled_hop(txs, ts, hp.H, float(hp.BETA), hp.M))
    noise_rng = np.random.default_rng(np.random.SeedSequence([seed_ns[1], K, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    if device != 'gpu':
        tr.record_graph(det)
    position = 0
    for symbol in range(len(txs)):
        count = det._geometry(symbol)[3]
        det.push(observed[position:position + count]); position += count
    stats = det.finish()
    t_detect = time.perf_counter() - start
    if device == 'gpu':
        gg = GraphGPU(det, hp.L, ngroups)
        app_fn = gg.word_app
        packed = pack_graph_gpu(det, hp.L, ngroups) if return_graph else None
        n_edges = int(sum(int(m.shape[0]) for m in gg.met))
    else:
        packed = tr.pack_graph(det, hp.L, ngroups)
        app_fn = lambda Lp: graph_word_app(packed, Lp)
        n_edges = int(len(packed['met']))
    # known zero bits of the last information symbol (code position code.info[28])
    mask = np.zeros((code.n, nb.Q))
    mask[code.info[28], (np.arange(nb.Q) & 31) != 0] = FLOOR
    inv = np.argsort(perm)
    Lp = mask[perm]                                          # detector word priors, channel order
    errors, accepted, symerr = [], [], []
    for it in range(iterations):
        app = normalize_log(app_fn(Lp))
        ext = normalize_log(app - Lp)[inv] + mask           # code order
        prior = np.exp(ext - ext.max(1, keepdims=True))
        prior /= prior.sum(1, keepdims=True)
        post, cext, hard, ok = code.decode(prior, bp_iters)
        bits = symbols_info(hard[code.info])
        errors.append(int(np.sum(bits[16:-16] != info[16:-16])))
        accepted.append(bool(ok and crc_ok(bits)))
        symerr.append(int(np.sum(hard != x)))
        if accepted[-1]:
            break
        Lp = normalize_log(np.log(np.maximum(cext, 1e-300)) + mask)[perm]
    out = dict(code=f'NB-LDPC GF(512) (43,29) dv={dv} seed={code_seed}', words=words, device=device, carrier_hz=carrier,
                amplitude_m_s2=amplitude, index=index, errors=errors, accepted=accepted, symbol_errors=symerr,
                pruned_states=stats['pruned_states'], graph_edges=n_edges,
                detect_s=t_detect, elapsed_s=time.perf_counter() - start)
    return (out, packed) if return_graph else out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='name:carrier:amplitude_factor_of_A0')
    p.add_argument('--packets', type=int, default=64)
    p.add_argument('--seed-start', type=int, default=7000000)
    p.add_argument('--seeds', default='20371103,20371104', help='payload,noise seed namespaces')
    p.add_argument('--iterations', type=int, default=16)
    p.add_argument('--bp', type=int, default=30)
    p.add_argument('--code-seed', type=int, default=1)
    p.add_argument('--dv', type=int, default=2)
    p.add_argument('--words', default='lex')
    p.add_argument('--workers', type=int, default=8, help='CPU detector workers')
    p.add_argument('--gpu-workers', type=int, default=0, help='GPU detector workers (CuPy)')
    a = p.parse_args()
    ns = tuple(int(x) for x in a.seeds.split(','))
    jobs = []
    for case in a.case:
        name, carrier, factor = case.split(':')
        folder = OUT/name; folder.mkdir(parents=True, exist_ok=True)
        for i in range(a.seed_start, a.seed_start + a.packets):
            if not (folder/f'{i}.json').exists():
                jobs.append(((float(carrier), float(factor)*mg.A0, i, a.iterations, ns, a.code_seed, a.bp, a.dv, a.words), folder))
    print(f'{len(jobs)} packets to run', flush=True)
    t0 = time.time(); done = 0
    for res, folder in run_hybrid(jobs, a.workers, a.gpu_workers):
        tmp = folder/f"{res['index']}.json.tmp"
        tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")
        done += 1
        if done % 8 == 0:
            print(f'  {done}/{len(jobs)} ({time.time() - t0:.0f}s)', flush=True)


def run_hybrid(jobs, cpu_workers, gpu_workers, fn=None):
    """Run (job, tag) pairs on a CPU pool and a GPU pool fed from one queue; yields (result, tag).
    The device is appended to each job tuple (index 9); jobs must have exactly 9 fields."""
    from concurrent.futures import wait, FIRST_COMPLETED
    fn = fn or run
    queue = list(jobs)[::-1]
    pools = []
    if cpu_workers:
        pools.append(('cpu', ProcessPoolExecutor(max_workers=cpu_workers), cpu_workers))
    if gpu_workers:
        pools.append(('gpu', ProcessPoolExecutor(max_workers=gpu_workers), gpu_workers))
    running = {}
    try:
        def fill():
            for dev, pool, cap in pools:
                while queue and sum(1 for v in running.values() if v[0] == dev) < cap:
                    job, tag = queue.pop()
                    running[pool.submit(fn, tuple(job[:9]) + (dev,))] = (dev, tag)
        fill()
        while running:
            finished, _ = wait(list(running), return_when=FIRST_COMPLETED)
            for f in finished:
                dev, tag = running.pop(f)
                yield f.result(), tag
            fill()
    finally:
        for _, pool, _ in pools:
            pool.shutdown(wait=True)

if __name__ == '__main__':
    main()
