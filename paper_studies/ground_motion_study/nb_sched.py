"""Compare detector/decoder schedules of the GF(512) iterative receiver on identical detector graphs.

For each packet the transmitter, channel and streaming detector of nb_packets.py are run once; every
schedule then iterates on the same recorded graph. A schedule is (bp, iters, damp, warm):
bp belief-propagation iterations per detector pass, iters detector passes, damp the weight of the
new decoder message in the detector prior (log domain), warm keeps the check-to-variable messages
between passes. Writes results/nb_sched_<tag>.json.
Usage: python nb_sched.py --carrier 24 --amp 0.50 --indices 7200003,7200014 --tag f24_a0.50
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
import hop_packets as hp
from hop_receiver import HopLimitedReceiver, sampled_hop
import turbo_rx as tr
import nb_ldpc as nb
import nb_packets as npk

OUT = Path(__file__).resolve().parent/'results'
SCHEDULES_OLD = {
    'S0 bp30 x16': (30, 16, 1.0, False),
    'S1 bp3 x50 warm': (3, 50, 1.0, True),
    'S2 bp3 x50 warm d.5': (3, 50, 0.5, True),
    'S3 bp1 x60 warm d.7': (1, 60, 0.7, True),
    'S4 bp10 x30 d.5': (10, 30, 0.5, False),
    'S5 bp30 x30 d.5': (30, 30, 0.5, False),
}
SCHEDULES = {
    'S0 bp30 x16': (30, 16, 1.0, False),
    'C1 chase S8 L2 app': ('chase', dict(S=8, Lalt=2, rank_by='app')),
    'C2 chase S5 L3 app': ('chase', dict(S=5, Lalt=3, rank_by='app')),
    'C3 chase S10 L2 app': ('chase', dict(S=10, Lalt=2, rank_by='app')),
}


def detect(carrier, amplitude, index, seed_ns=(20371103, 20371104), code_seed=1):
    code = npk.get_code(code_seed)
    rng = np.random.default_rng(np.random.SeedSequence([seed_ns[0], npk.K, index]))
    payload = rng.integers(0, 256, (npk.K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + binascii.crc_hqx(body, 0xffff).to_bytes(2, 'big'), dtype=np.uint8))
    x = code.encode(npk.info_symbols(info))
    perm = np.random.default_rng(4100 + code_seed).permutation(code.n)
    ranks = x[perm]; tones, prev = [], hp.M//2
    for w in ranks:
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    ts = 224/(len(tones) + 18); txs = np.r_[tones, np.full(6, hp.M//2)]
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    det = HopLimitedReceiver(r, ts, amplitude, len(ranks), history_symbols=4, beam_width=32768, lookahead_symbols=12,
                             h=hp.H, beta=hp.BETA, tone_count=hp.M, max_hop=hp.D, group_symbols=hp.L,
                             group_bits=hp.BITS, pad_bits=0)
    mean = r.stream(1).process(sampled_hop(txs, ts, hp.H, float(hp.BETA), hp.M))
    noise_rng = np.random.default_rng(np.random.SeedSequence([seed_ns[1], npk.K, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    tr.record_graph(det); position = 0
    for symbol in range(len(txs)):
        count = det._geometry(symbol)[3]
        det.push(observed[position:position + count]); position += count
    det.finish()
    return code, info, x, perm, tr.pack_graph(det, hp.L, len(ranks))


def iterate(code, info, x, perm, packed, bp, iters, damp, warm):
    mask = np.zeros((code.n, nb.Q)); mask[code.info[28], (np.arange(nb.Q) & 31) != 0] = npk.FLOOR
    inv = np.argsort(perm); Lp = mask[perm]; state = True if warm else None
    for it in range(iters):
        app = npk.normalize_log(npk.graph_word_app(packed, Lp))
        ext = npk.normalize_log(app - Lp)[inv] + mask
        prior = np.exp(ext - ext.max(1, keepdims=True)); prior /= prior.sum(1, keepdims=True)
        if warm:
            post, cext, hard, ok, state = code.decode(prior, bp, state)
        else:
            post, cext, hard, ok = code.decode(prior, bp)
        bits = npk.symbols_info(hard[code.info])
        if ok and npk.crc_ok(bits):
            return True, it + 1, int(np.sum(bits[16:-16] != info[16:-16]))
        new = npk.normalize_log(np.log(np.maximum(cext, 1e-300)) + mask)[perm]
        Lp = new if (damp >= 1 or it == 0) else npk.normalize_log(damp*new + (1 - damp)*Lp)
    return False, iters, int(np.sum(bits[16:-16] != info[16:-16]))


def chase(code, info, x, perm, packed, S=8, Lalt=2, iters=16, bp=30, rank_by='app'):
    """Turbo loop (S0); on failure branch over the top-Lalt values of the S least reliable symbols,
    clamping them in BP on the last detector extrinsic; pick the CRC-passing branch of best metric."""
    import itertools
    mask = np.zeros((code.n, nb.Q)); mask[code.info[28], (np.arange(nb.Q) & 31) != 0] = npk.FLOOR
    inv = np.argsort(perm); Lp = mask[perm]
    for it in range(iters):
        app = npk.normalize_log(npk.graph_word_app(packed, Lp))
        ext = npk.normalize_log(app - Lp)[inv] + mask
        prior = np.exp(ext - ext.max(1, keepdims=True)); prior /= prior.sum(1, keepdims=True)
        post, cext, hard, ok = code.decode(prior, bp)
        bits = npk.symbols_info(hard[code.info])
        if ok and npk.crc_ok(bits):
            return True, it + 1, 0, 0
        Lp = npk.normalize_log(np.log(np.maximum(cext, 1e-300)) + mask)[perm]
    # branch on the final detector extrinsic (prior) using full detector APP for reliability
    rel_src = np.exp(app[inv]) if rank_by == 'app' else prior
    order = np.argsort(rel_src.max(1))[:S]
    tops = [np.argsort(-rel_src[j])[:Lalt] for j in order]
    lp_app = app[inv]
    best = None; tried = 0
    for combo in itertools.product(range(Lalt), repeat=S):
        pr = prior.copy()
        for j, c, t in zip(order, combo, tops):
            pr[j] = 0.0; pr[j, t[c]] = 1.0
        post2, _, hard2, ok2 = code.decode(pr, bp); tried += 1
        if ok2:
            b2 = npk.symbols_info(hard2[code.info])
            if npk.crc_ok(b2):
                m = float(np.sum(lp_app[np.arange(code.n), hard2]))
                if best is None or m > best[0]:
                    best = (m, b2)
    if best is None:
        return False, iters, tried, -1
    return bool(np.array_equal(best[1], info)), iters, tried, int(np.sum(best[1][16:-16] != info[16:-16]))


def job(args):
    carrier, amp, index = args
    t0 = time.perf_counter()
    d = detect(carrier, amp*mg.A0, index)
    res = {}
    for name, sch in SCHEDULES.items():
        t1 = time.perf_counter()
        if sch[0] == 'chase':
            ok, its, tried, err = chase(*d, **sch[1])
        else:
            ok, its, err = iterate(*d, *sch)
        res[name] = dict(ok=ok, iters=its, errors=err, seconds=time.perf_counter() - t1)
    return index, res, time.perf_counter() - t0


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--carrier', type=float, required=True); p.add_argument('--amp', type=float, required=True)
    p.add_argument('--indices', required=True); p.add_argument('--tag', required=True)
    p.add_argument('--workers', type=int, default=18)
    a = p.parse_args()
    idx = [int(v) for v in a.indices.split(',')]
    out = {}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(job, (a.carrier, a.amp, i)) for i in idx]):
            i, res, el = f.result(); out[i] = res
    summary = {name: dict(ok=sum(out[i][name]['ok'] for i in out), n=len(out),
                          mean_iters=float(np.mean([out[i][name]['iters'] for i in out if out[i][name]['ok']] or [0])), wrong=sum(out[i][name]['errors'] > 0 for i in out),
                          mean_s=float(np.mean([out[i][name]['seconds'] for i in out])))
               for name in SCHEDULES}
    for k, v in summary.items():
        print(f'{k:22s} ok {v["ok"]}/{v["n"]} undetected {v["wrong"]}  iters {v["mean_iters"]:.1f}  {v["mean_s"]:.1f}s', flush=True)
    (OUT/f'nb_sched_{a.tag}.json').write_text(json.dumps(dict(summary=summary, packets=out), indent=1))


if __name__ == '__main__':
    main()
