"""Receiver timing of the NB-LDPC GF(512) joint design on one CPU core (run on an idle machine).

For each packet (qualification settings, fresh seeds) the CPU detector is driven block by block as in
streaming reception; every push is timed (block-processing time). After the final flush, one decoding
pass (word-level forward-backward over the detector graph + 30 BP iterations) is timed, and then the
full 16 passes are timed without early stopping (worst case). Writes results/timing_nb.json.
Usage: python timing_nb.py --packets 8
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, json, time
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
import hop_packets as hp
from hop_receiver import HopLimitedReceiver, sampled_hop
import turbo_rx as tr
import nb_ldpc as nb
import nb_packets as npk
import qualify_nb as q

RES = Path(__file__).resolve().parent/'results'
SEEDS = (20371541, 20371542)
START = 9200000


def one(fc, amp, index):
    code = npk.get_code(q.CODE_SEED, q.DV)
    rng = np.random.default_rng(np.random.SeedSequence([SEEDS[0], npk.K, index]))
    payload = rng.integers(0, 256, (npk.K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + binascii.crc_hqx(body, 0xffff).to_bytes(2, 'big'), dtype=np.uint8))
    x = code.encode(npk.info_symbols(info))
    perm = np.random.default_rng(4100 + q.CODE_SEED).permutation(code.n)
    tones, prev = [], hp.M//2
    for w in x[perm]:
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    ts = 224/(len(tones) + 18)
    txs = np.r_[tones, np.full(6, hp.M//2)]
    r = CausalReceiver(carrier_hz=fc, background=mg.background('baseline', fc))
    det = HopLimitedReceiver(r, ts, amp*mg.A0, code.n, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                             h=hp.H, beta=hp.BETA, tone_count=hp.M, max_hop=hp.D, group_symbols=hp.L,
                             group_bits=hp.BITS, pad_bits=0)
    mean = r.stream(1).process(sampled_hop(txs, ts, hp.H, float(hp.BETA), hp.M))
    noise_rng = np.random.default_rng(np.random.SeedSequence([SEEDS[1], npk.K, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amp*mg.A0)
    tr.record_graph(det)
    position = 0; blocks = []
    for symbol in range(len(txs)):
        count = det._geometry(symbol)[3]
        t0 = time.perf_counter(); det.push(observed[position:position + count]); blocks.append(time.perf_counter() - t0)
        position += count
    t0 = time.perf_counter(); det.finish(); flush = time.perf_counter() - t0
    t0 = time.perf_counter(); packed = tr.pack_graph(det, hp.L, code.n); pack = time.perf_counter() - t0
    mask = np.zeros((code.n, nb.Q)); mask[code.info[28], (np.arange(nb.Q) & 31) != 0] = npk.FLOOR
    inv = np.argsort(perm)

    def passes(n, stop):
        Lp = mask[perm]; t0 = time.perf_counter(); ok_at = None
        for it in range(n):
            app = npk.normalize_log(npk.graph_word_app(packed, Lp))
            ext = npk.normalize_log(app - Lp)[inv] + mask
            prior = np.exp(ext - ext.max(1, keepdims=True)); prior /= prior.sum(1, keepdims=True)
            post, cext, hard, ok = code.decode(prior, q.BP)
            if ok and npk.crc_ok(npk.symbols_info(hard[code.info])) and ok_at is None:
                ok_at = it + 1
                if stop:
                    break
            Lp = npk.normalize_log(np.log(np.maximum(cext, 1e-300)) + mask)[perm]
        return time.perf_counter() - t0, ok_at
    first, ok1 = passes(1, True)
    full, okf = passes(16, False)
    return dict(carrier_hz=fc, amplitude_A0=amp, index=index, max_block_s=max(blocks), mean_block_s=float(np.mean(blocks)),
                flush_s=flush, pack_s=pack, first_pass_s=first, sixteen_pass_s=full, decoded_first_pass=ok1 == 1,
                graph_edges=int(len(packed['met'])))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--packets', type=int, default=8)
    a = p.parse_args()
    one(24., 0.58, START - 1)                                  # warm-up (numba compilation)
    rows = []
    for fc, amp in ((24., 0.58), (18., 0.697)):
        for i in range(START, START + a.packets):
            r = one(fc, amp, i); rows.append(r); print(r, flush=True)
    summ = dict(max_block_s=max(r['max_block_s'] for r in rows),
                max_first_pass_s=max(r['first_pass_s'] + r['flush_s'] + r['pack_s'] for r in rows),
                max_sixteen_pass_s=max(r['sixteen_pass_s'] + r['flush_s'] + r['pack_s'] for r in rows),
                median_first_pass_s=float(np.median([r['first_pass_s'] for r in rows])),
                median_sixteen_pass_s=float(np.median([r['sixteen_pass_s'] for r in rows])))
    print(summ)
    (RES/'timing_nb.json').write_text(json.dumps(dict(summary=summ, rows=rows), indent=1) + '\n')


if __name__ == '__main__':
    main()
