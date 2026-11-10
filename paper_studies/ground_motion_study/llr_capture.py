"""Capture the released detector's soft bit outputs to split the coding gap.

Runs the released seven-tone link (baseline background) exactly as run_packets.py, but stores
the deinterleaved channel LLRs and the transmitted coded bits of every packet. From these,
gap_analysis.py computes the generalized mutual information (GMI) of the detector output,
i.e. the rate an ideal binary code could achieve on the actual soft metrics.
Seeds use a new namespace (20371105/6). Output: results/llr_capture/f<fc>_a<factor>.npz
Usage: python llr_capture.py --workers 8 [--packets 48]
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver, sampled_cffsk
from gravcomm.streaming_soft import StreamingSoftReceiver
from gravcomm import convolutional_candidates as tb

OUT = Path(__file__).resolve().parent/'results'/'llr_capture'
K, SCHEME, Q, BETA = 256, 'tb_2/3', .08, .6
CARRIERS = {24.: .8125, 18.: 1.0}
FACTORS = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


def crc(x):
    return binascii.crc_hqx(bytes(x), 0xffff)


def run(job):
    carrier, amplitude, index = job
    rng = np.random.default_rng(np.random.SeedSequence([20371105, K, index]))
    payload = rng.integers(0, 256, (K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + crc(body).to_bytes(2, 'big'), dtype=np.uint8))
    word = tb.encode(info, SCHEME[3:])
    permutation = np.random.default_rng(4000 + len(word)).permutation(len(word))
    packed = np.packbits(word[permutation]).astype(int)
    symbols = np.column_stack((packed//49, packed//7 % 7, packed % 7)).ravel() - 3
    ts = 224/(3*len(packed) + 18)
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    mean = r.stream(1).process(sampled_cffsk(np.r_[symbols, np.zeros(6, int)], ts, q=Q, transition_fraction=BETA))
    noise_rng = np.random.default_rng(np.random.SeedSequence([20371106, K, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    d = StreamingSoftReceiver(r, ts, amplitude, len(packed), history_symbols=4, beam_width=32768,
                              lookahead_symbols=12, q=Q, transition_fraction=BETA)
    outputs = []; position = 0
    for symbol in range(len(symbols) + 6):
        count = d._geometry(symbol)[3]
        outputs.extend(d.push(observed[position:position + count])); position += count
    stats = d.finish()
    channel = np.concatenate([x['llr'] for x in outputs])[:len(word)]
    llr = np.empty(len(word)); llr[permutation] = channel
    recovered = tb.decode(llr, SCHEME[3:], K); data = np.packbits(recovered)
    errors = int(np.unpackbits(data[2:-2] ^ payload).sum())
    return llr, word.astype(np.int8), errors, stats['pruned_states'], stats['missing_bit_hypotheses']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--packets', type=int, default=48)
    a = p.parse_args(); OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for fc, a_op in CARRIERS.items():
            for f in FACTORS:
                path = OUT/f'f{fc:g}_a{f:.2f}.npz'
                if path.exists():
                    continue
                t0 = time.time(); amp = f*a_op*mg.A0
                res = list(pool.map(run, [(fc, amp, i) for i in range(a.packets)]))
                np.savez(path, llr=np.array([r[0] for r in res]), word=np.array([r[1] for r in res]),
                         errors=np.array([r[2] for r in res]), pruned=np.array([r[3] for r in res]),
                         missing=np.array([r[4] for r in res]), amp_factor_of_op=f, amp_A0=f*a_op)
                print(f'{fc} Hz x{f:.2f}: errors {sum(r[2] for r in res)} in {224*a.packets} bits '
                      f'({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
