"""Packet simulation of the hop-limited CF-FSK link (requalification candidate).

Same frame, code, CRC, interleaver and receiver front end as the released link (run_packets.py);
the 384 coded bits plus 3 known pad bits are cut into 43 groups of 9 bits, each mapped to 4
hop-limited tones (HopCode), giving 172 data symbols, 6 centre-tone tail symbols and 12 budgeted
acquisition symbols in the 224 s frame (T = 224/190 s).
Each packet is written to results/<case>/<index>.json. Usage:
  python hop_packets.py --case hop7_f24:24:0.8125 --packets 64 --seed-start 5000000 --workers 10
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
from gravcomm import convolutional_candidates as tb
from hop_receiver import HopCode, HopLimitedReceiver, sampled_hop

OUT = Path(__file__).resolve().parent/'results'
K, SCHEME = 256, 'tb_2/3'
M, D, L, BITS, H, BETA = 7, 3, 4, 9, Fraction(3, 14), Fraction(1)
CODE = HopCode(M, D, L, BITS)


def frame(word):
    bits = np.r_[word, np.zeros((-len(word)) % BITS, int)]
    groups = bits.reshape(-1, BITS) @ (1 << np.arange(BITS - 1, -1, -1))
    tones, prev = [], M//2
    for g in groups:
        s = CODE.encode(int(g), prev); tones += s; prev = s[-1]
    return np.array(tones), len(groups), len(bits) - len(word)


def run(job):
    carrier, amplitude, index, history, beam, lookahead = job
    start = time.perf_counter()
    rng = np.random.default_rng(np.random.SeedSequence([20371103, K, index]))
    payload = rng.integers(0, 256, (K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + binascii.crc_hqx(body, 0xffff).to_bytes(2, 'big'), dtype=np.uint8))
    word = tb.encode(info, SCHEME[3:])
    permutation = np.random.default_rng(4000 + len(word)).permutation(len(word))
    tones, ngroups, pad = frame(word[permutation])
    ts = 224/(len(tones) + 18)
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    tx = np.r_[tones, np.full(6, M//2)]
    mean = r.stream(1).process(sampled_hop(tx, ts, H, float(BETA), M))
    noise_rng = np.random.default_rng(np.random.SeedSequence([20371104, K, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    d = HopLimitedReceiver(r, ts, amplitude, ngroups, history_symbols=history, beam_width=beam,
                           lookahead_symbols=lookahead, h=H, beta=BETA, tone_count=M, max_hop=D,
                           group_symbols=L, group_bits=BITS, pad_bits=pad)
    outputs = []; position = 0
    for symbol in range(len(tx)):
        count = d._geometry(symbol)[3]
        outputs.extend(d.push(observed[position:position + count])); position += count
    stats = d.finish()
    channel = np.concatenate([x['llr'] for x in outputs])[:len(word)]
    llr = np.empty(len(word)); llr[permutation] = channel
    recovered = tb.decode(llr, SCHEME[3:], K); data = np.packbits(recovered)
    accepted = bool(bytes(data[:2]) == b'\x00\x01' and
                    binascii.crc_hqx(bytes(data[:-2]), 0xffff) == int.from_bytes(bytes(data[-2:]), 'big'))
    errors = int(np.unpackbits(data[2:-2] ^ payload).sum())
    hard_tone_errors = int(sum(o['byte'] != CODE.decode(tones[L*o['byte_index']:L*o['byte_index'] + L],
                                                       M//2 if o['byte_index'] == 0 else tones[L*o['byte_index'] - 1])
                               for o in outputs))
    return dict(carrier_hz=carrier, amplitude_m_s2=amplitude, index=index, errors=errors, accepted=accepted,
                group_errors=hard_tone_errors, pruned_states=stats['pruned_states'], peak_states=stats['peak_states'],
                missing=stats['missing_bit_hypotheses'], history=history, beam=beam, lookahead=lookahead,
                whitening=d.whitener.spectral_error()['max_psd_error'], elapsed_s=time.perf_counter() - start)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='name:carrier:amplitude_factor_of_A0')
    p.add_argument('--packets', type=int, default=64)
    p.add_argument('--seed-start', type=int, default=5000000)
    p.add_argument('--workers', type=int, default=10)
    p.add_argument('--history', type=int, default=4)
    p.add_argument('--beam', type=int, default=32768)
    p.add_argument('--lookahead', type=int, default=12)
    a = p.parse_args()
    jobs = []
    for case in a.case:
        name, carrier, factor = case.split(':')
        folder = OUT/name; folder.mkdir(parents=True, exist_ok=True)
        for i in range(a.seed_start, a.seed_start + a.packets):
            if not (folder/f'{i}.json').exists():
                jobs.append(((float(carrier), float(factor)*mg.A0, i, a.history, a.beam, a.lookahead), folder))
    print(f'{len(jobs)} packets to run', flush=True)
    t0 = time.time(); done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures = {pool.submit(run, j): folder for j, folder in jobs}
        for f in as_completed(futures):
            res = f.result(); folder = futures[f]
            tmp = folder/f"{res['index']}.json.tmp"
            tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")
            done += 1
            if done % 32 == 0:
                print(f'  {done}/{len(jobs)} ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
