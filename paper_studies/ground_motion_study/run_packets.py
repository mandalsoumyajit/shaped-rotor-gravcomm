"""Packet screens of the released seven-tone rate-2/3 link under mitigation scenarios.

Follows streaming_dynamics_runtime.run (v0.5.0) with the background made a
parameter. Seeds use a new namespace; payload and noise streams are common to
all scenarios so that comparisons use common random numbers. These are
screening runs, not independent qualification.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver, sampled_cffsk
from gravcomm.streaming_soft import StreamingSoftReceiver
from gravcomm import convolutional_candidates as tb

OUT = Path(__file__).resolve().parent/'results'
K, SCHEME, Q, BETA = 256, 'tb_2/3', .08, .6


def crc(x):
    return binascii.crc_hqx(bytes(x), 0xffff)


def run(job):
    scenario, carrier, amplitude, index = job
    start = time.perf_counter()
    rng = np.random.default_rng(np.random.SeedSequence([20371103, K, index]))
    payload = rng.integers(0, 256, (K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + crc(body).to_bytes(2, 'big'), dtype=np.uint8))
    word = tb.encode(info, SCHEME[3:])
    permutation = np.random.default_rng(4000 + len(word)).permutation(len(word))
    packed = np.packbits(word[permutation]).astype(int)
    symbols = np.column_stack((packed//49, packed//7 % 7, packed % 7)).ravel() - 3
    ts = 224/(3*len(packed) + 18)
    r = CausalReceiver(carrier_hz=carrier, background=mg.background(scenario, carrier))
    mean = r.stream(1).process(sampled_cffsk(np.r_[symbols, np.zeros(6, int)], ts, q=Q, transition_fraction=BETA))
    noise_rng = np.random.default_rng(np.random.SeedSequence([20371104, K, index]))
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
    accepted = bool(bytes(data[:2]) == b'\x00\x01' and crc(data[:-2]) == int.from_bytes(bytes(data[-2:]), 'big'))
    errors = int(np.unpackbits(data[2:-2] ^ payload).sum())
    return dict(scenario=scenario, carrier_hz=carrier, amplitude_m_s2=amplitude, index=index, errors=errors,
                accepted=accepted, pruned_states=stats['pruned_states'],
                missing=stats['missing_bit_hypotheses'], whitening=d.whitener.spectral_error()['max_psd_error'],
                elapsed_s=time.perf_counter() - start)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='scenario:carrier:amplitude_factor_of_A0')
    p.add_argument('--packets', type=int, default=256)
    p.add_argument('--seed-start', type=int, default=3000000)
    p.add_argument('--workers', type=int, default=20)
    a = p.parse_args(); OUT.mkdir(exist_ok=True)
    jobs = []
    for case in a.case:
        name, carrier, factor = case.split(':')
        folder = OUT/f'{name}_f{float(carrier):g}_a{float(factor):.4f}'; folder.mkdir(exist_ok=True)
        for i in range(a.seed_start, a.seed_start + a.packets):
            if not (folder/f'{i}.json').exists():
                jobs.append((name, float(carrier), float(factor)*mg.A0, i, folder))
    print(f'{len(jobs)} packets to run', flush=True)
    done = 0; t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures = {pool.submit(run, j[:4]): j for j in jobs}
        for f in as_completed(futures):
            x = f.result(); folder = futures[f][4]
            (folder/f"{x['index']}.json").write_text(json.dumps(x) + '\n')
            done += 1
            if done % 50 == 0:
                print(f'{done}/{len(jobs)} {time.time() - t0:.0f}s', flush=True)
    print('DONE', flush=True)


if __name__ == '__main__':
    main()
