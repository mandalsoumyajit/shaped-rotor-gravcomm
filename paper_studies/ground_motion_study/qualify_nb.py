"""Independent qualification of the non-binary joint design (hop-limited CF-FSK + GF(512) LDPC (43,29),
symbol-level iterative receiver), following the released procedure (App. B / C of the paper).

Frozen configuration (recorded with hashes in results/qualnb_<case>/config.json before any packet runs):
  waveform M=7, max hop 3, beta=1, h=3/14, T=224/190 s, one GF(512) code symbol per 4-tone group
  (lexicographic enumerative words, hop_receiver.HopCode); code nb_ldpc.make_code(43, 29, dv=2, seed=1),
  random symbol interleaver (seed 4101), header and CRC as released; HopLimitedReceiver (history 4,
  beam 32768) + word-level forward-backward and FFT-BP (30 BP iterations per pass), 16 detector passes,
  early stop on syndrome + CRC, decisions from the last pass.
Fresh seeds: payload/noise namespaces 20371501/20371502, packet indices from 8500000.
Each packet is one bounded observation X_i = errors/224; a time-uniform betting confidence sequence
(gravcomm.packet_confidence, alpha = 0.025 per channel) is checked at 128-packet prefixes. Pass: upper
bound <= 1e-3; fail: lower bound > 1e-3; otherwise continue up to --max-packets.
Attempt 1 (0.55 / 0.66, seeds 20371501/02) stopped undecided with BER 1.25e-3 / 1.39e-3 after 6656 packets
(results/qualification_nb_attempt1.log); attempt 2 steps up with fresh seeds.
Usage: python qualify_nb.py --case 24=0.59 --case 18=0.71 --seeds 20371511,20371512 --start 8600000
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, hashlib, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import mitigation as mg
from gravcomm.packet_confidence import packet_ber_interval
import nb_packets as npk

HERE = Path(__file__).resolve().parent
RES = HERE/'results'
SEEDS = (20371501, 20371502)
START = 8500000
ITER, BP, CODE_SEED, DV, WORDS = 16, 30, 1, 2, 'lex'
FILES = ('hop_receiver.py', 'hop_packets.py', 'turbo_rx.py', 'nb_packets.py', 'nb_ldpc.py', 'word_tables.py',
         'qualify_nb.py', 'hop_receiver_gpu.py')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze(folder, fc, amp):
    cfg = dict(waveform='hop', carrier_hz=fc, amplitude_A0=amp, code='NB-LDPC GF(512) (43,29) dv=2 seed=1',
               words=WORDS, iterations=ITER, bp_iterations=BP, seeds=SEEDS, first_index=START, alpha=0.025,
               check_every=128, hashes={f: sha(HERE/f) for f in FILES})
    p = folder/'config.json'
    if p.exists():
        old = json.loads(p.read_text())
        if old['hashes'] != cfg['hashes']:
            raise RuntimeError(f'configuration changed since qualification of {folder.name} started')
    else:
        p.write_text(json.dumps(cfg, indent=1) + '\n')


def status(folder, max_packets):
    rows = {int(p.stem): json.loads(p.read_text()) for p in folder.glob('[0-9]*.json')}
    n = 0
    while START + n in rows:
        n += 1
    out = None
    for k in range(128, n + 1, 128):
        e = [rows[START + i]['errors'][-1] for i in range(k)]
        lo, hi = packet_ber_interval(e, 224, 0.025)
        out = dict(packets=k, errors=int(sum(e)), ber=sum(e)/(224*k), interval=[lo, hi],
                   rejected=sum(not rows[START + i]['accepted'][-1] for i in range(k)),
                   pruned=sum(rows[START + i]['pruned_states'] for i in range(k)))
        if hi <= 1e-3:
            return 'pass', out
        if lo > 1e-3:
            return 'fail', out
    return ('undecided' if n < max_packets else 'exhausted'), out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='carrier=amplitude_factor_of_A0')
    p.add_argument('--max-packets', type=int, default=16896)
    p.add_argument('--batch', type=int, default=256)
    p.add_argument('--workers', type=int, default=14, help='CPU detector workers')
    p.add_argument('--gpu-workers', type=int, default=4, help='GPU detector workers (CuPy)')
    p.add_argument('--seeds', default='20371501,20371502', help='payload,noise namespaces (fresh per attempt)')
    p.add_argument('--start', type=int, default=8500000)
    a = p.parse_args()
    global SEEDS, START
    SEEDS = tuple(int(x) for x in a.seeds.split(',')); START = a.start
    cases = []
    for c in a.case:
        fc, amp = (float(x) for x in c.split('='))
        folder = RES/f'qualnb_f{fc:g}_a{amp:.4f}'; folder.mkdir(parents=True, exist_ok=True)
        freeze(folder, fc, amp); cases.append((fc, amp, folder))
    t0 = time.time()
    while True:
        jobs = []
        for fc, amp, folder in cases:
            st, info = status(folder, a.max_packets)
            print(f'{folder.name}: {st} {info}', flush=True)
            if st != 'undecided':
                continue
            have = {int(q.stem) for q in folder.glob('[0-9]*.json')}
            todo = [i for i in range(START, START + a.max_packets) if i not in have][:a.batch]
            jobs += [((fc, amp*mg.A0, i, ITER, SEEDS, CODE_SEED, BP, DV, WORDS), folder) for i in todo]
        if not jobs:
            break
        # CPU and GPU detectors validated identical (gpu_check.py); the device is recorded per packet
        for res, folder in npk.run_hybrid(jobs, a.workers, a.gpu_workers):
            tmp = folder/f"{res['index']}.json.tmp"
            tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")
        print(f'  batch done ({time.time() - t0:.0f}s)', flush=True)
    summary = {f.name: dict(zip(('status', 'result'), status(f, a.max_packets))) for _, _, f in cases}
    (RES/'qualification_nb.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
