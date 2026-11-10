"""Independent qualification of the frozen joint design (hop-limited CF-FSK + K5 code + designed labeling,
iterative receiver), following the released procedure (App. B / C of the paper).

Frozen configuration (recorded with hashes in results/qual_<case>/config.json before any packet runs):
  waveform M=7, max hop 3, beta=1, h=3/14, T=224/190 s, 9 coded bits per 4 tones (hop_receiver.HopCode);
  labeling results/labels_hop.json ('optimized'); code K5 (23,35) tail-biting, punctured to 2/3 as released;
  released interleaver, header, CRC; HopLimitedReceiver (history 4, beam 32768) + turbo_rx, 16 iterations,
  decisions from the last iteration.
Fresh seeds: payload/noise namespaces 20371201/20371202, packet indices from 8000000.
Each independent packet is one bounded observation X_i = errors/224; a time-uniform betting confidence
sequence (gravcomm.packet_confidence, alpha = 0.025 per channel) is checked at 128-packet prefixes. Pass:
upper bound <= 1e-3; fail: lower bound > 1e-3; otherwise continue up to --max-packets.
Usage: python qualify.py --case 24=0.70 --case 18=0.84 [--max-packets 3072]
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, hashlib, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.packet_confidence import packet_ber_interval
import turbo_packets as tp

HERE = Path(__file__).resolve().parent
RES = HERE/'results'
SEEDS = (20371201, 20371202)
START = 8000000
CODE, LABEL, ITER = 'K5 (23,35)', 'optimized', 16


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze(folder, fc, amp):
    cfg = dict(waveform='hop', carrier_hz=fc, amplitude_A0=amp, code=CODE, labeling=LABEL, iterations=ITER,
               seeds=SEEDS, first_index=START, alpha=0.025, check_every=128,
               hashes={f: sha(HERE/f) for f in ('hop_receiver.py', 'hop_packets.py', 'turbo_rx.py', 'turbo_packets.py',
                                               'exit_chart.py', 'qualify.py', 'results/labels_hop.json')})
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
        e = [rows[START + i]['turbo_errors'][-1] for i in range(k)]
        lo, hi = packet_ber_interval(e, 224, 0.025)
        out = dict(packets=k, errors=int(sum(e)), ber=sum(e)/(224*k), interval=[lo, hi],
                   rejected=sum(not rows[START + i]['turbo_accepted'][-1] for i in range(k)),
                   pruned=sum(rows[START + i]['pruned_states'] for i in range(k)),
                   missing=sum(rows[START + i]['missing'] for i in range(k)))
        if hi <= 1e-3:
            return 'pass', out
        if lo > 1e-3:
            return 'fail', out
    return ('undecided' if n < max_packets else 'exhausted'), out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='carrier=amplitude_factor_of_A0')
    p.add_argument('--max-packets', type=int, default=3072)
    p.add_argument('--batch', type=int, default=256)
    p.add_argument('--workers', type=int, default=17)
    a = p.parse_args()
    cases = []
    for c in a.case:
        fc, amp = (float(x) for x in c.split('='))
        folder = RES/f'qual_hopk5opt_f{fc:g}_a{amp:.4f}'; folder.mkdir(parents=True, exist_ok=True)
        freeze(folder, fc, amp); cases.append((fc, amp, folder))
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as pool:
        while True:
            jobs = []
            for fc, amp, folder in cases:
                st, info = status(folder, a.max_packets)
                print(f'{folder.name}: {st} {info}', flush=True)
                if st != 'undecided':
                    continue
                have = {int(q.stem) for q in folder.glob('[0-9]*.json')}
                todo = [i for i in range(START, START + a.max_packets) if i not in have][:a.batch]
                jobs += [((('hop', fc, amp*mg.A0, i, ITER, SEEDS, CODE, LABEL)), folder) for i in todo]
            if not jobs:
                break
            futs = {pool.submit(tp.run, j): f for j, f in jobs}
            for fut in as_completed(futs):
                res = fut.result(); folder = futs[fut]
                tmp = folder/f"{res['index']}.json.tmp"
                tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")
            print(f'  batch done ({time.time() - t0:.0f}s)', flush=True)
    summary = {f.name: dict(zip(('status', 'result'), status(f, a.max_packets))) for _, _, f in cases}
    (RES/'qualification_joint.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
