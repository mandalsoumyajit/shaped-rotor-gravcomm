"""Packet checks of the frozen joint design under other isolation scenarios or amplitudes.

Uses turbo_packets.run unchanged (its hash is frozen by qualify.py); the isolation scenario is
switched by replacing mitigation.background inside the worker, so 'baseline' maps to the chosen
scenario. Same code, labeling, iterations as the qualification; fresh seed namespaces.
Usage: python validate_hop.py --case three_stage_hop_f24_a0.5956=three_stage:24:0.5956 --packets 3000
       -> results/<name>/<index>.json (resumable); summary via --summary
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import mitigation as mg
import turbo_packets as tp
import qualify as q

RES = Path(__file__).resolve().parent/'results'
SEEDS = (20371301, 20371302)
START = 9000000
_orig = mg.background


def run(args):
    scenario, fc, amp, index = args
    mg.background = lambda name, carrier: _orig(scenario if name == 'baseline' else name, carrier)
    tp.mg.background = mg.background
    return tp.run(('hop', fc, amp*mg.A0, index, q.ITER, SEEDS, q.CODE, q.LABEL))


def summary(names):
    from gravcomm.packet_confidence import packet_ber_interval
    out = {}
    for n in names:
        rows = [json.loads(p.read_text()) for p in (RES/n).glob('[0-9]*.json')]
        e = [r['turbo_errors'][-1] for r in rows]
        lo, hi = packet_ber_interval(e, 224, 0.05) if rows else (None, None)
        out[n] = dict(packets=len(rows), errors=int(sum(e)), ber=sum(e)/(224*max(len(rows), 1)), interval=[lo, hi],
                      failed_packets=sum(x > 0 for x in e), rejected=sum(not r['turbo_accepted'][-1] for r in rows),
                      wrong_accepted=sum(r['turbo_accepted'][-1] and r['turbo_errors'][-1] > 0 for r in rows),
                      pruned=sum(r['pruned_states'] for r in rows), missing=sum(r['missing'] for r in rows))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='name=scenario:carrier:amplitude_factor')
    p.add_argument('--packets', type=int, default=3000)
    p.add_argument('--workers', type=int, default=18)
    p.add_argument('--summary', action='store_true')
    a = p.parse_args()
    names = [c.split('=')[0] for c in a.case]
    if not a.summary:
        jobs = []
        for c in a.case:
            name, spec = c.split('=')
            scen, fc, amp = spec.split(':')
            folder = RES/name; folder.mkdir(parents=True, exist_ok=True)
            have = {int(x.stem) for x in folder.glob('[0-9]*.json')}
            jobs += [((scen, float(fc), float(amp), i), folder) for i in range(START, START + a.packets) if i not in have]
        print(f'{len(jobs)} packets to run', flush=True)
        t0 = time.time(); done = 0
        with ProcessPoolExecutor(a.workers) as pool:
            futs = {pool.submit(run, j): f for j, f in jobs}
            for fut in as_completed(futs):
                res = fut.result(); folder = futs[fut]
                tmp = folder/f"{res['index']}.json.tmp"
                tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")
                done += 1
                if done % 128 == 0:
                    print(f'  {done}/{len(jobs)} ({time.time() - t0:.0f}s)', flush=True)
    s = summary(names)
    (RES/'validation_hop_summary.json').write_text(json.dumps(s, indent=1) + '\n')
    print(json.dumps(s, indent=1))


if __name__ == '__main__':
    main()
