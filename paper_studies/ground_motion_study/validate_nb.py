"""Packet checks of the frozen NB-LDPC GF(512) joint design under other isolation scenarios or amplitudes.

Uses nb_packets.run unchanged (its hash is frozen by qualify_nb.py) with the qualification settings
(qualify_nb.ITER, BP, CODE_SEED, DV, WORDS); the isolation scenario is switched by replacing
mitigation.background inside the worker, so 'baseline' maps to the chosen scenario. Fresh seed
namespaces (20371531/32, indices from 9100000). One scenario per invocation; CPU + GPU workers.
Usage: python validate_nb.py --scenario three_stage --case vnb3_f24_a0.4935=24:0.4935 --packets 3000
       python validate_nb.py --summary --case vnb3_f24_a0.4935=24:0.4935
       -> results/<name>/<index>.json (resumable); summary in results/validation_nb_summary.json
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, functools, json, time
from pathlib import Path
import mitigation as mg
import nb_packets as npk
import qualify_nb as q

RES = Path(__file__).resolve().parent/'results'
SEEDS = (20371531, 20371532)
START = 9100000
if os.environ.get('VALIDATE_SEEDS'):                  # fresh namespaces for a repeated check
    SEEDS = tuple(int(x) for x in os.environ['VALIDATE_SEEDS'].split(','))
    START = int(os.environ.get('VALIDATE_START', START))


def run_scenario(scenario, job):
    if scenario != 'baseline' and not getattr(mg, '_validate_patched', False):
        orig = mg.background
        mg.background = lambda name, carrier: orig(scenario if name == 'baseline' else name, carrier)
        mg._validate_patched = True
    return npk.run(job)


def summary(names):
    from gravcomm.packet_confidence import packet_ber_interval
    out = {}
    for n in names:
        rows = [json.loads(p.read_text()) for p in (RES/n).glob('[0-9]*.json')]
        e = [r['errors'][-1] for r in rows]
        lo, hi = packet_ber_interval(e, 224, 0.05) if rows else (None, None)
        out[n] = dict(packets=len(rows), errors=int(sum(e)), ber=sum(e)/(224*max(len(rows), 1)), interval=[lo, hi],
                      failed_packets=sum(x > 0 for x in e), rejected=sum(not r['accepted'][-1] for r in rows),
                      wrong_accepted=sum(r['accepted'][-1] and r['errors'][-1] > 0 for r in rows),
                      pruned=sum(r['pruned_states'] for r in rows))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='name=carrier:amplitude_factor_of_A0')
    p.add_argument('--scenario', default='baseline')
    p.add_argument('--packets', type=int, default=3000)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--gpu-workers', type=int, default=4)
    p.add_argument('--summary', action='store_true')
    a = p.parse_args()
    names = [c.split('=')[0] for c in a.case]
    if not a.summary:
        jobs = []
        for c in a.case:
            name, spec = c.split('=')
            fc, amp = (float(x) for x in spec.split(':'))
            folder = RES/name; folder.mkdir(parents=True, exist_ok=True)
            (folder/'config.json').write_text(json.dumps(dict(scenario=a.scenario, carrier_hz=fc, amplitude_A0=amp,
                                                              seeds=SEEDS, start=START, iterations=q.ITER, bp=q.BP,
                                                              code_seed=q.CODE_SEED, dv=q.DV, words=q.WORDS)) + '\n')
            have = {int(x.stem) for x in folder.glob('[0-9]*.json')}
            jobs += [((fc, amp*mg.A0, i, q.ITER, SEEDS, q.CODE_SEED, q.BP, q.DV, q.WORDS), folder)
                     for i in range(START, START + a.packets) if i not in have]
        print(f'{len(jobs)} packets to run ({a.scenario})', flush=True)
        t0 = time.time(); done = 0
        fn = functools.partial(run_scenario, a.scenario)
        for res, folder in npk.run_hybrid(jobs, a.workers, a.gpu_workers, fn=fn):
            tmp = folder/f"{res['index']}.json.tmp"
            tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")
            done += 1
            if done % 128 == 0:
                print(f'  {done}/{len(jobs)} ({time.time() - t0:.0f}s)', flush=True)
    s = summary(names)
    path = RES/'validation_nb_summary.json'
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(s)
    path.write_text(json.dumps(old, indent=1) + '\n')
    print(json.dumps(s, indent=1))


if __name__ == '__main__':
    main()
