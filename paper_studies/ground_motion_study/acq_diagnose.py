"""Diagnose acquisition losses: rerun packets with true/estimated gain and timing in all combinations (GPU)."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import json, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import acq_packets as ap

RES = Path(__file__).resolve().parent/'results'


def one(args):
    fc, amp, idx = args
    b = ap.build(fc, amp, idx, ap.SEEDS)
    g_hat, tau_hat = ap.acquire(b, fc, amp)
    out = dict(index=idx, pe=float(np.degrees(np.angle(g_hat/b['g']))), te=tau_hat - b['tau'], ae=abs(g_hat) - 1)
    for name, g, tau in (('known', b['g'], b['tau']), ('acq', g_hat, tau_hat),
                         ('g_est_tau_true', g_hat, b['tau']), ('g_true_tau_est', b['g'], tau_hat)):
        r = ap.detect_decode(b, amp, g, tau, device='gpu')
        out[name] = (bool(r['accepted'][-1]), int(r['errors'][-1]), len(r['errors']))
    return out


if __name__ == '__main__':
    fc, amp, n = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    folder = RES/f'acq_f{fc:g}_a{amp:g}'
    rows = [json.loads(p.read_text()) for p in folder.glob('9*.json')]
    fails = sorted(r['index'] for r in rows if not r['accepted'][-1])[:n]
    with ProcessPoolExecutor(4) as pool:
        res = list(pool.map(one, [(fc, amp, i) for i in fails]))
    for r in res:
        print(r, flush=True)
    for k in ('known', 'acq', 'g_est_tau_true', 'g_true_tau_est'):
        print(k, 'failures', sum(not r[k][0] for r in res), 'of', len(res))
    (RES/f'acq_diagnose_f{fc:g}.json').write_text(json.dumps(res, indent=1))
