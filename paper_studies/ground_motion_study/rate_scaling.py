"""Information-rate scaling of the mitigation scenarios for the joint (hop-limited) waveform.

The matched-filter (equal whitened SNR) scaling of analysis.py was found optimistic for the joint
receiver by packet validation (three-stage isolation at the predicted amplitudes: FER 2e-3 / 7e-3
instead of ~5e-4). Here the amplitude factor of a scenario is the ratio of information-rate thresholds,
thr(scenario, carrier, instrument factor) / thr(baseline, reference carrier), computed with the
auxiliary-channel bound for the hop-limited waveform (M=7, h=3/14, beta=1, T=224/190 s, hop limit 3,
maxentropic input) at 256/172 bit/symbol, 15000 symbols, predictor order 18 (as hop_threshold.py).
The same seed is used for every threshold, so ratios are less noisy than the thresholds.
Usage: python rate_scaling.py [--cases scen:fc:k,...] [--workers 8]  -> results/rate_scaling.json (merged)
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path
from scipy.optimize import brentq
import mitigation as mg
import cpm_rate as cr

RES = Path(__file__).resolve().parent/'results'
OUT = RES/'rate_scaling.json'
BITS = 256/172


def threshold(args):
    scen, fc, k = args
    t0 = time.time()
    d = cr.Design(M=7, h=Fraction(3, 14), beta=1.0, T=224/190, max_hop=3)
    bg = mg.background(scen, fc)
    psd = lambda f: 2*mg.total_input_psd(bg, f, k)/mg.A0**2
    g = lambda a: cr.info_rate(d, psd, a, n_sym=15000, seed=1, p_max=18) - BITS
    lo, hi = 0.05, 1.6
    while g(hi) < 0:                                         # noisy scenarios: widen the bracket
        lo, hi = hi, 2*hi
    while g(lo) > 0:                                         # quiet scenarios: lower the bracket
        lo, hi = lo/2, lo
    a = brentq(g, lo, hi, xtol=1e-5, rtol=2e-3)
    return dict(scenario=scen, carrier_hz=fc, instrument=k, threshold_A0=a, seconds=time.time() - t0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cases', default='baseline:24:1,baseline:18:1,three_stage:24:1,three_stage:18:1')
    p.add_argument('--workers', type=int, default=8)
    a = p.parse_args()
    cases = []
    for c in a.cases.split(','):
        s, f, k = c.split(':'); cases.append((s, float(f), float(k)))
    old = json.loads(OUT.read_text()) if OUT.exists() else {}
    cases = [c for c in cases if f'{c[0]}|{c[1]:g}|{c[2]:g}' not in old]
    with ProcessPoolExecutor(a.workers) as pool:
        futs = {pool.submit(threshold, c): c for c in cases}
        for fut in as_completed(futs):
            try:
                r = fut.result()
            except Exception as e:                           # keep going; report the failed case
                print('FAILED', futs[fut], repr(e), flush=True)
                continue
            old[f"{r['scenario']}|{r['carrier_hz']:g}|{r['instrument']:g}"] = r
            print(r, flush=True)
            OUT.write_text(json.dumps(old, indent=1) + '\n')


if __name__ == '__main__':
    main()
