"""Torque-amplitude Pareto front of CF-FSK waveforms (formal design search).

Each candidate waveform theta = (M, excursion, T, beta, hop limit D) is scored by two numbers:
  * threshold amplitude A*(theta): the smallest amplitude at which its information rate
    (cpm_rate.info_rate, ideal coding) reaches R_req = 256 bits per 144 symbols of 224/162 s,
    reported in dB above the water-filling capacity threshold (a strict bound for any signal);
  * peak inertial torque at fixed rotor inertia, c * max_jump * df / (beta T), relative to the
    released CF-FSK waveform (quintic transitions throughout; c scales all designs equally).
Tones are equally spaced. (The water-filling quantile layout is close to equally spaced with a
smaller span, see waterfill_tones.py, so the excursion grid covers it.) With a hop limit D < M-1
the input is the maxentropic Markov chain on the allowed hops.
The non-dominated set of (A*, torque) is the Pareto front.

Results are appended line by line to results/pareto_f<fc>.jsonl, so the run can resume.
Usage: python cpm_optimize.py [--workers 18] [--carrier 24] [--carrier 18]
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import cpm_rate as cr
import mitigation as mg

OUT = Path(__file__).resolve().parent/'results'
T_CF = 224/162
REQ_BITS = 256/144                 # information bits per data symbol of the released frame
R_REQ = REQ_BITS/T_CF              # 1.286 bit/s during data symbols
CF = dict(M=7, h=Fraction(1, 5), beta=0.6, T=T_CF)
CARRIERS = {24.: .8125, 18.: 1.0}  # operating amplitude / A0
N_SYM, P_MAX, MAX_DEN = 5000, 12, 24
A_MAX = 2.0                        # thresholds above this (A0 units) are recorded as out of range


def noise_psd(fc):
    bg = mg.background('baseline', fc)
    return lambda f: 2*mg.total_input_psd(bg, f)/mg.A0**2


def rate(design, fc, amp, n_sym=N_SYM, seed=1, p_max=P_MAX):
    return cr.info_rate(design, noise_psd(fc), amp, n_sym=n_sym, seed=seed, p_max=p_max)


def threshold(design, fc, n_sym=N_SYM, p_max=P_MAX, lo=0.1, hi=A_MAX):
    g = lambda a: rate(design, fc, a, n_sym, p_max=p_max) - R_REQ*design.T
    if g(hi) < 0:
        return None
    return brentq(g, lo, hi, xtol=1e-3, rtol=2e-3)


def rational_h(M, h):
    two = 1 if M % 2 else 2
    return Fraction(h/two).limit_denominator(MAX_DEN)*two


def make(M, exc, T, beta, D):
    h = rational_h(M, 2*exc/(M - 1)*T)
    return cr.Design(M=M, h=h, beta=beta, T=T, max_hop=None if D >= M - 1 else D)


def evaluate(job):
    fc, M, exc, T, beta, D = job
    t0 = time.perf_counter()
    d = make(M, exc, T, beta, D)
    a = threshold(d, fc)
    ns = max(2, int(round(4.0*T))); mem_r, mem_p = cr.memory(d, ns, P_MAX, d.lattice()[2])
    return dict(predictor_order=mem_p, history=mem_r, fc=fc, M=M, excursion_req=exc, excursion=d.excursion, T=T, beta=beta, D=d.D, h=str(d.h),
                bits_per_symbol=R_REQ*T, entropy=d.entropy_rate(), fdot=d.fdot_max, threshold=a,
                elapsed=time.perf_counter() - t0)


def jobs_for(fc):
    out = []
    for M in range(2, 8):
        hops = sorted({M - 1, max(1, round((M - 1)/2)), 1})
        for D in hops:
            H = cr.Design(M=M, h=Fraction(1, 5), beta=1., T=1., max_hop=None if D >= M - 1 else D).entropy_rate()
            for bits in (0.5, 0.75, 1.0, 1.25, 1.5, 1.78, 2.0, 2.25, 2.5):
                if bits > 0.92*H:
                    continue
                for exc in (0.15, 0.2, 0.25, 0.3, 0.35, 0.434, 0.55):
                    for beta in (0.6, 0.8, 1.0):
                        out.append((fc, M, exc, bits/R_REQ, beta, D))
    return out


def pareto(rows, bound):
    pts = sorted((20*np.log10(r['threshold']/bound), r['torque_ratio'], i) for i, r in enumerate(rows)
                 if r['threshold'] is not None)
    front, best = [], np.inf
    for x, y, i in pts:
        if y < best:
            front.append(i); best = y
    return front


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--workers', type=int, default=18)
    p.add_argument('--carrier', type=float, action='append')
    a = p.parse_args()
    OUT.mkdir(exist_ok=True)
    cr.info_rate(cr.Design(**CF), noise_psd(24.), 0.5, n_sym=50)          # compile numba cache once
    for fc in (a.carrier or CARRIERS):
        path = OUT/f'pareto_f{fc:g}.jsonl'
        done = set()
        if path.exists():
            for line in path.read_text().splitlines():
                r = json.loads(line); done.add((r['fc'], r['M'], r['excursion_req'], round(r['T'], 9), r['beta'], r['D']))
        jobs = [j for j in jobs_for(fc) if (j[0], j[1], j[2], round(j[3], 9), j[4], j[5]) not in done]
        print(f'{fc} Hz: {len(jobs)} designs to evaluate ({len(done)} done)', flush=True)
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=a.workers) as pool, path.open('a') as fh:
            futs = [pool.submit(evaluate, j) for j in jobs]
            for i, f in enumerate(as_completed(futs)):
                fh.write(json.dumps(f.result()) + '\n'); fh.flush()
                if (i + 1) % 100 == 0:
                    print(f'  {i + 1}/{len(jobs)} ({time.time() - t0:.0f}s)', flush=True)
        print(f'{fc} Hz finished in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
