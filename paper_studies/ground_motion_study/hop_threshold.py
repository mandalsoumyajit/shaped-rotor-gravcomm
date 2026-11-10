"""Information-rate threshold of the joint design's waveform and its place in the gap budget.

Waveform: M=7 tones, h=3/14, beta=1 (no dwell), T=224/190 s, hop limit 3, maxentropic Markov tone
input (cpm_rate.Design). The frame carries 256 information bits in 172 data symbols, so the threshold
is the amplitude at which the auxiliary-channel information rate equals 256/172 bit/symbol.
Refined settings as for the Pareto front (predictor order 18, 15000 symbols).
Writes results/hop_threshold.json.
"""
import json, sys
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from scipy.optimize import brentq
import cpm_rate as cr
import cpm_optimize as co
import waterfill_tones as wt

RES = Path(__file__).resolve().parent/'results'
BITS = 256/172
import os
if os.environ.get("HOP_RATE") == "reference":       # Pareto convention: R_REQ during data symbols
    BITS = co.R_REQ*224/190


def one(fc):
    d = cr.Design(M=7, h=Fraction(3, 14), beta=1.0, T=224/190, max_hop=3)
    g = lambda a: co.rate(d, fc, a, n_sym=15000, p_max=18) - BITS
    a = brentq(g, 0.2, 1.2, xtol=1e-3, rtol=2e-3)
    return dict(carrier_hz=fc, threshold_A0=a, bound_A0=wt.bound_threshold(fc), entropy_rate=d.entropy_rate(),
                fdot=d.fdot_max, excursion=d.excursion)


def main():
    with ProcessPoolExecutor(2) as pool:
        res = list(pool.map(one, (24., 18.)))
    (RES/('hop_threshold_ref.json' if os.environ.get('HOP_RATE') == 'reference' else 'hop_threshold.json')).write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
