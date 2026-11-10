"""Water-filled CF-FSK: tone placement and tone probabilities from the colored-noise water-filling.

Upper bound. For complex-baseband input power A^2 in Gaussian noise of PSD N(f), no signal can
exceed the water-filling capacity C(A) = int log2(1 + [nu - N]^+/N) df, int [nu - N]^+ df = A^2.
A rotor's constant-envelope signal is one such input, so C(A) bounds every rotor modulation,
whatever its torque.

Design rule. For slowly varying frequency, the spectrum of a constant-envelope signal follows the
distribution of its instantaneous frequency (Woodward's theorem). Matching that distribution to the
water-filling spectrum S(f) = [nu - N(f)]^+ prescribes how long the rotor should spend at each
frequency. With a finite alphabet this is done either by
  (a) equally spaced tones used with probabilities p_k = int_cell S / A^2 ("probability shaping"), or
  (b) equally likely tones at the midpoint quantiles of S ("non-uniform spacing"),
and both are evaluated exactly with cpm_rate.info_rate.
The result is the threshold amplitude at which each waveform reaches the required information rate,
compared with the water-filling bound and with the released CF-FSK waveform.
"""
import json
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import cpm_rate as cr
import cpm_optimize as co

OUT = Path(__file__).resolve().parent/'results'
FGRID = np.linspace(-2, 2, 40001); DFG = FGRID[1] - FGRID[0]
T = 224/162
GRID_H = Fraction(1, 40)          # base spacing df0 = GRID_H/T = 0.018 Hz for non-uniform tones


def waterfill(fc, amp):
    """Water-filling spectrum S(f) on FGRID and its rate (bit/s) for input power amp^2."""
    n = co.noise_psd(fc)(FGRID)
    nu = brentq(lambda v: np.sum(np.maximum(v - n, 0))*DFG - amp**2, n.min(), 50*n.max())
    s = np.maximum(nu - n, 0)
    return s, float(np.sum(np.log2(1 + s/n))*DFG)


def bound_threshold(fc, rate=co.R_REQ):
    return brentq(lambda a: waterfill(fc, a)[1] - rate, 1e-3, 3)


def symmetric_density(s):
    return 0.5*(s + s[::-1])


def quantile_positions(s, M, h=GRID_H, T=T):
    """Equally likely tones at the midpoint quantiles of the (symmetrized) water-filling density."""
    d = symmetric_density(s); c = np.cumsum(d); c /= c[-1]
    q = (np.arange(M) + 0.5)/M
    f = np.interp(q, c, FGRID)
    df0 = float(h)/T
    pos = np.round(f/df0).astype(int)
    pos = (pos - pos[::-1])//2                         # enforce exact symmetry
    return tuple(int(x) for x in pos)


def cell_weights(s, tones_hz):
    """Probabilities of equally spaced tones from the water-filling mass in each tone's cell."""
    d = symmetric_density(s)
    edges = np.r_[-np.inf, 0.5*(tones_hz[1:] + tones_hz[:-1]), np.inf]
    w = np.array([d[(FGRID >= lo) & (FGRID < hi)].sum() for lo, hi in zip(edges[:-1], edges[1:])])
    return tuple(w/w.sum())


def threshold(design, fc, n_sym=6000, p_max=12, lo=0.15, hi=1.2):
    npsd = co.noise_psd(fc)
    return brentq(lambda a: cr.info_rate(design, npsd, a, n_sym=n_sym, p_max=p_max) - co.R_REQ*design.T,
                  lo, hi, xtol=2e-4)


def designs_for(fc, amp, M=7, beta=0.6):
    """CF-FSK and its water-filled variants, with the water-filling taken at amplitude amp."""
    s, _ = waterfill(fc, amp)
    cf = cr.Design(M=M, h=Fraction(1, 5), beta=beta, T=T)
    prob = cr.Design(M=M, h=Fraction(1, 5), beta=beta, T=T, weights=cell_weights(s, cf.tones()*cf.spacing))
    spaced = cr.Design(M=M, h=GRID_H, beta=beta, T=T, positions=quantile_positions(s, M))
    return dict(cffsk=cf, probability_shaped=prob, nonuniform_spacing=spaced)


def main():
    out = {}
    for fc, a_op in co.CARRIERS.items():
        a_bound = bound_threshold(fc)
        a_cf = threshold(cr.Design(**co.CF), fc)
        res = dict(amp_op=a_op, bound_threshold=a_bound, designs={})
        for name, d in designs_for(fc, a_cf).items():
            a = a_cf if name == 'cffsk' else threshold(d, fc)
            a_fine = None
            res['designs'][name] = dict(
                tones_hz=list(np.round(d.tones()*d.spacing, 4)), weights=list(np.round(d.transition()[0], 4)),
                fdot=d.fdot_max, entropy=d.entropy_rate(), threshold=a, threshold_over_op=a/a_op,
                db_above_bound=20*np.log10(a/a_bound), db_below_cffsk=20*np.log10(a_cf/a))
            print(fc, name, json.dumps(res['designs'][name]), flush=True)
        out[f'{fc:g}'] = res
        print(fc, 'bound threshold', a_bound, '= %.3f A_op' % (a_bound/a_op), flush=True)
    (OUT/'waterfill_tones.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
