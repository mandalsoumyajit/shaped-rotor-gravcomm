"""Torque-limited rate of rigid-rotor (angle) modulation: small-noise phase-channel estimate.

A rigid rotor's Newtonian field has fixed amplitude, so the received signal is
s(t) = A cos(2 pi f_c t + phi(t)), and only the phase trajectory can carry information.
The drive sets phi'' = u = 2 pi fdot(t), with |fdot| <= fdot_max = tau_max/(pi I_z), and the
receiver resources bound the frequency excursion, |phi'| <= 2 pi Delta.

At small noise the receiver observes phi(t) + w(t), where w is phase noise with one-sided
PSD N_w(f) = [N(f_c+f) + N(f_c-f)] / (2P), P = A^2/2 and N the total input PSD of the channel.
Since phi is the double integral of u, this is a linear Gaussian channel from u to the
observation with noise PSD (2 pi f)^4 N_w(f) referred to u. Under the mean-square limits
implied by the peak limits,
   int S_u df <= (2 pi fdot_max)^2,   int S_u/(2 pi f)^2 df <= (2 pi Delta)^2,
its capacity is a water-filling with level 1/(ln2 (a + b/(2 pi f)^2)); the multipliers
follow from the KKT conditions by bisection (this replaced an L-BFGS dual
minimization, for robustness; both give the same values). Gaussian u maximizes the rate of this linear channel, so the result bounds every
torque-limited modulation within the small-noise approximation; it neglects phase wrapping
and cycle slips and is therefore an estimate, not a strict bound.

Note: a spectrum-moment bound on the transmitted signal was tried first and rejected: it
cannot distinguish slow from fast frequency wander, so it does not constrain torque.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import mitigation as mg

OUT = Path(__file__).resolve().parent/'results'/'torque_bound.json'
TS = 224/162                     # symbol duration (s)
BETA = 0.6                       # transition fraction of the symbol
SPACING = 0.14464                # tone spacing (Hz)
N_TONES = 7
DELTA = (N_TONES - 1)/2*SPACING  # peak excursion (Hz)
FDOT = 15/8*(N_TONES - 1)*SPACING/(BETA*TS)   # peak slew of the quintic largest hop (Hz/s)
F = np.geomspace(1e-5, 8, 6000); DF = np.gradient(F)
W = 2*np.pi*F


def rate(nw, fdot, delta=DELTA, full=False, flo=0.):
    """Water-filling rate (bit/s) for phase-noise PSD nw(F) under slew and excursion limits.

    flo > 0 excludes information below flo (a proxy for limited long-term phase stability).
    KKT solution S_u = [1/(ln2 (a + b/w^2)) - K]^+ : a is found by bisection so that the slew
    (mean-square command) limit is met or a = 0; b = 0 unless the excursion limit is violated,
    in which case b is bisected to meet it. The returned rate is the primal value of this
    feasible, KKT-satisfying spectrum.
    """
    keep = F >= flo; wk, dfk, k = W[keep], DF[keep], (W**4*nw)[keep]
    kref = np.median(k); kk = k/kref
    ucap = (2*np.pi*fdot)**2/kref; ecap = (2*np.pi*delta)**2/kref
    iw2 = 1/wk**2

    def spec(a, b):
        return np.maximum(1/(np.log(2)*(a + b*iw2)) - kk, 0)

    def inner(b):
        if b > 0 and np.sum(spec(0., b)*dfk) <= ucap:
            return 0.
        lo, hi = -90., 90.                               # log a
        for _ in range(200):
            mid = .5*(lo + hi)
            if np.sum(spec(np.exp(mid), b)*dfk) > ucap:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-10:
                break
        return np.exp(hi)

    a, b = inner(0.), 0.
    if np.sum(spec(a, 0.)*iw2*dfk) > ecap:
        lo, hi = -90., 90.                               # log b
        for _ in range(200):
            mid = .5*(lo + hi); am = inner(np.exp(mid))
            if np.sum(spec(am, np.exp(mid))*iw2*dfk) > ecap:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-9:
                break
        b = np.exp(hi); a = inner(b)
    sv = spec(a, b)
    r = float(np.sum(np.log2(1 + sv/kk)*dfk))
    if not full:
        return r
    return dict(rate=r, slew_used=float(np.sum(sv*dfk)/ucap), excursion_used=float(np.sum(sv*iw2*dfk)/ecap),
                band_Hz=float(F[keep][sv > 0].max()) if np.any(sv > 0) else 0.)


def main():
    out = dict(assumptions=dict(Ts=TS, beta=BETA, spacing_Hz=SPACING, tones=N_TONES,
                                excursion_Hz=DELTA, peak_slew_Hz_per_s=FDOT), links={})
    for d, (fc, amp, *_r) in mg.DESIGNS.items():
        if d == 1.0:
            continue                                    # same channel as 0.5 m
        bg = mg.background('baseline', fc); p = amp**2/2
        nw = (mg.total_input_psd(bg, F) + mg.total_input_psd(bg, -F))/(2*p)
        at = rate(nw, FDOT, full=True)
        fmin = brentq(lambda fd: rate(nw, fd) - 1.0, 1e-3, FDOT, xtol=1e-6)
        scan = []
        for dl in (0.2, 0.3, DELTA, 0.6, 0.8, 1.0, 1.5):
            try:
                fm = brentq(lambda fd: rate(nw, fd, dl) - 1.0, 1e-4, 50, xtol=1e-7)
                scan.append(dict(excursion_Hz=dl, min_slew=fm, torque_ratio=FDOT/fm))
            except ValueError:
                scan.append(dict(excursion_Hz=dl, min_slew=None, torque_ratio=None))
        curve = [dict(slew_fraction=x, rate=rate(nw, x*FDOT)) for x in (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)]
        out['links'][f'{fc:g}Hz ({"0.5 and 1 m" if fc == 24 else "2 m"})'] = dict(
            rate_at_design=at, min_slew_for_1bit=fmin, torque_ratio=FDOT/fmin,
            excursion_scan=scan, rate_vs_slew=curve)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
