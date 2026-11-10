"""Short-block loss of the CF-FSK link: finite-blocklength normal approximation.

For a frame of n channel symbols carrying k information bits with frame error probability eps,
the normal approximation (Polyanskiy, Poor and Verdu, IEEE Trans. IT 56, 2307, 2010) gives
    k  ~  n I  -  sqrt(n V) Q^{-1}(eps)  +  (1/2) log2 n,
where I is the information rate (bit/symbol) and V the dispersion, the variance rate of the
information density. Both are estimated from one long trellis simulation (cpm_rate.info_rate):
I as the mean of the per-symbol information-density increments, V by batch means (batches much
longer than the correlation of the increments). Inputs are i.i.d. uniform tones, as in the
released waveform, so this is the finite-length limit for that waveform, not for the channel.

Output: amplitudes (in A0) at which the released frame (k = 256, n = 144 data symbols) is
possible with ideal coding at infinite length (threshold), at finite length (normal
approximation, eps = 1e-2 and 1e-3), and the operating amplitude of the qualified link.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm
import cpm_rate as cr
import cpm_optimize as co

OUT = Path(__file__).resolve().parent/'results'/'finite_length.json'
K_BITS, N_SYM_FRAME = 256, 144
import os
from fractions import Fraction
HOP = os.environ.get('FL_DESIGN') == 'hop'     # joint design: 256 bits in 172 hop-limited symbols
if HOP:
    N_SYM_FRAME = 172
    OUT = OUT.with_name('finite_length_hop.json')


def rate_and_dispersion(design, npsd, amp, n_sym=40000, batch=400, seed=7, p_max=18):
    _r, dens = cr.info_rate(design, npsd, amp, n_sym=n_sym, seed=seed, p_max=p_max, return_density=True)
    nb = len(dens)//batch
    sums = dens[:nb*batch].reshape(nb, batch).sum(axis=1)
    return float(dens.mean()), float(sums.var(ddof=1)/batch), float(sums.std(ddof=1)/batch/np.sqrt(nb))


def na_bits(I, V, n=N_SYM_FRAME, eps=1e-2):
    return n*I - np.sqrt(n*V)*norm.isf(eps) + 0.5*np.log2(n)


def grid_fit(design, npsd, amps, n_sym=60000, batch=500):
    """I(A) and V(A) on an amplitude grid; V smoothed by a quadratic fit in A (it varies slowly)."""
    rows = [rate_and_dispersion(design, npsd, a, n_sym=n_sym, batch=batch) for a in amps]
    I = np.array([r[0] for r in rows]); V = np.array([r[1] for r in rows]); se = np.array([r[2] for r in rows])
    vfit = np.polyfit(amps, V, 2)
    return I, V, se, vfit


def solve_curve(amps, I, vfit, eps=None):
    a = np.linspace(amps[0], amps[-1], 2001)
    Ii = np.interp(a, amps, I); Vi = np.polyval(vfit, a)
    bits = N_SYM_FRAME*Ii if eps is None else na_bits(Ii, Vi, eps=eps)
    k = np.argmax(bits >= K_BITS)
    if bits[k] < K_BITS:
        return None
    return float(np.interp(K_BITS, bits[k - 1:k + 1], a[k - 1:k + 1]))


def main():
    out = {}
    d = cr.Design(**co.CF)
    grids = {24.: np.linspace(0.30, 0.62, 9), 18.: np.linspace(0.36, 0.76, 9)}
    carriers = co.CARRIERS
    if HOP:
        d = cr.Design(M=7, h=Fraction(3, 14), beta=1.0, T=224/190, max_hop=3)
        carriers = {24.: 0.58, 18.: 0.697}          # NB-LDPC GF(512) design, qualified 2026-11-05 (K5: 0.70 / 0.84)
    for fc, a_op in carriers.items():
        npsd = co.noise_psd(fc)
        amps = grids[fc]
        I, V, se, vfit = grid_fit(d, npsd, amps)
        res = dict(amp_op=a_op, grid=dict(amp=list(amps), I=list(I), V=list(V), I_stderr=list(se), V_fit=list(vfit)))
        res['threshold'] = solve_curve(amps, I, vfit)
        for eps in (1e-2, 1e-3):
            res[f'normal_approx_eps{eps:g}'] = solve_curve(amps, I, vfit, eps)
        Iop, Vop, seop = rate_and_dispersion(d, npsd, a_op, n_sym=60000, batch=500)
        res['at_operating'] = {'I': Iop, 'V': Vop, 'I_stderr': seop, 'na_bits_eps0.01': na_bits(Iop, Vop)}
        thr = res['threshold']
        res['db_above_threshold'] = {k: 20*np.log10(res[k]/thr) for k in
                                     ('normal_approx_eps0.01', 'normal_approx_eps0.001', 'amp_op') if res[k]}
        out[f'{fc:g}'] = res
        print(fc, json.dumps({k: v for k, v in res.items() if k != 'grid'}, indent=1), flush=True)
        print('  I', np.round(I, 3), '\n  V', np.round(V, 3), flush=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
