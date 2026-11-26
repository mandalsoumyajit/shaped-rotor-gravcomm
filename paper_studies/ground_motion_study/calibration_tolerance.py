"""Calibration tolerance of the receiver: how accurately must f0 and Q be tracked?
(manuscript Sec. 2.4 and App. B.1; final revision of PHYSSCR-157981)

The receiver converts measured displacement to acceleration by dividing by the assumed response H'(f).
With true response H, every input-referred term (signal, thermal, support and Newtonian noise) and the
readout term s_x/H' are multiplied by the same ratio R(f) = H(f)/H'(f). Information is unchanged, but the
receiver's templates and noise model no longer match. For a decision statistic matched to a template
difference d under the assumed noise PSD S, the effective SNR relative to the ideal is

    eta = (Re <R>_w)^2 / <|R|^2>_w,     <g>_w = int g w dnu,   w = |d|^2/S / int |d|^2/S,

and the loss is -10 log10 eta. |d|^2 is taken as the PSD of the transmitted waveform (random coded data).
The resonance is tuned to the carrier (f0 = fc) with the reference Q, as in the noise model.

Run from paper_studies/ground_motion_study/. Writes results/calibration_tolerance.json
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.signal import welch
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mitigation as mg                      # noqa: E402  (puts the release package on the path)
import hop_packets as hp                     # noqa: E402
from hop_receiver import sampled_hop         # noqa: E402

Q0 = 6.37e5
FS = 16.
CARRIERS = (24., 18.)


def waveform_psd(n_words=20000, seed=11):
    rng = np.random.default_rng(seed); tones, prev = [], hp.M//2
    for w in rng.integers(0, 2**hp.BITS, n_words):
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    x = sampled_hop(tones, 224/190, hp.H, float(hp.BETA), hp.M, FS)
    f, p = welch(x, FS, nperseg=8192, return_onesided=False, detrend=False)
    o = np.argsort(f)
    return f[o], p[o]


def H(f, f0, Q):
    return 1/((2*np.pi)**2*(f0**2 - f**2 + 1j*f0**2/Q))


def grid(gamma):
    """Offsets within the 2 Hz band, dense within a few hundred linewidths of the resonance."""
    core = gamma*np.sinh(np.linspace(-np.arcsinh(0.05/gamma), np.arcsinh(0.05/gamma), 400001))
    outer = np.linspace(-1, 1, 20001)
    return np.unique(np.r_[core, outer[np.abs(outer) > 0.05]])


def loss_db(nu, w, fc, df0=0., q_ratio=1.):
    """True resonance fc + df0 with Q0*q_ratio; the receiver assumes fc and Q0."""
    f = fc + nu
    R = H(f, fc + df0, Q0*q_ratio)/H(f, fc, Q0)
    eta = np.trapezoid(R.real*w, nu)**2/np.trapezoid(np.abs(R)**2*w, nu)
    return float(-10*np.log10(eta))


def main():
    f, p = waveform_psd()
    i0 = np.argmin(np.abs(f)); line = float(p[i0]/np.median(p[i0-5:i0+6]))   # spectral line at the carrier?
    out = {'Q': Q0, 'line_ratio_at_carrier': line, 'carriers': {}}
    for fc in CARRIERS:
        gamma = fc/(2*Q0)                                   # half-width at half-maximum of |H|^2
        nu = grid(gamma)
        S = mg.total_input_psd(mg.background('baseline', fc), nu)
        wd = np.interp(nu, f, p)/S
        w = wd/np.trapezoid(wd, nu)
        res = {'hwhm_Hz': gamma, 'w0_per_Hz': float(np.interp(0, nu, w))}
        for target in (0.05, 0.1):
            df = brentq(lambda d: loss_db(nu, w, fc, df0=d) - target, 1e-9, 0.1)
            qhi = brentq(lambda r: loss_db(nu, w, fc, q_ratio=r) - target, 1.0001, 1e4)
            qlo = brentq(lambda r: loss_db(nu, w, fc, q_ratio=r) - target, 1e-4, 0.9999)
            res[f'{target}dB'] = {'df0_Hz': df, 'df0_over_hwhm': df/gamma, 'df0_over_fc': df/fc,
                                  'Q_ratio_range': [qlo, qhi]}
        # analytic small-error approximation: loss ~ 4.34 * pi * w(0) * delta^2 / gamma
        res['approx_df0_0.1dB'] = float(np.sqrt(0.1/4.343*gamma/(np.pi*res['w0_per_Hz'])))
        res['table'] = {f'{d:g}': loss_db(nu, w, fc, df0=d) for d in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3)}
        out['carriers'][f'{fc:g}'] = res
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results'/'calibration_tolerance.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
