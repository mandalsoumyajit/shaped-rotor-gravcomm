"""Compare the occupied-band weight model with the actual transmitted spectrum.

Equal-whitened-SNR scaling (mitigation.equivalent_amplitude_factor) weights 1/S(nu) by a model of
the occupied band: flat over the tone span with 0.15 Hz Gaussian shoulders. Here the weight is
replaced by the Welch PSD of the actual baseband waveform (random data) and the amplitude factors
of representative scenarios are compared. Writes results/weight_check.json.
"""
import json
from pathlib import Path
import numpy as np
from scipy.signal import welch
import mitigation as mg
from gravcomm.causal_receiver import sampled_cffsk
import hop_packets as hp
from hop_receiver import sampled_hop

RES = Path(__file__).resolve().parent/'results'
FS = 16.


def psd_hop(n_words=6000, seed=3):
    rng = np.random.default_rng(seed); tones, prev = [], hp.M//2
    for w in rng.integers(0, 2**hp.BITS, n_words):
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    x = sampled_hop(tones, 224/190, hp.H, float(hp.BETA), hp.M, FS)
    return welch(x, FS, nperseg=4096, return_onesided=False, detrend=False)


def psd_cf(n=24000, seed=3):
    rng = np.random.default_rng(seed)
    x = sampled_cffsk(rng.integers(-3, 4, n), 224/162, q=.08, transition_fraction=.6, input_rate_hz=FS)
    return welch(x, FS, nperseg=4096, return_onesided=False, detrend=False)


def factor(weight_fn, bg, bg_ref, k=1.):
    nu = np.linspace(-1, 1, 4001); w = weight_fn(nu)
    ref = np.trapezoid(w/mg.total_input_psd(bg_ref, nu), nu)
    new = np.trapezoid(w/mg.total_input_psd(bg, nu, k), nu)
    return float(np.sqrt(ref/new))


def main():
    out = {}
    for name, (f, p), edge in (('cf', psd_cf(), mg.RELEASED_TONE_EDGE_HZ), ('hop', psd_hop(), mg.HOP_TONE_EDGE_HZ)):
        order = np.argsort(f); f, p = f[order], p[order]
        actual = lambda nu: np.interp(nu, f, p)
        model = lambda nu: np.exp(-0.5*(np.maximum(np.abs(nu) - edge, 0)/0.15)**2)
        nu = np.linspace(-1, 1, 4001)
        rows = {'in_band_fraction_pm1Hz': float(np.trapezoid(actual(nu), nu)/np.trapezoid(p, f))}
        for fc_ref in (24., 18.):
            ref = mg.background('baseline', fc_ref)
            for scen, fc, k in (('three_stage', fc_ref, 1.), ('three_stage', 10., 1.), ('baseline', 10., 1.),
                                ('three_stage', fc_ref, .1), ('baseline', fc_ref, .1), ('gradiometer', fc_ref, 1.)):
                key = f'ref{fc_ref:g}|{scen}|{fc:g}Hz|inst{k:g}'
                rows[key] = dict(model=factor(model, mg.background(scen, fc), ref, k),
                                 actual=factor(actual, mg.background(scen, fc), ref, k))
        out[name] = rows
    (RES/'weight_check.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
