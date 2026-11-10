"""Spatial confinement and distance from capacity for the qualified released links.

Confinement: distance along the link axis at which the finite-volume A2 of each
qualified rotor falls to a given fraction of its value at the design distance.
Released-isolation packet screens at 0.80x the qualified amplitude
(baseline_f24_a0.6500, baseline_f18_a0.8000) give post-decoding BERs of about 8-10%.

Capacity: received power P = A2^2/2 against the one-sided total input PSD of the
released channel, harmonically weighted over the occupied band (same weight as the
equal-SNR sizing). Reports Eb/N0 at 1.00 bit/s, the gap to the Shannon limit
ln 2 (-1.59 dB), the wideband capacity P/(N0 ln 2) and the capacity restricted to
the occupied band.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import mitigation as mg

OUT = Path(__file__).resolve().parent/'results'/('link_limits_hop.json' if __import__('os').environ.get('GC_DESIGN') == 'hop' else 'link_limits.json')
RATE = 1.0          # nominal payload rate (bit/s)


def main():
    rotor = mg.ScaledRotor()
    out = {}
    for d, (fc, amp, s, *_rest) in mg.DESIGNS.items():
        a_design = abs(rotor.coefficient(s, (d, 0, 0)))
        amp_at = lambda x: abs(rotor.coefficient(s, (x, 0, 0)))
        conf = {}
        for frac in (0.8, 0.5, 0.25):
            de = brentq(lambda x: amp_at(x) - frac*a_design, d, 10*d, xtol=1e-5)
            conf[f'{frac:g}'] = dict(distance_m=de, ratio=de/d)
        twice = amp_at(2*d)/a_design
        # Noise: released two-stage channel at this carrier.
        bg = mg.background('baseline', fc)
        nu = np.linspace(-1, 1, 4001); w = mg.signal_weight(nu)
        S = mg.total_input_psd(bg, nu)
        n0_eff = np.trapezoid(w, nu)/np.trapezoid(w/S, nu)          # harmonic weighted mean
        n0_carrier = float(mg.total_input_psd(bg, np.array([0.]))[0])
        P = amp**2/2
        ebn0 = P/(n0_eff*RATE)
        band = 2*mg.TONE_EDGE_HZ
        out[f'{d:g}m'] = dict(
            carrier_hz=fc, A2_nGal=amp/1e-11, A2_check_nGal=a_design/1e-11,
            noise_asd_eff_nGal=np.sqrt(n0_eff)/1e-11, noise_asd_carrier_nGal=np.sqrt(n0_carrier)/1e-11,
            P_over_N0_Hz=P/n0_eff, EbN0_dB=10*np.log10(ebn0), gap_to_shannon_dB=10*np.log10(ebn0/np.log(2)),
            capacity_wideband_bit_s=P/(n0_eff*np.log(2)),
            occupied_span_Hz=band, capacity_in_span_bit_s=band*np.log2(1 + P/(n0_eff*band)),
            confinement=conf, amplitude_ratio_at_twice_distance=twice,
            integration_penalty_at_twice_distance=1/twice**2)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
