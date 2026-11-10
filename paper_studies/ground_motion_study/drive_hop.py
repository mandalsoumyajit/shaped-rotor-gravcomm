"""Drive torque, peak and mean absolute inertial power of the joint (hop-limited) design.

Commanded gravitational frequency f_g(t): quintic transitions between tone frequencies, full-symbol
transitions (beta = 1) for the hop design, beta = 0.6 with dwell for released CF-FSK (check).
Rigid-rotor inertial torque I*pi*fdot, mechanical power I*pi^2*f*fdot (omega = pi f_g).
Tone sequences come from random payload words through the enumerative hop code (hop_packets.CODE).
Writes results/drive_hop.json.
"""
import json
from pathlib import Path
import numpy as np
import mitigation as mg
import hop_packets as hp

RES = Path(__file__).resolve().parent/'results'


def trajectory(offs, T, beta, df, fc, n=200):
    """Sampled f_g and its derivative for tone offsets offs (index units), start at centre tone."""
    u = (np.arange(n) + .5)/n
    s, ds = 10*u**3 - 15*u**4 + 6*u**5, 30*u**2 - 60*u**3 + 30*u**4
    prev = np.r_[0, offs[:-1]]
    f, fd, w = [], [], []
    for a, b in zip(prev, offs):
        tr = beta*T
        f.append(fc + df*(a + (b - a)*s)); fd.append(df*(b - a)*ds/tr); w.append(np.full(n, tr/n))
        if beta < 1:
            f.append(np.array([fc + df*b])); fd.append(np.zeros(1)); w.append(np.array([(1 - beta)*T]))
    return np.concatenate(f), np.concatenate(fd), np.concatenate(w)


def stats(offs, T, beta, df, fc, inertia):
    f, fd, w = trajectory(offs, T, beta, df, fc)
    P = np.pi**2*inertia*f*fd
    hops = np.abs(np.diff(np.r_[0, offs]))
    return dict(peak_torque_nm=float(np.pi*inertia*np.max(np.abs(fd))), peak_power_w=float(np.max(np.abs(P))),
                mean_abs_power_w=float(np.sum(np.abs(P)*w)/np.sum(w)), mean_hop_spacings=float(hops.mean()),
                hop_distribution={int(k): float(np.mean(hops == k)) for k in range(int(hops.max()) + 1)})


def main(n_words=20000, seed=5):
    rng = np.random.default_rng(seed)
    tones, prev = [], hp.M//2
    for w in rng.integers(0, 2**hp.BITS, n_words):
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    hop_offs = np.array(tones) - hp.M//2
    T = 224/190; df = float(hp.H)/T
    cf_offs = rng.integers(-3, 4, 4*n_words)
    out = dict(hop={}, cf={})
    for d, (fc, a, s, m, I, t0, p0) in mg.HOP_DESIGNS.items():
        r = stats(hop_offs, T, 1.0, df, fc, I)
        r['ratio_mean_abs_to_peak'] = r['mean_abs_power_w']/r['peak_power_w']
        r['net_inertial_w_eta0.9'] = (1 - .9)*r['mean_abs_power_w']/2
        out['hop'][f'{d:g}m'] = r
    for d, (fc, a, s, m, I, t0, p0) in mg.RELEASED_DESIGNS.items():
        r = stats(cf_offs, 224/162, .6, .14464, fc, I)
        r['ratio_mean_abs_to_peak'] = r['mean_abs_power_w']/r['peak_power_w']
        out['cf'][f'{d:g}m'] = r
    (RES/'drive_hop.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
