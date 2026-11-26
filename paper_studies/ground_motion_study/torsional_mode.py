"""Does the 24.17 Hz drivetrain torsional mode of the 2 m rotor amplify disturbances near the 18 Hz carrier?
(manuscript App. A.8; final revision of PHYSSCR-157981)

Model: rotor inertia I on the combined carrier and stub-shaft torsional stiffness k, drive end held
(the screen of Appendix A, without motor or coupling inertia), f_t = 24.17 Hz.
1. Steady excitation (torque ripple, windage, drive jitter) across the carrier band 18 +/- 0.55 Hz:
   dynamic amplification 1/sqrt((1-r^2)^2 + (2 zeta r)^2).
2. Commanded modulation: single-mode response to quintic hop sequences (beta = 1, T_s = 224/190 s),
   peak relative to the quasistatic value, as in ground_motion_study/modal_sequences.py.
3. Rotor-angle jitter per newton-metre of steady torque at 18 Hz, and the resulting signal phase (2x).
Run from paper_studies/ground_motion_study/. Writes results/torsional_mode.json
"""
import json
from pathlib import Path
import numpy as np
from scipy.signal import cont2discrete, dlsim

HERE = Path(__file__).resolve().parent
I_Z = 59.893            # kg m^2, 2 m rotor (Sec. 4.2)
F_T = 24.17             # Hz, torsional screen (Appendix A)
FC, HALF_SPAN = 18.0, 0.545   # carrier and half the tone span (3 spacings of 0.1818 Hz)
ZETAS = (0.005, 0.02, 0.05)
T = 224/190
FS = 4000.0


def amplification(f, zeta):
    r = f/F_T
    return 1/np.sqrt((1 - r**2)**2 + (2*zeta*r)**2)


def pd_profile():
    tau = np.arange(0, 1, 1/(FS*T))
    return 30*tau**2*(1 - tau)**2


def modal_peak(hops, fn, zeta):
    pd = pd_profile(); a = np.concatenate([h*pd for h in hops]).astype(float)
    wn = 2*np.pi*fn
    A = np.array([[0, 1], [-wn**2, -2*zeta*wn]]); B = np.array([[0], [1]]); C = np.array([[1, 0]]); D = np.array([[0]])
    Ad, Bd, Cd, Dd, _ = cont2discrete((A, B, C, D), 1/FS, method='zoh')
    _, y, _ = dlsim((Ad, Bd, Cd, Dd, 1/FS), a[:, None])
    qs = a/wn**2
    return float(np.max(np.abs(y[:, 0]))/np.max(np.abs(qs)) - 1)


def random_hops(n, rng):
    s, out = 3, []
    for _ in range(n):
        nxt = rng.integers(max(0, s - 3), min(6, s + 3) + 1); out.append(int(nxt - s)); s = nxt
    return out


def main():
    k = I_Z*(2*np.pi*F_T)**2
    band = np.linspace(FC - HALF_SPAN, FC + HALF_SPAN, 201)
    rng = np.random.default_rng(1)
    cases = {'isolated': [3] + [0]*5, 'alternating +-3': [3, -3]*60, 'random D<=3': random_hops(400, rng)}
    out = {'I_z': I_Z, 'f_t': F_T, 'k_Nm_per_rad': k, 'steady': {}, 'modulation_overshoot': {}}
    for z in ZETAS:
        amp = amplification(band, z)
        out['steady'][f'{z:g}'] = {'band_min': float(amp.min()), 'band_max': float(amp.max()),
                                   'at_carrier': float(amplification(FC, z))}
        out['modulation_overshoot'][f'{z:g}'] = {name: modal_peak(h, F_T, z) for name, h in cases.items()}
    # rotor-angle jitter per N m of steady torque at the carrier, and the signal phase it causes
    theta = amplification(FC, 0.02)/k
    out['angle_per_Nm_rad'] = float(theta)
    out['signal_phase_per_Nm_deg'] = float(np.degrees(2*theta))
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results'/'torsional_mode.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
