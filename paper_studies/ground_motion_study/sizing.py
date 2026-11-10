"""Rotor sizes for a given qualified amplitude per carrier (Table 1 analogue).

For each released link distance the smallest geometric scale of the released rotor that produces the
required second-harmonic amplitude is found (mitigation.ScaledRotor.invert, finite-volume field), and
mass, inertia, peak inertial torque and peak power follow from geometric scaling of the released
designs (mitigation.sized). The released torques are those of CF-FSK; for another waveform they are
multiplied by its peak-slew ratio to CF-FSK at equal inertia (torque_ratio).
Usage: python sizing.py <amp24/A0> <amp18/A0> [torque_ratio] [excursion_hz]
"""
import json, sys
import mitigation as mg

ROTOR = None


def table(amp24, amp18, torque_ratio=1.0, excursion=mg.TONE_EDGE_HZ):
    global ROTOR
    ROTOR = ROTOR or mg.ScaledRotor()
    rows = {}
    for d, (fc, a_old, s_old, m_old, i_old, t_old, p_old) in mg.DESIGNS.items():
        a = (amp24 if fc == 24. else amp18)*mg.A0
        s = ROTOR.invert(a, d)
        z = mg.sized(d, s, fc)
        rows[f'{d:g} m'] = dict(carrier_hz=fc, amplitude_A0=a/mg.A0, scale=s, diameter_m=z['diameter_m'],
                                mass_kg=z['mass_kg'], inertia_kg_m2=z['inertia_kg_m2'],
                                peak_torque_nm=z['peak_torque_nm']*torque_ratio,
                                peak_power_w=z['peak_power_w']*torque_ratio,
                                rim_speed_ok=bool(s*(fc + excursion) <= 60),
                                vs_released=dict(mass=z['mass_kg']/m_old, torque=z['peak_torque_nm']*torque_ratio/t_old))
    return rows


if __name__ == '__main__':
    a24, a18 = float(sys.argv[1]), float(sys.argv[2])
    tr = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    exc = float(sys.argv[4]) if len(sys.argv) > 4 else mg.TONE_EDGE_HZ
    print(json.dumps(table(a24, a18, tr, exc), indent=1))
