"""A Brack-type bending-beam receiver at this paper's link geometry (Sec. 5.3).

Thermal floor of the beam of Brack et al. (Commun. Phys. 6, 270, 2023) in modal acceleration, compared with the
projected floor of the Carter-type receiver, and the fraction of the sized rotors' point acceleration that the
beam's first bending mode receives when the beam lies broadside to the rotor (axis perpendicular to the
rotor-receiver line, sensitive direction along that line). A free-free bending mode does not respond to uniform
acceleration, so the beam couples only to the variation of the field along its length.
Uses the finite-volume shaped rotor (mitigation.ROTOR_NPZ) scaled to each link design.
"""
import json, math
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.finite_rotor import FiniteRotor
from brack_comparison import mode_shape, gauss, BETA_L

HERE = Path(__file__).resolve().parent
KB, T = 1.380649e-23, 300.0
SCALES = {0.5: 0.9010437566657981, 1.0: 1.5761135310793883, 2.0: 2.854447700523453}


def thermal(mstar, f0, Q):
    return math.sqrt(4*KB*T*2*math.pi*f0/(mstar*Q))


def coupling(rotor, d, length, n=48, nphi=32):
    """Mode-weighted to point second-harmonic ratio (n beam Gauss points; checked against n=64, nphi=64)."""
    phi = 2*np.pi*np.arange(nphi)/nphi
    xb, w = gauss(n, 0.0, length)
    import brack_comparison as bc
    U = 0.5*(np.cosh(BETA_L/length*xb) + np.cos(BETA_L/length*xb) -
             (math.cosh(BETA_L) - math.cos(BETA_L))/(math.sinh(BETA_L) - math.sin(BETA_L)) *
             (np.sinh(BETA_L/length*xb) + np.sin(BETA_L/length*xb)))
    acc = np.zeros(nphi)
    for xi, wi, ui in zip(xb, w, U):
        yb = xi - length/2
        rho, alpha = math.hypot(d, yb), math.atan2(yb, d)
        ar = rotor.acceleration_component(rho, phi - alpha, 'x')
        at = rotor.acceleration_component(rho, phi - alpha, 'y')
        acc += wi*ui*(ar*math.cos(alpha) - at*math.sin(alpha))
    acc /= w@U**2
    point = rotor.acceleration_component(d, phi, 'x')
    c = lambda a: abs(2*np.mean(a*np.exp(-2j*phi)))
    return c(acc)/c(point)


def main():
    base = FiniteRotor.from_npz(mg.ROTOR_NPZ)
    bc = json.loads((HERE/'results'/'brack_comparison.json').read_text())
    mstar = bc['modal_mass_kg']
    out = {'T_K': T, 'brack_modal_mass_kg': mstar,
           'brack_thermal_nGal_rtHz': thermal(mstar, 42.6, 45000)/1e-11,
           'carter_thermal_nGal_rtHz': math.sqrt(8*math.pi*KB*T*50.3/(0.0031*6.37e5))/1e-11,
           'brack_linewidth_mHz': 42.6/45000*1e3, 'brack_ringdown_s': 45000/(math.pi*42.6),
           'coupling': {}, 'effective_floor_nGal_rtHz': {}}
    for name, (a, f0, Q) in {'brack': (out['brack_thermal_nGal_rtHz'], 42.6, 45000),
                             'carter': (out['carter_thermal_nGal_rtHz'], 50.3, 6.37e5)}.items():
        for df in (0.1, 0.5):
            out[f'{name}_thermal_displacement_fm_rtHz_at_{df}Hz'] = a*1e-11/((2*math.pi*f0)**2*2*df/f0)*1e15
    for d, s in SCALES.items():
        rotor = FiniteRotor(base.positions_m*s, base.masses_kg*s**3, base.bounding_radius_m*s)
        for L in (1.0, 0.5, 0.2):
            eta = coupling(rotor, d, L)
            out['coupling'][f'{d:g}m|L{L:g}'] = eta
            if L == 1.0 and d == 1.0:
                out['coupling_check_1m_L1_n64_nphi64'] = coupling(rotor, d, L, 64, 64)
            if L == 1.0:
                out['effective_floor_nGal_rtHz'][f'{d:g}m'] = out['brack_thermal_nGal_rtHz']/eta
            print(d, L, round(eta, 4), flush=True)
    (HERE/'results'/'brack_receiver.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
