"""Reproduce the two-rotating-bar experiment of Brack et al. (Commun. Phys. 6, 270, 2023; arXiv:2301.01644)
with this paper's finite-volume field code (gravcomm.FiniteRotor).

Geometry and parameters, all from the paper:
  transmitter: two tungsten bars 499.86 x 10.06 x 10.06 mm, measured masses 970.8 g and 971.8 g,
               rotating about z through their centres of mass; centres x0 = 800.0 mm apart, at y = d0.
  detector:    titanium beam 1000.00 x 16.97 x 8.49 mm along x, centred at the origin; first lateral
               (y) bending mode, f0 = 42.6 Hz, Q = 45000 (the value Fig. 4 is scaled to).
  The bars and the beam share the plane z = 0 (assumed; the paper does not state an offset).
Observable: complex end velocity of the first bending mode at resonance, normalised U(0)=1 at the beam end,
  v = Q * c_n / w0,  c_n = n-th rotation harmonic of  int a_y U dV / int U^2 dV   (their Eq. 2).
The beam density cancels from v. It enters only the power check, through rho = E h^2 (beta L)^4 / (12 L^4 w0^2)
from their fitted E = 107.4 GPa, which reproduces f0.

Targets: Table II power-law exponents (300-600 mm), |v| up to ~100 nm/s, gravitational power 3.75e-18 W at
d0 = 296 mm (w0/2, CW-CCW), ~8% lower amplitude for equal rotation directions, component angle ~23 deg.
Usage: python brack_comparison.py [--fine]
"""
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
import mitigation  # noqa: F401  (puts the release package on sys.path)
from gravcomm.constants import G
from gravcomm.finite_rotor import FiniteRotor

HERE = Path(__file__).resolve().parent
L_BAR, W_BAR = 0.49986, 0.01006
M_BARS = (0.9708, 0.9718)
X0 = 0.8000
L_B, B_Z, H_Y = 1.0, 0.01697, 0.00849          # beam length, z width, y thickness (bending direction)
F0, Q = 42.6, 45000.0
E_TI = 107.4e9
BETA_L = 4.730040744862704                      # first free-free root of cos(x)cosh(x) = 1
TABLE2 = {('n2', 'CCW-CW'): 4.35, ('n2', 'CW-CCW'): 4.35, ('n2', 'CCW-CCW'): 4.65, ('n2', 'CW-CW'): 4.65,
          ('n2', 'single'): 4.35, ('n4', 'CCW-CW'): 5.83, ('n4', 'CW-CCW'): 5.83, ('n4', 'CCW-CCW'): 5.78,
          ('n4', 'CW-CW'): 5.78, ('n4', 'single'): 5.83}


def mode_shape(x):
    """Free-free first bending mode on [0, L_B], normalised so U(0) = 1 (beam end)."""
    b = BETA_L/L_B*x
    sig = (math.cosh(BETA_L)-math.cos(BETA_L))/(math.sinh(BETA_L)-math.sin(BETA_L))
    return 0.5*(np.cosh(b)+np.cos(b)-sig*(np.sinh(b)+np.sin(b)))


def gauss(n, a, b):
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b-a)*t+0.5*(a+b), 0.5*(b-a)*w


def bar_rotor(mass, nl, nc, z_shift):
    """Tensor Gauss quadrature of a bar lying along body x, shifted by z_shift (receiver at z = -z_shift)."""
    xl, wl = gauss(nl, -L_BAR/2, L_BAR/2)
    xc, wc = gauss(nc, -W_BAR/2, W_BAR/2)
    X, Y, Z = np.meshgrid(xl, xc, xc, indexing='ij')
    W = np.einsum('i,j,k->ijk', wl, wc, wc)
    pos = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()+z_shift])
    m = mass*W.ravel()/W.sum()
    return FiniteRotor(pos, m, float(np.hypot(L_BAR/2, W_BAR/2))*(1+1e-9))


def modal_coefficients(d0, q):
    """Fourier coefficients a_{k,m} (m = -8..8) of the modal acceleration int a_y U dV / int U^2 dV
    produced by bar k as a function of its own angle phi (phi = 0: bar parallel to the beam)."""
    nphi = q['nphi']
    phi = 2*np.pi*np.arange(nphi)/nphi
    xb, wx = gauss(q['nx'], 0.0, L_B)                       # beam coordinate from the end
    yb, wy = gauss(q['ny'], -H_Y/2, H_Y/2)
    zb, wz = gauss(q['nz'], -B_Z/2, B_Z/2)
    U = mode_shape(xb)
    norm = (wx@U**2)*wy.sum()*wz.sum()                      # int U^2 dV
    coeffs = []
    for k, cx in enumerate((-X0/2, X0/2)):
        acc = np.zeros(nphi)
        for z, w_z in zip(zb, wz):
            rot = bar_rotor(M_BARS[k], q['nl'], q['nc'], -z)
            for xi, w_x, Ui in zip(xb, wx, U):
                for y, w_y in zip(yb, wy):
                    rx, ry = (xi-L_B/2)-cx, y-d0             # beam point relative to the bar centre
                    rho, alpha = math.hypot(rx, ry), math.atan2(ry, rx)
                    ar = rot.acceleration_component(rho, phi-alpha, 'x')   # along +r (away from bar)
                    at = rot.acceleration_component(rho, phi-alpha, 'y')   # along +alpha
                    ay = ar*math.sin(alpha)+at*math.cos(alpha)
                    acc += w_x*w_y*w_z*Ui*ay
        acc /= norm
        c = np.fft.fft(acc)/nphi                            # acc(phi) = sum_m c_m e^{i m phi}
        coeffs.append({m: complex(c[m % nphi]) for m in range(-8, 9)})
    return coeffs


CONFIGS = {'CW-CCW': (-1, +1), 'CCW-CW': (+1, -1), 'CW-CW': (-1, -1), 'CCW-CCW': (+1, +1), 'single': (-1, 0)}


def response(coeffs, n):
    """Complex modal-acceleration amplitude at e^{i n W t} for each configuration (peak = 2|c|)."""
    out = {}
    for name, (s1, s2) in CONFIGS.items():
        c = coeffs[0][n*s1]+(coeffs[1][n*s2] if s2 else 0)
        out[name] = 2*c
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fine', action='store_true', help='doubled quadrature, convergence check only')
    a = ap.parse_args()
    q = dict(nl=96, nc=4, nx=96, ny=3, nz=3, nphi=64)
    if a.fine:
        q = dict(nl=192, nc=6, nx=192, ny=4, nz=4, nphi=128)
    w0 = 2*np.pi*F0
    d0s = np.round(np.arange(0.290, 0.6001, 0.010), 3).tolist()
    if 0.296 not in d0s:
        d0s = sorted(d0s+[0.296])
    if a.fine:
        d0s = [0.296, 0.400, 0.600]
    res = {}
    for d0 in d0s:
        co = modal_coefficients(d0, q)
        res[d0] = {f'n{n}': {k: Q*v/w0 for k, v in response(co, n).items()} for n in (2, 4)}
        r = res[d0]['n2']
        print(f'd0={d0*1e3:.0f} mm  |v| n2 CW-CCW {abs(r["CW-CCW"])*1e9:7.2f} nm/s  CW-CW {abs(r["CW-CW"])*1e9:7.2f}'
              f'  n4 CW-CCW {abs(res[d0]["n4"]["CW-CCW"])*1e9:6.2f}', flush=True)
    rho = E_TI*H_Y**2*BETA_L**4/(12*L_B**4*w0**2)
    xb, wx = gauss(400, 0, L_B)
    mstar = rho*B_Z*H_Y*(wx@mode_shape(xb)**2)
    summary = {'quadrature': q, 'G': G, 'f0_Hz': F0, 'Q': Q, 'rho_from_E_and_f0': rho, 'modal_mass_kg': mstar,
               'beam_mass_kg': rho*L_B*B_Z*H_Y, 'velocity_m_s': {}, 'phase_deg': {}}
    for d0 in d0s:
        summary['velocity_m_s'][f'{d0:.3f}'] = {n: {k: abs(v) for k, v in r.items()} for n, r in res[d0].items()}
        summary['phase_deg'][f'{d0:.3f}'] = {n: {k: math.degrees(np.angle(v/r['CW-CW'])) for k, v in r.items()}
                                              for n, r in res[d0].items()}
    v296 = abs(res[0.296]['n2']['CW-CCW'])
    summary['at_296mm'] = {'v_n2_CWCCW_nm_s': v296*1e9, 'power_W': 0.5*mstar*w0/Q*v296**2,
                           'power_paper_W': 3.75e-18,
                           'same_vs_opposite_n2': abs(res[0.296]['n2']['CW-CW'])/v296,
                           'n4_vs_n2_CWCCW': abs(res[0.296]['n4']['CW-CCW'])/v296}
    if not a.fine:
        fit = {}
        sel = [d for d in d0s if 0.2995 <= d <= 0.6005]
        for n in ('n2', 'n4'):
            for k in CONFIGS:
                amp = np.array([abs(res[d][n][k]) for d in sel]); dd = np.array(sel)
                y = np.log(amp); x = np.log(dd)
                p = np.polyfit(x, y, 1)
                # Table II evidently fits amplitudes directly (linear residuals); the local exponent varies
                # from about 5 near 300 mm to about 4.2 at 450 mm, so the fitting choice matters.
                pl, _ = curve_fit(lambda xx, C, aa: C*xx**-aa, dd, amp, p0=(amp[0]*dd[0]**(-p[0]), -p[0]))
                r2 = 1-np.sum((amp-pl[0]*dd**-pl[1])**2)/np.sum((amp-amp.mean())**2)
                fit[f'{n}|{k}'] = {'a_amplitude_fit': pl[1], 'a_log_fit': -p[0], 'table2': TABLE2[(n, k)],
                                   'R2_amplitude_fit': r2}
        summary['power_law_300_600mm'] = fit
        # angle of each bar's component relative to their equal-direction sum (CW-CW), n = 2, 300 mm
        co = modal_coefficients(0.300, q)
        c1 = co[0][-2]; c2 = co[1][-2]
        summary['component_angle_deg_300mm'] = math.degrees(np.angle(c1/(c1+c2)))
    out = HERE/'results'/('brack_comparison_fine.json' if a.fine else 'brack_comparison.json')
    out.write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: summary[k] for k in summary if k in ('at_296mm', 'power_law_300_600mm',
                                                               'component_angle_deg_300mm', 'rho_from_E_and_f0')},
                     indent=1))


if __name__ == '__main__':
    main()
