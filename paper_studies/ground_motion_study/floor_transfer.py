"""Floor transfer from a force at the transmitter mount to the floor under the receiver (Sec. 5.6).

Infinite Kirchhoff plate (concrete slab) driven by a harmonic point force F
(Cremer, Heckl & Petersson, Structure-Borne Sound, 3rd ed.):
    w(r) = i F / (8 D k^2) * [H0^(1)(kr) + (2i/pi) K0(kr)],   k^4 = omega^2 rho h / D.
The equivalent free mass at distance r is M(r) = F / |omega^2 w(r)|; M(0) = 8 sqrt(D rho h) / omega.
The tolerable force follows from the coherent-acceleration budget at the receiver (10% of the
2 m signal, 3.6 nGal at 18 Hz) divided by the receiver isolation |H_iso(18 Hz)| = 2.9e-6.
"""
import numpy as np
from scipy.special import hankel1, k0

E, RHO, NU = 30e9, 2400.0, 0.2          # concrete
A_TOL = 3.6e-11                          # m/s^2 (3.6 nGal)
H_ISO = 2.9e-6                           # receiver isolation at 18 Hz
F_2F = 1.0                               # N, assumed twice-per-revolution force


def equivalent_mass(h, f, r):
    D = E*h**3/(12*(1 - NU**2))
    m = RHO*h
    w = 2*np.pi*f
    k = (w**2*m/D)**0.25
    g = 1.0 if r == 0 else abs(hankel1(0, k*r) + 2j/np.pi*k0(k*r))
    return 8*np.sqrt(D*m)/(w*g), 2*np.pi/k


if __name__ == '__main__':
    f, r = 18.0, 2.0
    for h in (0.15, 0.2, 0.3):
        M0, lam = equivalent_mass(h, f, 0.0)
        M, _ = equivalent_mass(h, f, r)
        F_tol = M*A_TOL/H_ISO
        print(f"h={h:.2f} m: bending wavelength {lam:4.1f} m, M(0)={M0:8.3g} kg, M({r:g} m)={M:8.3g} kg "
              f"(transfer {M0/M:.2f} of mount), F_tol={F_tol:.3f} N, extra isolation {20*np.log10(F_2F/F_tol):.0f} dB")
