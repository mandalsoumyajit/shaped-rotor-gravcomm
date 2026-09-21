# Rotor orientation audit

The source rotates about z. The existing link model places the receiver at (d,0,0) and measures x acceleration: in the rotor plane, radially toward the axis. This is not an arbitrary low-coupling direction.

An independent finite-volume calculation evaluated all three acceleration components at receiver elevations 0,5,15,30,45,60,75,90 degrees above the rotor plane, for equivalent unscaled distances 0.500,0.571,0.625 m. These represent the earlier three-distance rotor sizes. At each position it extracted the complex second-harmonic vector h and maximized |u dot h| over real unit sensor axes u by the largest eigenvalue of Re(h)Re(h)^T+Im(h)Im(h)^T. This allows the sensitive axis to change at every position rather than unfairly holding it radial away from the plane.

Every tested distance has its largest second-harmonic amplitude at zero elevation, with the optimal axis radial to numerical precision. For the 0.500 m equivalent distance, normalized best amplitudes are 1.0000,0.9887,0.9045,0.6829,0.4470,0.2555,0.1133,approximately zero. At 0.625 m, the corresponding 15/30-degree ratios are 0.9219/0.7282. Small axial residuals reflect quadrature asymmetry. This is a sampled orientation audit, not a continuous global optimization over arbitrary geometry.

Analytical cross-check: for the leading rotating quadrupole, with polar angle theta from the rotation axis, the radial and polar cosine-quadrature coefficients are proportional to 3 sin(theta)^2 and 2 sin(theta)cos(theta), while the azimuthal sine-quadrature coefficient is proportional to 2 sin(theta). The optimal fixed real-axis amplitude is therefore proportional to max[sqrt(9 sin(theta)^4+4 sin(theta)^2 cos(theta)^2),2 sin(theta)], maximized at theta=pi/2 with a radial sensor axis. At the pole the intended time-varying quadrupole signal vanishes. Receiver azimuth or rotor initial angle only changes phase for steady rotation.

Thus no signal-amplitude gain from tilting the rotor was found. A deployed orientation should nevertheless maximize SNR, which also depends on directional support/gravity noise and installation constraints. The current scalar background model cannot establish that SNR-optimal site orientation. This audit concerns the communication second harmonic, not the static monopole field or instantaneous total acceleration.

Reproduce with scripts/audit_rotor_orientation.py; raw vector phasors, optimizing axes and amplitudes are in orientation_audit.json. Model uses 64 mechanical phase samples and the saved finite-volume reference rotor.

Verification: equatorial radial results at all three distances agree with the existing 128-phase finite-volume field calculation within relative tolerance 1e-6.
