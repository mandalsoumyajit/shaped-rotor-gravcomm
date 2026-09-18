# Finite-volume gravitational waveforms

NPZ files contain SI quadrature positions and positive mass weights for the
complete aluminium carrier and steel inserts. Load with
`gravcomm.FiniteRotor.from_npz(path)`. Rotation is about z, receiver at (d,0,0),
and x acceleration is radial. AC traces subtract their spatial-phase mean;
the field method itself returns the full field including DC. No receiver noise,
filtering or deformation correction is included. Rotor speed is 1800 rpm.

The selected second-harmonic signal is at 60 Hz. Higher even harmonics are
reported separately; they must not be folded into its signal power. Tiny odd
harmonics can reflect mesh asymmetry. The material CSV gives COMPLEX component
coefficients: add these before taking magnitudes. The carrier is not generally
gravitationally silent. These files replace the point-insert approximation for
this shaped demonstrator; they do not redefine the separate deployment classes.

Quadrature: 4 points per quadratic tetrahedron, including curved Jacobians.
Waveforms: 256 mechanical phases per turn, distances 0.5, 1 and 2 m from axis.
The mesh-convergence tables and 128/256 phase check are in ../design_sweep/.
