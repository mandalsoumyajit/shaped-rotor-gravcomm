# Matched stresses, carrier modes and fixed-mass design study

## Scope and assumptions

Objective chosen by the user: maximize transmitted gravitational signal at fixed
total rotor mass. The study holds CAD mass at 7.91350949 kg, outer radius at
at most 250 mm, axial envelope at at most 90 mm, carrier thickness at 20 mm,
speed at 1800 rpm, bore diameter at 12 mm and nominal boss wall at 15 mm.
The signal objective is radial second-harmonic amplitude at 1 m from the axis;
0.5 m and 2 m are also recomputed. These are source fields, not receiver outputs.
The hub is fixed and inserts perfectly bonded. Contact/fasteners, shaft/bearings,
fatigue, damping, gyroscopic response and practical motor drive limits still
require separate work. This is a sampled trade study, not a global optimum.

## Location-matched analytical checks on the original shaped model

- Volume-weighted hoop stress in four mid-arc rim sectors (radius 220-245 mm;
  angle within 15 degrees of 45, 135, 225, 315 degrees) is 5.974-5.988 MPa.
  The existing thin-ring estimate using outer radius is 5.996 MPa. This validates
  the scale of nominal mid-arc hoop stress, not every point on the rim: sampled
  local stresses span about 3.2-8.7 MPa. A free-ring estimate evaluated at mean
  radius would differ; this loaded spoked ring is not an isolated thin ring.
- Each insert needs 14.1950 kN from m*omega^2*r. Consistent steel-element nodal
  force recovery gives 14.1945 and 14.1946 kN, within 0.004%. It converges from
  about 0.037% error on the coarse mesh. The recovery integrates B-transpose
  times elastic stress and subtracts consistent centrifugal nodal loading.
- A geometry-matched projected area is insert diameter times engagement
  thickness: 60*20=1200 mm^2. The associated NOMINAL load/area is 11.829 MPa.
  The older 20.279 MPa estimate used rim width*thickness=700 mm^2, a different
  gross-area proxy. Neither number is a local bonded/contact peak stress.
- Extrapolating elemental stresses onto the curved interface gave force errors
  of several percent and lacked satisfactory convergence. Those traction
  estimates are retained in JSON as rejected diagnostics, not evidence of
  validated local contact stress. The consistent-force method avoids that
  extrapolation, but still does not turn a bonded model into a contact model.

## Carrier modes

Modes below are for the original fixed-bore carrier at 1800 rpm, including
centrifugal prestress. They are not shaft-whirl critical speeds. Shapes are
mass-normalized and their amplitudes are arbitrary, not predicted vibration.

| Mode | Hz | Identification |
|---|---:|---|
| 1 | 78.58 | Predominantly axial bending: opposite inserts move in opposite axial directions |
| 2 | 93.73 | In-plane torsional motion of carrier relative to fixed hub |
| 3 | 102.98 | Predominantly axial bending: inserts move together |
| 4 | 123.94 | Axial bending dominated by the unweighted transverse spokes/rim |
| 5 | 232.66 | Higher axial bending with opposing weighted/unweighted sectors |
| 6 | 267.84 | Higher axial twisting/saddle-like carrier deformation |
| 7 | 467.79 | Higher in-plane bending with inserts moving mainly together tangentially |
| 8 | 569.40 | Higher in-plane rim/spoke deformation |

The second mode accounts for about 99.86% of effective polar inertia for a
distributed angular-acceleration excitation in this fixed-hub model. It is
therefore more directly relevant to speed modulation than the lowest axial
mode. This is a coupling calculation, not a forced-response or BER prediction.
Quadrature reconstruction of all eight generalized modal masses agrees with
the solver normalization to about 1e-7. See `results/stress_modes/carrier_modes.png`.

## Corners and local mesh limitations

Explicit 3 mm in-plane junction fillets were added. Their small mass change is
compensated by adjusting insert height, retaining fixed total CAD mass. The
through-bore and bonded pocket axial edges are not fully detailed hardware.
An initial curved quadratic mesh contained inverted elements and was rejected.
Final rounded cases use quadratic displacement elements with straight midside
geometry, local junction refinement and resolved circular boundaries. Every
accepted mesh passes positive-Jacobian checks. All bore nodes are constrained
by CAD surface membership, not by testing sampled node radius.

The earlier 38.9 MPa peak is NOT a converged stress validation: introducing local
refinement changed the peak substantially. Comparison with it must not be used
to claim that rounding increases stress. Global mesh-size agreement had missed
local stress sensitivity. Raw peaks are retained only as diagnostics in
`results/design_sweep/stress_regions.csv`.

The primary comparison uses volume-weighted average von Mises stress in fixed
5 mm and 10 mm radius circular neighbourhoods, spanning the full plate thickness
and including aluminium only. Windows are centered on nominal hub/spoke,
inner insert-boss/spoke and unweighted rim/spoke junctions, plus mid-arc rim
locations. Four symmetry-related windows per group are evaluated separately;
the table reports the largest regional average, not a pointwise maximum.
Centers remain fixed in physical coordinates between meshes. The CSV also
reports actual sampled volumes, sample counts and von Mises stress of the
averaged tensor. The latter can be lower because differently oriented stresses
cancel; it must not be substituted silently for the average of von Mises stress.
See `results/design_sweep/averaged_stress_regions.csv` for all windows.
Regional averages improve robustness for design comparisons but can conceal
local yielding or fatigue initiation; they are not alone a strength criterion.

## Fixed-mass rounded-design comparisons

All dimensions below are mm, insert mass is per insert. Signal includes the
complete finite aluminium and steel geometry. CAD mass is exact by construction;
quadrature/mesh mass error is separately recorded in the raw CSVs.

| Design | Spoke | Rim | Insert kg | A2 at 1 m, nGal | Gain | First mode Hz | Mean VM, 5 mm MPa | Mean VM, 10 mm MPa | 10 mm mean change, coarse to fine |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline_rounded | 25 | 35 | 1.997 | 2.4378 | +0.00% | 78.50 | 26.14 | 22.82 | +0.14% |
| medium_rounded | 30 | 20 | 2.394 | 2.9594 | +21.39% | 79.87 | 29.40 | 25.48 | +0.16% |
| wide_rounded | 40 | 20 | 2.214 | 2.6917 | +10.41% | 89.91 | 21.51 | 17.58 | -0.13% |

The controlling windows in this table are at the weighted hub/spoke junctions;
their 10 mm averages change by less than 0.2%. This does NOT imply that every
window has converged equally well. Across all individual windows, the largest
change is 5.59% for 10 mm windows and 10.15% for 5 mm windows.
Some small windows outside the locally refined junctions have sparse quadrature
sampling and appreciable sampled-volume changes. The 10 mm averages are the
primary comparison; `averaged_stress_convergence.csv` records each paired window
and its sampled-volume change so these limitations remain visible.

The final faceted meshes have mass errors from -0.138% to
-0.097% relative to exact CAD mass. Gravity weights are not
renormalized. This discretization error should be kept in mind when comparing
small signal changes; the larger design gains are well above that scale.

The 30 and 40 mm spoke alternatives use 70 mm diameter cylinders at 200 mm
radius, versus 60 mm cylinders in the original design. Reducing rim width
frees material budget for the signal-bearing steel; widening spokes can regain
stiffness. The 20 mm spoke, 20 mm rim sharp candidate gives about 32% more signal
but its first mode falls to about 68 Hz and its raw peak rises substantially.
It is not automatically the preferred design.

Sphere/cylinder isolation check: with the same 180 mm center radius, 20 mm spokes,
20 mm rim and total mass, spherical inserts change A2 at 1 m by -0.190%
relative to the cylinder control. This matched comparison shows essentially
the same signal, with finite-size corrections and carrier/pocket geometry
accounted for. Their larger radial footprint also limits how far outward
they fit with a 15 mm boss wall. A separately tested sphere near its radial
limit has center radius about 192.35 mm, not 200 mm.

For further optimization, vary continuous placement, cylinder aspect ratio,
spoke taper/thickness, rim width and fillet radius together. Preserve the full
finite-volume signal objective and include constraints on locally converged
stress, deformation, modal response and attachment design. Maximizing nominal
quadrupole alone pushes toward a light flexible carrier and is insufficient.

## Redone gravitational waveforms

`gravcomm.FiniteRotor` integrates the rotating full-body Newtonian field using
positive element quadrature masses. The carrier is not assumed silent. Its
fourfold structure contributes higher harmonics; bosses/pockets also affect
the second harmonic. The rotor's polar inertia and signal quadrupole anisotropy
are distinct quantities. Complex material harmonic coefficients must be summed
before taking magnitude.

For the original shaped geometry at 1 m, the full-body second-harmonic amplitude
is 2.4455 nGal, compared with
2.4690 nGal for the two-point
insert approximation (-0.95%). The aluminium contribution at the
second harmonic partly opposes the steel contribution; its fourth harmonic
adds to the steel fourth harmonic. The correction is small in this example,
but it is now calculated rather than assumed away.

`results/shaped_waveforms/` contains full phase traces and spectra at 0.5, 1 and
2 m, 256 mechanical phases, reusable NPZ source models, the separate material
harmonics and a comparison plot. Rotation is about z and the receiver lies on x.
These waveforms are for undeformed geometry at constant 1800 rpm; the same
FiniteRotor can be passed to the actuator trace to evaluate its achieved phase.
They do not include a receiver transfer function, noise or a demonstrated bit rate.

Validation includes the point-mass limit, phase/sign convention, axisymmetric
ring cancellation, actuator compatibility, C3D10 affine-strain reconstruction,
stress-component ordering, modal parsing, mass checks, mechanical mesh comparisons,
and a 128/256 phase comparison for the volume-integrated signal. The existing
model suite passes 51 tests and the FEA parsing/interpolation/averaging suite passes five.

Sources for implementation conventions: [CalculiX manual](https://www.dhondt.de/ccx_2.21.pdf),
[CalculiX Gauss order](https://github.com/Dhondtguido/CalculiX/blob/master/src/gauss.f),
[Gmsh API](https://gmsh.info/doc/texinfo/).
