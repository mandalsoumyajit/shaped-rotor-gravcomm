# CalculiX verification work

**Current paper default:** `medium_rounded`, a 7.9135 kg, 500 mm diameter
carrier with 30 mm spokes, 20 mm rim and 3 mm junction fillets. Start with
the repository reproduction README and `STRESS_MODES_AND_DESIGN.md`.
The next sections describe the original unfilleted carrier used for the
analytical checks; they are not the dimensions of the selected rounded design.

The intended rotor is a **shaped spoked carrier**, constructed by the Gmsh scripts,
not the paper's point-mass dumbbell or a literal tie rod. `shaped_rotor.py` models
that carrier and two finite steel inserts. `run_calculix.py` now runs only a
rotating-bar benchmark to check the toolchain against analytical expressions.

## Environment and commands

Installed in Ubuntu WSL: CalculiX 2.21, Gmsh 4.12.1, Ubuntu system Python,
NumPy and Matplotlib. No Windows-native solver or GUI is required.

From WSL, in this directory:

```sh
python3 -s run_calculix.py
python3 -s shaped_rotor.py
python3 -s summarize.py
python3 -s -m unittest discover -s . -p test_fea.py -v
```

Each simulation script accepts `--quick` for one mesh resolution. Full runs
use three resolutions and both stationary and rotating cases. Outputs include
the mesh, CalculiX input, solver log, `.dat`, `.frd`, metadata and summary CSV.
Gmsh can open `mesh.msh` or the shaped model's `geometry.brep`; CalculiX-compatible
postprocessors can inspect `.frd`. Inputs and results use SI units throughout.

## Shaped-carrier assumptions

- Aluminium carrier: 250 mm outer radius, 20 mm thickness, 35 mm rim width,
  four 25 mm wide spokes, 30 mm hub radius, and a 12 mm diameter through bore.
- Two steel cylinders: 60 mm diameter, 90 mm axial length, centered at opposite
  200 mm radii, with 90 mm diameter aluminium bosses. Each insert is 1.998 kg.
- Material properties: aluminium E=70 GPa, nu=0.33, rho=2700 kg/m^3;
  steel E=200 GPa, nu=0.30, rho=7850 kg/m^3. Homogeneous isotropic linear elasticity.
- 0 and 1800 rpm about the shaft axis (60 Hz gravitational second harmonic).
- The hub bore is fully fixed. Shaft, bearings, coupling and setscrew are omitted.
  In the original OpenSCAD boolean ordering, spokes can refill the hub bore;
  this FEA applies the intended through bore after uniting the carrier geometry.
- The CAD's 0.4 mm radial pocket clearance is removed for **ideal bonded,
  conformal insert interfaces**. No contact, fasteners, fit or insert-retention
  hardware is validated. The idealization must be revisited before hardware use.
- Sharp spoke/boss junctions have no fillets. Raw local maximum stresses may
  increase under refinement and must not be treated as converged strength margins.

All solids use quadratic tetrahedra (C3D10). Gmsh-to-CalculiX midside-node order
is determined from reference-element coordinates. Steel and aluminium use
separate material element sets and share interface nodes.

Linear centrifugal statics use `CENTRIF` with **omega squared**, in rad^2/s^2.
A subsequent `*STEP,PERTURBATION` frequency step includes the static preload.
These are fixed-bore structural modes, not a complete gyroscopic rotor-bearing
analysis or a Campbell diagram. Stress output from normalized mode shapes is
excluded from the physical static stress report.

The geometry's mass and polar inertia are computed by OpenCASCADE volume
integration. The carrier is not axisymmetric (especially the two bosses), so
its finite-body gravitational signal cannot simply be declared zero. This FEA
does not change the gravitational model or claim a communication rate.

## Analytical benchmark

Steel bar from radius r0=0.03 m to r1=0.30 m, 10 mm square section, fixed at
its inner face. At omega=250/0.3 rad/s, one-dimensional axial theory gives:

```
root force = rho*A*omega^2*(r1^2-r0^2)/2
tip extension = rho*omega^2/(2E) * [r1^2*(r1-r0) - (r1^3-r0^3)/3]
stationary first bending frequency = 1.875104^2/(2*pi*L^2)*sqrt(E*I/(rho*A))
```

Small differences are expected from three-dimensional elasticity, transverse
centrifugal loading, the fixed-face Poisson constraint, and mesh discretization.
The bar checks centrifugal units, node ordering, stiffness/mass assembly, and
the solver workflow; it is not the proposed transmitter geometry.

## Preliminary files

Folders immediately under `results/` named `rotor_*` are **superseded exploratory
arm-and-block calculations**, stopped after the user's geometry correction.
They are not the shaped design and must not be cited as its results. The earlier
`bar_h0.008_*` quick-run metadata predates a parser correction; use the
three-resolution benchmark results instead. Active results belong in
`results/benchmark/` and `results/shaped_rotor/`.

## Follow-up studies

`postprocess_carrier.py` compares rim hoop stress in specified sectors and
recovers consistent insert-interface nodal forces. It also identifies carrier
modes from mass-weighted displacement components and angular-acceleration
participation. Integration-point stresses are used; FRD nodal stress averaging
across materials is avoided. Extrapolated local interface tractions are retained
as diagnostics but fail a tight equilibrium/convergence check and are not used
as validated contact stresses.

`design_sweep.py` screens shapes at a fixed CAD mass of 7.91350949 kg, 250 mm
outer-radius limit, 90 mm axial envelope and 1800 rpm. It varies sphere/cylinder
inserts, placement, spoke and rim widths, adjusting insert size to maintain mass.
The objective is radial second-harmonic acceleration at a 1 m receiver radius;
0.5 m and 2 m are also evaluated. Exact CAD quadrupole anisotropy screens the
grid, then selected candidates use finite-volume gravity and structural FEA.
This is a bounded trade study, not a global optimum or a hardware qualification.

`refine_designs.py` refines selected designs and checks phase quadrature.
`filleted_designs.py` adds 3 mm in-plane junction fillets and maintains mass by
slightly resizing the inserts. Axial pocket-edge radii and insert-fastener/contact
details are still absent. The final `*_rounded` models use quadratic displacement
tetrahedra with midside nodes on straight edges (faceted geometric approximation)
to avoid inverted curved elements near small fillets. Targeted circular-boundary
resolution, local refinement, and positive-Jacobian checks are applied. Bore constraints use
CAD surface membership, including the midside nodes, not radius equality tests.
The earlier `*_filleted` curved-element attempts include a rejected invalid mesh;
they are diagnostic history and are excluded from final design comparisons.

`shaped_waveforms.py` regenerates full-carrier-plus-insert phase traces, spectra,
and reusable mass-quadrature NPZ files for `gravcomm.FiniteRotor`. In this model
rotation is about z, x is the radial receiver component, y transverse, z axial.
The earlier two-point `TwoMassRotor` uses x-z geometry; radial x comparisons are
equivalent, but transverse-component labels differ.

`static_regions.py` reports volume-weighted stress averages over fixed physical
5 mm and 10 mm radius windows across the full aluminium thickness. Both the
average von Mises stress and von Mises stress of the averaged tensor are saved,
along with sampled volumes. Raw peaks remain diagnostic. `design_report.py`
uses regional averages as the primary stress comparison.

Use Ubuntu's system packages consistently with `python3 -s` for these scripts.
The solver now runs in a temporary Linux directory and copies every output back
to the workspace, avoiding slow per-record writes through `/mnt/c`.

- [CalculiX 2.21 manual](https://www.dhondt.de/ccx_2.21.pdf): DLOAD, FREQUENCY,
  STEP/PERTURBATION and complex frequency analyses.
- [Gmsh manual](https://gmsh.info/doc/texinfo/): OpenCASCADE geometry and Python API.
