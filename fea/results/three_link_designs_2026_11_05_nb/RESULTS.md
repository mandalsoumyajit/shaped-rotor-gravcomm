# Structural reassessment of the hop-limited CPFSK 0.5, 1 and 2 m sources

Design points (2026-11-05): 5.79 kg at 0.5 m, 30.98 kg at 1 m, 184.05 kg at 2 m. Commanded gravitational frequency range carrier +- 0.545306 Hz; full-symbol (224/190 s) quintic hops of at most 3 tones; peak alpha = pi*0.867212 rad/s^2. All quantities below are calculations for the bare rotor with fixed bore and perfectly bonded inserts unless explicitly identified as a separate assembly screen.

## Fixed-bore linear FEA

Twenty-six static/prestressed-modal cases use two archived mesh resolutions, three actual operating speeds per design, and four independently solved gravity/angular-acceleration bases per mesh. No earlier cases were reused; all 26 cases were solved with CalculiX 2.21 (Ubuntu 24.04 package, user-local install, single-threaded). All solver jobs must finish successfully; static/modal outputs are separated. Five existing FEA regression tests pass. Positive Jacobians, exact consistent-load torque, force/torque equilibrium and similarity to the archived desktop solution are checked.

| Link (m) | Bare mass (kg) | Upper-speed centrifugal Al peak (MPa) | Horizontal combined Al bound (MPa) | Centrifugal displacement (um) | Horizontal combined displacement bound (um) |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 5.79 | 7.90 | 8.54 | 8.34 | 11.14 |
| 1 | 30.98 | 24.16 | 25.40 | 44.63 | 56.05 |
| 2 | 184.05 | 45.25 | 47.85 | 151.35 | 206.89 |

The combined columns are triangle-inequality bounds assembled from the linear static solutions, including maximum centrifugal load, either sign of maximum angular acceleration and any in-plane gravity direction. They are deliberately conservative; speed and angular-acceleration maxima need not occur simultaneously. Separate phase-sampled 3-tone full-symbol quintic hops use the actual trajectory. Neither calculation is a nonlinear contact simulation or a complete time-domain rotor-bearing response. Gravity/acceleration bases use unprestressed elastic stiffness, consistent with the original linear statics model.

Horizontal shaft orientation is the preliminary choice. All-attitude bounds (including axial gravity) are substantially larger and retained in compiled.json. Gravity direction rotates at the shaft frequency in rotor coordinates for a horizontal shaft. Loads at the 24/24/18 Hz gravitational carriers are not automatically mechanical excitation at those frequencies.

## Mesh convergence and strength interpretation

| Link (m) | Change in centrifugal raw Al peak (%) | Change in max displacement (%) | Change in maximum scaled-window combined average (%) |
|---:|---:|---:|---:|
| 0.5 | 15.85 | 0.089 | 0.226 |
| 1 | 15.85 | 0.089 | 0.229 |
| 2 | 15.85 | 0.089 | 0.225 |

Regional windows have reference radii 5/10 mm multiplied by the geometric scale; they span the carrier thickness. The reported average is the maximum of the 10 mm-reference windows. It is a comparison metric, not a replacement for local yield or fatigue analysis. Fine-mesh mass is about 0.11% below CAD because of the faceted approximation; density is not renormalized. See mesh_h*.json.

Third-mesh check (archived 2026-09-20 run at s=1, 1800 rpm; scale independent): a 10 mm reference mesh has 132043 nodes and 83555 elements. Its 1800 rpm aluminum peak is 58.492 MPa, a 0.628% increase from the 13 mm mesh. Displacement is 68.127 um and first frequency 79.825 Hz. This supports convergence of the centrifugal reference peak; it does not establish convergence of every combined-load or contact stress.

Preliminary material choices are 6061-T651 plate (240 MPa yield screen), specifically certified normalized 1045 inserts (300 MPa minimum procurement requirement), and 4140 Q&T shafts/flanges (650 MPa minimum study specification). A factor of two against yield is the initial static screening criterion. Material sources, product-form caveats and actual assembly choices are in ASSEMBLY_SPECIFICATION.md and materials.json. Raw mesh peaks and idealized attachments do not establish hardware safety factors. Fatigue, preload, real contact and local features remain unresolved.

All three idealized horizontal-shaft combined aluminum bounds are below the 120 MPa screening limit (240 MPa divided by two); the largest is about 47.85 MPa. The model therefore indicates useful static margin in the carrier, while insert retention and assembly dynamics still require explicit validation.

## Modes, coupling and response

| Link (m) | First three center-speed modes (Hz) | Dominant angular-acceleration mode | Polar-inertia fraction (%) | Shaft/carrier series torsion screen (Hz) |
|---:|:---|---:|---:|---:|
| 0.5 | 82.36, 100.90, 112.29 | 3 | 99.76 | 50.89 |
| 1 | 48.22, 58.57, 64.23 | 3 | 99.76 | 37.44 |
| 2 | 27.39, 32.94, 35.49 | 3 | 99.76 | 24.17 |

Mode matching uses the consistent mass matrix and modal assurance criterion (MAC), with node correspondence within each mesh. Mass normalization is independently recovered from degree-four exact element mass integrals. Component energies and mean insert motions identify axial/torsional motion. Angular and translational effective mass fractions distinguish modulation coupling from gravity coupling; low frequency alone does not establish strong excitation.

modal_response_h*.json contains damping sweeps at 0.5%, 2% and 5% for isolated quintic-transition scalar modal responses and rotation-frequency harmonic transfer factors. The transition simulation starts at rest and omits repeated-message resonance, motor control, gyroscopic terms and contact. It diagnoses time-scale separation, not arbitrary-message dynamic qualification. Shaft/carrier series values are a simple compliance screen using the chosen stub-shaft dimensions and dominant carrier mode; motor/coupling inertia and bearing/structure dynamics are absent. These values demonstrate why fixed-bore carrier frequencies must not be called assembly critical speeds.

## Assembly and next checks

ASSEMBLY_SPECIFICATION.md selects a horizontal straddled shaft, locating angular-contact pair plus floating cylindrical-roller support, piloted hub flanges and positive insert shoulder/collar retention. Plain-shaft and pocket-pressure estimates are saved separately in assembly_screen.json. These are preliminary calculations, not additional 3D assembly FEA. Before physical implementation, add flange bolts, contact/clearance, retention grooves/collars, preload, actual bearing and motor/coupling stiffness, and fatigue. Added hardware mass and off-axis retention features require a revised field/inertia audit before reusing the BER claims for that assembly.

The tungsten-heavy-alloy option in TUNGSTEN_OPTION.md retains insert mass and diameter while shortening its length; it is a promising windage/overhang reduction study. It has not been substituted into the qualified steel models.

## Reproduction

With GRAVCOMM_FEA_TAG=2026_11_05_nb: run fea/make_hop_frozen_designs.py, then fea/three_link_fea_hop.py --prepare --workers 8 in the WSL environment, then compile_three_link_fea_hop.py --mesh .018 and --mesh .013; three_link_modal_response_hop.py for each mesh; assembly_screen_hop.py; report_three_link_fea_hop.py. Inputs, complete solver outputs, fields, equilibrium checks and hashes are retained. No BER samples or source configuration was changed.

<!-- side-by-side comparison generated by summarize_hop_vs_archived.py -->

## Side-by-side: archived 2026-09-20 vs 2026_11_05_nb (fine 13 mm mesh unless stated)

Each cell is archived / new. Combined bounds are the horizontal-shaft triangle-inequality bounds of Appendix A.

| Quantity | 0.5 m | 1 m | 2 m |
|:---|---:|---:|---:|
| Scale s | 0.961185 / 0.901044 | 1.682515 / 1.576114 | 3.062574 / 2.854448 |
| Carrier f_g (Hz) | 24 / 24 | 24 / 24 | 18 / 18 |
| Max rotation (rpm) | 733.018 / 736.359 | 733.018 / 736.359 | 553.018 / 556.359 |
| Peak alpha (rad/s^2) | 6.1619 / 2.7244 | 6.1619 / 2.7244 | 6.1619 / 2.7244 |
| CAD mass (kg) | 7.0273 / 5.7890 | 37.6917 / 30.9836 | 227.3153 / 184.0498 |
| CAD polar inertia (kg m^2) | 0.25931 / 0.18772 | 4.26167 / 3.07413 | 85.15649 / 59.89568 |
| Peak rigid torque (N m) | 1.5979 / 0.5114 | 26.2600 / 8.3752 | 524.7264 / 163.1814 |
| Fine-mesh mass (kg) | 7.0196 / 5.7827 | 37.6505 / 30.9497 | 227.0667 / 183.8485 |
| Centrifugal Al peak, fine (MPa) | 8.91 / 7.90 | 27.29 / 24.16 | 51.46 / 45.25 |
| Centrifugal Al peak, coarse (MPa) | 7.69 / 6.82 | 23.55 / 20.86 | 44.42 / 39.05 |
| Horizontal combined Al bound (MPa) | 9.72 / 8.54 | 28.98 / 25.40 | 55.50 / 47.85 |
| Centrifugal displacement (um) | 10.03 / 8.34 | 53.80 / 44.63 | 184.69 / 151.35 |
| Horizontal combined displacement bound (um) | 15.26 / 11.14 | 77.77 / 56.05 | 314.95 / 206.89 |
| Angular-acceleration-only Al peak (MPa) | 0.7681 / 0.2984 | 2.3535 / 0.9131 | 7.7979 / 2.9951 |
| In-plane gravity Al bound (MPa) | 0.6040 / 0.5662 | 1.0572 / 0.9904 | 1.9244 / 1.7936 |
| All-attitude combined Al bound (MPa) | 12.33 / 11.00 | 33.50 / 29.65 | 63.69 / 55.52 |
| All-attitude combined displacement bound (um) | 51.68 / 43.14 | 189.35 / 153.96 | 684.64 / 528.02 |
| Phase-sampled transition Al bound, any attitude (MPa) | 12.18 / 10.93 | 33.03 / 29.44 | 62.30 / 54.90 |
| All-attitude steel-insert bound (MPa) | 4.97 / 4.40 | 15.12 / 13.38 | 28.54 / 25.06 |
| Regional 10 mm-ref. combined average, fine (MPa) | 4.762 / 4.222 | 13.569 / 11.962 | 25.851 / 22.455 |
| Regional 10 mm-ref. combined average, coarse (MPa) | 4.752 / 4.212 | 13.539 / 11.935 | 25.796 / 22.405 |
| Mesh change, raw centrifugal Al peak (%) | 15.85 / 15.85 | 15.85 / 15.85 | 15.85 / 15.85 |
| Mesh change, max displacement (%) | 0.089 / 0.089 | 0.089 / 0.089 | 0.089 / 0.089 |
| Mesh change, regional combined average (%) | 0.211 / 0.226 | 0.221 / 0.229 | 0.214 / 0.225 |
| Center-speed first three modes (Hz) | 77.33, 94.69, 105.27 / 82.36, 100.90, 112.29 | 45.39, 55.04, 60.18 / 48.22, 58.57, 64.23 | 25.76, 30.88, 33.09 / 27.39, 32.94, 35.49 |
| Torsional (max angular-coupling) mode no. | 3 / 3 | 3 / 3 | 3 / 3 |
| Torsional mode, center (Hz) | 105.27 / 112.29 | 60.18 / 64.23 | 33.09 / 35.49 |
| Torsional effective polar-inertia fraction (%) | 99.76 / 99.76 | 99.76 / 99.76 | 99.76 / 99.76 |
| Torsional mode, low speed, MAC-matched (Hz) | 105.270 / 112.292 | 60.175 / 64.230 | 33.084 / 35.488 |
| Torsional mode, high speed, MAC-matched (Hz) | 105.272 / 112.295 | 60.179 / 64.235 | 33.089 / 35.494 |
| Min MAC low/high vs center (all 12 modes) | 1.000000 / 1.000000 | 1.000000 / 1.000000 | 1.000000 / 1.000000 |
| Transition duration T (s) | 0.829630 / 1.178947 | 0.829630 / 1.178947 | 0.829630 / 1.178947 |
| First-mode period (ms) | 12.93 / 12.14 | 22.03 / 20.74 | 38.82 / 36.51 |
| Torsional period (ms) | 9.50 / 8.91 | 16.62 / 15.57 | 30.22 / 28.18 |
| T / first-mode period | 64.2 / 97.1 | 37.7 / 56.8 | 21.4 / 32.3 |
| T / torsional period | 87.3 / 132.4 | 49.9 / 75.7 | 27.4 / 41.8 |
| Isolated transition max modal peak/static, zeta=0.005 | 1.000164 / 1.000026 | 1.000515 / 1.000103 | 1.001205 / 1.000766 |
| Isolated transition max modal peak/static, zeta=0.02 | 1.000101 / 1.000043 | 1.000313 / 1.000119 | 1.000823 / 1.000459 |
| Isolated transition max modal peak/static, zeta=0.05 | 1.000097 / 1.000042 | 1.000284 / 1.000124 | 1.000856 / 1.000390 |
| Gravity-rotation amplification, mode 1, zeta=0.02 | 1.0247 / 1.0217 | 1.0751 / 1.0660 | 1.1389 / 1.1209 |
| Shaft-only torsion screen (Hz) | 48.58 / 57.09 | 39.13 / 46.08 | 27.69 / 33.01 |
| Shaft/carrier series torsion screen (Hz) | 44.11 / 50.89 | 32.81 / 37.44 | 21.23 / 24.17 |

### Same-mass tungsten inserts

| Quantity | 0.5 m | 1 m | 2 m |
|:---|---:|---:|---:|
| Horizontal Al bound, steel (MPa) | 9.72 / 8.54 | 28.98 / 25.40 | 55.50 / 47.85 |
| Horizontal Al bound, tungsten (MPa) | 9.72 / 8.55 | 28.99 / 25.41 | 55.52 / 47.87 |
| Horizontal displacement bound, steel (um) | 15.26 / 11.14 | 77.77 / 56.05 | 314.95 / 206.89 |
| Horizontal displacement bound, tungsten (um) | 15.19 / 11.08 | 77.39 / 55.74 | 313.67 / 205.82 |
| All-attitude insert bound, steel (MPa) | 4.97 / 4.40 | 15.12 / 13.38 | 28.54 / 25.06 |
| All-attitude insert bound, tungsten (MPa) | 3.80 / 3.37 | 11.58 / 10.24 | 21.85 / 19.19 |
| W centrifugal Al peak, coarse (MPa) | 7.58 / 6.72 | 23.21 / 20.56 | 43.77 / 38.49 |
| W centrifugal Al peak, fine (MPa) | 8.91 / 7.90 | 27.30 / 24.17 | 51.48 / 45.26 |
| W mesh change, displacement (%) | 0.092 / 0.092 | 0.092 / 0.092 | 0.092 / 0.092 |
| W mesh change, regional (%) | 0.203 / 0.218 | 0.223 / 0.230 | 0.215 / 0.226 |
| W center-speed first three modes (Hz) | 77.98, 95.50, 105.27 / 83.04, 101.77, 112.29 | 45.78, 55.52, 60.18 / 48.63, 59.08, 64.23 | 25.99, 31.16, 33.09 / 27.63, 33.23, 35.49 |
| W torsional mode (Hz) | 105.27 / 112.29 | 60.18 / 64.23 | 33.09 / 35.49 |
| W torsional inertia fraction (%) | 99.76 / 99.76 | 99.76 / 99.76 | 99.76 / 99.76 |
| Steel insert length (mm) | 76.16 / 71.40 | 133.32 / 124.89 | 242.68 / 226.19 |
| Tungsten insert length (mm) | 33.22 / 31.14 | 58.14 / 54.47 | 105.84 / 98.64 |
| Steel second harmonic (1e-10 m/s^2) | 4.20896 / 3.00455 | 4.20896 / 3.00455 | 5.18026 / 3.61064 |
| Tungsten second harmonic (1e-10 m/s^2) | 4.24586 / 3.02680 | 4.23517 / 3.02049 | 5.20591 / 3.62578 |
| Tungsten signal change (%) | 0.877 / 0.740 | 0.623 / 0.531 | 0.495 / 0.419 |

Tungsten 10 mm third mesh (archived, 0.5 m rotor, upper speed): Al peak 8.9622 MPa vs 8.9090 MPa on 13 mm (+0.597%); displacement 9.9532 um.

Tungsten 10 mm third mesh (2026_11_05_nb, 0.5 m rotor, upper speed): Al peak 7.9477 MPa vs 7.9005 MPa on 13 mm (+0.597%); displacement 8.2743 um.

### Similarity check (new / archived ratio vs exact geometric-similarity prediction)

| Link | Quantity | Ratio | Predicted | Relative error |
|---:|:---|---:|---:|---:|
| 0.5 | steel.centrifugal_alu_peak_mpa | 0.886804 | 0.886804 | -1.46e-07 |
| 0.5 | steel.centrifugal_u_um | 0.831317 | 0.831317 | +1.04e-07 |
| 0.5 | steel.angular_only_alu_mpa | 0.388542 | 0.388542 | +0.00e+00 |
| 0.5 | steel.angular_only_u_um | 0.364230 | 0.364230 | +2.22e-16 |
| 0.5 | steel.horizontal_gravity_alu_bound_mpa | 0.937430 | 0.937430 | +2.22e-16 |
| 0.5 | steel.modes.torsional_hz | 1.066707 | 1.066747 | -3.74e-05 |
| 0.5 | mass | 0.823789 | 0.823789 | +2.22e-16 |
| 0.5 | inertia | 0.723925 | 0.723925 | -3.33e-16 |
| 1 | steel.centrifugal_alu_peak_mpa | 0.885538 | 0.885538 | +7.29e-08 |
| 1 | steel.centrifugal_u_um | 0.829537 | 0.829537 | -7.96e-08 |
| 1 | steel.angular_only_alu_mpa | 0.387987 | 0.387987 | +2.22e-16 |
| 1 | steel.angular_only_u_um | 0.363451 | 0.363451 | +0.00e+00 |
| 1 | steel.horizontal_gravity_alu_bound_mpa | 0.936760 | 0.936760 | -2.22e-16 |
| 1 | steel.modes.torsional_hz | 1.067387 | 1.067509 | -1.14e-04 |
| 1 | mass | 0.822026 | 0.822026 | +2.22e-16 |
| 1 | inertia | 0.721344 | 0.721344 | -1.11e-16 |
| 2 | steel.centrifugal_alu_peak_mpa | 0.879232 | 0.879232 | -1.36e-07 |
| 2 | steel.centrifugal_u_um | 0.819481 | 0.819481 | -7.05e-08 |
| 2 | steel.angular_only_alu_mpa | 0.384088 | 0.384088 | -1.11e-16 |
| 2 | steel.angular_only_u_um | 0.357987 | 0.357987 | +2.22e-16 |
| 2 | steel.horizontal_gravity_alu_bound_mpa | 0.932042 | 0.932042 | +0.00e+00 |
| 2 | steel.modes.torsional_hz | 1.072670 | 1.072913 | -2.26e-04 |
| 2 | mass | 0.809667 | 0.809667 | +0.00e+00 |
| 2 | inertia | 0.703360 | 0.703360 | -3.33e-16 |

Predictions: centrifugal stress s^2 w^2, centrifugal displacement s^3 w^2, angular-acceleration stress s^2 alpha, displacement s^3 alpha, gravity stress s, mass s^3, inertia s^5; frequencies 1/s are exact only without prestress (prestress raises them by a speed-dependent amount).


## Run notes (2026-11-05, written by hand after the pipeline; not regenerated by the scripts)

- Design points: s = 0.9010437566657981 (0.5 m, 24 Hz), 1.5761135310793883 (1 m, 24 Hz), 2.854447700523453 (2 m, 18 Hz). Speeds fc +- 0.545306 Hz (rotation = f_g/2, i.e. 30 f_g rpm); peak alpha = pi*0.867212 rad/s^2; phase-sampled transitions are every 3-tone hop (k0 -> k0+3 and back, tone spacing 0.181760 Hz) with a quintic over the full 224/190 s symbol.
- CAD mass and inertia are the archived CAD values (7.913509 kg, 0.3160713 kg m^2 at s=1) scaled by s^3 and s^5. The 1 m and 2 m values agree with the link-sizing references (30.983 / 3.0745 / 8.376 and 184.05 / 59.893 / 163.17) to 1e-4. At 0.5 m the reference values (5.7912 kg, 0.18750 kg m^2, 0.5108 N m) are not self-consistent with s = 0.9010438: CAD scaling gives 5.7890 kg (-0.04%), 0.18772 kg m^2 (+0.12%) and 0.51143 N m (+0.12%). The FEA uses the CAD-scaled values (slightly conservative torque).
- Environment: WSL Ubuntu 24.04; CalculiX 2.21 from the Ubuntu calculix-ccx 2.21-1 deb and Gmsh 4.12.1 from the Ubuntu python3-gmsh 4.12.1+ds1-1.1build2 deb, both extracted user-locally (~/local/calculix_root, ~/local/gmsh_deb_root; env file ~/local/gravcomm_fea_env_deb.sh), with Ubuntu system python3 3.12.3, numpy 1.26.4, scipy 1.11.4. This is the archived package set. The PyPI gmsh 4.12.1 wheel (venv ~/venvs/gravcomm-fea) does NOT reproduce the archived tungsten meshes (13 mm: 82630 instead of 82565 nodes); the Ubuntu build reproduces both archived tungsten reference meshes byte-for-byte (SHA-256 identical) and the 10 mm mesh node count (100535).
- Verification before new cases: the unmodified three_link_fea.py regenerated 25/26 archived decks byte-identically (h0.018 basis_alpha differs only in the last printed digit of some CLOAD values under the venv numpy). Re-solving archived h0.018 design1_high and basis_gx gave max displacement, raw peak stress and all frequencies identical to the archived summary.json to every printed digit (7.687027794809792 MPa, 10.02244682 um, 77.51643/94.90422/105.4042 Hz). In the new run all 16 scale-independent basis decks (steel and tungsten) matched the archived input SHA-256 except h0.018 basis_alpha; all 16 basis results equal the archived ones (peak stress identical; displacement within 1 ulp).
- Differences from the archived procedure: ccx single-threaded (OMP_NUM_THREADS=1, 8 parallel jobs) instead of 2 threads; no case reuse (all 52 steel+tungsten cases solved); transition model is the full-symbol 3-tone hop instead of the 0.6-symbol extreme-to-extreme hop; peak alpha now 2.7244 rad/s^2 instead of 6.1619; horizontal support drawing not regenerated; the steel s=1 10 mm third-mesh check is taken from the archived record (scale independent). Five regression tests (test_fea.py) pass.
- Full solver outputs (decks, .dat/.frd, fields.npz, geometry npz; 9.4 GB steel + 7.1 GB tungsten) are kept only on the WSL side in ~/gravcomm_fea_nb_2026_11_05/public_release_v0_5_0/fea/results/*_2026_11_05_nb; this Windows folder mirrors the archived release layout (top-level records plus per-case summary.json). Log: pipeline_2026_11_05_nb.log. Command: bash fea/run_nb_pipeline.sh (GRAVCOMM_FEA_TAG=2026_11_05_nb).
