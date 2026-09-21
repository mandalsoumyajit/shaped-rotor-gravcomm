# Structural reassessment of the qualified 0.5, 1 and 2 m sources

Scope follows the user clarification: 7.03 kg at 0.5 m, 37.69 kg at 1 m and the smallest passing 227.32 kg source at 2 m. The larger 2 m alternatives are not the selected study. Communications timing is unchanged. All quantities below are calculations for the bare rotor with fixed bore and perfectly bonded inserts unless explicitly identified as a separate assembly screen.

## Fixed-bore linear FEA

Twenty-six static/prestressed-modal cases use two archived mesh resolutions, three actual operating speeds per design, and four independently solved gravity/angular-acceleration bases per mesh. Seven input-identical completed cases were reused after SHA-256 checks. All solver jobs must finish successfully; static/modal outputs are separated. Five existing FEA regression tests pass. Positive Jacobians, exact consistent-load torque, force/torque equilibrium and similarity to the archived desktop solution are checked.

| Link (m) | Bare mass (kg) | Upper-speed centrifugal Al peak (MPa) | Horizontal combined Al bound (MPa) | Centrifugal displacement (um) | Horizontal combined displacement bound (um) |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 7.03 | 8.91 | 9.72 | 10.03 | 15.26 |
| 1 | 37.69 | 27.29 | 28.98 | 53.80 | 77.77 |
| 2 | 227.32 | 51.46 | 55.50 | 184.69 | 314.95 |

The combined columns are triangle-inequality bounds assembled from the linear static solutions, including maximum centrifugal load, either sign of maximum angular acceleration and any in-plane gravity direction. They are deliberately conservative; speed and angular-acceleration maxima need not occur simultaneously. Separate phase-sampled extreme transitions use the actual quintic trajectory. Neither calculation is a nonlinear contact simulation or a complete time-domain rotor-bearing response. Gravity/acceleration bases use unprestressed elastic stiffness, consistent with the original linear statics model.

Horizontal shaft orientation is the preliminary choice. All-attitude bounds (including axial gravity) are substantially larger and retained in compiled.json. Gravity direction rotates at the shaft frequency in rotor coordinates for a horizontal shaft. Loads at the 24/24/18 Hz gravitational carriers are not automatically mechanical excitation at those frequencies.

## Mesh convergence and strength interpretation

| Link (m) | Change in centrifugal raw Al peak (%) | Change in max displacement (%) | Change in maximum scaled-window combined average (%) |
|---:|---:|---:|---:|
| 0.5 | 15.85 | 0.089 | 0.211 |
| 1 | 15.85 | 0.089 | 0.221 |
| 2 | 15.85 | 0.089 | 0.214 |

Regional windows have reference radii 5/10 mm multiplied by the geometric scale; they span the carrier thickness. The reported average is the maximum of the 10 mm-reference windows. It is a comparison metric, not a replacement for local yield or fatigue analysis. Fine-mesh mass is about 0.11% below CAD because of the faceted approximation; density is not renormalized. See mesh_h*.json.

Third-mesh check: a newly generated 10 mm reference mesh has 132043 nodes and 83555 elements. Its 1800 rpm aluminum peak is 58.492 MPa, a 0.628% increase from the 13 mm mesh. Displacement is 68.127 um and first frequency 79.825 Hz. This supports convergence of the centrifugal reference peak; it does not establish convergence of every combined-load or contact stress.

Preliminary material choices are 6061-T651 plate (240 MPa yield screen), specifically certified normalized 1045 inserts (300 MPa minimum procurement requirement), and 4140 Q&T shafts/flanges (650 MPa minimum study specification). A factor of two against yield is the initial static screening criterion. Material sources, product-form caveats and actual assembly choices are in ASSEMBLY_SPECIFICATION.md and materials.json. Raw mesh peaks and idealized attachments do not establish hardware safety factors. Fatigue, preload, real contact and local features remain unresolved.

All three idealized horizontal-shaft combined aluminum bounds are below the 120 MPa screening limit (240 MPa divided by two); the largest is about 55.50 MPa. The model therefore indicates useful static margin in the carrier, while insert retention and assembly dynamics still require explicit validation.

## Modes, coupling and response

| Link (m) | First three center-speed modes (Hz) | Dominant angular-acceleration mode | Polar-inertia fraction (%) | Shaft/carrier series torsion screen (Hz) |
|---:|:---|---:|---:|---:|
| 0.5 | 77.33, 94.69, 105.27 | 3 | 99.76 | 44.11 |
| 1 | 45.39, 55.04, 60.18 | 3 | 99.76 | 32.81 |
| 2 | 25.76, 30.88, 33.09 | 3 | 99.76 | 21.23 |

Mode matching uses the consistent mass matrix and modal assurance criterion (MAC), with node correspondence within each mesh. Mass normalization is independently recovered from degree-four exact element mass integrals. Component energies and mean insert motions identify axial/torsional motion. Angular and translational effective mass fractions distinguish modulation coupling from gravity coupling; low frequency alone does not establish strong excitation.

modal_response_h*.json contains damping sweeps at 0.5%, 2% and 5% for isolated quintic-transition scalar modal responses and rotation-frequency harmonic transfer factors. The transition simulation starts at rest and omits repeated-message resonance, motor control, gyroscopic terms and contact. It diagnoses time-scale separation, not arbitrary-message dynamic qualification. Shaft/carrier series values are a simple compliance screen using the chosen stub-shaft dimensions and dominant carrier mode; motor/coupling inertia and bearing/structure dynamics are absent. These values demonstrate why fixed-bore carrier frequencies must not be called assembly critical speeds.

## Assembly and next checks

ASSEMBLY_SPECIFICATION.md selects a horizontal straddled shaft, locating angular-contact pair plus floating cylindrical-roller support, piloted hub flanges and positive insert shoulder/collar retention. Plain-shaft and pocket-pressure estimates are saved separately in assembly_screen.json. These are preliminary calculations, not additional 3D assembly FEA. Before physical implementation, add flange bolts, contact/clearance, retention grooves/collars, preload, actual bearing and motor/coupling stiffness, and fatigue. Added hardware mass and off-axis retention features require a revised field/inertia audit before reusing the BER claims for that assembly.

The tungsten-heavy-alloy option in TUNGSTEN_OPTION.md retains insert mass and diameter while shortening its length; it is a promising windage/overhang reduction study. It has not been substituted into the qualified steel models.

## Reproduction

Run fea/three_link_fea.py --prepare --workers 4 in the documented WSL environment, then compile_three_link_fea.py --mesh .018 and --mesh .013; three_link_modal_response.py for each mesh; assembly_screen.py; report_three_link_fea.py. Completed input-identical cases are reusable. Inputs, complete solver outputs, fields, equilibrium checks and hashes are retained. No BER samples or source configuration was changed.
