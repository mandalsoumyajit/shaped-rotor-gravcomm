# Preliminary materials and assembly specification

This is an engineering study concept chosen after the user requested reasonable material and assembly assumptions. It is not a drawing release or a validated contact/rotordynamic assembly. The 3D FEA retains the original bare rotor and idealized boundary conditions; these assembly estimates are separate.

## Materials

- Carrier: unwelded 6061-T651 machined plate; retain E=70 GPa, nu=.33 and density 2700 kg/m^3. Use 240 MPa room-temperature yield for screening, rounded below the 35 ksi minimum listed for the applicable thickness range in [Kaiser sheet/plate data](https://online.kaiseraluminum.com/depot/PublicProductInformation/Document/1010/Kaiser_Aluminum_Sheet_Plate.pdf), p.20. Stock certification and product form must match.
- Inserts: normalized 1045 steel, E=200 GPa, nu=.30, density 7850 kg/m^3, with a specifically ordered and certified minimum yield of 300 MPa. [Atlas 1045](https://atlassteels.com.au/documents/Atlas1045.pdf) distinguishes typical values from guaranteed special-order properties; grade alone is insufficient.
- Shafts and hub flanges: 4140 quenched and tempered, with a conservative specified minimum yield of 650 MPa. [Atlas 4140](https://atlassteels.com.au/documents/Atlas4140.pdf) lists higher minimum proof stress for its relevant supplied condition/diameter. Retain the modeled steel elastic constants for continuity.
- Use a factor of two against the stated room-temperature yield values as an initial static design screen. This is a chosen study criterion, not a code-prescribed safety factor, fatigue allowable or certification.

## Arrangement

Choose a horizontal shaft and vertical rotor plane initially: gravity then acts in the plane and avoids the much larger static axial sag of the thin carrier. Retain both attitudes in the FEA bounds. Straddle the rotor with two bearing supports; use a matched back-to-back angular-contact locating pair on the drive side and an axially floating cylindrical-roller support opposite. This follows the general [SKF locating/non-locating arrangement](https://www.skf.com/binaries/pub12/Images/0901d196802809de-Rolling-bearings---17000_1-EN_tcm_12-121486.pdf). Bearing preload, exact part numbers, stiffness and foundation design are not yet validated.

Use two flanged stub shafts piloted into the existing scaled bore and clamped to opposite hub faces. This allows standard 20/40/80 mm bearing seats without enlarging the BER-qualified rotor bore. Six fitted flange bolts carry torque by an explicit positive load path; no friction-only torque credit is taken. Bolt holes and clamping pressure must be added to the assembly contact model. The plain-shaft stress calculation below omits stress concentrations at shoulders/fillets and bolt holes.

For each insert, use a machined positive shoulder on one side of the carrier and a removable split retaining collar engaged in a groove on the other. The cylindrical pocket carries radial centrifugal load by contact; shoulder/collar carries axial retention. Do not credit adhesive bonding or uncontrolled interference. Define running clearance, collar/groove dimensions, fastener preload, temperature and tolerances in a dedicated contact submodel before calling this attachment validated. Added off-axis metal/removal for grooves changes source mass and inertia and must be included in a revised finite-volume field calculation.

| Link (m) | Bearing-seat diameter (mm) | Bearing span (mm) | Hub flange thickness (mm) | Six flange bolts | Nominal shaft VM (MPa) | Insert radial force (kN each) | Nominal pocket pressure (MPa) |
|---:|---:|---:|---:|:---|---:|---:|---:|
| 0.5 | 20 | 120 | 5 | M4 | 3.32 | 2.41 | 1.86 |
| 1 | 40 | 200 | 8 | M6 | 4.81 | 22.61 | 5.70 |
| 2 | 80 | 360 | 12 | M10 | 10.02 | 141.26 | 10.76 |

The support calculation uses a central rotor weight on a simply supported solid circular shaft, M=W*L/4, nominal bending stress 32M/(pi*d^3), torsion 16T/(pi*d^3), and von Mises combination. Stiffness uses 48EI/L^3 in series with two radial bearing springs. The 1e7–1e9 N/m per-bearing sweep is a declared sensitivity range, not catalog stiffness. It gives a translational screen only, not a gyroscopic critical-speed prediction. Shaft torsional stiffness is GJ/L_drive; combining it with carrier torsion and motor/coupling dynamics is required before using fixed-bore modes as assembly modes.

The central shaft/flange mass inventory is illustrative and excludes bolts, collars, motor coupling and bearings. Its axisymmetric idealization contributes no rotating second harmonic, but its inertia still adds drive demand. The qualified masses refer to the bare source rotor. Real asymmetric features and the off-axis retention hardware require a new field/inertia audit.
