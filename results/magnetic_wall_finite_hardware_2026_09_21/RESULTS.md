# Finite-size magnetic link through a conducting wall

2026-09-21. Exploratory comparison; manuscript unchanged.

## Geometry and resource accounting

The barrier occupies the space remaining after the hardware and clearances are included. For gravitational axis-to-sensor distance d, rotor outer radius R, receiver half-depth 0.0127 m, and 0.01 m clearance at each face, wall thickness is t = d - R - 0.0127 - 0.02 m. Thus the nominal link distance cannot be treated as conducting thickness. The gravitational rotor orientation is the paper's configuration, with its radial extent toward the barrier. A full support/housing envelope could reduce available wall thickness further.

The magnetic alternative uses a circular coil parallel to the wall, coaxial with a circular receiving loop on the opposite side. Its outer radius equals the corresponding rotor radius. Copper mass is provisionally 10% of rotor mass; a square winding bundle with 60% copper packing determines its axial depth and mean radius. Each housing has 1 cm clearance from the wall. The receiving PCB is represented by an equivalent 0.3403 m radius loop with 5 mm half-depth. Its area is much larger than that of the gravitational receiver. These are explicit comparison assumptions, not an optimized magnetic design or equal receiver-volume allocation.

| Gravity distance (m) | Rotor diameter (m) | Wall thickness (m) | Magnetic coil mean radius (m) | Magnetic center separation (m) | Numerical power ceiling (W) |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.481 | 0.227 | 0.236 | 0.257 | 120.49 |
| 1.0 | 0.841 | 0.547 | 0.412 | 0.580 | 1980.24 |
| 2.0 | 1.531 | 1.202 | 0.751 | 1.242 | 29680.14 |

The ceilings equal the reported gravitational peak mechanical powers numerically. The magnetic calculation permits that much average dissipative source power. The gravitational average electrical demand is unvalidated, so these results do not demonstrate equal average electrical consumption. They deliberately give the magnetic alternative a generous power allowance.

## Model

The air/conductor/air slab is nonmagnetic, homogeneous, and laterally infinite. A magnetoquasistatic Hankel integral resolves finite source radius and finite receiver area, transmission through the slab, source-induced eddy-current dissipation, and receiver-averaged wall thermal noise at 300 K. With transverse wavenumber k, q = sqrt(k^2 + i omega mu0 sigma), r = (k-q)/(k+q), and slab transmission T = (1-r^2) exp(-qt)/(1-r^2 exp(-2qt)). The finite-loop Bessel factors are retained in the integral. This avoids a point-dipole or plane-wave skin-depth approximation at distances comparable to coil size.

The winding cross-section determines physical clearance and copper resistance; its electromagnetic response is approximated by a loop at the mean radius and axial plane. A resolved winding distribution remains a refinement. DC copper resistance is used. Driver voltage/current limits, skin/proximity losses, cooling, electronics power, and finite wall edges are omitted. Permeability and magnetic hysteresis require a separate model.

The ELF manuscript's 100-turn PCB example supplies summed turn area 36.388 m^2 and input-referred EMF noise 1.626 nV/sqrt(Hz) at 1 kHz. The present extrapolation uses a 20 Hz voltage-noise corner and divides by 2 pi f times summed area. This is an assumed low-frequency extension, not a measured sub-Hz noise curve. Ambient noise is initially 30 fT/sqrt(Hz), with 1 pT, 10 pT, and low-frequency ambient-knee sensitivity cases. Wall thermal noise is added separately. The scenario labelled `optimistic flat 1 fT total` in the raw output means a flat 1 fT input floor plus wall thermal noise; it is not a flat total noise spectrum.

## Results and interpretation

Two different boundaries are reported. The coherent BPSK benchmark has 224 payload bits plus 32 overhead bits per 224 s packet, raised-cosine roll-off 0.25, optimized carrier within 1–999 Hz, ideal pre-equalization/whitening, and analytical uncoded BER 1e-3. It is a restricted waveform benchmark, with no simulated finite-packet qualification. The water-filled Gaussian-channel calculation allows arbitrary spectral allocation across 0.0001–1000 Hz. Its required information rate is 1-h2(0.001) bit/s for balanced independent payload bits. This is an optimistic asymptotic feasibility bound; it does not impose the five-minute packet deadline or demonstrate a realizable code.

| Gravity distance (m) | BPSK power-ceiling crossing, sigma (S/m) | Ideal capacity-bound crossing, sigma (S/m) |
|---:|---:|---:|
| 0.5 | 7.91e8 | 8.06e9 |
| 1.0 | 1.91e8 | 2.32e9 |
| 2.0 | 4.47e7 | 6.19e8 |

At the illustrative conductivity 5.8e7 S/m, the restricted BPSK benchmark exceeds the ceiling for the 2 m gravitational geometry, whose wall is 1.20 m thick. However, the ideal capacity calculation remains feasible by a large margin. Thus BPSK failure does not establish that a magnetic link cannot deliver the target. No gravity-only operating region at the three reference geometries is established for this conductivity by the present calculation. Claims about all electromagnetic alternatives would additionally require broader source, receiver, material, and protocol optimization.

The very small calculated excitation powers for weakly attenuating barriers exclude electronics and other fixed power costs; they should not be quoted as complete transmitter consumption. Noise extrapolation, finite-block coding, and electrical hardware constraints need evaluation before a practical crossover can be claimed.

Ambient-level variations have little effect on these high-conductivity crossings because the assumed low-frequency instrumental noise dominates. Increasing quadrature resolution and comparing direct evaluation with the interpolated channel gives maximum relevant gain differences of 0.037, 0.046, and 0.056 dB. Changing each face clearance from 1 cm to 0.5 or 2 cm also changes wall thickness at fixed gravitational distance: ideal crossings span 7.55–9.29e9, 2.26–2.45e9, and 6.12–6.33e8 S/m, respectively. This is a combined clearance/thickness sensitivity.

`wall_boundary.png` and `.pdf` show the ideal bound across distance and conductivity for each fixed source envelope and power ceiling. Only the three reference gravitational distances have corresponding qualified gravitational designs; the extended map does not qualify gravity at other distances.

## Sources and reproducibility

- Finite-coil slab formulation: Vallecchi et al., *Coupling between coils in the presence of conducting medium* (2019), https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/iet-map.2018.5292 .
- Dissipation-based magnetic thermal noise: Lee and Romalis, *Calculation of Magnetic Field Noise from High-Permeability Magnetic Shields and Conducting Objects with Simple Geometry*, https://arxiv.org/abs/0709.2543 .
- Local noise evidence: `ELF_Paper/journal/ieee-access/IEEE-Access-ELF.tex` and `ELF_Paper/Coil_Design/pcb_loop_optimizer.py` (read only).
- Model: `AIP_Advances_submission/src/gravcomm/magnetic_wall.py`; runner: `AIP_Advances_submission/scripts/magnetic_wall_study.py`. Reusable noise/capacity functions are in `magnetic_link.py`; design constants are imported from `magnetic_conductivity_study.py`.
- Numerical outputs: `summary.json`, `maps.json`; tests: `test_magnetic_wall.py` and `test_magnetic_link.py`.

Earlier spherical-cavity and wall calculations that used the full nominal distance as available barrier thickness are superseded by this finite-hardware study.

## Background-noise qualification after review

The 30 fT/sqrt(Hz) case is a quiet illustrative scenario, and cannot serve as the realistic environmental baseline. The cited observations are frequency- and site-specific; the 1 kHz examples do not establish the sub-Hz spectrum used by optimized high-conductivity links. Flat 1 and 10 pT/sqrt(Hz) sensitivity cases also do not adequately test strong low-frequency environmental noise. The assumed PCB instrumental ASD itself rises to approximately 32 pT/sqrt(Hz) at 1 Hz and 1 nT/sqrt(Hz) at 0.1 Hz. This explains the weak sensitivity to those flat ambient sweeps; it does not establish that realistic ambient noise is unimportant. All reported boundaries remain conditional pending a justified frequency-dependent environmental model.

Neither link should receive an unmodeled gradiometric cancellation credit. The gravitational study includes a specified mechanical isolation chain and residual environmental noise, but no gradiometric subtraction. Magnetic noise must be specified at the receiver site: local interference on the receiving side does not automatically receive the transmitter-to-receiver wall attenuation. Narrow lines and nonstationary interference also require treatment beyond a broadband Gaussian floor.

A gradiometric extension must calculate both the differential source response h1 - w h2 and output noise PSD S11 + |w|^2 S22 - 2 Re(w* S12), where S12 = E[n1 n2*]. Include baseline geometry, independent sensor noise, environmental spatial coherence, and mismatch. Apply this framework to both magnetic and gravitational receivers; a rejection factor without the associated signal response and covariance is insufficient.

Next study step: use documented frequency-dependent receiver-site background scenarios across the optimized band; report occupied frequencies and instrumental/environmental/wall-noise contributions; recompute waveform and capacity boundaries. Present single-sensor results first, followed by any explicitly modeled gradiometric alternatives. No unconditional gravity-only conclusion is supported by the current results.
