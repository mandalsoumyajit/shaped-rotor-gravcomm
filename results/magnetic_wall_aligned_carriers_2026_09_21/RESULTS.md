# Aligned-carrier magnetic-wall comparison

This is the primary comparison requested by the user. Magnetic carrier centers equal the qualified gravitational centers: 24 Hz at 0.5 and 1 m, and 18 Hz at 2 m. Earlier independently optimized sub-Hz magnetic results are supplementary bounds and do not define this comparison.

The finite-hardware wall geometry is retained: thicknesses 0.227, 0.547 and 1.202 m. Source outer radii match the gravitational rotor radii; the magnetic source is parallel to the wall with its winding depth included. Noise comprises instrumental, wall thermal and ELF-paper background PSDs. The representative background is 1 pT/sqrt(Hz), bracketed by 0.38 and 3.14 pT/sqrt(Hz), with 10 pT/sqrt(Hz) as stress. These are broadband level assumptions. Neither link receives unmodeled gradiometric cancellation credit.

The analytic coherent BPSK benchmark carries 224 payload plus 32 overhead bits in 224 s, roll-off 0.25, and target uncoded BER 1e-3. Its support extends +/-0.7143 Hz about the fixed center. It assumes ideal pre-equalization and whitening. The optimistic Gaussian capacity bound is restricted to a 2 Hz interval about that same center, matching the stated nominal receive bandwidth. This is a hard allowed spectral window for the bound; it does not reproduce the gravitational receiver filter transition band. Water filling may favor one side within this window. Its information-rate target is 1-h2(0.001) bit/s; finite-length coding and acquisition are unqualified.

Power ceilings remain numerically equal to gravitational peak mechanical demand (120.49, 1980.24, 29680.14 W), while the magnetic model allows this as average source dissipation. Equal average electrical power remains unresolved. Model scope and omitted driver/thermal constraints follow the finite-hardware study.

| Distance (m) | Center (Hz) | Background (pT/sqrt(Hz)) | BPSK crossing (S/m) | Band-limited capacity crossing (S/m) |
|---:|---:|---:|---:|---:|
| 0.5 | 24 | 0.38 | 6.637e+07 | 7.404e+07 |
| 0.5 | 24 | 1 | 6.188e+07 | 6.93e+07 |
| 0.5 | 24 | 3.14 | 5.499e+07 | 6.193e+07 |
| 0.5 | 24 | 10 | 4.813e+07 | 5.453e+07 |
| 1 | 24 | 0.38 | 1.548e+07 | 1.714e+07 |
| 1 | 24 | 1 | 1.453e+07 | 1.614e+07 |
| 1 | 24 | 3.14 | 1.306e+07 | 1.458e+07 |
| 1 | 24 | 10 | 1.159e+07 | 1.3e+07 |
| 2 | 18 | 0.38 | 4.607e+06 | 5.122e+06 |
| 2 | 18 | 1 | 4.4e+06 | 4.908e+06 |
| 2 | 18 | 3.14 | 4.011e+06 | 4.493e+06 |
| 2 | 18 | 10 | 3.597e+06 | 4.046e+06 |

For the representative 1 pT/sqrt(Hz) background, at copper conductivity 5.8e7 S/m:

| Distance (m) | BPSK required source W | Ideal minimum source W |
|---:|---:|---:|
| 0.5 | 36.3 | 4.21 |
| 1 | 2.08e+21 | 1.18e+20 |
| 2 | 8.58e+55 | 3.07e+53 |

Large extrapolated powers diagnose infeasibility within the linear model; they are not physical drive designs. Failure of the capacity bound applies to the specified source, sensor, noise and allowed band. It does not exclude other electromagnetic frequencies, sensor types, geometries or communication pathways.

Maximum absolute change from doubling spectral integration resolution: 2.826e-07 dB. Earlier finite-wall direct quadrature checks and the fixed-frequency skin-depth audit also apply.

Reproduce: scripts/magnetic_wall_aligned_carriers.py. Raw data: summary.json.
