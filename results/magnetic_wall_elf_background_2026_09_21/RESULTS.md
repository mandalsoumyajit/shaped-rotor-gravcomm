# ELF-background magnetic-wall comparison

The representative background is 1 pT/sqrt(Hz), bracketed by the ELF paper's mine-surface examples of 0.38 and 3.14 pT/sqrt(Hz). The 10 pT/sqrt(Hz) tunnel example is retained as a stress case. The full 0.01–3 pT/sqrt(Hz) ELF sweep is also sampled. These are constant broadband ASD assumptions; the cited 1 kHz observations are not measurements across the entire optimized link band. No gradiometric cancellation is credited to either link. Ambient, instrumental, and wall-thermal PSDs add.

Geometry, finite-coil response, noise extrapolation, and power accounting follow ../magnetic_wall_finite_hardware_2026_09_21/RESULTS.md. In particular, the numerical magnetic average-power ceilings equal gravitational peak mechanical powers; equal average electrical consumption is not established. The Gaussian capacity boundary is asymptotic and does not demonstrate a finite-packet implementation.

| Distance (m) | Background (pT/sqrt(Hz)) | BPSK boundary (S/m) | Ideal capacity boundary (S/m) |
|---:|---:|---:|---:|
| 0.5 | 0.38 | 7.913e+08 | 8.055e+09 |
| 0.5 | 1 | 7.912e+08 | 8.055e+09 |
| 0.5 | 3.14 | 7.899e+08 | 8.055e+09 |
| 0.5 | 10 | 7.783e+08 | 8.055e+09 |
| 1 | 0.38 | 1.907e+08 | 2.316e+09 |
| 1 | 1 | 1.907e+08 | 2.316e+09 |
| 1 | 3.14 | 1.904e+08 | 2.316e+09 |
| 1 | 10 | 1.879e+08 | 2.316e+09 |
| 2 | 0.38 | 4.466e+07 | 6.187e+08 |
| 2 | 1 | 4.465e+07 | 6.187e+08 |
| 2 | 3.14 | 4.459e+07 | 6.187e+08 |
| 2 | 10 | 4.404e+07 | 6.187e+08 |

The boundaries change little across this environmental bracket because the assumed inductive receiver becomes noisy at low frequencies. This is a property of the specified receiver and broadband assumptions. It does not establish environmental insignificance for a different magnetometer or for an actual low-frequency site spectrum.

| Frequency (Hz) | Instrument ASD (pT/sqrt(Hz)) |
|---:|---:|
| 0.1 | 998.3 |
| 1 | 32.27 |
| 10 | 1.22 |
| 100 | 0.07714 |
| 1000 | 0.007112 |

At conductivity 5.8e7 S/m the 2 m geometry still exceeds the BPSK benchmark power budget while passing the optimistic capacity bound. This does not establish that magnetic communication is impossible. Only the three reference gravitational links are qualified.

Source: ELF_Paper/journal/ieee-access/IEEE-Access-ELF.tex, background discussion and PCB receiver example. Reproduce with scripts/magnetic_wall_elf_background.py. The earlier quiet-case map remains labelled as the earlier scenario; this directory contains the revised three-design comparison.

## Skin-depth audit

The user's 24 Hz check is correct. With sigma = 5.8e7 S/m, delta = 1/sqrt(pi f mu0 sigma) = 0.01349 m. Across the 1.20166 m wall, exp(-t/delta) = 2.06e-39 in amplitude. Direct finite-loop evaluation gives a field ratio relative to the same coils in free space of 2.84e-40. At 18 Hz the corresponding finite-loop ratio is 5.56e-35. The model therefore also predicts overwhelming attenuation at the gravitational carrier frequencies.

The favorable capacity results used independent magnetic frequency optimization. At sigma = 5.8e7 S/m and background 1 pT/sqrt(Hz), occupied positive-frequency bands are approximately 0.223–1.532 Hz, 0.0304–0.568 Hz, and 0.00241–0.241 Hz for the three geometries. For the 2 m geometry, the central 90% of information lies between 0.0146 and 0.192 Hz. At 0.1 Hz the skin depth is 0.209 m and bulk amplitude attenuation is 0.00318. Thus the favorable values do not describe an 18–24 Hz magnetic link.

An independent scalar water-level root solve reproduces the reported asymptotic capacity power. This arithmetic check does not qualify a realizable low-frequency modem: finite packet duration, acquisition, transients and receiver noise extrapolation remain unresolved. The representative ELF background is being extrapolated far below the frequencies of the cited measurements. Future reporting must show the selected magnetic spectrum next to every feasibility claim. Raw audit data are in skin_depth_audit.json; reproduce with scripts/audit_magnetic_wall_skin_depth.py.

## Primary comparison updated

Per the user's instruction, the primary magnetic comparison now uses the gravitational carrier centers (24, 24 and 18 Hz) and a 2 Hz allowed spectral window for the capacity bound. See ../magnetic_wall_aligned_carriers_2026_09_21/RESULTS.md. The independently optimized sub-Hz results in this directory are supplementary and do not represent the agreed comparison.
