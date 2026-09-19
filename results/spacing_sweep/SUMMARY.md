# Fixed-duration CF-FSK tone-spacing sweep

Desktop shaped rotor at 0.5 m; 1.75 s symbols; 0.7 s quintic transition; seven tones centered at 50.3 Hz. Eight payload bytes plus six known trailer symbols. Nominal payload rate 1.524 bit/s; packet rate including trailer 1.219 bit/s. Acquisition overhead is excluded. No added ECC.

Standard Viterbi decoding retains continuous phase and the full finite whitening-filter history, and enforces the existing byte-to-three-tone mapping. Phase and timing are acquired. Noise is generated from the physical receiver PSD. Primary runs use 8 complex samples/s and 5.25 s whitening memory. The source trajectory is ideal; actuator requirements are evaluated separately.

Each point stops after at least 100 bit errors in at least 30 independent packets, or 10,000 packets. Intervals are packet-aware, anytime-valid 95% intervals for each point individually. They are not simultaneous confidence intervals across the sweep. Sparse-error points are reported as upper bounds. Raw observed fractions remain in the CSV for auditing.

| Spacing (Hz) | Bits | Errors / errored packets | BER estimate or 95% upper bound | Worst RMS (N m) | Random RMS (N m) | Peak (N m) | Peak power (W) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.04762 | 16000 | 1026 / 130 | 0.0641 | 0.295 | 0.149 | 0.810 | 128.0 |
| 0.07143 | 48000 | 193 / 38 | 0.00402 | 0.438 | 0.216 | 1.190 | 188.0 |
| 0.09524 | 560000 | 152 / 30 | 0.000271 | 0.583 | 0.285 | 1.570 | 248.1 |
| 0.14286 | 640000 | 36 / 9 | < 0.000596 | 0.872 | 0.424 | 2.330 | 368.2 |
| 0.19048 | 240000 | 119 / 33 | 0.000496 | 1.162 | 0.564 | 3.090 | 488.3 |
| 0.28571 | 48000 | 160 / 43 | 0.00333 | 1.743 | 0.844 | 4.610 | 728.5 |
| 0.47619 | 16000 | 288 / 64 | 0.018 | 2.904 | 1.404 | 7.649 | 1209.2 |
| 0.71429 | 16000 | 1414 / 187 | 0.0884 | 4.356 | 2.106 | 11.449 | 1810.6 |
| 0.95238 | 16000 | 2776 / 233 | 0.173 | 5.807 | 2.807 | 15.248 | 2413.0 |

Torque includes 0.05 N m drag at the carrier. Worst-message RMS is the maximum cycle mean over all valid byte sequences; random-byte RMS assumes independent uniform bytes. The continuous rating is 0.5 N m, the peak rating 2 N m, and the bidirectional power rating 400 W. These duty metrics describe sustained payload transmission; trailer and acquisition duty depend on packet scheduling.

## Numerical checks

| Run | Packets | Bit errors | Errored packets |
|---|---:|---:|---:|
| q0p1_fs16_m3 | 1000 | 18 | 4 |
| q0p1_fs8_m4 | 1000 | 4 | 2 |
| q0p5_fs16_m3 | 250 | 364 | 72 |
| q0p5_fs8_m4 | 250 | 259 | 53 |
| q1_fs16_m3 | 250 | 2914 | 239 |
| q1_fs8_m4 | 250 | 2888 | 237 |

The wider-spacing deterioration persists with doubled sampling and longer whitening memory. At q=0.1, the first 1000 primary packets contain 15 bit errors, compared with 18 at doubled sampling and 4 with longer memory. These low counts leave uncertainty in the detailed low-BER curve. The checks support the broad trend; a precise optimum requires more sampling.

## Interpretation

The receiver exhibits a useful spacing window: narrow spacings produce poorly separated sequences, while wider spacings encounter increasing input-referred readout noise. The lowest observed error fraction occurs at 0.14286 Hz, with only 36 errors in nine packets; its supported statement is BER < 5.96e-4 at 95% confidence. The sweep does not resolve a precise optimal spacing.

At 0.09524 Hz the receiver meets the 1e-3 target with a per-point 95% upper bound of 8.98e-4. Peak torque (1.570 N m) and power (248 W) fit the drive ratings. Worst-message sustained RMS torque is 0.583 N m, exceeding the 0.5 N m continuous rating; independent random-byte RMS is 0.285 N m. At 0.14286 Hz, peak torque is 2.330 N m and worst-message RMS is 0.872 N m. Consequently, continuous transmission guaranteed for every valid message encounters motor torque constraints before reaching the apparent receiver optimum at this symbol duration. Traffic with specified statistical or burst duty can have different thermal requirements.

All 75 active regression tests pass. The manuscript and public release have not been changed by this sweep.
