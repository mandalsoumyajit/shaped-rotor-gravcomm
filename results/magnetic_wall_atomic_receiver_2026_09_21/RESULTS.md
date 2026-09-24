# Atomic magnetometer alternative

Aligned carrier centers remain 24, 24 and 18 Hz; the ideal capacity window remains +/-1 Hz. Source geometry and power ceilings are unchanged. The 0.3403 m radius PCB receiving loop is replaced by a compact 5 mm radius disk-average field sensor, 15 mm from the wall face. This is an explicit spatial surrogate for a vapor cell; wall thermal noise is recalculated for this smaller sensing region. It is not a model of a particular packaged instrument.

Two narrowband instrument scenarios are used: 35 fT/sqrt(Hz), motivated by a published SERF sensitivity of 33.5 fT/sqrt(Hz) at 10 Hz, and 5 pT/sqrt(Hz), from a NIST Mx instrument benchmark across 1–100 Hz. Constant ASD across the 17–25 Hz comparison bands is an assumption. The SERF case assumes a suitable compensated field and operating environment; it is not demonstrated unshielded performance. Any shielding or feedback that changes the wanted field requires its own transfer model. No gradiometric background suppression is credited. Neither scenario is extrapolated to sub-Hz in this study.

The external ELF backgrounds remain 0.38, 1, 3.14 and 10 pT/sqrt(Hz), added in PSD with instrument and wall thermal noise. Atomic sensing removes the induction conversion penalty but does not suppress actual environmental magnetic fluctuations.

| Distance (m) | Sensor | Background (pT/sqrt(Hz)) | BPSK boundary (S/m) | Capacity boundary (S/m) |
|---:|---|---:|---:|---:|
| 0.5 | optimistic compensated SERF | 0.38 | 6.867e+07 | 7.657e+07 |
| 0.5 | optimistic compensated SERF | 1 | 6.418e+07 | 7.179e+07 |
| 0.5 | optimistic compensated SERF | 3.14 | 5.724e+07 | 6.43e+07 |
| 0.5 | optimistic compensated SERF | 10 | 5.034e+07 | 5.68e+07 |
| 0.5 | Mx benchmark | 0.38 | 5.439e+07 | 6.121e+07 |
| 0.5 | Mx benchmark | 1 | 5.429e+07 | 6.11e+07 |
| 0.5 | Mx benchmark | 3.14 | 5.342e+07 | 6.015e+07 |
| 0.5 | Mx benchmark | 10 | 4.971e+07 | 5.61e+07 |
| 1 | optimistic compensated SERF | 0.38 | 1.435e+07 | 1.577e+07 |
| 1 | optimistic compensated SERF | 1 | 1.345e+07 | 1.482e+07 |
| 1 | optimistic compensated SERF | 3.14 | 1.22e+07 | 1.35e+07 |
| 1 | optimistic compensated SERF | 10 | 1.096e+07 | 1.218e+07 |
| 1 | Mx benchmark | 0.38 | 1.169e+07 | 1.296e+07 |
| 1 | Mx benchmark | 1 | 1.167e+07 | 1.294e+07 |
| 1 | Mx benchmark | 3.14 | 1.151e+07 | 1.277e+07 |
| 1 | Mx benchmark | 10 | 1.084e+07 | 1.205e+07 |
| 2 | optimistic compensated SERF | 0.38 | 4.698e+06 | 5.212e+06 |
| 2 | optimistic compensated SERF | 1 | 4.379e+06 | 4.872e+06 |
| 2 | optimistic compensated SERF | 3.14 | 3.975e+06 | 4.439e+06 |
| 2 | optimistic compensated SERF | 10 | 3.579e+06 | 4.012e+06 |
| 2 | Mx benchmark | 0.38 | 3.813e+06 | 4.263e+06 |
| 2 | Mx benchmark | 1 | 3.807e+06 | 4.257e+06 |
| 2 | Mx benchmark | 3.14 | 3.757e+06 | 4.203e+06 |
| 2 | Mx benchmark | 10 | 3.542e+06 | 3.971e+06 |

The receiver change includes both noise and sensing geometry. Boundary shifts therefore cannot be attributed to intrinsic sensitivity alone. The same unresolved average-electrical-power accounting and finite-packet limitations as the aligned-carrier study apply.

Sources:
- https://www.nature.com/articles/s41378-026-01226-z
- https://tf.nist.gov/ofm/smallclock/Mx_Magnetometer.html

Reproduce with scripts/magnetic_wall_atomic_receiver.py.
