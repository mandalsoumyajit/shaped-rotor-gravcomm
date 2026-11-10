# Magnetic-wall comparison for hop-limited CF-FSK designs (2026-11-05)

Released model code unchanged; only rotor radius (0.25 m x scale), rotor mass and power ceiling (new peak mechanical power) change. Coil copper mass 10% of rotor mass, coil outer radius = rotor radius, t = d - R - 0.0327 m, receivers, noise levels, BPSK benchmark and the fc +/- 1 Hz capacity window are as in the 2026-09-21 studies. The old columns are recomputed here and match the archived values.

| d (m) | R old/new (m) | t old/new (m) | ceiling old/new (W) | fc (Hz) | PCB 1 pT old/new (1e6 S/m) | SERF 1 pT old/new (1e6 S/m) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.2403/0.2253 | 0.2270/0.2420 | 120.49/38.95 | 24 | 69.3/56.2 | 71.8/61.6 |
| 1 | 0.4206/0.3940 | 0.5467/0.5733 | 1980.24/638.77 | 24 | 16.1/13.9 | 14.8/13.6 |
| 2 | 0.7656/0.7136 | 1.2017/1.2537 | 29680.14/9367.95 | 18 | 4.91/4.37 | 4.87/4.43 |

## New designs: aligned-carrier PCB receiver

| d (m) | fc | background (pT) | BPSK crossing (S/m) | capacity crossing (S/m) | copper BPSK W | copper capacity W |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 24 | 0.38 | 5.379e+07 | 6.016e+07 | 166 | 19.1 |
| 0.5 | 24 | 1 | 5.002e+07 | 5.618e+07 | 637 | 71.6 |
| 0.5 | 24 | 3.14 | 4.429e+07 | 5.004e+07 | 5.52e+03 | 615 |
| 0.5 | 24 | 10 | 3.857e+07 | 4.385e+07 | 5.52e+04 | 6.14e+03 |
| 1 | 24 | 0.38 | 1.337e+07 | 1.483e+07 | 2.45e+22 | 1.35e+21 |
| 1 | 24 | 1 | 1.252e+07 | 1.395e+07 | 9.46e+22 | 5.04e+21 |
| 1 | 24 | 3.14 | 1.123e+07 | 1.256e+07 | 8.2e+23 | 4.32e+22 |
| 1 | 24 | 10 | 9.919e+06 | 1.116e+07 | 8.21e+24 | 4.32e+23 |
| 2 | 18 | 0.38 | 4.099e+06 | 4.564e+06 | 6.86e+57 | 2.02e+55 |
| 2 | 18 | 1 | 3.91e+06 | 4.37e+06 | 1.99e+58 | 5.44e+55 |
| 2 | 18 | 3.14 | 3.554e+06 | 3.99e+06 | 1.54e+59 | 4.08e+56 |
| 2 | 18 | 10 | 3.176e+06 | 3.58e+06 | 1.52e+60 | 4.01e+57 |

## New designs: atomic receiver

| d (m) | sensor | background (pT) | BPSK (S/m) | capacity (S/m) | copper capacity W |
|---:|---|---:|---:|---:|---:|
| 0.5 | optimistic compensated SERF | 0.38 | 5.901e+07 | 6.567e+07 | 3.02 |
| 0.5 | optimistic compensated SERF | 1 | 5.524e+07 | 6.163e+07 | 11.4 |
| 0.5 | optimistic compensated SERF | 3.14 | 4.947e+07 | 5.538e+07 | 97.7 |
| 0.5 | optimistic compensated SERF | 10 | 4.375e+07 | 4.917e+07 | 976 |
| 0.5 | Mx benchmark | 0.38 | 4.711e+07 | 5.282e+07 | 247 |
| 0.5 | Mx benchmark | 1 | 4.703e+07 | 5.273e+07 | 255 |
| 0.5 | Mx benchmark | 3.14 | 4.63e+07 | 5.194e+07 | 341 |
| 0.5 | Mx benchmark | 10 | 4.322e+07 | 4.859e+07 | 1.22e+03 |
| 1 | optimistic compensated SERF | 0.38 | 1.315e+07 | 1.448e+07 | 5.84e+22 |
| 1 | optimistic compensated SERF | 1 | 1.228e+07 | 1.357e+07 | 2.19e+23 |
| 1 | optimistic compensated SERF | 3.14 | 1.109e+07 | 1.231e+07 | 1.88e+24 |
| 1 | optimistic compensated SERF | 10 | 9.905e+06 | 1.105e+07 | 1.88e+25 |
| 1 | Mx benchmark | 0.38 | 1.06e+07 | 1.179e+07 | 4.76e+24 |
| 1 | Mx benchmark | 1 | 1.059e+07 | 1.177e+07 | 4.92e+24 |
| 1 | Mx benchmark | 3.14 | 1.044e+07 | 1.161e+07 | 6.58e+24 |
| 1 | Mx benchmark | 10 | 9.793e+06 | 1.093e+07 | 2.35e+25 |
| 2 | optimistic compensated SERF | 0.38 | 4.277e+06 | 4.757e+06 | 9.17e+55 |
| 2 | optimistic compensated SERF | 1 | 3.972e+06 | 4.431e+06 | 3.16e+56 |
| 2 | optimistic compensated SERF | 3.14 | 3.589e+06 | 4.02e+06 | 2.63e+57 |
| 2 | optimistic compensated SERF | 10 | 3.215e+06 | 3.616e+06 | 2.62e+58 |
| 2 | Mx benchmark | 0.38 | 3.436e+06 | 3.854e+06 | 6.64e+57 |
| 2 | Mx benchmark | 1 | 3.43e+06 | 3.848e+06 | 6.86e+57 |
| 2 | Mx benchmark | 3.14 | 3.383e+06 | 3.797e+06 | 9.18e+57 |
| 2 | Mx benchmark | 10 | 3.18e+06 | 3.578e+06 | 3.28e+58 |

## Copper (5.8e7 S/m) checks

| d (m) | t (m) old/new | bulk amp 24 Hz old/new | bulk amp 18 Hz old/new | finite-loop ratio at fc old/new | critical Cu thickness PCB 1 pT old/new (m) | critical Cu thickness SERF 1 pT old/new (m) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.2270/0.2420 | 4.92e-08/1.61e-08 | 4.69e-07/1.78e-07 | 2.38e-08/7.03e-09 | 0.2503/0.2378 | 0.2552/0.2508 |
| 1 | 0.5467/0.5733 | 2.51e-18/3.5e-19 | 5.73e-16/1.04e-16 | 4.26e-18/7.78e-19 | 0.2833/0.2779 | 0.2547/0.2492 |
| 2 | 1.2017/1.2537 | 2.06e-39/4.34e-41 | 3.13e-34/1.11e-35 | 5.56e-35/4.71e-36 | 0.3030/0.2978 | 0.2968/0.2904 |

Critical copper thickness: wall thickness at which the capacity-bound minimum source power at 5.8e7 S/m equals the design ceiling, holding source coil, receiver, clearances, carrier and 2 Hz window fixed (new diagnostic; the manuscript quotes the evaluated wall thicknesses).

Reproduce: PYTHONPATH=src python scripts/magnetic_wall_nb_2026_11_05.py. Raw data: summary.json here and in the aligned_carriers / atomic_receiver folders with the same tag.
