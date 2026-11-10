# Magnetic-wall comparison for hop-limited CF-FSK designs (2026-11-05)

Released model code unchanged; only rotor radius (0.25 m x scale), rotor mass and power ceiling (new peak mechanical power) change. Coil copper mass 10% of rotor mass, coil outer radius = rotor radius, t = d - R - 0.0327 m, receivers, noise levels, BPSK benchmark and the fc +/- 1 Hz capacity window are as in the 2026-09-21 studies. The old columns are recomputed here and match the archived values.

| d (m) | R old/new (m) | t old/new (m) | ceiling old/new (W) | fc (Hz) | PCB 1 pT old/new (1e6 S/m) | SERF 1 pT old/new (1e6 S/m) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.2403/0.2335 | 0.2270/0.2338 | 120.49/46.66 | 24 | 69.3/61.6 | 71.8/65.4 |
| 1 | 0.4206/0.4087 | 0.5467/0.5586 | 1980.24/766.59 | 24 | 16.1/14.8 | 14.8/14 |
| 2 | 0.7656/0.7401 | 1.2017/1.2272 | 29680.14/11238.96 | 18 | 4.91/4.56 | 4.87/4.58 |

## New designs: aligned-carrier PCB receiver

| d (m) | fc | background (pT) | BPSK crossing (S/m) | capacity crossing (S/m) | copper BPSK W | copper capacity W |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 24 | 0.38 | 5.897e+07 | 6.591e+07 | 34 | 3.99 |
| 0.5 | 24 | 1 | 5.487e+07 | 6.159e+07 | 131 | 15 |
| 0.5 | 24 | 3.14 | 4.862e+07 | 5.489e+07 | 1.13e+03 | 129 |
| 0.5 | 24 | 10 | 4.238e+07 | 4.815e+07 | 1.13e+04 | 1.28e+03 |
| 1 | 24 | 0.38 | 1.42e+07 | 1.575e+07 | 2.94e+21 | 1.68e+20 |
| 1 | 24 | 1 | 1.33e+07 | 1.481e+07 | 1.13e+22 | 6.27e+20 |
| 1 | 24 | 3.14 | 1.193e+07 | 1.334e+07 | 9.84e+22 | 5.38e+21 |
| 1 | 24 | 10 | 1.055e+07 | 1.186e+07 | 9.84e+23 | 5.37e+22 |
| 2 | 18 | 0.38 | 4.276e+06 | 4.76e+06 | 4.65e+56 | 1.54e+54 |
| 2 | 18 | 1 | 4.079e+06 | 4.557e+06 | 1.34e+57 | 4.13e+54 |
| 2 | 18 | 3.14 | 3.71e+06 | 4.163e+06 | 1.05e+58 | 3.1e+55 |
| 2 | 18 | 10 | 3.318e+06 | 3.739e+06 | 1.03e+59 | 3.04e+56 |

## New designs: atomic receiver

| d (m) | sensor | background (pT) | BPSK (S/m) | capacity (S/m) | copper capacity W |
|---:|---|---:|---:|---:|---:|
| 0.5 | optimistic compensated SERF | 0.38 | 6.254e+07 | 6.976e+07 | 1.21 |
| 0.5 | optimistic compensated SERF | 1 | 5.843e+07 | 6.536e+07 | 4.56 |
| 0.5 | optimistic compensated SERF | 3.14 | 5.214e+07 | 5.854e+07 | 39.2 |
| 0.5 | optimistic compensated SERF | 10 | 4.594e+07 | 5.176e+07 | 392 |
| 0.5 | Mx benchmark | 0.38 | 4.958e+07 | 5.574e+07 | 98.9 |
| 0.5 | Mx benchmark | 1 | 4.949e+07 | 5.564e+07 | 102 |
| 0.5 | Mx benchmark | 3.14 | 4.87e+07 | 5.478e+07 | 137 |
| 0.5 | Mx benchmark | 10 | 4.537e+07 | 5.114e+07 | 489 |
| 1 | optimistic compensated SERF | 0.38 | 1.356e+07 | 1.492e+07 | 7.95e+21 |
| 1 | optimistic compensated SERF | 1 | 1.268e+07 | 1.4e+07 | 2.98e+22 |
| 1 | optimistic compensated SERF | 3.14 | 1.147e+07 | 1.271e+07 | 2.57e+23 |
| 1 | optimistic compensated SERF | 10 | 1.026e+07 | 1.143e+07 | 2.56e+24 |
| 1 | Mx benchmark | 0.38 | 1.097e+07 | 1.219e+07 | 6.48e+23 |
| 1 | Mx benchmark | 1 | 1.096e+07 | 1.217e+07 | 6.69e+23 |
| 1 | Mx benchmark | 3.14 | 1.08e+07 | 1.201e+07 | 8.96e+23 |
| 1 | Mx benchmark | 10 | 1.015e+07 | 1.131e+07 | 3.2e+24 |
| 2 | optimistic compensated SERF | 0.38 | 4.419e+06 | 4.91e+06 | 4.06e+54 |
| 2 | optimistic compensated SERF | 1 | 4.108e+06 | 4.579e+06 | 1.4e+55 |
| 2 | optimistic compensated SERF | 3.14 | 3.717e+06 | 4.159e+06 | 1.17e+56 |
| 2 | optimistic compensated SERF | 10 | 3.334e+06 | 3.747e+06 | 1.16e+57 |
| 2 | Mx benchmark | 0.38 | 3.56e+06 | 3.99e+06 | 2.94e+56 |
| 2 | Mx benchmark | 1 | 3.555e+06 | 3.984e+06 | 3.04e+56 |
| 2 | Mx benchmark | 3.14 | 3.506e+06 | 3.932e+06 | 4.07e+56 |
| 2 | Mx benchmark | 10 | 3.298e+06 | 3.708e+06 | 1.45e+57 |

## Copper (5.8e7 S/m) checks

| d (m) | t (m) old/new | bulk amp 24 Hz old/new | bulk amp 18 Hz old/new | finite-loop ratio at fc old/new | critical Cu thickness PCB 1 pT old/new (m) | critical Cu thickness SERF 1 pT old/new (m) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.2270/0.2338 | 4.92e-08/2.98e-08 | 4.69e-07/3.04e-07 | 2.38e-08/1.38e-08 | 0.2503/0.2417 | 0.2552/0.2500 |
| 1 | 0.5467/0.5586 | 2.51e-18/1.04e-18 | 5.73e-16/2.66e-16 | 4.26e-18/2e-18 | 0.2833/0.2780 | 0.2547/0.2492 |
| 2 | 1.2017/1.2272 | 2.06e-39/3.09e-40 | 3.13e-34/6.07e-35 | 5.56e-35/1.59e-35 | 0.3030/0.2972 | 0.2968/0.2904 |

Critical copper thickness: wall thickness at which the capacity-bound minimum source power at 5.8e7 S/m equals the design ceiling, holding source coil, receiver, clearances, carrier and 2 Hz window fixed (new diagnostic; the manuscript quotes the evaluated wall thicknesses).

Reproduce: PYTHONPATH=src python scripts/magnetic_wall_hop_2026_11_05.py. Raw data: summary.json here and in the aligned_carriers / atomic_receiver folders with the same tag.
