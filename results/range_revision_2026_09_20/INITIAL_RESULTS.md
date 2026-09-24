# Initial range-revision calculations

These are new calculations for the coauthor revision, separate from the original manuscript and baseline results. The manuscript has not been rewritten.

## Methods and status

- Screen: 112 parameter combinations, 32 packets per mechanically feasible point. Screening zeros are not treated as BER evidence.
- Current seven-tone, byte-constrained, full-waveform whitened Viterbi decoder; eight-byte payload and six-symbol trailer; acquired clock/phase and stationary receiver noise.
- Accepted finite-volume source quadrature; 128/256 phase agreement and reproduction of the original 0.5/1/2 m fields. No new FEA.
- Four regression tests cover source similarity, physical carrier/noise translation, baseline message-duty torque, and agreement with the original decoder at its baseline.
- Validation uses separate seeds and packet-aware time-uniform intervals. When two candidates are tested for one scenario, the table uses 97.5% bounds for each (95% simultaneous within that scenario). Bounds are not simultaneous across scenarios.
- Mechanics columns are analytic trajectory/duty screens. New achieved-drive trajectories, environmental noise, acquisition, and hardware phase stability remain to be evaluated.

## Desktop source: unchanged receiver and drive

| Distance (m) | Source A2 (nGal) | Gaussian benchmark (bit/s) | Tested mapped rate (bit/s) | Payload incl. trailer (bit/s) | BER interval | Status |
|---:|---:|---:|---:|---:|---|---|
| 0.5 | 51.803 | 3.90452 | 1.33333 | 1.06667 | [9.19e-07, 0.00096] | BER upper bound below target |
| 0.75 | 9.5625 | 0.760687 | 0.190172 | 0.152138 | [5.86e-06, 0.000988] | BER upper bound below target |
| 1 | 2.9594 | 0.154948 | 0.038737 | 0.0309896 | [5.14e-06, 0.000979] | BER upper bound below target |
| 1.5 | 0.57565 | 0.00846606 | 0.00211652 | 0.00169322 | [0.000463, 0.00099] | BER upper bound below target |
| 2 | 0.18118 | 0.000880663 | 0.000220166 | 0.000176132 | [0.000643, 0.000999] | BER upper bound below target |

The tested 2 m desktop candidate requires 3.36 hours per symbol and 100.9 hours per 64-bit packet including its trailer. Stationary noise and accurate acquired phase over this duration are major untested assumptions; this is not a practical 2 m desktop demonstration.

## Enlarged source at 2 m

Both enlarged scenarios use 63.31 kg moving mass, 1.00 m diameter, and 32 times the desktop polar inertia. Assumed drive ratings are 16 N m peak, 4 N m continuous, and 1.6 kW bidirectional mechanical power; center drag is assumed to be 0.4 N m. These are different drive specifications from the desktop.

| Scenario | Carrier (Hz) | Receiver resonance (Hz) | Gaussian benchmark (bit/s) | Tested payload rate (bit/s) | BER upper bound | Status |
|---|---:|---:|---:|---:|---:|---|
| same receiver | 50.300 | 50.300 | 0.420224 | 0.084045 | 0.0009997 | BER upper bound below target |
| slow same receiver | 25.150 | 50.300 | 4.74997e-05 | — | — | benchmark only |
| slow retuned receiver | 25.150 | 25.150 | 0.742796 | 0.14856 | 0.00098946 | BER upper bound below target |

Maintaining the original rotation speed increases the centrifugal stress scale by four relative to the desktop at that speed. Halving the speed preserves that scale but moves the carrier to 25.15 Hz. The unchanged resonant receiver performs poorly there; retuning it is a separate design assumption, retaining the assumed proof mass, Q, temperature, and readout noise. Neither case is mechanically qualified.

## One large rotor versus smaller coherent rotors

Equal total moving mass (63.31 kg) and equal nominal centrifugal stress scale; actual nonintersecting rotor locations are included. All rotors in an array carry the same waveform and phases are ideally aligned at the receiver. Footprints and drive/support requirements differ.

| Rotors | Individual mass (kg) | Carrier (Hz) | Total A2 (nGal) | Benchmark: fixed receiver (bit/s) | Benchmark: retuned receiver (bit/s) |
|---:|---:|---:|---:|---:|---:|
| 1 | 63.31 | 25.15 | 5.919 | 4.75e-05 | 0.7428 |
| 2 | 31.65 | 31.69 | 3.342 | 2.4172e-05 | 0.28 |
| 4 | 15.83 | 39.92 | 2.142 | 2.896e-05 | 0.11358 |
| 8 | 7.91 | 50.30 | 1.282 | 0.03775 | 0.03775 |

These are source/noise benchmarks, not array BER results. Smaller rotors reduce total inertia but also reduce coherent gravitational coupling at fixed total mass. The receiver frequency response can reverse the ranking; the retuned column does not hold receiver design fixed. Parallel independent subchannels remain to be simulated.

## Numerical diagnostics

- desktop_2m, Ts=12112.1 s, q=0.5: 1000 packets; 3/4-symbol filters give 51/51 bit errors with 0 disagreeing decoded bits. Independent doubled-sampling run: 63 bit errors. Diagnostic counts are not precise rare-event BER estimates.
- scale2_2m_slow_retuned_receiver, Ts=7.18008 s, q=0.25: 1000 packets; 3/4-symbol filters give 68/68 bit errors with 0 disagreeing decoded bits. Independent doubled-sampling run: 27 bit errors. Diagnostic counts are not precise rare-event BER estimates.
- scale2_2m_same_receiver, Ts=25.3833 s, q=0.5: 1000 packets; 3/4-symbol filters give 5/5 bit errors with 0 disagreeing decoded bits. Independent doubled-sampling run: 7 bit errors. Diagnostic counts are not precise rare-event BER estimates.
- scale2_2m_slow_retuned_receiver, Ts=14.3602 s, q=0.5: 1000 packets; 3/4-symbol filters give 6/6 bit errors with 0 disagreeing decoded bits. Independent doubled-sampling run: 2 bit errors. Diagnostic counts are not precise rare-event BER estimates.

## Next calculations

1. Refine duration/spacing around the qualifying candidates; the coarse grid does not identify maximum rates. Include the known 1.75 s / q=0.075 desktop candidate in the refined search.
2. Complete sampling/filter checks for the final chosen candidates, and evaluate achieved rather than prescribed actuator phase.
3. Add receiver environmental-noise and long-term phase/frequency-stability sensitivity; multi-hour symbols make these essential.
4. Reassess enlarged-source mechanical loads and receiver assumptions, then evaluate an explicit parallel-subchannel array.
5. Restructure the main text/appendix and update abstract/conclusions after settling those results.

## Reproduction

Use Python with NumPy, SciPy, Matplotlib, and Numba, with one BLAS thread. From the submission folder:

```text
python scripts/range_revision.py prepare
python scripts/range_revision.py screen --packets 32 --workers 4
python scripts/range_revision.py validate --max-packets 25000 --workers 4
python scripts/range_array_comparison.py
python scripts/range_revision_checks.py --packets 1000 --workers 4
python tests/test_range_revision.py
python scripts/summarize_range_revision.py
```

Candidate-specific follow-up commands and seeds are recorded in RUN_COMMANDS.md and individual JSON files. Validation can resume with a larger packet limit. Do not rerun scripts concurrently against the same output file.
