# Communications checks before manuscript revision

All random-message runs use the selected shaped rotor at 0.5 m, 2 s symbols, 40% transitions, spacing–dwell product 0.1, eight-byte packets and six known trailer symbols. Initial timing, phase, frequency, tone and prehistory are acquired. Nominal and trailer-inclusive rates are 1.333 and 1.067 bit/s. No added ECC.

## Sampling, whitening and technical-noise checks

| Sampling (Hz) | FIR memory (s) | Added noise ASD (SI/√Hz) | Bits | Errors / errored packets | BER estimate or upper bound | Packet-aware 95% interval |
|---:|---:|---:|---:|---:|---:|---|
| 4 | 6 | 0 | 2816000 | 132 / 30 | 4.688e-05 | [1.51e-05, 0.000171] |
| 4 | 6 | 1e-10 | 32000 | 688 / 135 | 0.0215 | [0.0157, 0.0328] |
| 4 | 6 | 2e-10 | 32000 | 5805 / 480 | 0.1814 | [0.166, 0.196] |
| 4 | 8 | 0 | 384000 | 13 / 3 | < 0.0009307 | [1.22e-06, 0.000931] |
| 8 | 6 | 0 | 2912000 | 144 / 31 | 4.945e-05 | [1.59e-05, 0.00017] |

The error-count stopping criterion is at least 100 bit errors in at least 30 errored packets. The longer-memory diagnostic has a 6,000-packet cap; when that criterion is unmet it is reported as an upper bound. Confidence sequences operate on independent packet error fractions and account for arbitrary within-packet error clustering. Intervals are individual 95% statements, with no pooled or simultaneous coverage claim. Configurations use common seed indices and can be correlated with one another.

Added technical noise is independent flat Gaussian acceleration noise added in quadrature to the physical PSD. The whitening filter is redesigned using that total PSD. This assumes the additional spectrum is known to the receiver. Symbol timing and tone spacing are held fixed; the waveform has not been reoptimized for either added-noise case. Coherent interference is outside this test. The corresponding Gaussian benchmarks are 2.9202 and 1.8848 bit/s, compared with 3.9045 bit/s without added noise. Poor BER at the fixed waveform therefore does not establish a fundamental rate ceiling under those noise levels.

## Achieved-drive sensitivity

| Symbol duration (s) | Message | Max envelope error vs nominal | Change under 200→400 Hz drive-grid refinement | Noisy packet decisions changed / paired trials |
|---:|---|---:|---:|---:|
| 2.0 | worst_duty | 3.98e-06 | 1.19e-05 | 0 / 500 |
| 2.0 | random | 3.97e-06 | 1.19e-05 | 0 / 500 |
| 1.75 | worst_duty | 5.01e-06 | 1.5e-05 | 0 / 500 |
| 1.75 | random | 5e-06 | 1.5e-05 | 0 / 500 |

The tests include repeated bytes 251,42 (the worst sustained-duty cycle) and a fixed random message at both 2 s and 1.75 s. Each pair uses identical physical noise and the same nominal decoder templates. These are fixed-message sensitivity checks; their error counts are not estimates of random-message BER. Trailer duty lowers packet RMS relative to the indefinitely repeated payload. At 1.75 s the drive can follow the waveform while the worst-message continuous-duty rating remains exceeded.

## Energy and centrifugal references

The closed all-transition cycle lasts 98 s. Positive motor work is 1390.879749 J, absorbed work 616.563642 J, and drag heat 774.316106 J. The energy-balance residual including kinetic-energy change is -1.1e-09 J.

At the selected maximum speed of 1516.5 rpm, stored energy is 3985.632 J and the scaled 10 mm regional centrifugal stress average is 18.0828 MPa. This remains a centrifugal-only calculation.

## Capacity convergence

| Distance (m) | Capacity in 2 Hz (bit/s) |
|---:|---:|
| 0.5 | 3.90452497 |
| 1.0 | 0.154947973 |
| 2.0 | 0.000880663433 |

Largest relative change on doubling the frequency grid: 1.07e-09. The 0.5 m value is unchanged to displayed precision between 2, 4 and 8 Hz bands.

## Reproduction and scope

Scripts: `check_communications_revision.py`, `check_revision_drive.py`, `check_revision_capacity.py`, and `summarize_communications_checks.py`. Run from the submission directory with PYTHONPATH=src in the configured WSL Python environment, one BLAS thread, and the sequence/test dependencies. BER configurations are stored in each JSON file; use --fs and --memory accordingly, with --technical-asd for the noise checks. The 8 s FIR diagnostic uses --max-packets 6000. The JSON records retain all packet error counts and are resumable.

See `test_report.txt` for the regression result and `manifest.json` for SHA-256 hashes of the checked scripts, core modules and numerical records. The manuscript and public release remain unchanged during these checks.
