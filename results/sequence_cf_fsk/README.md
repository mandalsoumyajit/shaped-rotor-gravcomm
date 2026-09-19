# CF-FSK sequence-decoder search

The source is the selected 7.9135 kg shaped desktop rotor, evaluated at 0.5 m.
The carrier remains 50.3 Hz; the receiver has the paper's 3.1 g proof mass,
Q = 637000, temperature 300 K and 100 fm/sqrt(Hz) displacement readout.
The drive retains 2 N m peak, 0.5 N m continuous, 400 W motoring/absorption,
20 ms compensated demand response and 0.05 N m center-speed drag.

## Receiver and waveform

`src/gravcomm/whitened_viterbi.py` fits a causal minimum-phase FIR to the inverse
square root of the receiver acceleration-noise PSD. The complex-baseband sample
PSD is `2 fs S_a(fc + f) / A2^2`. Both signal and observations pass through the
same whitening filter. Simulation noise is synthesized independently from the
original receiver spectrum, so residual whitening error remains in the trials.

Standard Viterbi uses the complete quintic transition and dwell waveform, keeps
accumulated phase and all finite whitening-filter memory, performs
add-compare-select updates, and traces back the surviving path. An algebraic
reachability constraint removes unreachable phase states. There is no beam
pruning. The final BER runs also restrict each three-tone word to the 256 valid
byte encodings. The transmitted data are never supplied to the decoder.

Initial tone, phase, frequency and timing are acquired inputs. Six known zero-tone
trailer symbols supply the trailing observations for sequence detection. Each
reported packet carries eight unknown bytes (24 payload symbols). Packet rates
include this trailer; acquisition and synchronization overhead remain unspecified.
Independent packet trials begin from the same known prehistory, with independent
random bytes and stationary colored noise. They are conditional acquired-packet
results.

## Search and validation

- `screen.json`: original coarse grid, before the trailer was added.
- `refinement.json`: refined spacing and duration grid with a known trailer.
- `fast_refinement.json`: periods below 2 s, including explicit actuator and
  trellis-size exclusions. Some short-period cases exceeded the one-million-state
  computation budget; the search does not establish a global rate maximum.
- `boundary_refinement.json`: refinement near the motor duty constraint.
- `selected_ber.json`: error-count-driven bit-level trials for the 2 s candidate.
- `boundary_ber.json`: the faster 1.75 s comparison, using the same final decoder.
- `selected_actuator.json`: constrained-drive trajectory, energy and decoding checks.
- `whitening_convergence.json`: paired six/eight-second filter checks on identical
  received packets. The 13-tap filter has 1.05% RMS PSD error; the 17-tap filter
  reduces that to 0.37%. These figures refer to the entire sampled band.
- `sampling_convergence.json`: doubled-sampling diagnostic. Low counts in the
  convergence runs supply checks and bounds, not precise BER estimates.

The principal sampling rate is 2 complex samples/s, spanning a 2 Hz baseband.
The Gaussian benchmark for that band is 3.9045 bit/s, versus 3.9031 bit/s for the
earlier 1.5625 Hz band. Sample-rate checks use 4 samples/s. These are discretized
calibrated-baseband models; physical analog filtering and acquisition remain
implementation tasks.

The final error-count target is at least 100 bit errors in at least 30 independent
errored packets. All packet error counts are retained. Bits within a packet are
correlated and are not treated as independent Bernoulli trials. The script gives
an empirical-Bernstein interval for independent packet error fractions and a
separate bound from one independently sampled bit per packet. Summable confidence
spending makes each bound valid across scheduled checks and adaptive stopping.
The two bounds are separately 95% statements. A low-count run is labeled as a
bound rather than a precise BER estimate.

The empirical-Bernstein formula follows Maurer and Pontil (2009),
https://arxiv.org/abs/0907.3740 . The two-sided form uses log(4/delta); scheduled
look j receives delta = 0.05/[j(j+1)].

The final report uses a separate packet-aware confidence sequence implemented in
`packet_confidence.py`. For each hypothesized mean m, it mixes fixed nonnegative
products of `1 + lambda(m) * (X - m)`, where X is a packet's bit-error fraction.
Positive and negative bets give the two sides, each with a 0.025 error budget.
Ville's inequality supplies time-uniform coverage. This uses all packet fractions
and avoids the independence assumption between bits in a burst. The construction
follows the betting framework of Waudby-Smith and Ramdas,
https://arxiv.org/abs/2010.09686 . Each reported method has its own coverage
statement; the report does not intersect several 95% intervals and call the
intersection 95%.

## Reproduction

Install the package with optional `sequence` and `test` dependencies. Use one BLAS
thread; the BER driver can distribute independent packets over four workers.

```
python -m pytest tests -q
python scripts/search_sequence_cf_fsk.py
python scripts/search_sequence_cf_fsk.py --refine
python scripts/refine_fast_sequence.py
python scripts/refine_sequence_boundary.py
python scripts/measure_sequence_ber.py --legacy-prefix 7500 --max-packets 100000
python scripts/measure_sequence_ber.py --boundary --legacy-prefix 750 --max-packets 3000
python scripts/check_selected_whitening.py
python scripts/verify_sequence_actuator.py --selected
python scripts/summarize_sequence_search.py
```

The legacy-prefix options reproduce the sequential random-number streams used
before parallel execution was added. Subsequent packets use independent seeded
streams indexed by packet number. The public v0.2.0 code release predates these
new files; these results currently belong to the manuscript workspace.
