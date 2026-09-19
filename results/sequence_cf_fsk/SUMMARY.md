# Whitened-Viterbi CF-FSK search: verified results

The selected tested desktop operating point uses 2 s symbols and carries
1.333 bit/s before packet overhead, four times the earlier 0.333 bit/s point.
Eight-byte packets with six known trailer symbols carry 1.067 bit/s.
The source remains the selected shaped rotor at a 0.5 m receiver distance.

## Error-count validation

Every final estimate has at least 100 observed bit errors in at least 30
independent errored packets. Errors within a packet are correlated. The 95%
intervals below operate on independent packet bit-error fractions and remain
valid under the error-count stopping rule; they do not assume independent bits.

| Symbol period (s) | Nominal bit/s | Packet bit/s | Bits tested | Bit errors | Errored packets | Measured BER | 95% BER interval |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1.333 | 1.067 | 4,032,000 | 134 | 30 | 3.3e-05 | 7.2e-06 to 0.00012 |
| 1.75 | 1.524 | 1.219 | 48,000 | 181 | 33 | 0.0038 | 0.0019 to 0.011 |

The 1.75 s case uses a 40% transition fraction and spacing-dwell product 0.075.
Its BER interval lies above the 0.001 target. The 2 s case uses a 40% transition
fraction and spacing-dwell product 0.1; its interval lies below that target.
This is the best verified point in the evaluated grid. The retained search
records include computation-budget exclusions at some shorter periods.

## Selected operating point

- Carrier: 50.3 Hz; seven tones from 50.05 to 50.55 Hz.
- Tone spacing: 0.083333 Hz; transition: 0.8 s; dwell: 1.2 s.
- Rotor speed: 1501.5 to 1516.5 rpm.
- Peak total torque: 1.2134 N m.
- All-transition cycle RMS torque: 0.2267 N m.
- Maximum individual-symbol RMS torque: 0.4974 N m,
  a conservative bound on long-message duty.
- Peak motoring / absorbed mechanical power:
  191.75 / 175.94 W.
- Maximum sampled carrier tracking error: 4.73e-05 Hz.

These values satisfy the assumed 2 N m peak, 0.5 N m continuous and 400 W
bidirectional drive specifications. A separate noiseless packet constructed from
the achieved drive phase decodes correctly; its maximum complex-envelope
deviation from the nominal waveform is 4.9e-05.
The Monte Carlo trials use the nominal sequence waveform. Assembly strength,
fatigue, containment and environmental coupling
remain the separate engineering validation described in the manuscript.

## Why throughput improved

The old detector discarded transitions and required orthogonal dwell tones.
The new decoder uses the entire continuous-phase waveform and its sequence
memory. Tone spacing times dwell duration is now 0.1, allowing closer tones and
shorter symbols. The known byte mapping is enforced directly in the trellis.
The earlier operating point used 8 s symbols and a spacing-dwell product of 1.

The Gaussian benchmark is 3.9045 bit/s over the new 2 Hz sampled baseband;
the earlier 1.5625 Hz calculation was 3.9031 bit/s. The nominal mapped rate
uses about 34% of this benchmark, and the short-packet rate uses about 27%.
The revised bandwidth therefore contributes very little to the rate gain.

The old error calculation was a conditional symbol-by-symbol test over 49
contexts. The new results are bit-error measurements on complete random-byte
packets with acquired initial phase, carrier and timing. Packet rates include
the trailer and exclude acquisition overhead.

## Whitening and numerical checks

The causal 13-tap filter retains 6 s of memory at 2 complex samples/s.
Whitened noise has 1.05% RMS PSD deviation over the full sampled band; an
8 s, 17-tap filter reduces this to 0.37%. Both filters gave identical payload
decisions on 500 paired noisy packets. A separate 500-packet check doubled
the sampling rate. These small checks are convergence diagnostics, and their
low error counts are not presented as BER estimates.

Noise is synthesized from the original receiver spectrum independently of
the whitening approximation. Viterbi keeps phase, tone and filter-memory
states, with exact add-compare-select updates and traceback. There is no beam
pruning. Initial phase and timing are assumed acquired. The receiver and drive
parameters retain the paper's stated modeling assumptions.

The full regression suite is recorded in `test_report.txt`. Reproduction,
search grids, complete packet counts and statistical details are in `README.md`
and the adjacent JSON files. The manuscript and public release have not yet
been updated to this new operating point.
