# Bounded-state streaming soft receiver

## Implemented architecture

CompactInnovations uses finite-startup Cholesky coefficients for the first 65 samples, then a fixed 64th-order complex linear predictor. Signal and noise use the same transform. The exact marginal stationary covariance is retained at startup; the later whitening approximation is tested against the aliased physical noise spectrum. Memory and filtering work per block are independent of packet length.

StreamingSoftReceiver consumes fixed-clock samples incrementally. Each trellis state labels the discrete modulation phase and recent tone history, carrying its own causal SOS filter mean state and 64-sample mean history. At a state merge only the best forward survivor carries its filter history; all incoming retained edges are kept in a fixed-lag max-log graph. This is reduced-state sequence estimation, not exact merging of equal continuous states. That approximation is assessed by increasing recent history and comparing to wider full-history paths. No packet reruns, dense packet covariance, or full-packet graph are used in the streaming decoder.

The working candidate retains four recent symbols, permits 16,384 states (observed peak 10,878), and uses 12 symbols of lookahead. The graph holds at most 15 symbol layers before emission/removal. No state-budget pruning or missing bit alternatives occurred on the independent diagnostic records. Byte timing, carrier, phase and initial deterministic mean are still acquired inputs. End-of-frame tail information is used and remaining outputs are flushed causally when the last sample arrives.

## Processing deadlines and bounded storage

Tests feed each block only when it is available, measuring the complete push call including soft output emission. This is an accelerated input replay, not a hardware real-time scheduler test. The 4 Hz output clock yields blocks of seven or eight samples, separated by 1.75 or 2 s. The nominal symbol duration is 224/114 = 1.964912 s. One Python process and one numerical-library thread were used. Startup model preparation takes about 0.1 s in the primary diagnostic. All measured push times are well below the minimum inter-block interval; the replay therefore predicts no accumulating processing backlog under this workload.

| Record | History | Peak states | Max block processing (s) | Total compute (s) | Persistent arrays (MB) |
|---|---:|---:|---:|---:|---:|
| 32 bytes, 20 Hz, reference amplitude | 4 | 10878 | 0.0735 | 4.655 | 41.09 |
| 32 bytes, 50.3 Hz, amplitude scale 0.5 | 3 | 1764 | 0.0097 | 0.677 | 6.33 |
| 32 bytes, 50.3 Hz, amplitude scale 0.5 | 4 | 10878 | 0.0772 | 4.683 | 41.09 |
| 32 bytes, 50.3 Hz, amplitude scale 0.5 | 5 | 75264 | 0.4416 | 33.504 | 267.70 |
| 128 bytes, 50.3 Hz, amplitude scale 1 | 4 | 10878 | 0.0893 | 18.615 | 41.09 |

The 128-byte record is deliberately a long-stream stress test, not a proposed latency-compliant frame. The four-symbol receiver retains the same persistent array storage for 32 and 128 bytes. Reported memory excludes Python object overhead, temporary work arrays, imported libraries and caller-owned input/output storage; it is not peak process RAM. Those other working arrays are also bounded by the configured state budget, history order and block size, not packet length.

Twelve-symbol lookahead corresponds to 23.58 s after a byte finishes, with sample-grid rounding of at most 0.25 s. This delay overlaps continued reception; it is not an extra complete decoding pass after the frame. The last push includes final traceback/flush and is included in the maximum timing above. No FEC decoder or acquisition algorithm is yet included in the end-to-end latency claim.

## Whitening validation

| Predictor order | Maximum whitened PSD error | RMS PSD error |
|---:|---:|---:|
| 16 | 0.009575 | 0.003087 |
| 32 | 6.872e-05 | 1.9e-05 |
| 64 | 2.812e-09 | 7.873e-10 |
| 128 | 1.676e-14 | 4.563e-15 |

These figures use the 20 Hz carrier with the specified NHNM/isolation background. The weak-record output additionally saves its own spectral-error check. The 64-sample history corresponds to 16 s of sampled noise memory. Coefficients remain fixed during a record; carrier/background changes require a newly validated model.

## Soft-output validation and correction of the old reference

The earlier width-512/1024 packet reference was not converged on four bit alternatives, despite those widths agreeing. A 16,384-path forced-opposite search finds lower-cost physical candidate sequences at all four. The streaming result agrees closely with those widened full-history checks:

| Bit | Streaming history 4 | Widened forced reference | Absolute difference |
|---:|---:|---:|---:|
| 40 | 65.745878 | 65.745764 | 0.0001142 |
| 56 | -37.947676 | -37.947481 | 0.0001946 |
| 112 | -45.341762 | -45.342183 | 0.0004209 |
| 136 | -50.068779 | -50.068595 | 0.0001842 |

These wider paths are stronger references, not a proof of global full-packet optimality. Agreement between two beam widths alone did not certify the previous soft decoder. Do not reuse the old confidence vectors as ground truth.

History convergence on independent records:

- strong_history4_vs5: max absolute change 0.0007847; relative L2 change 9.359e-06; 0 sign changes across 256 finite bit metrics.
- weak_history4_vs5: max absolute change 0.01955; relative L2 change 0.0004034; 0 sign changes across 256 finite bit metrics.

Doubling the lookahead from 12 to 24 symbols in the strong history-4 case also produced the saved comparison; no missing alternatives occurred. Four streaming unit tests cover exact one-byte likelihoods, exhaustive 65,536-message two-byte graph likelihoods when no state merging is possible, causal chunk/future independence, bounded graph size, and compact whitening/startup. The 20 prior receiver/sequence/environmental tests remain relevant.

## Outcome and remaining work

Computational feasibility is demonstrated for these streaming replay cases: fixed persistent storage, bounded measured work per block with substantial deadline margin, and decisions emitted during reception. The expensive full-packet rerun architecture is no longer required for the online prototype. This is not yet post-decoding BER qualification, a hard real-time guarantee, or acquisition validation. Next integrate the emitted max-log outputs with error correction, test more signal levels and modulation settings, and validate post-decoding performance against independent trials. A C++ port is optional optimization at this stage, not a prerequisite for keeping pace in these tests.
