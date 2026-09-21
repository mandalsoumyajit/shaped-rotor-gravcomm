# Causal sampled receiver model

The implemented receiver processes data causally using a stable sixth-order SOS low-pass at 16 complex samples/s and decimates to 4 complex samples/s. Signal and noise undergo the identical streaming operations. Cholesky innovations use the full finite-record covariance and past samples only. No hard spectral projection is used by the receiver.

## Fixed receiver specification

- Full -3 dB bandwidth: 2 Hz; two-sided noise-equivalent bandwidth: 2.020583 Hz.
- Maximum folded/main noise PSD within +-1 Hz: 9.289e-06. Order 6 is the lowest tested even order satisfying the 1e-4 design budget; order 4 fails.
- Group delay: 0.6070 s at center; up to 1.0559 s within the -3 dB band.
- Step error remains below 1e-4 after 5.4375 s. Remaining impulse-response energy falls below 1e-6 after 4.5625 s. These are settling criteria, not a hard delay or per-packet warm-up requirement.

Bandwidth, noise spectrum and filter coefficients are held fixed versus frequency offset while center tuning remains allowed. The 4 Hz Nyquist span and the 2 Hz -3 dB bandwidth are different quantities. Transition-band information can remain useful; the previous ideal-band capacity bound is not automatically applicable.

## Noise and likelihood checks

- Exact whitened finite-record covariance differs from identity by at most 1.35e-14.
- Independent causal-noise simulation: covariance-lag error 0.797% of variance; innovation power 1.00340; largest measured lag correlation 0.007658.
- 1,024-sample covariance factorization: 0.1239 s; streaming 1,024 innovations: 0.0028 s in this local run.
- Common-record initialization: 16 s burn-in differs from the longer-prehistory reference by 3.5e-11 noise RMS; 32 s by 0.
- Seven tests pass: stream/chunk invariance and no future dependence; alias sum versus lag selection; alias budget; covariance/Cholesky prefix consistency; exhaustive 256-message metric comparison; continuous waveform sampling on a fixed clock against independent phase quadrature; filter continuity across packet boundaries.

The metric is (y-mu)^H C^-1 (y-mu). It is exact for the specified finite stationary Gaussian record, to numerical covariance accuracy. Both candidate means and observations use the same lower-triangular transform. Dense innovations do not make a finite-memory sequence trellis exact; a sequence detector must carry the required history or pass a separate approximation check.

## Scope and remaining integration

The model starts at an explicitly sampled, calibrated acceleration-envelope boundary. Upstream transducer calibration and analog acquisition/anti-aliasing are not silently assumed to have been validated. Stationary noise state is continuous across packets; the offline simulation burn-in is not per-packet latency. Known deterministic mean state is separately specified. Acquisition, complete receiver startup, implementation timing and BER remain unvalidated. Existing transmitter/noise assumptions are preserved; no receiver sensitivity improvement is credited.

Next: connect the matched innovations metric to sequence detection, using the exact short-message likelihood as the reference and the saved memory diagnostics to control approximations. Only then resume coded comparisons.
