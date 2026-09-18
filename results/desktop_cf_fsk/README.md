# Desktop CF-FSK numerical operating point

Selected source: `medium_rounded`, the highest-signal member of the three
completed rounded designs at fixed CAD mass. This is a bounded design choice,
not a global optimum. It is distinct from the 50 kg Bench deployment class.

| Parameter | Value |
|---|---|
| Carrier / steel inserts | 30 mm spokes, 20 mm rim, 20 mm plate, 3 mm junction fillets; two 70 mm diameter x 79.24 mm cylinders at 200 mm radius |
| CAD mass / polar inertia | 7.913509 kg / 0.31607131 kg m^2 |
| Second-harmonic peak at 0.5 m from axis | 5.1802594e-10 m/s^2; full aluminium and steel mass distribution |
| Alphabet | Seven received-signal tones, 49.83125 to 50.76875 Hz |
| Adjacent spacing | 0.15625 Hz |
| Symbol / transition / dwell | 8 s / 1.6 s / 6.4 s |
| Frequency interpolation | 10u^3 - 15u^4 + 6u^5, u=t/1.6 |
| Actual byte mapping | 8 bits to three base-7 digits; 256 of 343 words |
| Mapped payload rate | 1/3 bit/s = 20 bit/min, excluding framing, acquisition, pilots and FEC |
| Speed range | 1494.9375 to 1523.0625 rpm |
| Peak inertial / total torque | 1.09091 / 1.14085 N m |
| RMS total torque, all-transition pattern | 0.15489 N m |
| RMS total torque, repeated extreme hops | Approximately 0.315 N m |
| Peak motoring / absorbed mechanical power | 180.29 / 164.48 W |
| Assumed drive | 2 N m peak, 0.5 N m continuous, 400 W motoring and absorption, 20 ms first-order response |
| Assumed drag | Viscous; 0.05 N m at the 50.3 Hz signal centre |

The 0.351 bit/s alphabet information rate log2(7)/8 is NOT substituted for
the mapped payload rate. Each codeword carries one byte. Some invalid received
words can be detected; no error correction or packet BER is claimed.

The source field uses finite-element volume quadrature. Its faceted mesh mass
is about 0.11% below the exact CAD mass. Drive sizing uses exact CAD inertia.
The second harmonic alone is decoded; the stronger higher harmonics near the
source are not folded into this channel's signal power.

## Receiver and detector

Parameters of the illustrative oscillator model: 3.1 g, 50.3 Hz, Q=637000,
300 K, flat displacement readout ASD 100 fm/sqrt(Hz). Structural thermal and
input-referred readout noise are both included. The input acceleration ASD is
5.15e-11 m/s^2/sqrt(Hz) at the centre and about 1.94e-10 at the outer tones.
These are model predictions, not a measured noise specification across this
band. No cooling or altered resonance is assumed for this operating point.

The receiver is an ideal calibrated analytic acceleration channel with a
rectangular bandpass 50.3 +/- 0.78125 Hz. It discards transition samples and
compares noncoherent complex-bin energies over the 6.4 s dwell. Tone spacing
times dwell length is one. Symbol timing and carrier frequency are assumed
already acquired; instantaneous transmitter phase is NOT supplied to the
decoder. Oscillator calibration and rejection of unknown startup motion are
idealized; raw resonator displacement is not treated as an instantaneous
acceleration measurement.

The full achieved-phase waveform is bandpass filtered before bin integration,
including transition leakage. Bin noise covariance is integrated from the
frequency-dependent PSD: K_kl=2 integral S_a(fc+nu) W_k(nu) conj(W_l(nu)) dnu.
S_a is the real one-sided PSD; the factor two supplies the complex analytic
baseband PSD. W is the normalized dwell-window response. Noise bins can be
correlated and have different variances; no independent equal-noise-bin
approximation or resonance-only sensitivity is used for the reported result.

The ideal rectangular filter is a reproducible offline reference. A physical
causal filter, acquisition scheme, oscillator-state estimator and their delay
and calibration errors remain implementation work.

## Selection and checks

The discrete symbol-period screen used 4, 6 and 8 s, with 5000 noise draws for
each of 49 cyclic transition contexts. The first two cases gave worst empirical
symbol-error fractions 0.16 and 0.005. Eight seconds was the first screened
case with no observed errors. This is not a continuous rate optimization.

The selected point was then simulated with the constrained drive and a separate
random seed: 100000 independent noise draws per fixed context, 4.9 million
decisions total. There were five errors. Mean empirical symbol-error fraction
is 1.02e-6; the maximum per-context observed fraction is 1e-5. A Bonferroni
correction across the 49 one-sided exact binomial confidence intervals gives
a simultaneous 95% upper limit of 9.21e-5 for the largest tested-context error
probability. This confidence limit is conditional on the numerical model.

The trials are repeated decisions on fixed signal contexts, not millions of
successive symbols of a physical packet. Temporal packet-noise correlations
are not simulated, and the contexts do not exhaust arbitrary longer histories
of the ideal bandpass. Accordingly, these results are not a universal packet
BER bound. A separate padded 256-byte random packet passes noiseless modulation,
bandpass filtering, and decoding.

Additional independent flat acceleration-noise ASD of 1e-10 m/s^2/sqrt(Hz)
gives mean SER 1.43e-5; 2e-10 gives 2.14e-3; 5e-10 gives 0.288 in separate
980000-trial tests. These are sensitivity calculations, not site measurements.
Coherent environmental lines and mechanical coupling from the source require
separate measurements and controls.

Checks include all 256 byte mappings, all 49 ordered tone transitions, exact
quintic phase/impulse/torque-squared integrals, detector normalization against
the exact white-noise noncoherent orthogonal M-FSK formula, and a receiver
sampling refinement from 100 to 200 samples/s (largest bin change 1.27e-5 of
source amplitude). The actuator work-energy residual is about 5e-8 J and its
largest sampled frequency error is 1.99e-5 Hz with nominal lag compensation.
These small tracking errors describe the assumed mathematical actuator.

The continuous-torque estimate 0.155 N m applies to the test pattern, which
visits every transition once. A real byte stream need not have uniform symbol
statistics; the 0.315 N m alternating-extremes estimate supplies a separate
duty check. Both fit the illustrative 0.5 N m continuous envelope.

## Reproduction and files

From the submission directory under Ubuntu WSL:

```sh
OPENBLAS_NUM_THREADS=2 python3 -s scripts/desktop_cf_fsk.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=2 python3 -s -m unittest discover -s tests -v
```

- `summary.json`: all parameters, drive diagnostics and confidence calculation.
- `screen.csv`, `technical_noise_sensitivity.csv`: selection and sensitivity.
- `transition_errors.csv`: errors and sample counts by transition.
- `trace.csv`: target/achieved frequency, torque and power at 10 samples/s.
- `detector_inputs.npz`: decision means, noise covariance and achieved phase.
- `operating_point.pdf` / `.png`: column-sized manuscript figure.
- `grav_comm_before_cf_fsk.tex`: source snapshot before this addition.

Manuscript text is in `../../sections/desktop_cf_fsk.tex`, included by the main
paper. No actual safe speed has been established: fixed-bore, bonded-insert
FEA omits real shaft/attachment/retention, fatigue and rotor-bearing dynamics.
The present data rate and actuator demand are a numerical design target.

Receiver physics reference: Carter et al., Scientific Reports 14, 17775 (2024),
https://doi.org/10.1038/s41598-024-68623-0 (structural damping and oscillator
parameters; the flat readout approximation is an explicit modeling assumption).
