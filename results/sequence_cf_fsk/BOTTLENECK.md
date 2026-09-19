# Rate-limiting constraints

The current design has a coupled continuous-torque/noise constraint. Peak motor
torque, peak mechanical power and centrifugal speed have headroom at the tested
points. The 1.333 bit/s point is a verified operating point; its selection used
a conservative duty screen, and the ultimate rate remains to be refined.

## Controlled comparison at the same 0.5 m range and receiver noise

All cases use seven tones, the same rotor, a 40% transition fraction, full-waveform
whitened Viterbi decoding and the same byte mapping. Rates below exclude trailers
and acquisition. Torque columns refer to long streams of valid encoded bytes.

| Symbol period | Tone spacing | Mapped bit/s | BER result | Worst-message RMS torque | Random-byte RMS torque |
|---:|---:|---:|---|---:|---:|
| 2 s | 0.08333 Hz | 1.333 | Measured 3.3e-5; 134 errors | 0.447 N m | 0.221 N m |
| 1.75 s | 0.07143 Hz | 1.524 | Measured 0.00377; 181 errors | 0.438 N m | 0.216 N m |
| 1.75 s | 0.09524 Hz | 1.524 | 95% upper bound 0.000873 | 0.583 N m | 0.285 N m |

The final row has 112 bit errors in 528,000 bits from
8,250 independent packets. Its packet-aware interval is
[8.62e-05, 0.000873]. This row is used to test the BER threshold;
it is not promoted to a precise BER estimate unless the full error-count
requirement is met.

Widening tones at 1.75 s allows the modeled receiver to meet BER < 0.001.
It raises guaranteed worst-message RMS torque above the assumed 0.5 N m
continuous rating. The narrower tones fit that rating but miss the BER target.
This isolates the interaction between receiver discrimination and motor duty.

## Which transmitter constraint matters

The wider-tone 1.75 s waveform needs 1.57 N m peak torque and 248 W peak
mechanical power, below the assumed 2 N m and 400 W limits. Its maximum speed
is 1517.57 rpm, versus the 1800 rpm design ceiling; centrifugal stress changes
very little relative to the 2 s point. Continuous torque is the relevant
transmitter constraint in this comparison. Combined-load stress, fatigue and
thermal time constants still require assembly-level validation.

## Duty policy changes the conclusion

The original screen bounded each individual symbol at 0.5 N m RMS. At the
selected 2 s point this gave 0.497 N m. Exact maximum-mean-cycle analysis over
valid byte sequences gives 0.447 N m for the worst indefinitely repeated data
pattern, bytes 251 and 42. The independent waveform calculation confirms this
value. Uniform random bytes average 0.221 N m.

At the wider-tone 1.75 s point, uniform random bytes average only 0.285 N m.
Thus the continuous-torque restriction depends on whether the design guarantees
arbitrary sustained data patterns or accepts a specified statistical duty and
thermal protection policy. A 0.5 N m continuous rating does not impose a hard
per-symbol torque ceiling. The next optimization should use this exact
valid-message duty calculation or an explicit motor thermal model.

The receiver contribution is noise-driven sequence discrimination. These tones
remain inside the processed band; the 2 Hz Gaussian benchmark is approximately
3.9045 bit/s. Both improved receiver sensitivity and increased allowable motor
duty can move the feasible frontier. Assigning one universal hardware bottleneck
would go beyond the evaluated design and duty assumptions.

## Reproduction

```
python scripts/exact_sequence_duty.py
python scripts/measure_sequence_ber.py --symbol-s 1.75 --transition-fraction 0.4 --spacing-product 0.1 --output-name torque_relaxed_ber --max-packets 10000 --workers 4
python scripts/summarize_bottleneck.py
```

The torque-relaxed test may be stopped once its anytime-valid upper bound is
below the target; its retained packet counts specify the actual stopping point.
