# Parallel qualification of larger 2 m designs

Three larger replacement candidates were frozen before this run; the failed 218.91 kg design remains unchanged. All qualification packets are independent of the earlier screening samples. The streaming receiver, code, packet length, carrier, amplitude, instrumental profile and declared environmental spectrum were held fixed. These conclusions concern the simulated, acquired-packet link; they do not qualify physical transmitter motion, acquisition, site noise or hardware.

All three candidates are at 2 m and preserve the 18 Hz carrier, rate-2/3 code, r=.6, q=.08, fixed noise/bandwidth and 224 s frame. Only source size and received amplitude change. Confidence applies simultaneously to this new family of three candidates; it is not a combined confidence statement across historical studies.

## Decision rule

The error fraction in each 224-bit tentative payload is one independent bounded observation. Errors in CRC-rejected frames count. Two-sided time-uniform confidence sequences use alpha=0.05/3 per candidate, giving at least 95% simultaneous confidence across all three candidates by a union bound. Pass requires the upper bound <=1e-3; fail requires the lower bound >1e-3. Decisions are tested only after deterministic 128-packet prefixes. Bit independence is not assumed. No screening data or diagnostic reruns enter the sample counts.

| Distance (m) | Mass (kg) | Carrier (Hz) | Amplitude / reference | Packets | Payload errors / bits | BER | Simultaneous-confidence interval | Decision |
|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| 2 | 227.32 | 18 | 1 | 10624 | 833 / 2379776 | 0.000350033 | [0.000168351, 0.000992984] | PASS |
| 2 | 235.50 | 18 | 1.0625 | 7040 | 92 / 1576960 | 5.83401e-05 | [6.50239e-06, 0.00099634] | PASS |
| 2 | 243.48 | 18 | 1.125 | 6656 | 10 / 1490944 | 6.70716e-06 | [4.77416e-11, 0.000996526] | PASS |

## Packet acceptance and timing

| Channel | CRC/header erasures | Packet failures | Wrong accepted payloads | Empirical accepted payload rate (bit/s) | Max block time (s) | Max budgeted latency (s) |
|:---|---:|---:|---:|---:|---:|---:|
| tb_2_3_f18_r0.6_q0.08_a1 | 58 | 58 | 0 | 0.994541 | 0.329 | 224.038 |
| tb_2_3_f18_r0.6_q0.08_a1.0625 | 9 | 9 | 0 | 0.998722 | 0.302 | 224.042 |
| tb_2_3_f18_r0.6_q0.08_a1.125 | 1 | 1 | 0 | 0.999850 | 0.305 | 224.041 |

Accepted rate is nominal offered payload rate (1 bit/s) times observed acceptance fraction, without retransmissions or extra inter-frame processing overhead. It is not a qualified throughput guarantee. A BER pass does not imply 1 bit/s accepted goodput, acquisition success or a bounded undetected-error probability. Zero observed wrong accepted packets is not a zero-probability claim.

Execution uses the implemented bounded streaming receiver, with packet processing concurrent with reception. Timing under many-worker CPU load is diagnostic rather than worst-case real-time certification. The modeled frame budget is 224 s, including assumed acquisition overhead and simulated tail; known timing/phase and calibrated input are prerequisites.

## Scope and reproduction

All individual packet records, frozen design parameters, seed starts, source hashes and final bounds are retained. Run scripts/qualify_replacement_ber.py --workers 64 --batch 128 on the same runtime; completed deterministic prefixes resume without reusing earlier screening samples. Run scripts/compile_replacement_ber.py to verify source hashes, packet identity/quality, first stopping checkpoint and final decisions. The protocol is PLAN.md. Frozen failures must be redesigned and independently requalified; their amplitude or code was not changed during this run.

## Error-packet convergence checks

Replayed the first three error-containing packets per candidate (or all if fewer): 7 checks passed with wider history and longer lookahead. These reruns do not add BER samples.

## Physical comparison

| Amplitude / reference | Diameter (m) | Mass (kg) | Peak torque (N m) | Peak mechanical power (kW) | Stress-similarity ratio |
|---:|---:|---:|---:|---:|---:|
| 1 | 1.531 | 227.32 | 524.73 | 29.68 | 0.9409 |
| 1.0625 | 1.549 | 235.50 | 556.59 | 31.48 | 0.9521 |
| 1.125 | 1.567 | 243.48 | 588.39 | 33.28 | 0.9627 |

Mechanical estimates use uniform finite-volume scaling and achieved commanded motion with zero drag. The stress-similarity ratio is the existing screen, not a new yield-stress calculation; peak mechanical power is not average electrical power. New structure and actuator validation remain necessary.

Smallest passing tested candidate: 227.32 kg, diameter 1.531 m. This is not a continuous optimum. Delivered-rate compensation still requires changed timing and independent testing.
