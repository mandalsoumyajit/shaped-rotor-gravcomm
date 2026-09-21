# Independent BER qualification of frozen designs

The minimum-mass candidates were frozen before this run. All qualification packets are independent of the earlier screening samples. The streaming receiver, code, packet length, carrier, amplitude, instrumental profile and declared environmental spectrum were held fixed. These conclusions concern the simulated, acquired-packet link; they do not qualify physical transmitter motion, acquisition, site noise or hardware.

Two distinct normalized channels cover the three distances. The 0.5 and 1 m transmitters were sized to the same received harmonic amplitude and use the same carrier/receiver/noise model. Their shared BER result is mapped to both source sizes; field-inversion discretization tolerances are retained in the frozen source records.

## Decision rule

The error fraction in each 224-bit tentative payload is one independent bounded observation. Errors in CRC-rejected frames count. Two-sided time-uniform confidence sequences use alpha=0.025 per channel, giving at least 95% simultaneous confidence across both channels by a union bound. Pass requires the upper bound <=1e-3; fail requires the lower bound >1e-3. Decisions are tested only after deterministic 128-packet prefixes. Bit independence is not assumed. No screening data or diagnostic reruns enter the sample counts.

| Distance (m) | Mass (kg) | Carrier (Hz) | Amplitude / reference | Packets | Payload errors / bits | BER | Simultaneous-confidence interval | Decision |
|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| 0.5 | 7.03 | 24 | 0.8125 | 10112 | 835 / 2265088 | 0.000368639 | [0.000205877, 0.000993232] | PASS |
| 1 | 37.69 | 24 | 0.8125 | 10112 | 835 / 2265088 | 0.000368639 | [0.000205877, 0.000993232] | PASS |
| 2 | 218.91 | 18 | 0.9375 | 2560 | 1106 / 573440 | 0.00192871 | [0.00107831, 0.00444024] | FAIL |

## Packet acceptance and timing

| Channel | CRC/header erasures | Packet failures | Wrong accepted payloads | Empirical accepted payload rate (bit/s) | Max block time (s) | Max budgeted latency (s) |
|:---|---:|---:|---:|---:|---:|---:|
| tb_2_3_f24_r0.6_q0.08_a0.8125 | 80 | 80 | 0 | 0.992089 | 0.109 | 224.019 |
| tb_2_3_f18_r0.6_q0.08_a0.9375 | 82 | 82 | 0 | 0.967969 | 0.081 | 224.019 |

Accepted rate is nominal offered payload rate (1 bit/s) times observed acceptance fraction, without retransmissions or extra inter-frame processing overhead. It is not a qualified throughput guarantee. A BER pass does not imply 1 bit/s accepted goodput, acquisition success or a bounded undetected-error probability. Zero observed wrong accepted packets is not a zero-probability claim.

Execution uses the implemented bounded streaming receiver, with packet processing concurrent with reception. Timing under many-worker CPU load is diagnostic rather than worst-case real-time certification. The modeled frame budget is 224 s, including assumed acquisition overhead and simulated tail; known timing/phase and calibrated input are prerequisites.

## Scope and reproduction

All individual packet records, frozen design parameters, seed starts, source hashes and final bounds are retained. Run scripts/qualify_streaming_ber.py --workers 40 --batch 128 on the same runtime; completed deterministic prefixes resume without reusing earlier screening samples. Run scripts/compile_ber_qualification.py to verify source hashes, packet identity/quality, first stopping checkpoint and final decisions. The protocol is PROTOCOL.md. Frozen failures must be redesigned and independently requalified; their amplitude or code was not changed during this run.
