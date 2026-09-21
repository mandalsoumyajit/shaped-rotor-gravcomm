# Fine amplitude thresholds

This pass repeats endpoints and adds intermediate amplitudes using 64 fresh packets per gate-approved configuration. No earlier adaptive pilot samples are pooled into the fresh statistics. Fixed instrumental floor, 2 Hz bandwidth, declared environmental noise and the 224-bit payload/224 s nominal frame budget are unchanged. Errors in rejected frames count toward BER.

## Fresh error measurements

| Configuration | Errors / payload bits | Erasures | Wrong accepted | BER | Packetwise 95% upper bound | Accepted payload rate (bit/s) |
|:---|---:|---:|---:|---:|---:|---:|
| tb_2_3_f18_r0.6_q0.08_a0.9375 | 3 / 14336 | 2 / 64 | 0 | 0.00020926 | 0.0805 | 0.9688 |
| tb_2_3_f18_r0.6_q0.08_a1 | 3 / 14336 | 2 / 64 | 0 | 0.00020926 | 0.0805 | 0.9688 |
| tb_2_3_f24_r0.6_q0.08_a0.8125 | 3 / 14336 | 2 / 64 | 0 | 0.00020926 | 0.0805 | 0.9688 |
| tb_2_3_f24_r0.6_q0.08_a0.875 | 0 / 14336 | 0 / 64 | 0 | 0 | 0.0803 | 1.0000 |
| tb_3_4_f18_r0.5_q0.125_a0.8125 | 114 / 14336 | 4 / 64 | 0 | 0.007952 | 0.0882 | 0.9375 |
| tb_3_4_f18_r0.5_q0.125_a0.875 | 32 / 14336 | 2 / 64 | 0 | 0.0022321 | 0.0825 | 0.9688 |
| tb_3_4_f24_r0.6_q0.08_a0.8125 | 27 / 14336 | 1 / 64 | 0 | 0.0018834 | 0.0822 | 0.9844 |
| tb_3_4_f24_r0.6_q0.08_a0.875 | 0 / 14336 | 0 / 64 | 0 | 0 | 0.0803 | 1.0000 |
| uncoded_f18_r0.6_q0.08_a1.125 | 20 / 14336 | 4 / 64 | 0 | 0.0013951 | 0.0816 | 0.9375 |
| uncoded_f18_r0.6_q0.08_a1.1875 | 11 / 14336 | 2 / 64 | 0 | 0.0007673 | 0.081 | 0.9688 |
| uncoded_f24_r0.5_q0.1_a0.875 | 22 / 14336 | 4 / 64 | 0 | 0.0015346 | 0.0818 | 0.9375 |
| uncoded_f24_r0.5_q0.1_a0.9375 | 11 / 14336 | 2 / 64 | 0 | 0.0007673 | 0.081 | 0.9688 |

The nominal offered payload rate is 1 bit/s. Accepted rate above reflects CRC erasures with no retransmissions or additional inter-frame processing overhead; it is an empirical sample quantity, not qualified goodput. The confidence bounds allow within-packet correlation and are per fixed configuration, not simultaneous after selection. They do not establish BER <=1e-3.

## Mass and torque alternatives

Each row passes the empirical BER, latency, clearance, receiver convergence and centrifugal stress-similarity screens. This is not drive, structure or statistical reliability qualification.

| Distance (m) | Objective | Configuration | Diameter (m) | Mass (kg) | Peak torque (N m) | RMS torque (N m) | Peak mechanical power (kW) |
|---:|:---|:---|---:|---:|---:|---:|---:|
| 0.5 | mass | tb_2_3_f24_r0.6_q0.08_a0.8125 | 0.481 | 7.03 | 1.60 | 0.75 | 0.120 |
| 0.5 | torque | uncoded_f24_r0.5_q0.1_a0.9375 | 0.494 | 7.63 | 1.09 | 0.46 | 0.082 |
| 1 | mass | tb_2_3_f24_r0.6_q0.08_a0.8125 | 0.841 | 37.69 | 26.26 | 12.29 | 1.980 |
| 1 | torque | uncoded_f24_r0.5_q0.1_a0.9375 | 0.865 | 40.95 | 17.92 | 7.65 | 1.351 |
| 2 | mass | tb_2_3_f18_r0.6_q0.08_a0.9375 | 1.512 | 218.91 | 492.78 | 230.60 | 27.873 |
| 2 | torque | uncoded_f18_r0.6_q0.08_a1.1875 | 1.583 | 251.28 | 307.08 | 143.70 | 17.367 |

Maximum measured streaming block time: 0.072 s. Maximum budgeted latency: 224.015 s. Timings are diagnostic measurements under parallel host load.


Mechanical values exclude drag and assume achieved commanded motion. Peak absolute mechanical power is not average electrical draw. Uniform finite-volume scaling and nominal stress similarity do not substitute for new FEA, modal/shaft analysis or drive selection.

## Interpretation and qualification

Observed errors at the lower amplitude and zero/few errors at the higher amplitude provide an empirical bracket only. A rare decoded error burst can change the apparent threshold between 32- and 64-packet runs. Choose a reliability margin before freezing a design, then use disjoint packets for qualification. The present short trials are intended to reject clearly poor candidates and compare physical costs.

With the existing time-uniform packet-fraction bound, 5374 zero-error independent packets are needed for a per-configuration 95% upper BER bound <=1e-3. Any errors or simultaneous coverage requirements change that budget; see qualification_budget.json.

The manuscript now explains why the source uses an equatorial receiver with radial sensitivity: it maximizes the leading rotating-quadrupole signal at fixed distance, and sampled finite-volume checks support that choice. Directional site noise could change the SNR-optimal orientation. See ../streaming_dynamics_v2/ORIENTATION.md. The revised manuscript builds successfully; these unqualified communications results have not been inserted.

## Reproduction

Run scripts/streaming_thresholds.py --stage gates --workers 8, then --stage trials --workers 12. Source hashes, full configuration lists, separate payload/noise seed construction and seed ranges are saved in the manifests and runtime. Amplitude maps reuse earlier finite-volume results and add streaming_size_map.one at .8125/.9375/1.1875 for all distances. Run scripts/compile_streaming_thresholds.py to verify every record and source hash and generate this report.


## Targeted error-packet checks

Replayed 6 error-containing packets of the mass/torque candidates with history 5 and longer lookahead. 6 passed the convergence criterion. These are reruns of existing records and add no independent BER samples. Raw comparisons and source hashes are saved in error_checks/.
