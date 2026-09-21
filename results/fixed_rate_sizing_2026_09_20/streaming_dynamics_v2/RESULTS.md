# Modulation dynamics and amplitude refinement

Completed a bounded refinement with fixed instrumental noise, fixed 2 Hz full -3 dB bandwidth, tunable center frequency and the declared environmental background. Nominal offered payload rate is 1 bit/s; packets carry 224 payload bits within 256 information bits. CRC overhead, 12 assumed acquisition symbols and six simulated tail symbols are included in the 224 s frame budget. Acquisition itself remains unvalidated.

Execution: 48 waveform convergence checks, 576 initial pilot packets, 104 lower-bracket pilot packets, 12 candidate-amplitude convergence checks and 384 fresh confirmation packets. All source hashes and individual outcomes are preserved.

## Conditional three-distance candidates

The following candidates pass only the empirical screen (observed payload BER <=1e-3, latency <300 s, sampled detector convergence, clearance and nominal stress similarity). They are not statistically BER-qualified or hardware-qualified. Ranking minimizes mass, breaking ties by peak torque; the minimum-torque alternatives are also retained in selected.json.

| Distance (m) | Diameter (m) | Mass (kg) | Carrier (Hz) | Code | Transition fraction | q | Amplitude / reference | Errors / 7168 | Peak torque (N m) | Peak mechanical power (kW) |
|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.487 | 7.33 | 24 | uncoded | 0.5 | 0.1 | 0.875 | 6 | 1.02 | 0.077 |
| 1 | 0.853 | 39.35 | 24 | uncoded | 0.5 | 0.1 | 0.875 | 6 | 16.76 | 1.264 |
| 2 | 1.492 | 210.26 | 18 | tb_3/4 | 0.5 | 0.125 | 0.875 | 0 | 569.07 | 32.191 |

## Mass versus torque

At 0.5 and 1 m, rate-3/4 coding at r=0.6, q=0.08 and amplitude 0.875 uses the same rotor size as the tabulated uncoded case and produced zero errors/erasures in 32 fresh packets. Peak torque is 1.41/23.23 N m versus 1.02/16.76 N m uncoded. The uncoded sample is close to the BER cutoff, so its torque advantage remains provisional.

At 2 m, the mass-minimum rate-3/4 candidate uses 210.26 kg, 569.07 N m and 32.19 kW peak mechanical power. The lower-torque uncoded candidate uses 243.48 kg, 291.37 N m and 16.48 kW, with six errors and two erasures in 32 packets. Rate-2/3 at amplitude 1 gives 227.32 kg, 524.73 N m and zero observed errors. No candidate dominates all objectives.

The previous smallest sampled masses were 7.91/42.51/258.89 kg; this refinement gives 7.33/39.35/210.26 kg. At 2 m, that mass reduction increases peak torque relative to the earlier uncoded case (403.44 N m). These are mechanical comparisons of empirically screened candidates, not statistically validated reliability improvements.

At 24 Hz with r=0.6 and q=0.08, fresh trials reject amplitude 0.75 and retain 0.875: rate-3/4 has 38 versus zero payload errors and rate-2/3 has ten versus zero, each over 7168 bits. This is an empirical threshold bracket, not a confidence-certified BER crossing.


Peak power is absolute torque times angular speed, not average electrical consumption. Values here exclude drag and assume commanded motion is achieved. Larger drives, regenerative braking, bearings, modes and mechanical supports have not been designed. Geometric scaling is uniform, with mass proportional to scale cubed and inertia to scale to the fifth power.

## Fresh confirmation outcomes

| Configuration | Payload errors | Erasures / 32 | Wrong accepted | BER | Packetwise 95% upper bound |
|:---|---:|---:|---:|---:|---:|
| tb_2_3_f18_r0.6_q0.08_a0.875 | 63 | 2 | 0 | 0.0087891 | 0.162 |
| tb_2_3_f18_r0.6_q0.08_a1 | 0 | 0 | 0 | 0 | 0.154 |
| tb_2_3_f24_r0.6_q0.08_a0.75 | 10 | 1 | 0 | 0.0013951 | 0.155 |
| tb_2_3_f24_r0.6_q0.08_a0.875 | 0 | 0 | 0 | 0 | 0.154 |
| tb_3_4_f18_r0.5_q0.0625_a1.25 | 195 | 11 | 0 | 0.027204 | 0.178 |
| tb_3_4_f18_r0.5_q0.125_a0.875 | 0 | 0 | 0 | 0 | 0.154 |
| tb_3_4_f18_r0.5_q0.125_a1 | 0 | 0 | 0 | 0 | 0.154 |
| tb_3_4_f24_r0.6_q0.08_a0.75 | 38 | 4 | 0 | 0.0053013 | 0.158 |
| tb_3_4_f24_r0.6_q0.08_a0.875 | 0 | 0 | 0 | 0 | 0.154 |
| uncoded_f18_r0.6_q0.08_a1.125 | 6 | 2 | 0 | 0.00083705 | 0.154 |
| uncoded_f24_r0.5_q0.1_a0.875 | 6 | 2 | 0 | 0.00083705 | 0.154 |
| uncoded_f24_r0.6_q0.05_a1 | 71 | 13 | 0 | 0.0099051 | 0.162 |

All tentative payload errors count, including CRC-rejected packets. Packetwise intervals allow within-packet error clustering; they are per-configuration bounds, not simultaneous confidence after adaptive selection. Zero errors in 32 packets is insufficient for a BER <=1e-3 claim. CRC erasures reduce delivered goodput below the nominal offered rate; no retransmission protocol has been optimized.

## Amplitude brackets

The pilot uses four shared-seed packets at each grid amplitude: 0.75/0.875/1 at 24 Hz, and 1/1.125/1.25 at 18 Hz. Where the lowest grid point passes, four-packet probes at two lower amplitudes extend the bracket (0.5/0.625 at 24 Hz and 0.75/0.875 at 18 Hz). thresholds.json records every waveform/code/carrier group, including failures and unbracketed boundaries. These are empirical brackets, not measured BER-target crossing confidence intervals. Candidate-specific convergence and fresh confirmation can overturn a pilot bracket.

## Why these waveform variables matter

Let r be transition duration / symbol duration and q be adjacent-tone spacing times dwell duration. Tone spacing is q/[(1-r)T]. For the quintic smooth transition, peak angular acceleration scales as q/[r(1-r)T^2]. Changing r at fixed q changes both spacing and acceleration, so matched-spacing cases are included. The grid spans spacing*T of 1/8,1/6,1/5,1/4 and r from 1/4 to 3/5. Wider separation can improve sequence discrimination but increases excursion and off-center receiver noise; gentler transitions reduce acceleration but change waveform distinguishability. A lower amplitude threshold reduces rotor inertia as well as mass.

## Receiver and orientation checks

All seven non-default waveforms agree with exhaustive one-byte soft likelihoods in the added regression test. Existing four streaming tests plus 15 causal and five environmental tests also pass (25 tests total). Initial waveform gates compare history 4 versus 5 and lookahead 12 versus 24 with relative LLR change <1%, identical signs/decoded information, no missing hypotheses and no pruning. State caps are 32768/200000; phase lattices are exact rational, not rounded.

Maximum main-search block time: 0.096 s. Maximum FEC time: 0.002 s. Maximum budgeted latency: 224.021 s. Timing is measured under parallel CPU load, not a deployment guarantee.

The user-requested finite-volume orientation audit found no larger second-harmonic amplitude than the assumed equatorial radial arrangement among eight elevations at three representative scaled separations. At the desktop equivalent distance, 15/30-degree tilts retain 90.4%/68.3% of the best signal; the axial signal vanishes. The optimizing real sensor axis was allowed to change at every position. See ORIENTATION.md, orientation_audit.json and orientation.pdf. Directional environmental noise could change the SNR optimum and is outside the scalar noise scenario.

## Scope and reproduction

This is a bounded waveform/threshold refinement, not a global optimum over packet length, code family, geometry or receiver hardware. The 224-bit payload, two carriers, finite waveform grid and scalar environmental spectrum remain restrictions. In particular, more severe low-spacing cases may require amplitudes above the searched range. The previously screened LDPC and rate-1/2 convolutional codes are not reoptimized here.

Run streaming_dynamics_search.py --stage gates --workers 8, then --stage screen --workers 12. Run extend_streaming_dynamics.py to probe below passing lower boundaries. Build amplitude_size_map.json from streaming_size_map.one for the grid (the previous maps are reused). Run compile_streaming_dynamics.py --select to write shortlist.json; then run search stages confirm_gates and confirm, and compile without --select. Run audit_streaming_dynamics.py and report_streaming_dynamics.py. Remote CPU runs use the saved runtime and scripts; all job manifests include seeds/configurations/source hashes. Old study results and manuscript numerical claims are preserved.
