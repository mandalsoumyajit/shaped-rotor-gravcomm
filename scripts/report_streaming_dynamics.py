"""Report conditional waveform refinement and empirical amplitude brackets."""
import json,hashlib
from collections import defaultdict
from pathlib import Path
from streaming_dynamics_runtime import OUT,ROOT,save
rows=json.loads((OUT/'compiled.json').read_text());fresh=[x for x in rows if x['confirmation'] and x['packets']==32 and x['candidate_gate']]
selected=[];torque_best=[]
for d in (.5,1.,2.):
 choices=[x for x in fresh if x['distance_m']==d and x['eligible']]
 if choices:
  selected.append(min(choices,key=lambda x:(x['mass_kg'],x['peak_torque_nm'])))
  torque_best.append(min(choices,key=lambda x:(x['peak_torque_nm'],x['mass_kg'])))
save(OUT/'selected.json',dict(minimum_mass=selected,minimum_peak_torque=torque_best))
counts={stage:len(list((OUT/stage).glob('*.json'))) for stage in ('gates','screen','extra','confirm_gates','confirm')}
records=[json.loads(p.read_text()) for stage in ('screen','extra','confirm') for p in (OUT/stage).glob('*.json')]
lines=['# Modulation dynamics and amplitude refinement','',
'Completed a bounded refinement with fixed instrumental noise, fixed 2 Hz full -3 dB bandwidth, tunable center frequency and the declared environmental background. Nominal offered payload rate is 1 bit/s; packets carry 224 payload bits within 256 information bits. CRC overhead, 12 assumed acquisition symbols and six simulated tail symbols are included in the 224 s frame budget. Acquisition itself remains unvalidated.','',
f'Execution: {counts["gates"]} waveform convergence checks, {counts["screen"]} initial pilot packets, {counts["extra"]} lower-bracket pilot packets, {counts["confirm_gates"]} candidate-amplitude convergence checks and {counts["confirm"]} fresh confirmation packets. All source hashes and individual outcomes are preserved.','',
'## Conditional three-distance candidates','',
'The following candidates pass only the empirical screen (observed payload BER <=1e-3, latency <300 s, sampled detector convergence, clearance and nominal stress similarity). They are not statistically BER-qualified or hardware-qualified. Ranking minimizes mass, breaking ties by peak torque; the minimum-torque alternatives are also retained in selected.json.','',
'| Distance (m) | Diameter (m) | Mass (kg) | Carrier (Hz) | Code | Transition fraction | q | Amplitude / reference | Errors / 7168 | Peak torque (N m) | Peak mechanical power (kW) |',
'|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|---:|']
for x in selected:lines.append(f"| {x['distance_m']:g} | {x['diameter_m']:.3f} | {x['mass_kg']:.2f} | {x['carrier_hz']:g} | {x['scheme']} | {x['transition_fraction']:g} | {x['q']:g} | {x['amplitude_scale']:g} | {x['errors']} | {x['peak_torque_nm']:.2f} | {x['peak_power_w']/1000:.3f} |")
lines+=['','## Mass versus torque','',
'At 0.5 and 1 m, rate-3/4 coding at r=0.6, q=0.08 and amplitude 0.875 uses the same rotor size as the tabulated uncoded case and produced zero errors/erasures in 32 fresh packets. Peak torque is 1.41/23.23 N m versus 1.02/16.76 N m uncoded. The uncoded sample is close to the BER cutoff, so its torque advantage remains provisional.', '',
'At 2 m, the mass-minimum rate-3/4 candidate uses 210.26 kg, 569.07 N m and 32.19 kW peak mechanical power. The lower-torque uncoded candidate uses 243.48 kg, 291.37 N m and 16.48 kW, with six errors and two erasures in 32 packets. Rate-2/3 at amplitude 1 gives 227.32 kg, 524.73 N m and zero observed errors. No candidate dominates all objectives.', '',
'The previous smallest sampled masses were 7.91/42.51/258.89 kg; this refinement gives 7.33/39.35/210.26 kg. At 2 m, that mass reduction increases peak torque relative to the earlier uncoded case (403.44 N m). These are mechanical comparisons of empirically screened candidates, not statistically validated reliability improvements.', '',
'At 24 Hz with r=0.6 and q=0.08, fresh trials reject amplitude 0.75 and retain 0.875: rate-3/4 has 38 versus zero payload errors and rate-2/3 has ten versus zero, each over 7168 bits. This is an empirical threshold bracket, not a confidence-certified BER crossing.', '']
lines+=['','Peak power is absolute torque times angular speed, not average electrical consumption. Values here exclude drag and assume commanded motion is achieved. Larger drives, regenerative braking, bearings, modes and mechanical supports have not been designed. Geometric scaling is uniform, with mass proportional to scale cubed and inertia to scale to the fifth power.','',
'## Fresh confirmation outcomes','',
'| Configuration | Payload errors | Erasures / 32 | Wrong accepted | BER | Packetwise 95% upper bound |',
'|:---|---:|---:|---:|---:|---:|']
for x in sorted([x for x in fresh if x['distance_m']==.5],key=lambda x:x['id']):lines.append(f"| {x['id']} | {x['errors']} | {x['erasures']} | {x['wrong_accepted']} | {x['ber']:.5g} | {x['ber_ci'][1]:.3g} |")
lines+=['','All tentative payload errors count, including CRC-rejected packets. Packetwise intervals allow within-packet error clustering; they are per-configuration bounds, not simultaneous confidence after adaptive selection. Zero errors in 32 packets is insufficient for a BER <=1e-3 claim. CRC erasures reduce delivered goodput below the nominal offered rate; no retransmission protocol has been optimized.','',
'## Amplitude brackets','',
'The pilot uses four shared-seed packets at each grid amplitude: 0.75/0.875/1 at 24 Hz, and 1/1.125/1.25 at 18 Hz. Where the lowest grid point passes, four-packet probes at two lower amplitudes extend the bracket (0.5/0.625 at 24 Hz and 0.75/0.875 at 18 Hz). thresholds.json records every waveform/code/carrier group, including failures and unbracketed boundaries. These are empirical brackets, not measured BER-target crossing confidence intervals. Candidate-specific convergence and fresh confirmation can overturn a pilot bracket.','']
groups=defaultdict(list)
for p in [p for stage in ('screen','extra') for p in (OUT/stage).glob('*.json')]:
 x=json.loads(p.read_text());c=x['config'];groups[(c['scheme'],c['carrier_hz'],c['transition_fraction'],c['q'],c['amplitude_scale'])].append(x)
waveforms=defaultdict(list)
for (scheme,fc,r,q,a),rs in groups.items():waveforms[(scheme,fc,r,q)].append(dict(amplitude_scale=a,ber=sum(x['errors'] for x in rs)/(224*len(rs)),errors=sum(x['errors'] for x in rs),packets=len(rs)))
thresholds=[]
for (scheme,fc,r,q),points in sorted(waveforms.items()):
 points=sorted(points,key=lambda x:x['amplitude_scale']);passing=[x['amplitude_scale'] for x in points if x['ber']<=.001];high=min(passing) if passing else None;lower=[x['amplitude_scale'] for x in points if x['ber']>.001 and (high is None or x['amplitude_scale']<high)]
 thresholds.append(dict(scheme=scheme,carrier_hz=fc,transition_fraction=r,q=q,points=points,empirical_lower_failing_amplitude=max(lower) if lower else None,empirical_upper_passing_amplitude=high,status='bracketed_in_pilot' if lower and high else 'above_grid' if high is None else 'at_or_below_grid'))
save(OUT/'thresholds.json',thresholds)
lines+=['## Why these waveform variables matter','',
'Let r be transition duration / symbol duration and q be adjacent-tone spacing times dwell duration. Tone spacing is q/[(1-r)T]. For the quintic smooth transition, peak angular acceleration scales as q/[r(1-r)T^2]. Changing r at fixed q changes both spacing and acceleration, so matched-spacing cases are included. The grid spans spacing*T of 1/8,1/6,1/5,1/4 and r from 1/4 to 3/5. Wider separation can improve sequence discrimination but increases excursion and off-center receiver noise; gentler transitions reduce acceleration but change waveform distinguishability. A lower amplitude threshold reduces rotor inertia as well as mass.','',
'## Receiver and orientation checks','',
'All seven non-default waveforms agree with exhaustive one-byte soft likelihoods in the added regression test. Existing four streaming tests plus 15 causal and five environmental tests also pass (25 tests total). Initial waveform gates compare history 4 versus 5 and lookahead 12 versus 24 with relative LLR change <1%, identical signs/decoded information, no missing hypotheses and no pruning. State caps are 32768/200000; phase lattices are exact rational, not rounded.', '',
f'Maximum main-search block time: {max(x["max_block_s"] for x in records):.3f} s. Maximum FEC time: {max(x["fec_s"] for x in records):.3f} s. Maximum budgeted latency: {max(x["latency_budget_s"] for x in records):.3f} s. Timing is measured under parallel CPU load, not a deployment guarantee.', '',
'The user-requested finite-volume orientation audit found no larger second-harmonic amplitude than the assumed equatorial radial arrangement among eight elevations at three representative scaled separations. At the desktop equivalent distance, 15/30-degree tilts retain 90.4%/68.3% of the best signal; the axial signal vanishes. The optimizing real sensor axis was allowed to change at every position. See ORIENTATION.md, orientation_audit.json and orientation.pdf. Directional environmental noise could change the SNR optimum and is outside the scalar noise scenario.', '',
'## Scope and reproduction','',
'This is a bounded waveform/threshold refinement, not a global optimum over packet length, code family, geometry or receiver hardware. The 224-bit payload, two carriers, finite waveform grid and scalar environmental spectrum remain restrictions. In particular, more severe low-spacing cases may require amplitudes above the searched range. The previously screened LDPC and rate-1/2 convolutional codes are not reoptimized here.', '',
'Run streaming_dynamics_search.py --stage gates --workers 8, then --stage screen --workers 12. Run extend_streaming_dynamics.py to probe below passing lower boundaries. Build amplitude_size_map.json from streaming_size_map.one for the grid (the previous maps are reused). Run compile_streaming_dynamics.py --select to write shortlist.json; then run search stages confirm_gates and confirm, and compile without --select. Run audit_streaming_dynamics.py and report_streaming_dynamics.py. Remote CPU runs use the saved runtime and scripts; all job manifests include seeds/configurations/source hashes. Old study results and manuscript numerical claims are preserved.','']
(OUT/'RESULTS.md').write_text('\n'.join(lines))
save(OUT/'report_manifest.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/compile_streaming_dynamics.py',ROOT/'scripts/audit_rotor_orientation.py',ROOT/'tests/test_streaming_dynamics.py',ROOT/'fea/results/shaped_waveforms/medium_rounded.npz']})
print(json.dumps(dict(selected=selected,minimum_torque=torque_best,counts=counts),indent=2))
