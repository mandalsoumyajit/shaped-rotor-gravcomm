"""Generate a reviewable bounded-search report from saved packet records."""
import json,datetime,hashlib
from pathlib import Path
from streaming_link_search import OUT,ROOT,save
from range_revision import duty
rows=json.loads((OUT/'compiled_candidates.json').read_text())
fresh=[x for x in rows if x['statistical_sample']=='fresh_refinement' and x['packets']==32]
selected=[];pareto=[]
for d in (.5,1.,2.):
 feasible=[x for x in fresh if x['distance_m']==d and x['screen_eligible'] and x['candidate_convergence_pass']]
 selected+=sorted(feasible,key=lambda x:(x['mass_kg'],x['peak_torque_nm'],x['peak_mechanical_power_w']))[:1]
 keys=('mass_kg','peak_torque_nm','peak_mechanical_power_w')
 for x in feasible:
  if not any(all(y[k]<=x[k] for k in keys) and any(y[k]<x[k] for k in keys) for y in feasible):pareto.append(x)
save(OUT/'selected_screen_candidates.json',selected);save(OUT/'pareto_screen_candidates.json',pareto)
lines=['# Streaming receiver: three-distance link screen','',
'Completed a bounded nominal 1 payload bit/s search at 0.5, 1 and 2 m. Receiver instrumental noise versus offset and the full 2 Hz -3 dB bandwidth remain fixed; center frequency varies. These are empirical screening results under the declared environmental spectrum, not qualified BER-compliant hardware designs.','',
'## Search and validation','',
'- 10 initial protocol/length convergence checks; 8 passed. 256 pilot packets cover 128 configurations: carrier 12/16/24/40 Hz, amplitude 0.5/1/2/4 times the reference, 96/224 payload bits, uncoded and three convolutional rates plus short LDPC where the detector gate passed.',
'- Candidate-specific convergence checks add 18 Hz and amplitudes 0.75/1.25. The first shortlist had 12 checks (11 passed) and 352 fresh packets. A second four-candidate shortlist addresses the 2 m failures: all four gates passed and 128 fresh packets were completed. Total: 736 packet trials and 26 convergence records. Payload/noise seeds are separate; common random numbers pair configurations. Refinement seeds are disjoint from selection and gate seeds.',
'- History 4 versus 5 and lookahead 12 versus 24 must agree on hard decisions/signs and within 1% in relative LLR norm; no state pruning/missing hypotheses and whitening spectral error <0.1%. These are sampled convergence diagnostics, not an exact likelihood guarantee.',
'- All 24 causal receiver, streaming detector and environmental-noise regression tests pass.','',
'## Smallest rotor passing the empirical refinement screen','',
'Ranking first minimizes mass within the sampled uniform-scaling family; ties use lower peak inertial torque, then lower peak mechanical power. Other mass/torque/power tradeoffs are retained in pareto_screen_candidates.json. Eligibility here means observed BER <=1e-3, latency <300 s, convergence checks and clearance/stress-similarity screens; it does not mean the BER confidence bound meets the target.','',
'| Distance (m) | Diameter (m) | Mass (kg) | Carrier (Hz) | Protocol | Payload errors / bits | Peak torque (N m) | Peak mechanical power (kW) |',
'|---:|---:|---:|---:|:---|---:|---:|---:|']
for x in selected:lines.append(f"| {x['distance_m']:g} | {x['diameter_m']:.3f} | {x['mass_kg']:.2f} | {x['carrier_hz']:g} | {x['scheme']} | {x['payload_errors']} / {x['packets']*x['payload_bits']} | {x['peak_torque_nm']:.2f} | {x['peak_mechanical_power_w']/1000:.3f} |")
lines+=['','The frequency is the gravitational signal carrier; rotor speed is 30 times this value in rpm. Torque and power are worst-transition inertial requirements. Power is peak absolute torque-times-angular-speed, not average electrical consumption. Regeneration can recover braking energy but does not remove motor/inverter peak torque and power requirements.','',
'## Independent refinement results','',
'Each row uses only the 32 fresh packets (224 useful bits each); pilot and diagnostic gate packets are saved separately. BER counts errors in all tentative payloads, including CRC-rejected packets. Erasures and undetected errors are separate outputs. The intervals below are packet-cluster-aware 95% confidence sequences for each fixed configuration, not simultaneous coverage of the adaptively selected winner.','',
'| Protocol | Carrier (Hz) | Amplitude / reference | Errors | Erasures / 32 | Wrong accepted | BER | Packetwise BER upper bound | Max latency (s) |',
'|:---|---:|---:|---:|---:|---:|---:|---:|---:|']
for x in sorted([v for v in fresh if v['distance_m']==.5],key=lambda x:x['id']):
 lines.append(f"| {x['scheme']} | {x['carrier_hz']:g} | {x['amplitude_scale']:g} | {x['payload_errors']} | {x['erasures']} | {x['undetected_wrong_packets']} | {x['ber']:.5g} | {x['ber_interval95_packetwise'][1]:.3g} | {x['latency_budget_s']:.3f} |")
lines+=['','No candidate is statistically qualified at BER <=1e-3 by these short runs. Zero observed errors are not a zero upper bound. Freeze a candidate and validate it independently before making a reliability claim. CRC-rejected frames reduce delivered goodput below the nominal 1 payload bit/s offered rate; no retransmission budget is included.','',
'## Mechanics and interpretation','',
'Uniform scaling uses a_s(d)=s a_1(d/s), mass proportional to s^3 and inertia to s^5. Fields use the finite-volume medium rounded rotor, 64 angular samples with 128-sample refinement. Enforced clearance is >=50 mm beyond a 12.7 mm receiver radius. The stress proxy s*(carrier+3*tone_spacing)<=60 Hz preserves the reference centrifugal stress scale. This is not new FEA or modal/rotordynamic validation.',
'',
'At greater distance, inertia rises much faster than diameter. Lower-rate coding requires faster modulation for the same useful throughput and increases acceleration torque; a coding gain must offset that mechanical penalty. Lower carriers reduce torque-speed power and stress but encounter more assumed environmental noise. Thus code rate and carrier cannot be selected independently of Tx size.',
'',
'The main table assumes zero drag to expose unavoidable inertial requirements. The following illustrative sensitivity uses viscous drag b=b_ref*s^3 with b_ref=0.05/(pi*50.3) N m s. This is a declared extrapolation, not a bearing/windage prediction. Stored kinetic energy and sustained worst-valid-byte RMS torque are also reported.',
'',
'| Distance (m) | Inertial RMS torque (N m) | Max stored energy (kJ) | With illustrative drag: peak torque (N m) | Peak power (kW) |',
'|---:|---:|---:|---:|---:|']
for x in selected:
 m=duty(x['symbol_s'],.4,.1,inertia=x['inertia_kg_m2'],fc=x['carrier_hz'],drag=.05*x['scale']**3*x['carrier_hz']/50.3)
 lines.append(f"| {x['distance_m']:g} | {x['worst_valid_message_rms_nm']:.2f} | {x['maximum_stored_energy_j']/1000:.2f} | {m['peak_torque_nm']:.2f} | {m['peak_mechanical_power_w']/1000:.3f} |")
lines+=['','## Coding and runtime conclusion','', 'The mass/torque ranking selects uncoded transmission because its measured BER is below 1e-3 in this small sample. This is not a recommendation to omit FEC: two of 32 frames were erased at 24 Hz/reference amplitude, and one of 32 at 18 Hz/1.25 amplitude. Rate-3/4 convolutional coding produced zero payload errors and zero erasures at both settings, using the same rotor sizes but about 1.66 times the peak inertial torque. Rate-2/3 coding also had zero errors there with still higher torque. LDPC at 24 Hz/reference amplitude had zero errors, but its lower code rate costs more torque and its 0.75-amplitude detector gate changed a decoded bit with history refinement. These trials do not establish a superior code-specific amplitude threshold.', '', 'Across all 736 packet trials, the maximum streaming block computation was 0.107 s and maximum FEC time 0.031 s. Maximum budgeted latency was 224.045 s. All computations finished; no background workers or scheduled jobs remain.', '', '## Scope and remaining work','',
'- Current CF-FSK q=0.1 and transition fraction 0.4 are fixed in this pass. Uniform geometric scaling and a coarse amplitude grid do not establish a global optimum. The shorter frames were screened but not statistically ruled out.',
'- NHNM-derived ground noise, Rayleigh coupling, the two-stage isolation transfer and the above-10-Hz extension are explicit assumptions (see ../environmental_background_v1/MODEL.md). Site measurements and isolation sensitivity may change the preferred frequency and size.',
'- Acquisition, timing/phase estimation and amplitude estimation are assumed. Twelve acquisition symbols are budgeted, but not simulated. Six tail symbols are simulated. The 224 s frame budget includes both; online lookahead overlaps reception, and only measured residual processing is added. Setup/JIT runs before reception. CPU timings under parallel workstation load are diagnostic, not a real-time certification.',
'- A single acquired complex envelope harmonic is modeled. Higher source harmonics, drive tracking, mechanical resonances, mounting and analog acquisition remain outside this link model. No enlarged motor or structure has been qualified.',
'- Next search: bracket amplitude near the promising code-specific thresholds, vary modulation excursion/transition duration with matching receiver convergence checks, and compare scalable rotor/drive alternatives before spending thousands of packets on final BER qualification.','',
'## Reproduction','',
'Run scripts/streaming_link_search.py --stage gates --workers 6, then --stage screen --workers 12 --packets 2. Run scripts/streaming_size_map.py. The saved shortlist.json specifies the refinement; scripts/refine_streaming_links.py --gates --workers 6, then --workers 12 --packets 32. The second refinement uses scripts/refine_streaming_links2.py and shortlist2.json with the same command options and separate packet seeds. All 18 size entries are saved in amplitude_size_map.json. Finally run scripts/compile_streaming_links.py and scripts/report_streaming_links.py. Source hashes, configuration IDs, seeds and individual packet records are retained. No manuscript performance claims were changed.','']
(OUT/'RESULTS.md').write_text('\n'.join(lines))
save(OUT/'report_manifest.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'scripts/compile_streaming_links.py',Path(__file__),ROOT/'scripts/streaming_size_map.py',ROOT/'fea/results/shaped_waveforms/medium_rounded.npz']})
print(json.dumps(selected,indent=2))
