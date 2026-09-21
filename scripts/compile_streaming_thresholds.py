"""Audit and compile fine thresholds, retaining packet clustering and drive tradeoffs."""
import json,hashlib
from pathlib import Path
from collections import defaultdict
from streaming_thresholds import OUT,ROOT,save
from compile_streaming_dynamics import mechanical
from gravcomm.packet_confidence import packet_ber_interval
for stage in ('gates','trials'):
 manifest=json.loads((OUT/f'{stage}_manifest.json').read_text());paths=list((OUT/stage).glob('*.json'));assert len(paths)==len(manifest['jobs'])
 for c,i,g in manifest['jobs']:
  x=json.loads((OUT/stage/f'{c["id"]}_{i}.json').read_text());assert x['config']==c and x['index']==i and x['gate']==g
  assert x['bytes_emitted']==x['mapped_bytes'] and not x['pruned_states'] and not x['missing_bit_hypotheses']
 for p,h in manifest['sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 print(stage,len(paths),'audited')
error_checks=[]
if (OUT/'error_checks_manifest.json').exists():
 em=json.loads((OUT/'error_checks_manifest.json').read_text())
 for c,i,g in em['jobs']:
  x=json.loads((OUT/'error_checks'/f'{c["id"]}_{i}.json').read_text());assert x['config']==c and x['index']==i;error_checks.append(x)
 for p,h in em['sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 print('error_checks',len(error_checks),'audited')
groups=defaultdict(list)
for p in (OUT/'trials').glob('*.json'):
 x=json.loads(p.read_text());groups[x['config']['id']].append(x)
sizes=json.loads((OUT/'amplitude_size_map.json').read_text());rows=[]
for ident,records in sorted(groups.items()):
 x=records[0];c=x['config'];ts=x['symbol_s'];errors=[r['errors'] for r in records];n=len(records);m=mechanical(ts,c['transition_fraction'],c['q'],c['carrier_hz']);gate=json.loads((OUT/'gates'/f'{ident}_50000.json').read_text())['passes_detector_gate']
 for size in sizes:
  if size['amplitude_scale']!=c['amplitude_scale']:continue
  stress=size['scale']*(c['carrier_hz']+3*c['q']/((1-c['transition_fraction'])*ts))/60
  row=dict(c|size,packets=n,payload_bits=x['payload_bits'],symbol_s=ts,errors=sum(errors),raw_coded_ber=sum(r['raw_errors'] for r in records)/(n*x['coded_bits']),ber=sum(errors)/(n*x['payload_bits']),ber_ci=packet_ber_interval(errors,x['payload_bits']),erasures=sum(not r['accepted'] for r in records),wrong_accepted=sum(r['wrong_accepted'] for r in records),accepted_payload_rate_bps=sum(r['accepted'] for r in records)/n,latency_s=max(r['latency_budget_s'] for r in records),max_block_s=max(r['max_block_s'] for r in records),stress_ratio=stress,peak_torque_nm=m['peak_torque_nm']*size['inertia_kg_m2'],rms_torque_nm=m['worst_valid_message_rms_nm']*size['inertia_kg_m2'],peak_power_w=m['peak_mechanical_power_w']*size['inertia_kg_m2'],gate=gate)
  row['error_packet_convergence']=all(v['passes_detector_gate'] for v in error_checks if v['config']['id']==ident) if error_checks else None
  row['ber_qualified95']=row['ber_ci'][1]<=.001
  row['ber_rejected95']=row['ber_ci'][0]>.001
  row['screen_eligible']=bool(row['ber']<=.001 and stress<=1 and size['radial_gap_m']>=.05 and row['latency_s']<300 and gate and row['error_packet_convergence'] is not False and all(r['passes_detector_gate'] for r in records));rows.append(row)
save(OUT/'compiled.json',rows)
selected={}
for d in (.5,1.,2.):
 choices=[x for x in rows if x['distance_m']==d and x['screen_eligible']]
 selected[str(d)]={}
 for label,key in [('mass',lambda x:(x['mass_kg'],x['peak_torque_nm'])),('torque',lambda x:(x['peak_torque_nm'],x['mass_kg']))]:
  if choices:selected[str(d)][label]=min(choices,key=key)
save(OUT/'selected.json',selected)
lines=['# Fine amplitude thresholds','',
'This pass repeats endpoints and adds intermediate amplitudes using 64 fresh packets per gate-approved configuration. No earlier adaptive pilot samples are pooled into the fresh statistics. Fixed instrumental floor, 2 Hz bandwidth, declared environmental noise and the 224-bit payload/224 s nominal frame budget are unchanged. Errors in rejected frames count toward BER.', '',
'## Fresh error measurements','',
'| Configuration | Errors / payload bits | Erasures | Wrong accepted | BER | Packetwise 95% upper bound | Accepted payload rate (bit/s) |',
'|:---|---:|---:|---:|---:|---:|---:|']
for x in rows:
 if x['distance_m']==.5:
  lines.append(f"| {x['id']} | {x['errors']} / {x['packets']*x['payload_bits']} | {x['erasures']} / {x['packets']} | {x['wrong_accepted']} | {x['ber']:.5g} | {x['ber_ci'][1]:.3g} | {x['accepted_payload_rate_bps']:.4f} |")
lines+=['','The nominal offered payload rate is 1 bit/s. Accepted rate above reflects CRC erasures with no retransmissions or additional inter-frame processing overhead; it is an empirical sample quantity, not qualified goodput. The confidence bounds allow within-packet correlation and are per fixed configuration, not simultaneous after selection. They do not establish BER <=1e-3.','',
'## Mass and torque alternatives','',
'Each row passes the empirical BER, latency, clearance, receiver convergence and centrifugal stress-similarity screens. This is not drive, structure or statistical reliability qualification.','',
'| Distance (m) | Objective | Configuration | Diameter (m) | Mass (kg) | Peak torque (N m) | RMS torque (N m) | Peak mechanical power (kW) |',
'|---:|:---|:---|---:|---:|---:|---:|---:|']
for d,choices in selected.items():
 for label,x in choices.items():lines.append(f"| {float(d):g} | {label} | {x['id']} | {x['diameter_m']:.3f} | {x['mass_kg']:.2f} | {x['peak_torque_nm']:.2f} | {x['rms_torque_nm']:.2f} | {x['peak_power_w']/1000:.3f} |")
lines += ['', f'Maximum measured streaming block time: {max(x["max_block_s"] for x in rows):.3f} s. Maximum budgeted latency: {max(x["latency_s"] for x in rows):.3f} s. Timings are diagnostic measurements under parallel host load.', '']
lines+=['','Mechanical values exclude drag and assume achieved commanded motion. Peak absolute mechanical power is not average electrical draw. Uniform finite-volume scaling and nominal stress similarity do not substitute for new FEA, modal/shaft analysis or drive selection.','',
'## Interpretation and qualification','',
'Observed errors at the lower amplitude and zero/few errors at the higher amplitude provide an empirical bracket only. A rare decoded error burst can change the apparent threshold between 32- and 64-packet runs. Choose a reliability margin before freezing a design, then use disjoint packets for qualification. The present short trials are intended to reject clearly poor candidates and compare physical costs.', '',
'With the existing time-uniform packet-fraction bound, 5374 zero-error independent packets are needed for a per-configuration 95% upper BER bound <=1e-3. Any errors or simultaneous coverage requirements change that budget; see qualification_budget.json.', '',
'The manuscript now explains why the source uses an equatorial receiver with radial sensitivity: it maximizes the leading rotating-quadrupole signal at fixed distance, and sampled finite-volume checks support that choice. Directional site noise could change the SNR-optimal orientation. See ../streaming_dynamics_v2/ORIENTATION.md. The revised manuscript builds successfully; these unqualified communications results have not been inserted.', '',
'## Reproduction','',
'Run scripts/streaming_thresholds.py --stage gates --workers 8, then --stage trials --workers 12. Source hashes, full configuration lists, separate payload/noise seed construction and seed ranges are saved in the manifests and runtime. Amplitude maps reuse earlier finite-volume results and add streaming_size_map.one at .8125/.9375/1.1875 for all distances. Run scripts/compile_streaming_thresholds.py to verify every record and source hash and generate this report.','']
if error_checks:
 lines += ['','## Targeted error-packet checks','',f'Replayed {len(error_checks)} error-containing packets of the mass/torque candidates with history 5 and longer lookahead. {sum(x["passes_detector_gate"] for x in error_checks)} passed the convergence criterion. These are reruns of existing records and add no independent BER samples. Raw comparisons and source hashes are saved in error_checks/.', '']
(OUT/'RESULTS.md').write_text('\n'.join(lines))
save(OUT/'report_manifest.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/compile_streaming_dynamics.py',ROOT/'scripts/range_revision.py',ROOT/'fea/results/shaped_waveforms/medium_rounded.npz']})
for x in rows:
 if x['distance_m']==.5:print(x['id'],x['errors'],x['ber'],x['erasures'])
print('SELECTED',json.dumps(selected))
