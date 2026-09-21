"""Verify every frozen qualification record and reproduce final confidence decisions."""
import json,hashlib
import numpy as np
from pathlib import Path
from qualify_replacement_ber import OUT,ROOT,qualify
from streaming_dynamics_runtime import save
manifest=json.loads((OUT/'manifest.json').read_text());frozen=json.loads((OUT/'frozen_designs.json').read_text());final=json.loads((OUT/'final.json').read_text());assert manifest['frozen']==frozen
for relative,digest in manifest['source_sha256'].items():assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest,relative
lines=['# Parallel qualification of larger 2 m designs','',
'Three larger replacement candidates were frozen before this run; the failed 218.91 kg design remains unchanged. All qualification packets are independent of the earlier screening samples. The streaming receiver, code, packet length, carrier, amplitude, instrumental profile and declared environmental spectrum were held fixed. These conclusions concern the simulated, acquired-packet link; they do not qualify physical transmitter motion, acquisition, site noise or hardware.', '',
'All three candidates are at 2 m and preserve the 18 Hz carrier, rate-2/3 code, r=.6, q=.08, fixed noise/bandwidth and 224 s frame. Only source size and received amplitude change. Confidence applies simultaneously to this new family of three candidates; it is not a combined confidence statement across historical studies.', '',
'## Decision rule','',
'The error fraction in each 224-bit tentative payload is one independent bounded observation. Errors in CRC-rejected frames count. Two-sided time-uniform confidence sequences use alpha=0.05/3 per candidate, giving at least 95% simultaneous confidence across all three candidates by a union bound. Pass requires the upper bound <=1e-3; fail requires the lower bound >1e-3. Decisions are tested only after deterministic 128-packet prefixes. Bit independence is not assumed. No screening data or diagnostic reruns enter the sample counts.', '',
'| Distance (m) | Mass (kg) | Carrier (Hz) | Amplitude / reference | Packets | Payload errors / bits | BER | Simultaneous-confidence interval | Decision |',
'|---:|---:|---:|---:|---:|---:|---:|:---|:---|']
verified={}
for channel in frozen['channels']:
 c=channel['config'];ident=c['id'];rs=sorted([json.loads(p.read_text()) for p in (OUT/ident).glob('*.json')],key=lambda x:x['index']);expected=final['channels'][ident];assert len(rs)==expected['packets'];assert len(rs)%manifest['batch']==0
 for i,x in enumerate(rs):
  assert x['config']==c and x['index']==channel['seed_start']+i and not x['gate'];assert x['payload_bits']==224 and 0<=x['errors']<=224
  assert x['passes_detector_gate'] and not x['missing_bit_hypotheses'] and not x['pruned_states'] and x['bytes_emitted']==x['mapped_bytes']
 result=qualify(rs,manifest['per_channel_alpha']);np.testing.assert_allclose(result['confidence_interval'],expected['confidence_interval'],rtol=1e-10,atol=1e-12);assert {k:v for k,v in result.items() if k!='confidence_interval'}=={k:v for k,v in expected.items() if k!='confidence_interval'};assert result['status'] in ('pass','fail')
 # Verify this is the first decisive checkpoint, not a selected later prefix.
 for n in range(manifest['batch'],len(rs),manifest['batch']):assert qualify(rs[:n],manifest['per_channel_alpha'])['status']=='running'
 for i in range(2):
  gate=json.loads((OUT/'gates'/f"{ident}_{channel['seed_start']-10+i}.json").read_text());assert gate['config']==c and gate['passes_detector_gate'] and gate['gate']
 verified[ident]=result
 for source in channel['source_designs']:
  lo,hi=result['confidence_interval'];lines.append(f"| {source['distance_m']:g} | {source['mass_kg']:.2f} | {c['carrier_hz']:g} | {c['amplitude_scale']:g} | {result['packets']} | {result['errors']} / {result['payload_bits']} | {result['ber']:.6g} | [{lo:.6g}, {hi:.6g}] | {result['status'].upper()} |")
lines+=['','## Packet acceptance and timing','',
'| Channel | CRC/header erasures | Packet failures | Wrong accepted payloads | Empirical accepted payload rate (bit/s) | Max block time (s) | Max budgeted latency (s) |',
'|:---|---:|---:|---:|---:|---:|---:|']
for ident,r in verified.items():lines.append(f"| {ident} | {r['erasures']} | {r['packet_failures']} | {r['wrong_accepted']} | {r['accepted_payload_rate_bps']:.6f} | {r['max_block_s']:.3f} | {r['max_latency_s']:.3f} |")
lines+=['','Accepted rate is nominal offered payload rate (1 bit/s) times observed acceptance fraction, without retransmissions or extra inter-frame processing overhead. It is not a qualified throughput guarantee. A BER pass does not imply 1 bit/s accepted goodput, acquisition success or a bounded undetected-error probability. Zero observed wrong accepted packets is not a zero-probability claim.', '',
'Execution uses the implemented bounded streaming receiver, with packet processing concurrent with reception. Timing under many-worker CPU load is diagnostic rather than worst-case real-time certification. The modeled frame budget is 224 s, including assumed acquisition overhead and simulated tail; known timing/phase and calibrated input are prerequisites.', '',
'## Scope and reproduction','',
'All individual packet records, frozen design parameters, seed starts, source hashes and final bounds are retained. Run scripts/qualify_replacement_ber.py --workers 64 --batch 128 on the same runtime; completed deterministic prefixes resume without reusing earlier screening samples. Run scripts/compile_replacement_ber.py to verify source hashes, packet identity/quality, first stopping checkpoint and final decisions. The protocol is PLAN.md. Frozen failures must be redesigned and independently requalified; their amplitude or code was not changed during this run.', '']
em=json.loads((OUT/'error_checks_manifest.json').read_text())
for relative,digest in em['sha256'].items():assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest
for c,i,g in em['jobs']:
 x=json.loads((OUT/'error_checks'/f"{c['id']}_{i}.json").read_text());original=json.loads((OUT/c['id']/f'{i}.json').read_text())
 assert x['config']==c and x['index']==i and x['gate'] and x['passes_detector_gate'] and x['errors']==original['errors']
for ch in frozen['channels']:
 ident=ch['config']['id'];rs=sorted([json.loads(p.read_text()) for p in (OUT/ident).glob('*.json')],key=lambda x:x['index'])
 expected=[r['index'] for r in rs if r['errors']][:3]
 assert sorted(i for c,i,g in em['jobs'] if c['id']==ident)==expected
lines+=['## Error-packet convergence checks','',f"Replayed the first three error-containing packets per candidate (or all if fewer): {len(em['jobs'])} checks passed with wider history and longer lookahead. These reruns do not add BER samples.",'']
passing=[ch for ch in frozen['channels'] if verified[ch['config']['id']]['status']=='pass']
selected=min(passing,key=lambda ch:ch['source_designs'][0]['mass_kg']) if passing else None
save(OUT/'selected.json',selected)
lines+=['## Physical comparison','', '| Amplitude / reference | Diameter (m) | Mass (kg) | Peak torque (N m) | Peak mechanical power (kW) | Stress-similarity ratio |','|---:|---:|---:|---:|---:|---:|']
for ch in frozen['channels']:
 x=ch['source_designs'][0];lines.append(f"| {x['amplitude_scale']:g} | {x['diameter_m']:.3f} | {x['mass_kg']:.2f} | {x['peak_torque_nm']:.2f} | {x['peak_power_w']/1000:.2f} | {x['stress_ratio']:.4f} |")
lines+=['','Mechanical estimates use uniform finite-volume scaling and achieved commanded motion with zero drag. The stress-similarity ratio is the existing screen, not a new yield-stress calculation; peak mechanical power is not average electrical power. New structure and actuator validation remain necessary.','']
if selected:
 x=selected['source_designs'][0];lines+= [f"Smallest passing tested candidate: {x['mass_kg']:.2f} kg, diameter {x['diameter_m']:.3f} m. This is not a continuous optimum. Delivered-rate compensation still requires changed timing and independent testing.",'']
(OUT/'RESULTS.md').write_text('\n'.join(lines));save(OUT/'audit.json',dict(source_hashes_verified=True,record_and_prefix_decisions_verified=True,channels=verified));print(json.dumps(verified,indent=2))
