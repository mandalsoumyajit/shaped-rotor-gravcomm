"""Verify every frozen qualification record and reproduce final confidence decisions."""
import json,hashlib
import numpy as np
from pathlib import Path
from qualify_streaming_ber import OUT,ROOT,qualify
from streaming_dynamics_runtime import save
manifest=json.loads((OUT/'manifest.json').read_text());frozen=json.loads((OUT/'frozen_designs.json').read_text());final=json.loads((OUT/'final.json').read_text());assert manifest['frozen']==frozen
selection_path=OUT.parent/'streaming_thresholds_v3/selected.json'
assert hashlib.sha256(selection_path.read_bytes()).hexdigest()==frozen['selected_source_sha256']
selection=json.loads(selection_path.read_text())
for channel in frozen['channels']:
 for source in channel['source_designs']:assert source==selection[str(source['distance_m'])]['mass']
for relative,digest in manifest['source_sha256'].items():assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest,relative
lines=['# Independent BER qualification of frozen designs','',
'The minimum-mass candidates were frozen before this run. All qualification packets are independent of the earlier screening samples. The streaming receiver, code, packet length, carrier, amplitude, instrumental profile and declared environmental spectrum were held fixed. These conclusions concern the simulated, acquired-packet link; they do not qualify physical transmitter motion, acquisition, site noise or hardware.', '',
'Two distinct normalized channels cover the three distances. The 0.5 and 1 m transmitters were sized to the same received harmonic amplitude and use the same carrier/receiver/noise model. Their shared BER result is mapped to both source sizes; field-inversion discretization tolerances are retained in the frozen source records.', '',
'## Decision rule','',
'The error fraction in each 224-bit tentative payload is one independent bounded observation. Errors in CRC-rejected frames count. Two-sided time-uniform confidence sequences use alpha=0.025 per channel, giving at least 95% simultaneous confidence across both channels by a union bound. Pass requires the upper bound <=1e-3; fail requires the lower bound >1e-3. Decisions are tested only after deterministic 128-packet prefixes. Bit independence is not assumed. No screening data or diagnostic reruns enter the sample counts.', '',
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
'All individual packet records, frozen design parameters, seed starts, source hashes and final bounds are retained. Run scripts/qualify_streaming_ber.py --workers 40 --batch 128 on the same runtime; completed deterministic prefixes resume without reusing earlier screening samples. Run scripts/compile_ber_qualification.py to verify source hashes, packet identity/quality, first stopping checkpoint and final decisions. The protocol is PROTOCOL.md. Frozen failures must be redesigned and independently requalified; their amplitude or code was not changed during this run.', '']
(OUT/'RESULTS.md').write_text('\n'.join(lines));save(OUT/'audit.json',dict(source_hashes_verified=True,record_and_prefix_decisions_verified=True,channels=verified));print(json.dumps(verified,indent=2))
