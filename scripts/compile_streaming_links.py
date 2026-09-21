"""Compile streaming pilot results; no statistical qualification is inferred."""
import json,sys
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from streaming_link_search import ROOT,OUT,save
from gravcomm.packet_confidence import packet_ber_interval
from range_revision import duty
@lru_cache(None)
def unit(ts,fc):
 return duty(ts,.4,.1,inertia=1.,fc=fc,drag=0.)
def main():
 groups=defaultdict(list)
 for folder in ('screen','refinement','refinement2'):
  for p in (OUT/folder).glob('*.json'):
   x=json.loads(p.read_text());groups[x['config']['id']].append(x)
 sizes=json.loads((OUT/'amplitude_size_map.json').read_text());rows=[]
 for ident,records in sorted(groups.items()):
  pilot_records=[v for v in records if v['index']<20000];fresh_records=[v for v in records if v['index']>=20000];records=fresh_records or pilot_records
  x=records[0];c=x['config'];errors=[v['errors'] for v in records];n=len(records);bits=x['payload_bits'];ts=x['symbol_s'];fc=c['carrier_hz']
  g=OUT/'candidate_gates'/f'{ident}_10000.json'
  if not g.exists():g=OUT/'candidate_gates2'/f'{ident}_10000.json'
  gate=json.loads(g.read_text())['passes_detector_gate'] if g.exists() else None
  for size in sizes:
   if size['amplitude_scale']!=c['amplitude_scale'] or 'excluded' in size:continue
   stress=size['scale']*(fc+3*.1/(.6*ts))/60
   mechanics={k:v*size['inertia_kg_m2'] for k,v in unit(ts,fc).items() if k in ('peak_torque_nm','worst_valid_message_rms_nm','peak_mechanical_power_w','maximum_stored_energy_j')}
   latency=max(v['latency_budget_s'] for v in records)
   row=dict(c | size,symbol_s=ts,payload_bits=bits,packets=n,payload_errors=sum(errors),ber=sum(errors)/(n*bits),ber_interval95_packetwise=packet_ber_interval(errors,bits),packet_failures=sum(v['failure'] for v in records),erasures=sum(not v['accepted'] for v in records),undetected_wrong_packets=sum(v['wrong_accepted'] for v in records),latency_budget_s=latency,max_block_s=max(v['max_block_s'] for v in records),max_fec_s=max(v['fec_s'] for v in records),stress_similarity_ratio=stress,geometry_stress_screen=stress<=1 and size['radial_gap_m']>=.05,**mechanics,carrier_rpm=30*fc,record_quality_pass=all(v['passes_detector_gate'] for v in records),candidate_convergence_pass=gate,statistical_sample='fresh_refinement' if fresh_records else 'pilot',pilot_packets=len(pilot_records),pilot_payload_errors=sum(v['errors'] for v in pilot_records),candidate_gate_payload_errors=json.loads(g.read_text())['errors'] if g.exists() else None)
   row['screen_eligible']=bool(row['geometry_stress_screen'] and row['record_quality_pass'] and latency<300 and row['ber']<=.001 and gate is not False)
   rows.append(row)
 save(OUT/'compiled_candidates.json',rows)
 for d in (.5,1.,2.):
  choices=sorted([x for x in rows if x['distance_m']==d and x['screen_eligible']],key=lambda x:(x['mass_kg'],x['peak_torque_nm'],x['peak_mechanical_power_w']))
  print('DISTANCE',d)
  for x in choices[:8]:print(json.dumps({k:x[k] for k in ('id','mass_kg','diameter_m','packets','payload_errors','peak_torque_nm','peak_mechanical_power_w','candidate_convergence_pass')}))
if __name__=='__main__':main()


