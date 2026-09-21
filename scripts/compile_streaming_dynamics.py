"""Map modulation trials to rotor and drive tradeoffs; no BER certification."""
import argparse,json
from collections import defaultdict
from functools import lru_cache
from streaming_dynamics_runtime import OUT,save
from range_revision import duty
from gravcomm.packet_confidence import packet_ber_interval
@lru_cache(None)
def mechanical(ts,r,q,fc):return duty(ts,r,q,inertia=1,fc=fc,drag=0)
def main():
 p=argparse.ArgumentParser();p.add_argument('--select',action='store_true');a=p.parse_args();groups=defaultdict(list)
 for stage in ('screen','extra','confirm'):
  for path in (OUT/stage).glob('*.json'):
   x=json.loads(path.read_text());groups[x['config']['id']].append(x)
 sizes=json.loads((OUT/'amplitude_size_map.json').read_text());rows=[]
 for ident,allrecords in groups.items():
  fresh=[x for x in allrecords if 43000<=x['index']<44000];records=fresh or allrecords;c=records[0]['config'];ts=records[0]['symbol_s'];errors=[x['errors'] for x in records];bits=records[0]['payload_bits'];m=mechanical(ts,c['transition_fraction'],c['q'],c['carrier_hz']);g=OUT/'confirm_gates'/f'{ident}_42000.json';gate=json.loads(g.read_text())['passes_detector_gate'] if g.exists() else None
  for size in sizes:
   if size['amplitude_scale']!=c['amplitude_scale']:continue
   stress=size['scale']*(c['carrier_hz']+3*c['q']/((1-c['transition_fraction'])*ts))/60
   row=dict(c|size,symbol_s=ts,packets=len(records),errors=sum(errors),ber=sum(errors)/(bits*len(records)),ber_ci=packet_ber_interval(errors,bits),erasures=sum(not x['accepted'] for x in records),wrong_accepted=sum(x['wrong_accepted'] for x in records),latency_s=max(x['latency_budget_s'] for x in records),max_block_s=max(x['max_block_s'] for x in records),stress_ratio=stress,confirmation=bool(fresh),candidate_gate=gate,peak_torque_nm=m['peak_torque_nm']*size['inertia_kg_m2'],rms_torque_nm=m['worst_valid_message_rms_nm']*size['inertia_kg_m2'],peak_power_w=m['peak_mechanical_power_w']*size['inertia_kg_m2'],stored_energy_j=m['maximum_stored_energy_j']*size['inertia_kg_m2'])
   row['eligible']=bool(row['ber']<=.001 and stress<=1 and size['radial_gap_m']>=.05 and row['latency_s']<300 and all(x['passes_detector_gate'] for x in records) and gate is not False);rows.append(row)
 save(OUT/'compiled.json',rows)
 if a.select:
  configs={}
  for fc,d in ((24.,.5),(18.,2.)):
   for scheme in ('uncoded','tb_3/4','tb_2/3'):
    choices=[x for x in rows if x['carrier_hz']==fc and x['distance_m']==d and x['scheme']==scheme and x['eligible']]
    for key in (lambda x:(x['mass_kg'],x['peak_torque_nm']),lambda x:(x['peak_torque_nm'],x['mass_kg'])):
     if choices:
      row=min(choices,key=key);configs[row['id']]=groups[row['id']][0]['config']
  # Include one higher-amplitude neighbor of each coded mass-minimum candidate,
  # so confirmation can distinguish a marginal threshold from a safer point.
  for fc,d in ((24.,.5),(18.,2.)):
   for scheme in ('tb_3/4','tb_2/3'):
    choices=[x for x in rows if x['carrier_hz']==fc and x['distance_m']==d and x['scheme']==scheme and x['eligible']]
    if not choices:continue
    best=min(choices,key=lambda x:(x['mass_kg'],x['peak_torque_nm']))
    neighbors=[x for x in choices if x['transition_fraction']==best['transition_fraction'] and x['q']==best['q'] and x['amplitude_scale']>best['amplitude_scale']]
    if neighbors:
     row=min(neighbors,key=lambda x:x['amplitude_scale']);configs[row['id']]=groups[row['id']][0]['config']
  save(OUT/'shortlist.json',list(configs.values()));print('SHORTLIST',len(configs))
  for c in configs.values():print(c)
 for d in (.5,1.,2.):
  choices=[x for x in rows if x['distance_m']==d and x['eligible'] and (x['confirmation'] if not a.select else True)]
  print('DISTANCE',d)
  for x in sorted(choices,key=lambda x:(x['mass_kg'],x['peak_torque_nm']))[:4]:print({k:x[k] for k in ('id','mass_kg','diameter_m','peak_torque_nm','peak_power_w','errors','packets')})
if __name__=='__main__':main()
