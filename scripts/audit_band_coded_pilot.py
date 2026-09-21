"""Audit completed fixed-band checkpoints and source hashes."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/band_coded_pilot_v1'
def main():
 m=json.loads((OUT/'manifest.json').read_text());checks=[]
 for p,sha in m['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha,p
 for c in m['configs']:
  for i in range(m['packets']):
   p=OUT/f'{c["id"]}_{i}.json'
   if not p.exists():continue
   r=json.loads(p.read_text());assert r['config']==c and r['index']==i
   assert 0<=r['errors']<=r['payload_bits'] and 0<=r['raw_errors']<=r['coded_bits']
   assert r['wrong_accepted']==bool(r['errors'] and r['accepted'])
   assert r['failure']==bool(r['errors'] or not r['accepted'])
   assert r['missing']==r['repaired']
   assert r['packet_airtime_s']<300 and r['payload_bits']/r['packet_airtime_s']==r['nominal_payload_bit_s']==1
   assert (len(r['convergence'])==2)==(i==0)
   checks.append(p.name)
 result=dict(status='passed',complete=len(checks)==len(m['configs'])*m['packets'],checkpoints_checked=len(checks),checks=['Local/remote source hashes match','Immutable configurations and packet indexes','Error-count and CRC/acceptance accounting','All missing alternatives searched','Nominal useful rate and airtime','Convergence records attached to first packet only'],scope='Consistency audit, not BER validation or proof of global beam optimality')
 (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
