"""Cost/convergence probe for narrower full-packet beam searches; no BER claim."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import json,time
from pathlib import Path
import numpy as np
import band_coded_pilot as pilot
Original=pilot.BeamByteDetector
root=Path(__file__).resolve().parents[1]
out=root/'results/fixed_rate_sizing_2026_09_20/bandwidth_audit_v1'
manifest=json.loads((root/'results/fixed_rate_sizing_2026_09_20/band_coded_pilot_v1/manifest.json').read_text())
c=next(c for c in manifest['configs'] if c['id']=='s0.95_tb_3_4')
records=[]
for width in (512,2048,8192):
 captured=[]
 class Probe(Original):
  def __init__(self,model,taps,beam_width):super().__init__(model,taps,width)
  def decode(self,*args,**kwargs):
   r=super().decode(*args,**kwargs)
   if len(args)<3 and not kwargs.get('forced_bits'):captured.append(r['llr'].copy())
   return r
 pilot.BeamByteDetector=Probe
 before=time.monotonic();result=pilot.run((c,1))
 records.append(dict(width=width,llr=captured[-1].tolist(),result=result,seconds=time.monotonic()-before))
 pilot.save(out/'full_packet_width_probe.json',records)
 print(json.dumps(dict(width=width,errors=result['errors'],missing=result['missing'],seconds=records[-1]['seconds'])),flush=True)
ref=np.array(records[-1]['llr'])
for row in records:
 row['relative_llr_change']=float(np.linalg.norm(np.array(row['llr'])-ref)/np.linalg.norm(ref))
 row['llr_sign_changes']=int(np.count_nonzero((np.array(row['llr'])<0)!=(ref<0)))
pilot.save(out/'full_packet_width_probe.json',records)
print('WIDTH_PROBE_COMPLETE',flush=True)
