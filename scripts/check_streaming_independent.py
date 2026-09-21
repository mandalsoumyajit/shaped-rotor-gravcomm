import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.streaming_soft import StreamingSoftReceiver
from gravcomm.compact_innovations import CompactInnovations
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/streaming_soft_v1';rows=[];T=224/114
for n,scale,seed,histories in ((32,.5,862,(3,4,5)),(128,1.,863,(4,))):
 A=5.180259430854157e-10*scale;r=CausalReceiver(carrier_hz=50.3,background=SeismicBackground());rng=np.random.default_rng(seed)
 v=rng.integers(0,256,n,dtype=np.uint8).astype(int);symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
 mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],T));y=mean+r.stationary_noise(len(mean),rng,A)
 for history in histories:
  d=StreamingSoftReceiver(r,T,A,n,history_symbols=history,beam_width=100000 if history==5 else 16384,lookahead_symbols=12)
  out=[];times=[];pos=0
  for symbol in range(3*n+6):
   count=d._geometry(symbol)[3];t=time.perf_counter();out.extend(d.push(y[pos:pos+count]));times.append(time.perf_counter()-t);pos+=count
  llr=np.concatenate([x['llr'] for x in out]);row=dict(nbytes=n,scale=scale,seed=seed,carrier_hz=50.3,history_symbols=history,**d.finish(),max_symbol_processing_s=max(times),total_processing_s=sum(times),llr=llr.tolist(),whitening_error=d.whitener.spectral_error())
  rows.append(row);(OUT/'independent_records.json').write_text(json.dumps(rows,indent=2)+'\n');print({k:v for k,v in row.items() if k!='llr'},flush=True)
