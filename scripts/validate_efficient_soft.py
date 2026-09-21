"""Soft-search convergence audit; diagnostic records, not BER qualification."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.causal_sequence import CausalBeamByteDetector
from gravcomm.environmental_noise import SeismicBackground
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/efficient_soft_v1';OUT.mkdir(exist_ok=True)
A=5.180259430854157e-10;rows=[];checks=[]
def save():
 (OUT/'validation.json').write_text(json.dumps(dict(status='diagnostic_in_progress',runs=rows,reference_checks=checks),indent=2,allow_nan=False)+'\n')
def measure(d,y,n,mode,quota=0):
 t=time.perf_counter()
 got=d.decode_diverse(y,n,bit_hypotheses=quota) if mode=='diverse' else d.decode(y,n)
 elapsed=time.perf_counter()-t
 row=dict(mode=mode,nbytes=n,carrier_hz=d.receiver.carrier_hz,amplitude=d.amplitude,global_width=d.width,quota=quota,runtime_s=elapsed,missing=got['missing_bit_hypotheses'],peak_survivors=got['peak_survivors'],cost=got['path_cost'],exact=got['exact_search'])
 rows.append(row);save();print(json.dumps(row),flush=True)
 return got,row
def record(n,duration,fc,amplitude,seed):
 r=CausalReceiver(carrier_hz=fc,background=SeismicBackground());rng=np.random.default_rng(seed)
 data=rng.integers(0,256,n,dtype=np.uint8);v=data.astype(int)
 symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
 mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],duration))
 return r,mean+r.stationary_noise(len(mean),rng,amplitude)
def compare(a,b):
 good=np.isfinite(a)&np.isfinite(b)
 return dict(max_llr_difference=float(np.max(abs(a[good]-b[good]))) if good.any() else None,relative_llr_l2=float(np.linalg.norm(a[good]-b[good])/max(np.linalg.norm(b[good]),1e-30)),sign_changes=int(np.count_nonzero(np.signbit(a[good])!=np.signbit(b[good]))),compared_bits=int(good.sum()))
for fc,scale,seed in ((20.,1.,618),(50.3,.5,619),(20.,.5,620)):
 r,y=record(2,.5,fc,A*scale,seed)
 exact,refrow=measure(CausalBeamByteDetector(r,.5,A*scale,65536),y,2,'plain')
 for quota in (1,4,16):
  got,row=measure(CausalBeamByteDetector(r,.5,A*scale,64),y,2,'diverse',quota)
  row['versus_exhaustive']=compare(got['llr'],exact['llr']);save()
# All records use a complete causal likelihood with the new seismic background.
for fc,scale,seed in ((20.,1.,761),(50.3,.5,762)):
 n=40;duration=256/126;r,y=record(n,duration,fc,A*scale,seed)
 outputs=[]
 for width,quota in ((128,0),(128,1),(128,4),(256,8)):
  got,row=measure(CausalBeamByteDetector(r,duration,A*scale,width),y,n,'diverse' if quota else 'plain',quota)
  outputs.append((got,row))
 best=outputs[-1][0]
 for got,row in outputs:
  row['versus_largest_diverse']=compare(got['llr'],best['llr']);save()
 # Select weak bits plus bits distributed over packet; forced reference retains
 # unrestricted global beam within each constraint. This is not exhaustive.
 indexes=sorted(set(np.linspace(0,8*n-1,8,dtype=int).tolist()+np.argsort(abs(best['llr']))[:4].tolist()))
 d=CausalBeamByteDetector(r,duration,A*scale,512)
 for bit in indexes:
  value=int(np.unpackbits(best['bytes'])[bit]);t=time.perf_counter()
  alt=d.decode(y,n,forced_bits={int(bit):1-value})
  labels=np.unpackbits(best['survivor_bytes'],axis=1)[:,bit]
  costs=best['survivor_costs'];pool_min=[float(costs[labels==v].min()) for v in (0,1)]
  altlabels=np.unpackbits(alt['survivor_bytes'],axis=1)[:,bit]
  for v in (0,1):
   relevant=alt['survivor_costs'][altlabels==v]
   if len(relevant):pool_min[v]=min(pool_min[v],float(relevant.min()))
  pooled=pool_min[1]-pool_min[0]
  check=dict(nbytes=n,carrier_hz=fc,scale=scale,bit=int(bit),diverse_llr=float(best['llr'][bit]),pooled_forced_llr=pooled,improvement_in_opposite_cost=float((costs[labels==1-value].min())-alt['path_cost']),runtime_s=time.perf_counter()-t)
  checks.append(check);save();print(json.dumps(check),flush=True)
result=dict(status='complete_diagnostic_not_BER_qualification',runs=rows,reference_checks=checks,source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('src/gravcomm/causal_sequence.py','src/gravcomm/causal_receiver.py','src/gravcomm/environmental_noise.py','scripts/validate_efficient_soft.py','tests/test_causal_sequence.py')})
(OUT/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
