import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.causal_sequence import CausalBeamByteDetector
from gravcomm.streaming_soft import StreamingSoftReceiver
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/streaming_soft_v1';A=5.180259430854157e-10;T=224/114;n=32
r=CausalReceiver(carrier_hz=20,background=SeismicBackground());rng=np.random.default_rng(861)
v=rng.integers(0,256,n,dtype=np.uint8).astype(int);symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],T));y=mean+r.stationary_noise(len(mean),rng,A)
ref=json.loads((ROOT/'results/fixed_rate_sizing_2026_09_20/decoder_latency_v2/frame224.json').read_text())[-1];bits=np.unpackbits(np.array(ref['bytes'],dtype=np.uint8));rows=[]
for width in (4096,16384):
 d=CausalBeamByteDetector(r,T,A,width)
 for bit in (40,56,112,136):
  t=time.perf_counter();got=d.decode(y,n,forced_bits={bit:1-int(bits[bit])})
  llr=(got['path_cost']-ref['cost'])*(1-2*int(bits[bit]))
  row=dict(width=width,bit=bit,llr=llr,path_cost=got['path_cost'],runtime_s=time.perf_counter()-t,bytes=got['bytes'].tolist());rows.append(row)
  (OUT/'wide_reference.json').write_text(json.dumps(rows,indent=2)+'\n');print({k:v for k,v in row.items() if k!='bytes'},flush=True)
# Longer recent-history model, still bounded independently of packet length.
t=time.perf_counter();d=StreamingSoftReceiver(r,T,A,n,history_symbols=5,beam_width=100000,lookahead_symbols=12);out=[];times=[];pos=0
for symbol in range(3*n+6):
 count=d._geometry(symbol)[3];tic=time.perf_counter();out.extend(d.push(y[pos:pos+count]));times.append(time.perf_counter()-tic);pos+=count
result=dict(**d.finish(),history_symbols=5,total_s=sum(times),max_symbol_s=max(times),llr=np.concatenate([x['llr'] for x in out]).tolist())
(OUT/'history5.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k!='llr'},flush=True)
