"""Online block-timing and reduced-state convergence audit."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.compact_innovations import CompactInnovations
from gravcomm.streaming_soft import StreamingSoftReceiver
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/streaming_soft_v1';OUT.mkdir(exist_ok=True)
A=5.180259430854157e-10;T=224/114;n=32
r=CausalReceiver(carrier_hz=20,background=SeismicBackground());rng=np.random.default_rng(861)
v=rng.integers(0,256,n,dtype=np.uint8).astype(int);symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],T));y=mean+r.stationary_noise(len(mean),rng,A)
reference=json.loads((ROOT/'results/fixed_rate_sizing_2026_09_20/decoder_latency_v2/frame224.json').read_text())[-1]
ref=np.array(reference['llr']);whitening=[];rows=[]
for order in (16,32,64,128):
 c=CompactInnovations(r,A,order);whitening.append(dict(order=order,**c.spectral_error()))
for history,width,lag in ((2,4096,12),(3,4096,12),(4,16384,12),(4,16384,24)):
 setup=time.perf_counter();d=StreamingSoftReceiver(r,T,A,n,noise_order=64,history_symbols=history,beam_width=width,lookahead_symbols=lag);setup=time.perf_counter()-setup
 outputs=[];times=[];position=0
 for symbol in range(3*n+6):
  count=d._geometry(symbol)[3];t=time.perf_counter();outputs.extend(d.push(y[position:position+count]));times.append(time.perf_counter()-t);position+=count
 stats=d.finish();llr=np.concatenate([x['llr'] for x in outputs]);valid=np.isfinite(llr)
 row=dict(history_symbols=history,beam_width=width,lookahead_symbols=lag,setup_s=setup,total_processing_s=sum(times),max_symbol_processing_s=max(times),p99_symbol_processing_s=float(np.quantile(times,.99)),symbol_duration_s=T,real_time_factor=sum(times)/(len(times)*T),max_output_delay_after_byte_s=max((x['available_after_symbol']-3*(x['byte_index']+1))*T for x in outputs),**stats,max_llr_difference=float(max(abs(llr[valid]-ref[valid]))) if valid.any() else None,relative_llr_l2=float(np.linalg.norm(llr[valid]-ref[valid])/np.linalg.norm(ref[valid])),sign_changes=int(np.count_nonzero(np.signbit(llr[valid])!=np.signbit(ref[valid]))),llr=[float(x) if np.isfinite(x) else None for x in llr],decoded_bytes=[x['byte'] for x in outputs],symbol_processing_s=times)
 rows.append(row)
 (OUT/'validation.json').write_text(json.dumps(dict(status='in_progress',whitening=whitening,runs=rows),indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:v for k,v in row.items() if k not in ('llr','decoded_bytes','symbol_processing_s')}),flush=True)
result=dict(status='complete_diagnostic_not_BER_qualification',whitening=whitening,runs=rows,source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('src/gravcomm/compact_innovations.py','src/gravcomm/streaming_soft.py','src/gravcomm/causal_receiver.py','src/gravcomm/environmental_noise.py','tests/test_streaming_soft.py','scripts/validate_streaming_soft.py')})
(OUT/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
