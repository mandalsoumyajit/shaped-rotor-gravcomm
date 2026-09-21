import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,cProfile,pstats,io,json
from pathlib import Path
import numpy as np
root=Path('AIP_Advances_Submission');sys.path.insert(0,str(root/'src'))
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.causal_sequence import CausalBeamByteDetector
from gravcomm.environmental_noise import SeismicBackground
out=root/'results/fixed_rate_sizing_2026_09_20/decoder_latency_v2';out.mkdir(exist_ok=True)
n=32;duration=224/(3*n+18);A=5.180259430854157e-10
r=CausalReceiver(carrier_hz=20,background=SeismicBackground());rng=np.random.default_rng(861);v=rng.integers(0,256,n)
symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],duration));y=mean+r.stationary_noise(len(mean),rng,A)
d=CausalBeamByteDetector(r,duration,A,512)
pr=cProfile.Profile();pr.enable();base=d.decode(y,n,_save_checkpoints=True)
for bit in (0,32,64,96,128,160,192,224):
 d.decode(y,n,forced_bits={bit:1-int(np.unpackbits(base['bytes'])[bit])},_resume=base['_checkpoints'][bit//8])
pr.disable();pr.dump_stats(str(out/'receiver_profile_optimized.prof'));stream=io.StringIO();pstats.Stats(pr,stream=stream).strip_dirs().sort_stats('tottime').print_stats(30)
(out/'receiver_profile_optimized.txt').write_text(stream.getvalue());print(stream.getvalue())
