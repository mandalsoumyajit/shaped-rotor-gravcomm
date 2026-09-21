import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,argparse,hashlib
from pathlib import Path
import numpy as np
root=Path('AIP_Advances_Submission');sys.path.insert(0,str(root/'src'))
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.causal_sequence import CausalBeamByteDetector
from gravcomm.environmental_noise import SeismicBackground
def main():
 out=root/'results/fixed_rate_sizing_2026_09_20/efficient_soft_v1';rows=[]
 parser=argparse.ArgumentParser()
 parser.add_argument('--nbytes',type=int,default=40)
 parser.add_argument('--payload-bits',type=int,default=None)
 parser.add_argument('--backend',choices=['thread','process'],default='thread')
 parser.add_argument('--widths',type=int,nargs='+',default=[32,128])
 parser.add_argument('--carrier',type=float,default=20.)
 parser.add_argument('--scale',type=float,default=1.)
 parser.add_argument('--seed',type=int,default=761)
 parser.add_argument('--output',default='checkpoint_packet.json')
 args=parser.parse_args()
 A=5.180259430854157e-10*args.scale
 r=CausalReceiver(carrier_hz=args.carrier,background=SeismicBackground());rng=np.random.default_rng(args.seed)
 n=args.nbytes
 if args.payload_bits is not None and args.payload_bits!=8*n-32:raise ValueError('Payload plus 16-bit header and 16-bit CRC must fill mapped bytes')
 duration=256/126 if args.payload_bits is None else args.payload_bits/(3*n+18)
 data=rng.integers(0,256,n,dtype=np.uint8);v=data.astype(int)
 symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
 mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],duration));y=mean+r.stationary_noise(len(mean),rng,A)
 for width in args.widths:
  d=CausalBeamByteDetector(r,duration,A,width);t=time.perf_counter();got=d.decode_checkpoint_soft(y,n,workers=4,parallel_backend=args.backend);elapsed=time.perf_counter()-t
  row=dict(nbytes=n,symbol_s=duration,payload_bits=args.payload_bits,airtime_with_acquisition_budget_s=duration*(3*n+18),width=width,workers=4,backend=args.backend,runtime_s=elapsed,cost=got['path_cost'],llr=got['llr'].tolist(),bytes=got['bytes'].tolist(),method=got['method'],carrier_hz=args.carrier,scale=args.scale,seed=args.seed,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'src/gravcomm/causal_sequence.py',root/'scripts/benchmark_checkpoint_soft.py')});rows.append(row)
  (out/args.output).write_text(json.dumps(rows,indent=2)+'\n');print(width,elapsed,got['path_cost'],flush=True)
if __name__=='__main__':main()
