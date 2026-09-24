"""Resumeable controlled alphabet/spacing screen; outputs never overwrite prior studies."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse,json,math,hashlib,time,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from tone_count_runtime import ROOT,OUT,run,save

def configurations():
 configs=[]
 for M in range(2,10):
  options=[]
  for L in range(1,5):
   bits=(M**L).bit_length()-1
   if bits>12 or (M&(M-1) and L==1):continue
   options.append((bits/L,-L,L,bits))
  _,_,L,bits=max(options)
  mappings=[(L,bits)]+([(3,8)] if M==7 else [])
  for L,bits in mappings:
   for spanT in [.4,.8,1.2,1.6]:
    for fc,amp in [(24.,.8125),(18.,1.)]:
     q=.4*spanT/(M-1)
     configs.append(dict(id=f'M{M}_g{L}b{bits}_spanT{spanT:g}_f{fc:g}',scheme='tb_2/3',information_bits=256,carrier_hz=fc,transition_fraction=.6,q=q,amplitude_scale=amp,tone_count=M,group_symbols=L,group_bits=bits,span_symbol_product=spanT,beam_width=32768))
 return configs

def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['smoke','screen','gates','confirm'],default='screen');p.add_argument('--workers',type=int,default=8);args=p.parse_args()
 OUT.mkdir(parents=True,exist_ok=True);configs=configurations()
 if args.stage=='smoke':configs=[c for c in configs if c['tone_count'] in (2,7,8) and c['span_symbol_product']==1.2 and c['carrier_hz']==24]
 if args.stage in ('gates','confirm'):configs=json.loads((OUT/'shortlist.json').read_text())
 seeds={'smoke':range(51000,51001),'screen':range(52000,52004),'gates':range(53000,53001),'confirm':range(54000,54032)}[args.stage]
 jobs=[(c,i,args.stage=='gates') for c in configs for i in seeds];folder=OUT/args.stage;folder.mkdir(exist_ok=True)
 paths=[Path(__file__),ROOT/'scripts/tone_count_runtime.py',ROOT/'src/gravcomm/alphabet_streaming.py']
 manifest=dict(jobs=jobs,sha256={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in paths},scope='Exploratory screen; acquired packets; fixed qualified source amplitudes; no BER qualification.')
 dest=OUT/f'{args.stage}_manifest.json'
 if dest.exists():assert json.loads(dest.read_text())==json.loads(json.dumps(manifest))
 else:save(dest,manifest)
 pending=[j for j in jobs if not (folder/f'{j[0]["id"]}_{j[1]}.json').exists()];done=len(jobs)-len(pending);start=time.time()
 save(OUT/'progress.json',dict(stage=args.stage,completed=done,total=len(jobs),status='running'))
 with ProcessPoolExecutor(max_workers=args.workers) as pool:
  futures={pool.submit(run,j):j for j in pending}
  for f in as_completed(futures):
   c,index,gate=futures[f]
   try:x=f.result()
   except Exception as exc:x=dict(config=c,index=index,error=repr(exc),traceback=traceback.format_exc())
   save(folder/f'{c["id"]}_{index}.json',x);done+=1
   save(OUT/'progress.json',dict(stage=args.stage,completed=done,total=len(jobs),elapsed_s=time.time()-start,status='complete' if done==len(jobs) else 'running'))
   print(json.dumps(dict(done=done,total=len(jobs),id=c['id'],errors=x.get('errors'),pruned=x.get('pruned_states'),seconds=x.get('elapsed_s'),error=x.get('error'))),flush=True)
if __name__=='__main__':main()
