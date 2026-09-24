"""Matched adjacent-spacing control: q=0.08, across tone counts."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import json,hashlib,time,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from tone_count_runtime import ROOT,OUT,run,save
from tone_count_study import configurations

def main():
 folder=OUT/'extra';folder.mkdir(exist_ok=True);configs=[]
 for M in [2,4,6,8]:
  for fc in [24.,18.]:
   c=next(dict(c) for c in configurations() if c['tone_count']==M and c['carrier_hz']==fc)
   c['q']=.08;c['span_symbol_product']=.2*(M-1);c['id']=f'M{M}_g{c["group_symbols"]}b{c["group_bits"]}_spanT{c["span_symbol_product"]:g}_f{fc:g}';configs.append(c)
 jobs=[(c,i,False) for c in configs for i in range(52000,52004)]
 save(OUT/'extra_manifest.json',dict(jobs=jobs,reason='Matched adjacent spacing q=0.08 control; avoids excluding eight-tone settings solely because the equal-span grid produces costly phase lattices.',sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/tone_count_runtime.py',ROOT/'src/gravcomm/alphabet_streaming.py']}))
 with ProcessPoolExecutor(max_workers=4) as pool:
  fs={pool.submit(run,j):j for j in jobs if not (folder/f'{j[0]["id"]}_{j[1]}.json').exists()}
  for n,f in enumerate(as_completed(fs),1):
   c,i,_=fs[f]
   try:x=f.result()
   except Exception as e:x=dict(config=c,index=i,error=repr(e),traceback=traceback.format_exc())
   save(folder/f'{c["id"]}_{i}.json',x);save(OUT/'extra_progress.json',dict(completed=n,total=len(fs),status='complete' if n==len(fs) else 'running'))
   print(json.dumps(dict(done=n,id=c['id'],errors=x.get('errors'),pruned=x.get('pruned_states'))),flush=True)
if __name__=='__main__':main()
