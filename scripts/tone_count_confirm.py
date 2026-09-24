"""Fresh 32-packet confirmation of candidates passing detector diagnostics."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import json,time,hashlib,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from tone_count_runtime import OUT,ROOT,run,save
from tone_count_study import configurations

def main():
 configs={}
 for p in (OUT/'gates').glob('*_53000.json'):
  x=json.loads(p.read_text())
  if x.get('passes_detector_gate'):configs[x['config']['id']]=x['config']
 # The existing qualified waveform is always retained as the paired control.
 for c in configurations():
  if c['tone_count']==7 and c['group_symbols']==3 and c['span_symbol_product']==1.2:configs[c['id']]=c
 folder=OUT/'confirm';folder.mkdir(exist_ok=True);jobs=[]
 for c in configs.values():
  manifest=dict(config=c,seeds=list(range(54000,54032)),sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/tone_count_runtime.py',ROOT/'src/gravcomm/alphabet_streaming.py']})
  dest=folder/f'{c["id"]}.manifest.json'
  if dest.exists():assert json.loads(dest.read_text())==manifest
  else:save(dest,manifest)
  jobs.extend((c,i,False) for i in range(54000,54032) if not (folder/f'{c["id"]}_{i}.json').exists())
 done=0;start=time.time()
 with ProcessPoolExecutor(max_workers=8) as pool:
  fs={pool.submit(run,j):j for j in jobs}
  for f in as_completed(fs):
   c,i,_=fs[f]
   try:x=f.result()
   except Exception as e:x=dict(config=c,index=i,error=repr(e),traceback=traceback.format_exc())
   save(folder/f'{c["id"]}_{i}.json',x);done+=1
   save(OUT/'confirm_progress.json',dict(completed=done,total_new_jobs=len(jobs),candidate_count=len(configs),elapsed_s=time.time()-start,status='complete' if done==len(jobs) else 'running'))
   if done%8==0 or done==len(jobs):print(json.dumps(dict(done=done,total=len(jobs),last_id=c['id'],last_errors=x.get('errors'))),flush=True)
if __name__=='__main__':main()
