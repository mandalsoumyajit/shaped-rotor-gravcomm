"""Fresh-seed detector convergence for screened alphabet candidates."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import json,time,hashlib,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import tone_count_runtime as runtime
from gravcomm.alphabet_streaming import AlphabetStreamingReceiver
OUT=runtime.OUT

def wider(*args,**kwargs):
 if kwargs.get('history_symbols',4)>4:kwargs['beam_width']=200000
 return AlphabetStreamingReceiver(*args,**kwargs)

def run_one(c):
 runtime.StreamingSoftReceiver=wider
 return runtime.run((c,53000,True))

def main():
 configs=json.loads((OUT/'shortlist.json').read_text());folder=OUT/'gates';folder.mkdir(exist_ok=True)
 # Each candidate has its own immutable job manifest, allowing additional candidates later.
 pending=[]
 for c in configs:
  path=folder/f'{c["id"]}_53000.json'
  if path.exists():continue
  deps=[Path(__file__),Path(runtime.__file__),runtime.ROOT/'src/gravcomm/alphabet_streaming.py']
  manifest=dict(config=c,seed=53000,wide_history_budget=200000,sha256={str(p.relative_to(runtime.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
  runtime.save(folder/f'{c["id"]}.manifest.json',manifest);pending.append(c)
 done=len(configs)-len(pending);start=time.time()
 with ProcessPoolExecutor(max_workers=4) as pool:
  fs={pool.submit(run_one,c):c for c in pending}
  for f in as_completed(fs):
   c=fs[f]
   try:x=f.result()
   except Exception as e:x=dict(config=c,error=repr(e),traceback=traceback.format_exc())
   runtime.save(folder/f'{c["id"]}_53000.json',x);done+=1
   runtime.save(OUT/'gate_progress.json',dict(completed=done,total=len(configs),elapsed_s=time.time()-start,status='complete' if done==len(configs) else 'running'))
   print(json.dumps(dict(id=c['id'],passes=x.get('passes_detector_gate'),convergence=x.get('convergence'),error=x.get('error'))),flush=True)
if __name__=='__main__':main()
