"""Extend pilot brackets below the initial grid where its lowest amplitude passes."""
import json,hashlib
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
from streaming_dynamics_runtime import OUT,ROOT,run,save
from streaming_dynamics_search import config
if __name__=='__main__':
 groups=defaultdict(list)
 for p in (OUT/'screen').glob('*.json'):
  x=json.loads(p.read_text());groups[x['config']['id']].append(x)
 jobs=[]
 for rs in groups.values():
  c=rs[0]['config'];lower=.75 if c['carrier_hz']==24 else 1.
  if c['amplitude_scale']!=lower or sum(x['errors'] for x in rs)/(224*len(rs))>.001:continue
  for amp in ((.5,.625) if c['carrier_hz']==24 else (.75,.875)):
   d=config(c['scheme'],c['carrier_hz'],c['transition_fraction'],c['q'],amp);jobs.extend((d,i,False) for i in range(45000,45004))
 manifest=dict(jobs=jobs,sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/streaming_dynamics_search.py',ROOT/'scripts/streaming_dynamics_runtime.py',*(ROOT/'src/gravcomm').glob('*.py')]});mp=OUT/'extra_manifest.json'
 if mp.exists():assert json.loads(mp.read_text())==json.loads(json.dumps(manifest))
 else:save(mp,manifest)
 folder=OUT/'extra';folder.mkdir(exist_ok=True);pending=[j for j in jobs if not (folder/f'{j[0]["id"]}_{j[1]}.json').exists()];done=len(jobs)-len(pending)
 with ProcessPoolExecutor(max_workers=12) as pool:
  fs={pool.submit(run,j):j for j in pending}
  for f in as_completed(fs):
   j=fs[f];x=f.result();save(folder/f'{j[0]["id"]}_{j[1]}.json',x);done+=1;save(OUT/'progress.json',dict(stage='extra',completed=done,total=len(jobs)));print(json.dumps(dict(done=done,total=len(jobs),id=x['config']['id'],errors=x['errors'])),flush=True)
