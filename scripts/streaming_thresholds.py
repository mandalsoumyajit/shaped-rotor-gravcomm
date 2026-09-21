"""Fine amplitude brackets with disjoint seeds and unchanged streaming runtime."""
import json,argparse,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from streaming_dynamics_runtime import ROOT,run,save
from streaming_dynamics_search import config
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/streaming_thresholds_v3'
def configs():
 rows=[]
 for scheme,fc,r,q,amps in [('tb_3/4',24.,.6,.08,(.8125,.875)),('tb_2/3',24.,.6,.08,(.8125,.875)),('tb_3/4',18.,.5,.125,(.8125,.875)),('tb_2/3',18.,.6,.08,(.9375,1.)),('uncoded',24.,.5,.1,(.875,.9375)),('uncoded',18.,.6,.08,(1.125,1.1875))]:
  rows.extend(config(scheme,fc,r,q,a) for a in amps)
 return rows
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['gates','trials'],required=True);p.add_argument('--workers',type=int,default=12);a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True);cs=configs();gate=a.stage=='gates'
 if not gate:cs=[c for c in cs if json.loads((OUT/'gates'/f'{c["id"]}_50000.json').read_text())['passes_detector_gate']]
 jobs=[(c,i,gate) for c in cs for i in ([50000] if gate else range(51000,51064))];mp=OUT/f'{a.stage}_manifest.json';manifest=dict(jobs=jobs,sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/streaming_dynamics_runtime.py',ROOT/'scripts/streaming_dynamics_search.py',*(ROOT/'src/gravcomm').glob('*.py')]})
 if mp.exists():assert json.loads(mp.read_text())==json.loads(json.dumps(manifest))
 else:save(mp,manifest)
 folder=OUT/a.stage;folder.mkdir(exist_ok=True);pending=[j for j in jobs if not (folder/f'{j[0]["id"]}_{j[1]}.json').exists()];done=len(jobs)-len(pending)
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  fs={pool.submit(run,j):j for j in pending}
  for f in as_completed(fs):
   j=fs[f];x=f.result();save(folder/f'{j[0]["id"]}_{j[1]}.json',x);done+=1;save(OUT/'progress.json',dict(stage=a.stage,completed=done,total=len(jobs),workers_running=done<len(jobs)));print(json.dumps(dict(done=done,total=len(jobs),id=c['id'] if False else x['config']['id'],errors=x['errors'],gate=x['passes_detector_gate'])),flush=True)
