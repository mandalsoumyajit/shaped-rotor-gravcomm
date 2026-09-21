"""Independent shortlist checks and additional packets, preserving pilot provenance."""
import argparse,json,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
from streaming_link_search import OUT,ROOT,run,save
from pathlib import Path
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--gates',action='store_true');p.add_argument('--workers',type=int,default=8);p.add_argument('--packets',type=int,default=32);a=p.parse_args()
 configs=json.loads((OUT/'shortlist2.json').read_text());stage='candidate_gates2' if a.gates else 'refinement2';folder=OUT/stage;folder.mkdir(exist_ok=True)
 if not a.gates:
  configs=[c for c in configs if json.loads((OUT/'candidate_gates2'/f'{c["id"]}_10000.json').read_text())['passes_detector_gate']]
 jobs=[(c,i,a.gates) for c in configs for i in ([10000] if a.gates else range(30000,30000+a.packets))]
 manifest=dict(configs=configs,jobs=len(jobs),packet_start=10000 if a.gates else 30000,sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/streaming_link_search.py',*(ROOT/'src/gravcomm').glob('*.py')]})
 mp=OUT/f'{stage}_manifest.json'
 if mp.exists():assert json.loads(mp.read_text())==manifest
 else:save(mp,manifest)
 pending=[j for j in jobs if not (folder/f'{j[0]["id"]}_{j[1]}.json').exists()];done=len(jobs)-len(pending)
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  fs={pool.submit(run,j):j for j in pending}
  for f in as_completed(fs):
   j=fs[f];x=f.result();save(folder/f'{j[0]["id"]}_{j[1]}.json',x);done+=1
   save(OUT/'progress.json',dict(stage=stage,completed=done,total=len(jobs),status='complete' if done==len(jobs) else 'running'))
   print(json.dumps(dict(done=done,total=len(jobs),id=j[0]['id'],errors=x['errors'],gate=x['passes_detector_gate'])),flush=True)
