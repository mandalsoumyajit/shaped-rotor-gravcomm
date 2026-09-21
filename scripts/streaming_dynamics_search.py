"""Staged waveform search; immutable manifests and fresh confirmation seeds."""
import argparse,json,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from streaming_dynamics_runtime import OUT,ROOT,run,save
WAVES=[(.4,.1),(.25,.125),(.5,.1),(.6,.08),(.6,.05),(.5,.0625),(.4,.075),(.5,.125)]
def config(scheme,fc,r,q,a):return dict(id=f'{scheme.replace("/","_")}_f{fc:g}_r{r:g}_q{q:g}_a{a:g}',scheme=scheme,information_bits=256,carrier_hz=fc,transition_fraction=r,q=q,amplitude_scale=a)
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['gates','screen','confirm_gates','confirm'],required=True);p.add_argument('--workers',type=int,default=12);a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
 gate=a.stage in ('gates','confirm_gates');jobs=[]
 if a.stage=='gates':
  jobs=[(config(s,fc,r,q,1.25 if fc==18 else 1.),40000,True) for s in ('uncoded','tb_3/4','tb_2/3') for fc in (18.,24.) for r,q in WAVES]
 elif a.stage=='screen':
  for p in (OUT/'gates').glob('*.json'):
   x=json.loads(p.read_text());c=x['config']
   if x['passes_detector_gate']:
    for amp in ((1.,1.125,1.25) if c['carrier_hz']==18 else (.75,.875,1.)):
     d=config(c['scheme'],c['carrier_hz'],c['transition_fraction'],c['q'],amp);jobs.extend((d,i,False) for i in range(41000,41004))
 else:
  configs=json.loads((OUT/'shortlist.json').read_text())
  if a.stage=='confirm':configs=[c for c in configs if json.loads((OUT/'confirm_gates'/f'{c["id"]}_42000.json').read_text())['passes_detector_gate']]
  jobs=[(c,i,gate) for c in configs for i in ([42000] if gate else range(43000,43032))]
 manifest=dict(jobs=jobs,sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/streaming_dynamics_runtime.py',*(ROOT/'src/gravcomm').glob('*.py')]})
 mp=OUT/f'{a.stage}_manifest.json'
 if mp.exists():assert json.loads(mp.read_text())==json.loads(json.dumps(manifest))
 else:save(mp,manifest)
 folder=OUT/a.stage;folder.mkdir(exist_ok=True);pending=[j for j in jobs if not (folder/f'{j[0]["id"]}_{j[1]}.json').exists()];done=len(jobs)-len(pending)
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  fs={pool.submit(run,j):j for j in pending}
  for f in as_completed(fs):
   j=fs[f];x=f.result();save(folder/f'{j[0]["id"]}_{j[1]}.json',x);done+=1;save(OUT/'progress.json',dict(stage=a.stage,completed=done,total=len(jobs),status='complete' if done==len(jobs) else 'running'))
   print(json.dumps(dict(done=done,total=len(jobs),id=x['config']['id'],errors=x['errors'],gate=x['passes_detector_gate'])),flush=True)
if __name__=='__main__':main()
