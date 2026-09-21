"""Recheck failing packets without counting them as additional independent trials."""
import json,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from streaming_thresholds import OUT,ROOT,run,save
if __name__=='__main__':
 jobs=json.loads((OUT/'error_checks_jobs.json').read_text());folder=OUT/'error_checks';folder.mkdir(exist_ok=True);save(OUT/'error_checks_manifest.json',dict(jobs=jobs,sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/streaming_thresholds.py',ROOT/'scripts/streaming_dynamics_runtime.py',*(ROOT/'src/gravcomm').glob('*.py')]}))
 with ProcessPoolExecutor(max_workers=8) as pool:
  fs={pool.submit(run,j):j for j in jobs}
  for f in as_completed(fs):
   j=fs[f];x=f.result();save(folder/f'{j[0]["id"]}_{j[1]}.json',x);print(json.dumps(dict(id=x['config']['id'],index=x['index'],errors=x['errors'],gate=x['passes_detector_gate'],convergence=x['convergence'])),flush=True)
