"""Audit waveform search outputs against immutable jobs and source hashes."""
import json,hashlib
from streaming_dynamics_runtime import OUT,ROOT
for stage in ('gates','screen','extra','confirm_gates','confirm'):
 m=json.loads((OUT/f'{stage}_manifest.json').read_text());jobs=m['jobs'];paths=list((OUT/stage).glob('*.json'));assert len(paths)==len(jobs),(stage,len(paths),len(jobs))
 for c,i,g in jobs:
  x=json.loads((OUT/stage/f'{c["id"]}_{i}.json').read_text());assert x['config']==c and x['index']==i and x['gate']==g
  assert x['bytes_emitted']==x['mapped_bytes'] and 0<=x['errors']<=x['payload_bits']
  assert abs(x['symbol_s']*(3*x['mapped_bytes']+18)-x['airtime_s'])<1e-9
  assert x['missing_bit_hypotheses']==0 and x['pruned_states']==0
 for path,digest in m['sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,(stage,path)
 print(stage,len(paths),'verified')
print('PASS: job/configuration, framing, metrics and runtime provenance.')
