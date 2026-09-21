"""Audit saved streaming search records and runtime provenance."""
import json,hashlib
from pathlib import Path
from streaming_link_search import OUT,ROOT
counts={'gates':10,'screen':256,'candidate_gates':12,'refinement':352,'candidate_gates2':4,'refinement2':128}
for folder,expected in counts.items():
 files=list((OUT/folder).glob('*.json'));assert len(files)==expected,(folder,len(files))
 for p in files:
  x=json.loads(p.read_text());assert 0<=x['errors']<=x['payload_bits'];assert x['bytes_emitted']==x['mapped_bytes'];assert x['missing_bit_hypotheses']==0;assert x['pruned_states']==0
  assert x['payload_bits']/x['airtime_s']==1
  assert abs(x['symbol_s']*(3*x['mapped_bytes']+18)-x['airtime_s'])<1e-9
for name in ('gates_manifest.json','screen_manifest.json','candidate_gates_manifest.json','refinement_manifest.json','candidate_gates2_manifest.json','refinement2_manifest.json'):
 x=json.loads((OUT/name).read_text())
 for relative,expected in x.get('source_sha256',x.get('sha256',{})).items():assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,(name,relative)
sizes=json.loads((OUT/'amplitude_size_map.json').read_text());assert len(sizes)==18
assert len({(x['distance_m'],x['amplitude_scale']) for x in sizes})==18
for x in sizes:assert x['field_refinement_relative']<1e-6 and x['radial_gap_m']>=.05
print('PASS: 736 packet trials, 26 convergence records, 18 field maps, framing and runtime hashes.')
