"""Verify current release summaries, code provenance and optional data archive."""
from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
def check_manifest(path):
 if not path.exists():return
 for line in path.read_text().splitlines():
  digest,rel=line.split('  ',1);assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
for group in ['ber_qualification_v1','ber_replacement_v2']:
 folder=ROOT/'results/fixed_rate_sizing_2026_09_20'/group
 m=json.loads((folder/'manifest.json').read_text())
 for rel,digest in m['source_sha256'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
 f=json.loads((folder/'final.json').read_text())
 for name,row in f['channels'].items():
  assert row['payload_bits']==224*row['packets']
  assert abs(row['ber']-row['errors']/row['payload_bits'])<1e-14
  assert (row['confidence_interval'][1]<=.001) if row['status']=='pass' else (row['confidence_interval'][0]>.001)
check_manifest(ROOT/'MANIFEST.sha256');check_manifest(ROOT/'DATA_MANIFEST.sha256')
print('Release code hashes, data manifests and qualification summaries verified.')
