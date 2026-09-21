"""Audit returned Exxact results and repeat two noisy configurations locally."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import json,collections
from pathlib import Path
from joint_link_search import run_config,source_hashes
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/joint_link_study_2026_09_20'
run=OUT/'exxact_crossed_v1'
m=json.loads((run/'manifest.json').read_text())
normalize=lambda hashes:{k.replace(chr(92),'/'):v for k,v in hashes.items()}
assert normalize(m['source_sha256'])==normalize(source_hashes())
rows=[];totals=collections.defaultdict(lambda:dict(configurations=0,packets=0,zero_failure_configurations=0))
for c in m['configs']:
    r=json.loads((run/('case_'+c['id']+'.json')).read_text())
    assert r['config']==c and r['packets']==len(r['records'])==50
    assert r['mechanics']['feasible']
    assert r['coded_bits']+r['padding_bits']==8*((r['coded_bits']+7)//8)
    assert all(0<=x['errors']<=c['information_bits'] and 0<=x['raw_errors']<=r['coded_bits'] for x in r['records'])
    assert abs(r['packet_failure_rate']-sum(x['failure'] for x in r['records'])/50)<1e-14
    t=totals[c['case']['distance_m']];t['configurations']+=1;t['packets']+=50
    t['zero_failure_configurations']+=r['packet_failure_rate']==0
    rows.append(r)
checks=[]
for scheme in ('ldpc','tb_3/4'):
    candidates=[r for r in rows if r['config']['scheme']==scheme and 0<r['packet_failure_rate']<1]
    assert candidates
    original=candidates[0];again=run_config(original['config'],4,{})
    assert original['records']==again['records'], 'Local/remote packet outcomes differ'
    checks.append(dict(id=original['config']['id'],scheme=scheme,packets=50,packet_records_identical=True,
                       packet_failure_rate=original['packet_failure_rate']))
report=dict(status='passed',source_hashes_match=True,configurations=len(rows),total_packets=sum(r['packets'] for r in rows),
            by_distance=totals,local_remote_checks=checks,
            caveat='Agreement of these two configurations is not a cross-platform equivalence proof for all future runs. No rare-event reliability qualification.')
(OUT/'exxact_results_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))