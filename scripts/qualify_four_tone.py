"""Independent four-tone qualification; frozen inputs and resumable packet records."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import argparse,json,hashlib,time,sys,platform
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import tone_count_runtime as runtime
from tone_count_gates import wider
from resume_tone_count_screen import resilient_save as save
from gravcomm.packet_confidence import packet_ber_interval
ROOT=runtime.ROOT
OUT=ROOT/'results/tone_count_study_2026_09_21/four_tone_qualification_v1'

def execute(job):
 runtime.StreamingSoftReceiver=wider
 return runtime.run(job)

def summarize(rs):
 counts=[x['errors'] for x in rs];n=len(rs)
 ci=packet_ber_interval(counts,224,.025) if n else [0.,1.]
 quality=all(x['passes_detector_gate'] for x in rs)
 latency=max((x['latency_budget_s'] for x in rs),default=0)
 return dict(packets=n,payload_bits=224*n,errors=sum(counts),ber=sum(counts)/(224*n) if n else None,confidence_interval=ci,alpha=.025,status='invalid' if not quality else 'latency_failed' if latency>=300 else 'pass' if ci[1]<=.001 else 'fail' if ci[0]>.001 else 'running',crc_rejections=sum(not x['accepted'] for x in rs),wrong_accepted=sum(x['wrong_accepted'] for x in rs),accepted_payload_rate_bps=sum(x['accepted'] for x in rs)/n if n else None,max_latency_s=latency,max_block_s=max((x['max_block_s'] for x in rs),default=0))

def main():
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=12);p.add_argument('--max-packets',type=int,default=32768);a=p.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 configs=[json.loads((runtime.OUT/'confirm'/f'M4_g1b2_spanT0.6_f{f}.manifest.json').read_text())['config'] for f in (24,18)]
 deps=[Path(__file__),ROOT/'scripts/tone_count_runtime.py',ROOT/'scripts/tone_count_gates.py',ROOT/'scripts/resume_tone_count_screen.py',*(ROOT/'src/gravcomm').glob('*.py')]
 manifest=dict(configs=configs,seed_start=10000000,gate_seeds=[9999990,9999991],batch=128,max_packets=a.max_packets,target_ber=.001,latency_limit_s=300,familywise_alpha=.05,per_channel_alpha=.025,confidence='time-uniform packet-fraction mixture; all tentative payload errors including rejected packets',sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
 mp=OUT/'manifest.json'
 if mp.exists():assert json.loads(mp.read_text())==manifest,'Frozen inputs changed'
 else:
  save(mp,manifest)
  for dep in deps:
   dest=OUT/'code_snapshot'/dep.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(dep.read_bytes())
 save(OUT/'environment.json',dict(python=sys.version,platform=platform.platform(),workers=a.workers,logical_cpus=os.cpu_count()))
 start=time.time()
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  gate_jobs=[(c,i,True) for c in configs for i in manifest['gate_seeds']]
  for job in gate_jobs:
   dest=OUT/'gates'/f'{job[0]["id"]}_{job[1]}.json'
   x=json.loads(dest.read_text()) if dest.exists() else pool.submit(execute,job).result()
   save(dest,x)
   if not x['passes_detector_gate']:
    save(OUT/'progress.json',dict(status='blocked_detector_convergence',config=job[0],seed=job[1]));raise RuntimeError('Fresh convergence gate failed')
  print('FRESH_CONVERGENCE_GATES_PASS',flush=True)
  records={}
  for c in configs:
   folder=OUT/c['id'];folder.mkdir(exist_ok=True)
   rs=[json.loads(p.read_text()) for p in folder.glob('*.json')]
   assert all(x['config']==c and 10000000<=x['index']<10000000+a.max_packets for x in rs)
   records[c['id']]={x['index']:x for x in rs}
  while True:
   states={};jobs=[]
   for c in configs:
    rows=records[c['id']];n=0
    while 10000000+n in rows:n+=1
    # Decisions use deterministic complete batches, including on restart.
    prefix=n//128*128
    states[c['id']]=summarize([rows[10000000+i] for i in range(prefix)])
    if states[c['id']]['status']=='running':
     if prefix>=a.max_packets:states[c['id']]['status']='inconclusive_budget'
     else:jobs.extend((c,10000000+i,False) for i in range(prefix,min(prefix+128,a.max_packets)) if 10000000+i not in rows)
   state=dict(status='running' if jobs else 'complete',channels=states,packets_saved=sum(len(x) for x in records.values()),elapsed_this_run_s=time.time()-start,workers=a.workers)
   save(OUT/'progress.json',state);print(json.dumps(state),flush=True)
   if not jobs:
    save(OUT/'final.json',state);break
   fs={pool.submit(execute,j):j for j in jobs}
   for future in as_completed(fs):
    job=fs[future];x=future.result();save(OUT/job[0]['id']/f'{job[1]}.json',x);records[job[0]['id']][job[1]]=x
    if not x['passes_detector_gate']:raise RuntimeError('Packet detector quality failed')
if __name__=='__main__':main()
