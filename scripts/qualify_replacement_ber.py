"""Frozen-design sequential BER qualification with familywise confidence control."""
import argparse,json,hashlib,time,platform
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from streaming_dynamics_runtime import ROOT,run,save
from gravcomm.packet_confidence import packet_ber_interval
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/ber_replacement_v2'
def qualify(records,alpha):
 errors=[x['errors'] for x in records];n=len(records);lo,hi=packet_ber_interval(errors,224,alpha) if n else (0.,1.)
 return dict(packets=n,payload_bits=224*n,errors=sum(errors),ber=sum(errors)/(224*n) if n else None,confidence_interval=[lo,hi],alpha=alpha,status='pass' if hi<=.001 else 'fail' if lo>.001 else 'running',erasures=sum(not x['accepted'] for x in records),packet_failures=sum(x['failure'] for x in records),wrong_accepted=sum(x['wrong_accepted'] for x in records),accepted_payload_rate_bps=sum(x['accepted'] for x in records)/n if n else None,max_block_s=max((x['max_block_s'] for x in records),default=0),max_latency_s=max((x['latency_budget_s'] for x in records),default=0))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=40);p.add_argument('--batch',type=int,default=128);a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
 frozen=json.loads((OUT/'frozen_designs.json').read_text());cs=frozen['channels'];alpha=.05/len(cs);start=time.time()
 manifest=dict(frozen=frozen,confidence='existing two-sided time-uniform packet-fraction mixture CS',familywise_alpha=.05,per_channel_alpha=alpha,target=.001,batch=a.batch,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/streaming_dynamics_runtime.py',*(ROOT/'src/gravcomm').glob('*.py')]})
 mp=OUT/'manifest.json'
 if mp.exists():assert json.loads(mp.read_text())==manifest,'Inputs changed; use a fresh qualification directory'
 else:save(mp,manifest)
 # Independent convergence diagnostics do not enter BER sample counts.
 with ProcessPoolExecutor(max_workers=min(a.workers,6)) as gate_pool:
  jobs=[(ch['config'],ch['seed_start']-10+i,True) for ch in cs for i in range(2)]
  fs={gate_pool.submit(run,j):j for j in jobs}
  for future in as_completed(fs):
   job=fs[future];result=future.result();save(OUT/'gates'/f'{job[0]["id"]}_{job[1]}.json',result)
   assert result['passes_detector_gate'],'Convergence gate failed'
 print('ALL_CONVERGENCE_GATES_PASS',flush=True)
 records={};statuses={} 
 for channel in cs:
  ident=channel['config']['id'];folder=OUT/ident;folder.mkdir(exist_ok=True);rs=sorted([json.loads(p.read_text()) for p in folder.glob('*.json')],key=lambda x:x['index']);assert all(x['index']==channel['seed_start']+i and x['config']==channel['config'] for i,x in enumerate(rs));records[ident]=rs;statuses[ident]=qualify(rs,alpha)
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  while any(s['status']=='running' for s in statuses.values()):
   jobs=[]
   for channel in cs:
    c=channel['config'];ident=c['id'];n=len(records[ident])
    if statuses[ident]['status']=='running':jobs.extend((c,channel['seed_start']+i,False) for i in range(n,((n//a.batch)+1)*a.batch))
   jobs.sort(key=lambda j:(j[1]%1000000,j[0]['id']))
   futures={pool.submit(run,j):j for j in jobs}
   for f in as_completed(futures):
    j=futures[f];x=f.result()
    if not x['passes_detector_gate']:raise RuntimeError('Receiver quality condition failed; qualification invalid')
    save(OUT/j[0]['id']/f'{j[1]}.json',x);records[j[0]['id']].append(x)
   for c in cs:
    ident=c['config']['id'];records[ident].sort(key=lambda x:x['index']);statuses[ident]=qualify(records[ident],alpha)
   state=dict(channels=statuses,workers_running=any(x['status']=='running' for x in statuses.values()),elapsed_this_run_s=time.time()-start,packets_total=sum(len(x) for x in records.values()))
   save(OUT/'progress.json',state);print(json.dumps(state),flush=True)
 save(OUT/'final.json',dict(channels=statuses,workers_running=False,elapsed_this_run_s=time.time()-start,packets_total=sum(len(x) for x in records.values())));print('QUALIFICATION_COMPLETE',flush=True)
