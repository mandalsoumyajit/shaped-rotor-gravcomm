"""Streaming framed-link screening. Pilot results are not BER qualification."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,binascii,argparse,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.streaming_soft import StreamingSoftReceiver
from gravcomm.coded_sequence import ShortLDPC
from gravcomm import convolutional_candidates as tb
from gravcomm.packet_confidence import packet_ber_interval
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/streaming_link_search_v1'
A0=5.180259430854157e-10
SCHEMES=('uncoded','tb_3/4','tb_2/3','tb_1/2','ldpc')
def save(path,x):
 path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def crc(x):return binascii.crc_hqx(bytes(x),0xffff)
def layout(scheme,k):
 coded=k if scheme=='uncoded' else 2*k if scheme=='ldpc' else tb.metadata(k,scheme[3:])['n']
 nbytes=(coded+7)//8;payload=k-32;ts=payload/(3*nbytes+18)
 return dict(coded_bits=coded,mapped_bytes=nbytes,payload_bits=payload,symbol_s=ts,airtime_s=payload)
def run(job):
 c,index,gate=job;start=time.perf_counter();k=c['information_bits'];scheme=c['scheme'];spec=layout(scheme,k)
 rng=np.random.default_rng(np.random.SeedSequence([20370930,k,index]))
 payload=rng.integers(0,256,(k-32)//8,dtype=np.uint8);body=b'\x00\x01'+bytes(payload)
 info=np.unpackbits(np.frombuffer(body+crc(body).to_bytes(2,'big'),dtype=np.uint8))
 ldpc=ShortLDPC(k) if scheme=='ldpc' else None
 word=info if scheme=='uncoded' else ldpc.encode(info) if ldpc else tb.encode(info,scheme[3:])
 permutation=np.random.default_rng(4000+len(word)).permutation(len(word)) if scheme!='uncoded' else np.arange(len(word))
 packed=np.packbits(word[permutation]).astype(int)
 symbols=np.column_stack((packed//49,packed//7%7,packed%7)).ravel()-3
 amplitude=A0*c['amplitude_scale'];r=CausalReceiver(carrier_hz=c['carrier_hz'],background=SeismicBackground())
 ts=spec['symbol_s'];mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],ts,q=c.get('q',.1)))
 # Separate payload and noise seeds, with common standardized noise across amplitudes/codes.
 noise_rng=np.random.default_rng(np.random.SeedSequence([20370931,k,index]))
 observed=mean+r.stationary_noise(len(mean),noise_rng,amplitude)
 if ldpc:ldpc.decode(np.ones(len(word)))
 elif scheme!='uncoded':tb.decode(np.ones(len(word)),scheme[3:],k)
 def detect(history,lag):
  setup=time.perf_counter();d=StreamingSoftReceiver(r,ts,amplitude,len(packed),history_symbols=history,beam_width=100000 if history==5 else 16384,lookahead_symbols=lag,q=c.get('q',.1))
  setup=time.perf_counter()-setup;outputs=[];position=0;times=[];budget=[]
  for symbol in range(len(symbols)+6):
   count=d._geometry(symbol)[3];tick=time.perf_counter();outputs.extend(d.push(observed[position:position+count]));times.append(time.perf_counter()-tick);budget.append(count/4);position+=count
  stats=d.finish();channel=np.concatenate([x['llr'] for x in outputs])[:len(word)]
  if not np.all(np.isfinite(channel)):raise RuntimeError('Missing soft hypotheses')
  llr=np.empty(len(word));llr[permutation]=channel
  tick=time.perf_counter();syndrome=True;iterations=0
  if scheme=='uncoded':recovered=(llr<0).astype(np.uint8)
  elif ldpc:recovered,syndrome,iterations=ldpc.decode(llr);recovered=recovered[:k]
  else:recovered=tb.decode(llr,scheme[3:],k)
  fec_s=time.perf_counter()-tick;data=np.packbits(recovered)
  accepted=bool(syndrome and bytes(data[:2])==b'\x00\x01' and crc(data[:-2])==int.from_bytes(bytes(data[-2:]),'big'))
  errors=int(np.unpackbits(data[2:-2]^payload).sum())
  # Replay backlog under the measured execution durations. Parallel host load is
  # reported; this is not a deployment worst-case timing certification.
  finish=0.;arrival=0.;maxbacklog=0.
  for elapsed,interval in zip(times,budget):
   arrival+=interval;finish=max(finish,arrival)+elapsed;maxbacklog=max(maxbacklog,finish-arrival)
  result=dict(errors=errors,raw_errors=int(np.count_nonzero((llr<0)!=word)),accepted=accepted,wrong_accepted=bool(accepted and errors),failure=bool(errors or not accepted),fec_iterations=int(iterations),fec_s=fec_s,setup_s=setup,total_detector_s=sum(times),max_block_s=max(times),max_backlog_s=maxbacklog,final_processing_s=finish-arrival+fec_s,latency_budget_s=spec['airtime_s']+finish-arrival+fec_s,whitening=d.whitener.spectral_error(),**stats)
  return result,llr,recovered
 result,llr,recovered=detect(4,12);convergence=[]
 if gate:
  for history,lag in ((5,12),(4,24)):
   other,metrics,decoded=detect(history,lag)
   comparison=dict(history_symbols=history,lookahead_symbols=lag,relative_llr_change=float(np.linalg.norm(metrics-llr)/max(np.linalg.norm(metrics),1e-30)),max_llr_change=float(max(abs(metrics-llr))),sign_changes=int(np.count_nonzero(np.signbit(metrics)!=np.signbit(llr))),decoded_information_changes=int(np.count_nonzero(decoded!=recovered)),errors=other['errors'],pruned_states=other['pruned_states'],missing=other['missing_bit_hypotheses'])
   comparison['passes']=bool(comparison['relative_llr_change']<.01 and comparison['decoded_information_changes']==0 and comparison['sign_changes']==0 and comparison['pruned_states']==0 and comparison['missing']==0)
   convergence.append(comparison)
 return dict(config=c,index=index,gate=gate,**spec,**result,convergence=convergence,passes_detector_gate=bool((not gate or all(x['passes'] for x in convergence)) and result['whitening']['max_psd_error']<.001 and result['pruned_states']==0 and result['missing_bit_hypotheses']==0),elapsed_s=time.perf_counter()-start)
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['gates','screen'],required=True);p.add_argument('--workers',type=int,default=12);p.add_argument('--packets',type=int,default=2);args=p.parse_args();OUT.mkdir(exist_ok=True,parents=True)
 jobs=[]
 if args.stage=='gates':
  for k in (128,256):
   for scheme in SCHEMES:
    c=dict(id=f'g_{k}_{scheme.replace("/","_")}',scheme=scheme,information_bits=k,carrier_hz=20.,amplitude_scale=1.,q=.1);jobs.append((c,0,True))
 else:
  gates=[json.loads(p.read_text()) for p in (OUT/'gates').glob('*.json')]
  eligible={(x['config']['scheme'],x['config']['information_bits']) for x in gates if x['passes_detector_gate']}
  if not eligible:raise RuntimeError('No gate-approved protocol')
  for scheme,k in sorted(eligible):
   for fc in (12.,16.,24.,40.):
    for scale in (.5,1.,2.,4.):
     c=dict(id=f'k{k}_{scheme.replace("/","_")}_f{fc:g}_a{scale:g}',scheme=scheme,information_bits=k,carrier_hz=fc,amplitude_scale=scale,q=.1)
     jobs.extend((c,index,False) for index in range(args.packets))
 folder=OUT/args.stage;folder.mkdir(exist_ok=True);alljobs=jobs;pending=[j for j in jobs if not (folder/f'{j[0]["id"]}_{j[1]}.json').exists()];done=len(jobs)-len(pending)
 provenance={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*(ROOT/'src/gravcomm').glob('*.py')]}
 manifest=dict(stage=args.stage,jobs=len(jobs),packets=args.packets,source_sha256=provenance,scope='Pilot only. Assumed acquisition, fixed instrumental profile and declared environmental scenario; no BER qualification.')
 path=OUT/f'{args.stage}_manifest.json'
 if path.exists():assert json.loads(path.read_text())==manifest,'Changed inputs require fresh directory'
 else:save(path,manifest)
 with ProcessPoolExecutor(max_workers=args.workers) as pool:
  futures={pool.submit(run,j):j for j in pending}
  for f in as_completed(futures):
   j=futures[f]
   try:record=f.result()
   except Exception as e:
    save(OUT/f'failure_{j[0]["id"]}_{j[1]}.json',dict(error=repr(e),job=j));raise
   save(folder/f'{j[0]["id"]}_{j[1]}.json',record);done+=1
   save(OUT/'progress.json',dict(stage=args.stage,completed=done,total=len(alljobs),status='complete' if done==len(alljobs) else 'running'))
   print(json.dumps(dict(done=done,total=len(alljobs),id=record['config']['id'],errors=record['errors'],gate=record['passes_detector_gate'],elapsed_s=record['elapsed_s'])),flush=True)
 print('STAGE_COMPLETE',args.stage,flush=True)
if __name__=='__main__':main()
