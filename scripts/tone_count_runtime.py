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
from gravcomm.alphabet_streaming import AlphabetStreamingReceiver as StreamingSoftReceiver
from gravcomm.coded_sequence import ShortLDPC
from gravcomm import convolutional_candidates as tb
from gravcomm.packet_confidence import packet_ber_interval
OUT=ROOT/'results/tone_count_study_2026_09_21/bounded_v2'
A0=5.180259430854157e-10
SCHEMES=('uncoded','tb_3/4','tb_2/3','tb_1/2','ldpc')
def save(path,x):
 path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def crc(x):return binascii.crc_hqx(bytes(x),0xffff)
def layout(scheme,k,M,L,bits):
 coded=tb.metadata(k,scheme[3:])['n'];groups=(coded+bits-1)//bits;payload=k-32
 ts=payload/(L*groups+18)
 return dict(coded_bits=coded,mapped_bytes=groups,mapped_groups=groups,payload_bits=payload,symbol_s=ts,airtime_s=payload,padding_bits=groups*bits-coded,modulated_symbols=L*groups,tone_span_hz=None)
def run(job):
 c,index,gate=job;start=time.perf_counter();k=c['information_bits'];scheme=c['scheme'];M=c['tone_count'];L=c['group_symbols'];bits=c['group_bits'];spec=layout(scheme,k,M,L,bits);spec['tone_span_hz']=(M-1)*c['q']/((1-c['transition_fraction'])*spec['symbol_s'])
 rng=np.random.default_rng(np.random.SeedSequence([20370930,k,index]))
 payload=rng.integers(0,256,(k-32)//8,dtype=np.uint8);body=b'\x00\x01'+bytes(payload)
 info=np.unpackbits(np.frombuffer(body+crc(body).to_bytes(2,'big'),dtype=np.uint8))
 ldpc=ShortLDPC(k) if scheme=='ldpc' else None
 word=info if scheme=='uncoded' else ldpc.encode(info) if ldpc else tb.encode(info,scheme[3:])
 permutation=np.random.default_rng(4000+len(word)).permutation(len(word)) if scheme!='uncoded' else np.arange(len(word))
 padded=np.pad(word[permutation],(0,spec['padding_bits']))
 packed=(padded.reshape(-1,bits)@(2**np.arange(bits-1,-1,-1))).astype(int)
 digits=np.column_stack([packed//M**j%M for j in range(L-1,-1,-1)]).ravel()
 symbols=2*digits-(M-1)
 amplitude=A0*c['amplitude_scale'];r=CausalReceiver(carrier_hz=c['carrier_hz'],background=SeismicBackground())
 ts=spec['symbol_s'];mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],ts,q=c['q']/2,transition_fraction=c.get('transition_fraction',.4)))
 # Separate payload and noise seeds, with common standardized noise across amplitudes/codes.
 noise_rng=np.random.default_rng(np.random.SeedSequence([20370931,k,index]))
 observed=mean+r.stationary_noise(len(mean),noise_rng,amplitude)
 if ldpc:ldpc.decode(np.ones(len(word)))
 elif scheme!='uncoded':tb.decode(np.ones(len(word)),scheme[3:],k)
 def detect(history,lag):
  setup=time.perf_counter();d=StreamingSoftReceiver(r,ts,amplitude,len(packed),history_symbols=history,beam_width=c.get("beam_width",200000),lookahead_symbols=lag,q=c.get('q',.1),transition_fraction=c.get('transition_fraction',.4),tone_count=M,group_symbols=L,group_bits=bits,pad_bits=spec['padding_bits'])
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
 baseT=224/162
 history=max(L-1,4)
 lag=int(np.ceil(12*baseT/spec['symbol_s']-1e-10))
 result,llr,recovered=detect(history,lag);convergence=[]
 spec.update(history_symbols=history,lookahead_symbols=lag)
 if gate:
  for history,lag in ((history+1,lag),(history,2*lag)):
   other,metrics,decoded=detect(history,lag)
   comparison=dict(history_symbols=history,lookahead_symbols=lag,relative_llr_change=float(np.linalg.norm(metrics-llr)/max(np.linalg.norm(metrics),1e-30)),max_llr_change=float(max(abs(metrics-llr))),sign_changes=int(np.count_nonzero(np.signbit(metrics)!=np.signbit(llr))),decoded_information_changes=int(np.count_nonzero(decoded!=recovered)),errors=other['errors'],pruned_states=other['pruned_states'],missing=other['missing_bit_hypotheses'])
   comparison['passes']=bool(comparison['relative_llr_change']<.01 and comparison['decoded_information_changes']==0 and comparison['sign_changes']==0 and comparison['pruned_states']==0 and comparison['missing']==0)
   convergence.append(comparison)
 return dict(config=c,index=index,gate=gate,peak_torque_per_inertia=float(15*np.pi*(M-1)*c['q']/(8*c['transition_fraction']*(1-c['transition_fraction'])*spec['symbol_s']**2)),**spec,**result,convergence=convergence,passes_detector_gate=bool((not gate or all(x['passes'] for x in convergence)) and result['whitening']['max_psd_error']<.001 and result['pruned_states']==0 and result['missing_bit_hypotheses']==0),elapsed_s=time.perf_counter()-start)