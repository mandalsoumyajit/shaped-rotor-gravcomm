"""General-alphabet study fork of the frozen streaming reduced-state detector.

Channel/filter history is retained per survivor; paths with matching recent
symbols and phase merge approximately. A bounded graph supplies fixed-lag soft
outputs, without packet reruns. This is a prototype requiring convergence tests.
"""
from collections import deque
import numpy as np
from scipy.signal import sosfilt,sosfilt_zi
from .sequence_decoder import CFFSK
from .compact_innovations import CompactInnovations

class AlphabetStreamingReceiver:
 def __init__(self,receiver,symbol_s,amplitude,nbytes,tail_symbols=6,*,noise_order=64,history_symbols=3,beam_width=4096,lookahead_symbols=12,q=.1,transition_fraction=.4,tone_count=7,group_symbols=3,group_bits=8,pad_bits=0):
  if nbytes<1 or history_symbols<2 or beam_width<1 or lookahead_symbols<0:raise ValueError('Invalid stream configuration')
  from fractions import Fraction
  from math import lcm
  M=int(tone_count);L=int(group_symbols);bits=int(group_bits)
  if M<2 or L<1 or not 0<bits<=12 or 2**bits>M**L or history_symbols<L-1 or not 0<=pad_bits<bits:raise ValueError('Invalid alphabet mapping/history')
  self.M=M;self.L=L;self.bits=bits;self.pad=pad_bits
  # Even alphabets have half-integer offsets; zero is a separate acquired/tail state.
  offsets=[Fraction(2*i-(M-1),2) for i in range(M)]
  if M%2:self.zero=M//2
  else:self.zero=M;offsets.append(Fraction(0))
  self.B=len(offsets);self.offsets=np.array([float(x) for x in offsets])
  qr=Fraction(str(q)).limit_denominator(100000);rr=Fraction(str(transition_fraction)).limit_denominator(100000)
  if abs(float(qr)-q)>1e-12 or abs(float(rr)-transition_fraction)>1e-12:raise ValueError('Nonrational waveform setting')
  increments=[[qr/(1-rr)*(rr*a/2+(1-rr/2)*b) for b in offsets] for a in offsets]
  self.phases=lcm(*(v.denominator for row in increments for v in row))
  if self.phases>100000:raise ValueError('Excessive phase lattice')
  self.inc=np.array([[int(v*self.phases)%self.phases for v in row] for row in increments])
  self.receiver=receiver;self.T=symbol_s;self.nbytes=nbytes;self.tail=tail_symbols
  self.total=L*nbytes+tail_symbols;self.order=noise_order;self.h=history_symbols;self.width=beam_width;self.lag=lookahead_symbols
  self.q=q;self.r=transition_fraction;self.whitener=CompactInnovations(receiver,amplitude,noise_order)
  self.sos=receiver.sos;self.state=sosfilt_zi(self.sos).astype(complex)[:,None,:]
  self.hist=np.array([sum(self.zero*self.B**j for j in range(history_symbols))],np.int64);self.modulus=self.B**history_symbols
  self.phase=np.array([0]);self.cost=np.zeros(1);self.mean_past=np.empty((1,0),complex)
  self.obs_past=np.empty(0,complex);self.graph=deque();self.buffer=np.empty(0,complex)
  self.symbol=0;self.sample=0;self.next_byte=0;self.peak_states=1;self.peak_graph=0;self.peak_state_bytes=0;self.missing=0;self.pruned_states=0
 def _geometry(self,t):
  fs=self.receiver.input_rate_hz;D=self.receiver.decimation
  start=int(np.ceil(t*self.T*fs-1e-12));end=int(np.ceil((t+1)*self.T*fs-1e-12))
  first=(-start)%D
  return start,end,first,len(range(first,end-start,D))
 def push(self,samples):
  x=np.asarray(samples,complex)
  if x.ndim!=1:raise ValueError('One sample stream required')
  outputs=[];position=0
  while position<len(x) or len(self.buffer):
   if self.symbol>=self.total:raise ValueError('Samples beyond configured frame')
   count=self._geometry(self.symbol)[3]
   need=count-len(self.buffer);take=min(need,len(x)-position)
   self.buffer=np.r_[self.buffer,x[position:position+take]];position+=take
   if len(self.buffer)<count:break
   y=self.buffer;self.buffer=np.empty(0,complex)
   outputs.extend(self._step(y))
  return outputs
 def _step(self,y):
  t=self.symbol;start,end,first,count=self._geometry(t);fs=self.receiver.input_rate_hz;D=self.receiver.decimation
  # Fixed rational symbol clocks repeat, but keeping only this block's templates
  # bounds storage even for an irrational/incommensurate clock ratio.
  local=np.arange(start,end)/fs-t*self.T;tr=self.r*self.T;u=np.minimum(local/tr,1)
  integral=tr*(2.5*u**4-3*u**5+u**6)+np.maximum(local-tr,0)
  spacing=self.q/((1-self.r)*self.T)
  before=np.repeat(self.offsets,self.B);after=np.tile(self.offsets,self.B)
  raw=np.exp(2j*np.pi*spacing*(before[:,None]*local+(after-before)[:,None]*integral))
  branch,bstate=sosfilt(self.sos,raw,axis=1,zi=np.zeros((len(self.sos),self.B**2,2),complex))
  dim=2*len(self.sos);initial=np.eye(dim,dtype=complex).reshape(dim,len(self.sos),2).transpose(1,0,2)
  response,transition=sosfilt(self.sos,np.zeros((dim,len(local)),complex),axis=1,zi=initial)
  branch=branch[:,first::D];response=response[:,first::D]
  bstate=bstate.transpose(1,0,2).reshape(self.B**2,dim);transition=transition.transpose(1,0,2).reshape(dim,dim)
  choices=np.array([self.zero]) if t>=self.L*self.nbytes else np.arange(self.M)
  parent=np.repeat(np.arange(len(self.cost)),len(choices));choice=np.tile(choices,len(self.cost));old=self.hist%self.B
  values=choice.copy()
  if t<self.L*self.nbytes:
   position=t%self.L
   # Radix-M prefix reconstructed from the retained base-B history.
   for j in range(position):values+=((self.hist[parent]//self.B**j)%self.B)*self.M**(j+1)
   remaining=self.L-position-1
   # A prefix is valid iff some suffix reaches a permitted group value.
   low=values*self.M**remaining;high=np.minimum((values+1)*self.M**remaining-1,2**self.bits-1)
   stride=2**self.pad if t//self.L==self.nbytes-1 else 1
   valid=((low+stride-1)//stride)*stride<=high
   parent=parent[valid];choice=choice[valid];values=values[valid]
  index=old[parent]*self.B+choice;factor=np.exp(2j*np.pi*self.phase[parent]/self.phases)
  flat=self.state.transpose(1,0,2).reshape(len(self.cost),dim)
  means=(flat@response)[parent]+factor[:,None]*branch[index]
  Wp,Wc=self.whitener.block_matrices(self.sample,count)
  wy=Wp@self.obs_past+Wc@y
  prediction=(self.mean_past@Wp.T)[parent]+means@Wc.T
  edge=np.sum(abs(prediction-wy)**2,axis=1);scores=self.cost[parent]+edge
  next_hist=(self.hist[parent]*self.B+choice)%self.modulus
  next_phase=(self.phase[parent]+self.inc[old[parent],choice])%self.phases
  keys=next_hist*self.phases+next_phase
  order=np.argsort(scores,kind='stable');_,first_key=np.unique(keys[order],return_index=True)
  selection=order[np.sort(first_key)[:self.width]];selected_keys=keys[selection]
  self.pruned_states+=max(0,len(first_key)-self.width)
  sorting=np.argsort(selected_keys);sorted_keys=selected_keys[sorting]
  where=np.searchsorted(sorted_keys,keys);valid=where<len(selection)
  valid &= sorted_keys[np.minimum(where,len(selection)-1)]==keys
  destinations=sorting[np.minimum(where,len(selection)-1)]
  self.graph.append((t,self.cost.copy(),parent[valid],destinations[valid],edge[valid],values[valid]))
  self.cost=scores[selection];self.cost-=min(self.cost)
  updated=(flat@transition)[parent[selection]]+factor[selection,None]*bstate[index[selection]]
  self.state=updated.reshape(-1,len(self.sos),2).transpose(1,0,2)
  self.mean_past=np.concatenate((self.mean_past[parent[selection]],means[selection]),axis=1)[:,-self.order:]
  self.obs_past=np.r_[self.obs_past,y][-self.order:]
  self.hist=next_hist[selection];self.phase=next_phase[selection]
  self.symbol+=1;self.sample+=count;self.peak_states=max(self.peak_states,len(selection));self.peak_graph=max(self.peak_graph,len(self.graph))
  self.peak_state_bytes=max(self.peak_state_bytes,self.storage_bytes())
  outputs=[]
  while self.next_byte<self.nbytes and (t>=self.L*self.next_byte+self.L-1+self.lag or self.symbol==self.total):
   target=self.L*self.next_byte+self.L-1;back=np.zeros(len(self.cost));llr=np.full(self.bits,np.nan);hard=None
   for step,forward,src,dst,metric,words in reversed(self.graph):
    if step==target:
     total=forward[src]+metric+back[dst];hard=int(words[np.argmin(total)])
     for bit in range(self.bits):
      if self.next_byte==self.nbytes-1 and bit>=self.bits-self.pad:
       llr[bit]=np.inf;continue
      labels=(words>>(self.bits-1-bit))&1;zero=total[labels==0];one=total[labels==1]
      if len(zero) and len(one) and np.isfinite(zero.min()) and np.isfinite(one.min()):llr[bit]=one.min()-zero.min()
     break
    prior=np.full(len(forward),np.inf);np.minimum.at(prior,src,metric+back[dst]);back=prior-np.min(prior)
   missing=int(np.count_nonzero(~np.isfinite(llr[:self.bits-self.pad] if self.next_byte==self.nbytes-1 and self.pad else llr)));self.missing+=missing
   outputs.append(dict(byte_index=self.next_byte,byte=hard,llr=llr,available_after_symbol=self.symbol,missing_bit_hypotheses=missing))
   self.next_byte+=1
   while self.graph and self.graph[0][0]<self.L*self.next_byte:self.graph.popleft()
  return outputs
 def storage_bytes(self):
  """Persistent NumPy buffers only; excludes Python overhead and work arrays."""
  arrays=(self.state,self.hist,self.phase,self.cost,self.mean_past,self.obs_past,self.buffer,self.whitener.startup,self.whitener.taps)
  return sum(x.nbytes for x in arrays)+sum(sum(x.nbytes for x in row[1:]) for row in self.graph)
 def finish(self):
  if self.symbol!=self.total or len(self.buffer):raise ValueError('Incomplete frame')
  return dict(symbols=self.symbol,bytes_emitted=self.next_byte,peak_states=self.peak_states,peak_graph_symbols=self.peak_graph,peak_persistent_array_bytes=self.peak_state_bytes,missing_bit_hypotheses=self.missing,pruned_states=self.pruned_states)
