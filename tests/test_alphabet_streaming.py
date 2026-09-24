"""Mapping and independent exhaustive likelihood checks for tone-count study."""
from pathlib import Path
import sys
import numpy as np
import pytest
from numpy.testing import assert_allclose
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.alphabet_streaming import AlphabetStreamingReceiver
from gravcomm.streaming_soft import StreamingSoftReceiver
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
A=5.180259430854157e-10

def wave(value,M,L,T,tail=0):
 digits=[value//M**i%M for i in range(L-1,-1,-1)]
 offsets=2*np.array(digits)-(M-1)
 return sampled_cffsk(np.r_[offsets,np.zeros(tail,int)],T,q=.04,transition_fraction=.6)

@pytest.mark.parametrize('M,L,bits,pad',[(2,4,4,0),(3,2,3,0),(4,2,4,0),(5,3,6,0),(6,2,5,0),(7,3,8,0),(8,2,6,0),(9,2,6,1),(5,4,9,3),(7,4,11,1)])
def test_exhaustive_group(M,L,bits,pad):
 r=CausalReceiver();T=1.4;tail=2;values=np.arange(0,2**bits,2**pad)
 means=np.array([r.stream(1).process(wave(int(v),M,L,T,tail)) for v in values])
 y=means[len(values)//3]+r.stationary_noise(means.shape[1],np.random.default_rng(M),A)
 white=r.finite_record(len(y),A)
 scores=np.array([np.sum(abs(white.whiten(y-m))**2) for m in means])
 expected=[]
 for k in range(bits-pad):
  one=((values>>(bits-1-k))&1).astype(bool);expected.append(scores[one].min()-scores[~one].min())
 d=AlphabetStreamingReceiver(r,T,A,1,tail,history_symbols=max(L,4),beam_width=10000,lookahead_symbols=20,q=.08,transition_fraction=.6,tone_count=M,group_symbols=L,group_bits=bits,pad_bits=pad)
 got=d.push(y);st=d.finish();assert st['pruned_states']==0 and st['missing_bit_hypotheses']==0
 assert_allclose(got[0]['llr'][:bits-pad],expected,atol=1e-8)
 assert got[0]['byte']==values[np.argmin(scores)]

def test_seven_tone_equivalence():
 r=CausalReceiver();T=224/162;v=np.array([31,203,17,249]);symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
 mean=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],T,q=.08,transition_fraction=.6))
 y=mean+r.stationary_noise(len(mean),np.random.default_rng(921),A)
 options=dict(history_symbols=4,beam_width=32768,lookahead_symbols=12,q=.08,transition_fraction=.6)
 outputs=[]
 for cls in [StreamingSoftReceiver,AlphabetStreamingReceiver]:
  d=cls(r,T,A,len(v),**options);got=[]
  for j in range(0,len(y),7):got.extend(d.push(y[j:j+7]))
  stats=d.finish();outputs.append(np.concatenate([x['llr'] for x in got]));assert stats['pruned_states']==0
 assert_allclose(*outputs,atol=1e-9)


def test_even_alphabet_across_group_boundaries():
 r=CausalReceiver();T=1.4;tail=2
 labels=np.arange(16);digits=np.column_stack([labels//4,labels%4])
 means=np.array([r.stream(1).process(sampled_cffsk(np.r_[2*v-3,np.zeros(tail,int)],T,q=.04,transition_fraction=.6)) for v in digits])
 y=means[11]+r.stationary_noise(means.shape[1],np.random.default_rng(419),A)
 white=r.finite_record(len(y),A);scores=np.array([np.sum(abs(white.whiten(y-m))**2) for m in means])
 expected=[]
 for bit in range(4):
  one=((labels>>(3-bit))&1).astype(bool);expected.append(scores[one].min()-scores[~one].min())
 d=AlphabetStreamingReceiver(r,T,A,2,tail,history_symbols=4,beam_width=1000,lookahead_symbols=20,q=.08,transition_fraction=.6,tone_count=4,group_symbols=1,group_bits=2)
 got=[]
 for j in range(0,len(y),3):got.extend(d.push(y[j:j+3]))
 assert d.finish()['missing_bit_hypotheses']==0
 assert_allclose(np.concatenate([x['llr'] for x in got]),expected,atol=1e-8)


@pytest.mark.parametrize('M,L,bits,q',[(6,2,5,.064),(8,2,6,.4*.8/7),(9,2,6,.02)])
def test_nondefault_spacing_likelihood(M,L,bits,q):
 r=CausalReceiver();T=1.4;tail=2;values=np.arange(2**bits)
 def mean(v):
  digits=np.array([int(v)//M**j%M for j in range(L-1,-1,-1)])
  return r.stream(1).process(sampled_cffsk(np.r_[2*digits-(M-1),np.zeros(tail,int)],T,q=q/2,transition_fraction=.6))
 means=np.array([mean(v) for v in values]);y=means[11]+r.stationary_noise(means.shape[1],np.random.default_rng(198+M),A)
 white=r.finite_record(len(y),A);scores=np.array([sum(abs(white.whiten(y-m))**2) for m in means]);expected=[]
 for k in range(bits):
  one=((values>>(bits-1-k))&1).astype(bool);expected.append(scores[one].min()-scores[~one].min())
 d=AlphabetStreamingReceiver(r,T,A,1,tail,history_symbols=4,beam_width=1000,lookahead_symbols=20,q=q,transition_fraction=.6,tone_count=M,group_symbols=L,group_bits=bits)
 got=d.push(y);d.finish();assert_allclose(got[0]['llr'],expected,atol=1e-8)
