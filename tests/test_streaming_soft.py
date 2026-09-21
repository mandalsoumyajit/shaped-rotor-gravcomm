import unittest
import numpy as np
from numpy.testing import assert_allclose,assert_array_equal
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.causal_sequence import CausalBeamByteDetector
from gravcomm.compact_innovations import CompactInnovations
from gravcomm.streaming_soft import StreamingSoftReceiver
from gravcomm.environmental_noise import SeismicBackground
A=5.180259430854157e-10
class StreamingTests(unittest.TestCase):
 def test_startup_and_stationary_whitening(self):
  r=CausalReceiver(carrier_hz=20,background=SeismicBackground());c=CompactInnovations(r,A,64)
  Wp,Wc=c.block_matrices(0,60)
  x=np.random.default_rng(3).normal(size=60)
  assert_allclose(Wc@x,r.finite_record(60,A).whiten(x),atol=1e-10)
  self.assertLess(c.spectral_error()['max_psd_error'],1e-3)
 def test_exhaustive_one_byte(self):
  r=CausalReceiver();v=141;T=224/114
  symbols=np.array([v//49,v//7%7,v%7])-3
  mean=r.stream(1).process(sampled_cffsk(symbols,T));y=mean+r.stationary_noise(len(mean),np.random.default_rng(18),A)
  exact=CausalBeamByteDetector(r,T,A,256).decode(y,1,0)
  d=StreamingSoftReceiver(r,T,A,1,0,history_symbols=3,beam_width=256,lookahead_symbols=3)
  got=d.push(y);d.finish()
  self.assertEqual(got[0]['byte'],exact['bytes'][0]);assert_allclose(got[0]['llr'],exact['llr'],atol=1e-9)
 def test_chunks_future_independence_and_bounded_graph(self):
  r=CausalReceiver();n=20;T=2.;rng=np.random.default_rng(14);v=rng.integers(0,256,n)
  symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
  y=r.stream(1).process(sampled_cffsk(np.r_[symbols,np.zeros(6,int)],T))
  def run(x,chunk):
   d=StreamingSoftReceiver(r,T,A,n,history_symbols=2,beam_width=256,lookahead_symbols=6);out=[]
   for i in range(0,len(x),chunk):out.extend(d.push(x[i:i+chunk]))
   return out,d.finish()
  a,stats=run(y,17);b,_=run(y,10000)
  assert_allclose(np.concatenate([x['llr'] for x in a]),np.concatenate([x['llr'] for x in b]),equal_nan=True)
  changed=y.copy();changed[240:]+=100;future,_=run(changed,13)
  for x,z in zip(a,future):
   if x['available_after_symbol']*8<=240:assert_allclose(x['llr'],z['llr'],equal_nan=True)
  self.assertEqual(len(a),n);self.assertLessEqual(stats['peak_graph_symbols'],9)
  self.assertLessEqual(stats['peak_states'],256)
 def test_two_byte_unmerged_graph_matches_full_reference(self):
  r=CausalReceiver();T=.5;v=np.array([31,203]);symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
  mean=r.stream(1).process(sampled_cffsk(symbols,T));y=mean+r.stationary_noise(len(mean),np.random.default_rng(915),A)
  reference=CausalBeamByteDetector(r,T,A,65536).decode(y,2,0)
  d=StreamingSoftReceiver(r,T,A,2,0,history_symbols=6,beam_width=65536,lookahead_symbols=6)
  got=d.push(y);d.finish()
  assert_allclose(np.concatenate([x['llr'] for x in got]),reference['llr'],atol=1e-9)
if __name__=='__main__':unittest.main()
