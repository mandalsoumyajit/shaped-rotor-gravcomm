import unittest
import numpy as np
from numpy.testing import assert_allclose,assert_array_equal
from scipy.linalg import solve
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.causal_sequence import CausalBeamByteDetector
A=5.180259430854157e-10
class CausalSequenceTests(unittest.TestCase):
 def test_exhaustive_hard_soft_and_cost(self):
  r=CausalReceiver()
  for duration in (2.,256/126):
   means=[]
   for b in range(256):
    symbols=np.r_[np.array([b//49,b//7%7,b%7])-3,np.zeros(2,int)]
    means.append(r.stream(1).process(sampled_cffsk(symbols,duration)))
   means=np.asarray(means).T
   y=means[:,137]+r.stationary_noise(len(means),np.random.default_rng(721),A)
   model=r.finite_record(len(y),A);residual=y[:,None]-means
   costs=np.real(np.sum(residual.conj()*solve(model.covariance,residual,assume_a='her'),axis=0))
   got=CausalBeamByteDetector(r,duration,A,256).decode(y,1,2)
   self.assertTrue(got['exact_search']);self.assertEqual(got['bytes'][0],np.argmin(costs))
   assert_allclose(got['survivor_costs'],costs[got['survivor_bytes'][:,0]],atol=1e-9)
   for bit in range(8):
    mask=(np.arange(256)>>(7-bit))&1
    self.assertAlmostEqual(got['llr'][bit],costs[mask==1].min()-costs[mask==0].min(),places=8)
 def test_pruning_and_constraints_are_explicit(self):
  r=CausalReceiver();symbols=np.r_[np.array([91//49,91//7%7,91%7])-3,np.zeros(2,int)]
  y=r.stream(1).process(sampled_cffsk(symbols,2.))
  d=CausalBeamByteDetector(r,2.,A,1);got=d.decode(y,1,2)
  self.assertFalse(got['exact_search']);self.assertEqual(got['missing_bit_hypotheses'],8)
  self.assertTrue(np.all(np.isnan(got['llr'])))
  other=d.decode(y,1,2,forced_bits={0:1})
  self.assertEqual(int(other['bytes'][0])>>7,1)
 def test_explicit_prepacket_state_and_sample_phase(self):
  r=CausalReceiver();duration=256/126
  # Build a known zero-tone prehistory with nontrivial filter transient and clock phase.
  stream=r.stream(0);stream.process(np.ones(19,complex));zi=stream.state.copy()
  b=183;symbols=np.r_[np.array([b//49,b//7%7,b%7])-3,np.zeros(2,int)]
  y=stream.process(sampled_cffsk(symbols,duration))
  got=CausalBeamByteDetector(r,duration,A,256).decode(y,1,2,initial_filter_state=zi,sample_offset=19)
  assert_array_equal(got['bytes'],[b]);self.assertLess(got['path_cost'],1e-20)
  assert_allclose(got['final_mean_filter_state'],stream.state,atol=1e-12)
 def test_forced_soft_reference(self):
  r=CausalReceiver();b=73
  symbols=np.r_[np.array([b//49,b//7%7,b%7])-3,np.zeros(2,int)]
  y=r.stream(1).process(sampled_cffsk(symbols,2.))
  exact=CausalBeamByteDetector(r,2.,A,256)
  expected=exact.decode(y,1,2);got=exact.decode_soft(y,1,2)
  assert_allclose(got['llr'],expected['llr'],atol=1e-10)
  self.assertTrue(got['exact_search'])
  approximate=CausalBeamByteDetector(r,2.,A,1).decode_soft(y,1,2)
  self.assertTrue(np.all(np.isfinite(approximate['llr'])))
  self.assertFalse(approximate['exact_search']);self.assertEqual(approximate['search_runs'],9)
 def test_diverse_search_matches_exact_and_preserves_alternatives(self):
  r=CausalReceiver();b=137;duration=256/126
  symbols=np.r_[np.array([b//49,b//7%7,b%7])-3,np.zeros(2,int)]
  y=r.stream(1).process(sampled_cffsk(symbols,duration))
  y+=r.stationary_noise(len(y),np.random.default_rng(737),A)
  exact=CausalBeamByteDetector(r,duration,A,256)
  reference=exact.decode(y,1,2)
  got=exact.decode_diverse(y,1,2,bit_hypotheses=1)
  assert_allclose(got['llr'],reference['llr'],atol=1e-10)
  assert_allclose(got['path_cost'],reference['path_cost'],atol=1e-10)
  sparse=CausalBeamByteDetector(r,duration,A,1).decode_diverse(y,1,2,bit_hypotheses=1)
  self.assertEqual(sparse['missing_bit_hypotheses'],0)
  self.assertFalse(sparse['exact_search']);self.assertGreater(sparse['peak_survivors'],1)
 def test_diverse_multibyte_background_metric(self):
  from gravcomm.environmental_noise import SeismicBackground
  r=CausalReceiver(carrier_hz=20,background=SeismicBackground())
  data=np.array([41,208,17]);symbols=np.column_stack((data//49,data//7%7,data%7)).ravel()-3
  duration=.5;signal=sampled_cffsk(np.r_[symbols,np.zeros(6,int)],duration)
  y=r.stream(1).process(signal)+r.stationary_noise(len(signal[::4]),np.random.default_rng(152),A)
  got=CausalBeamByteDetector(r,duration,A,8).decode_diverse(y,3,bit_hypotheses=2)
  self.assertEqual(got['missing_bit_hypotheses'],0)
  means=r.stream(1).process(sampled_cffsk(np.r_[got['symbols'],np.zeros(6,int)],duration))
  self.assertAlmostEqual(got['path_cost'],r.finite_record(len(y),A).metric(y,means),places=9)
 def test_checkpoint_replay_matches_full_forced_search(self):
  r=CausalReceiver();data=np.array([73,189]);duration=.5
  symbols=np.column_stack((data//49,data//7%7,data%7)).ravel()-3
  signal=sampled_cffsk(np.r_[symbols,np.zeros(2,int)],duration)
  y=r.stream(1).process(signal)+r.stationary_noise(len(signal[::4]),np.random.default_rng(381),A)
  d=CausalBeamByteDetector(r,duration,A,16)
  reference=d.decode_soft(y,2,2)
  for workers in (1,2):
   got=d.decode_checkpoint_soft(y,2,2,workers=workers)
   assert_allclose(got['llr'],reference['llr'],atol=1e-9)
   assert_allclose(got['path_cost'],reference['path_cost'],atol=1e-9)
   assert_array_equal(got['bytes'],reference['bytes'])
 def test_cached_linear_branches_match_samplewise_filtering(self):
  r=CausalReceiver();duration=224/114;data=np.array([51,206,27])
  v=np.column_stack((data//49,data//7%7,data%7)).ravel()-3
  stream=r.stream(0);stream.process(np.ones(19,complex));zi=stream.state.copy()
  signal=sampled_cffsk(np.r_[v,np.zeros(2,int)],duration)
  y=stream.process(signal)+r.stationary_noise(len(signal[1::4]),np.random.default_rng(814),A)
  kwargs=dict(initial_filter_state=zi,sample_offset=19)
  direct=CausalBeamByteDetector(r,duration,A,64,cached_branches=False).decode_checkpoint_soft(y,3,2,**kwargs)
  cached=CausalBeamByteDetector(r,duration,A,64,cached_branches=True).decode_checkpoint_soft(y,3,2,**kwargs)
  assert_array_equal(cached['bytes'],direct['bytes'])
  assert_allclose(cached['llr'],direct['llr'],atol=1e-8,rtol=1e-9)
  self.assertAlmostEqual(cached['path_cost'],direct['path_cost'],places=8)
if __name__=='__main__':unittest.main()
