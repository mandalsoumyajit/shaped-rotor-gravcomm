import unittest
import numpy as np
from numpy.testing import assert_allclose
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.causal_sequence import CausalBeamByteDetector
from gravcomm.streaming_soft import StreamingSoftReceiver
class WaveformDynamicsTests(unittest.TestCase):
 def test_nondefault_dynamics_exhaustive_byte(self):
  receiver=CausalReceiver();T=224/147;A=5.180259430854157e-10;symbols=np.array([2,-1,3]);rng=np.random.default_rng(515)
  for r,q in ((.25,.125),(.5,.1),(.6,.08),(.6,.05),(.5,.0625),(.4,.075),(.5,.125)):
   with self.subTest(r=r,q=q):
    mean=receiver.stream(1).process(sampled_cffsk(symbols,T,transition_fraction=r,q=q));y=mean+receiver.stationary_noise(len(mean),rng,A)
    reference=CausalBeamByteDetector(receiver,T,A,256,transition_fraction=r,q=q).decode(y,1,0)
    detector=StreamingSoftReceiver(receiver,T,A,1,0,history_symbols=3,beam_width=256,lookahead_symbols=3,transition_fraction=r,q=q)
    got=detector.push(y);detector.finish();assert_allclose(got[0]['llr'],reference['llr'],atol=1e-9)
if __name__=='__main__':unittest.main()
