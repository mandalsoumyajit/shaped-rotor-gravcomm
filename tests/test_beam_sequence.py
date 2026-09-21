import unittest
import numpy as np
from scipy.signal import lfilter
from gravcomm.beam_sequence import BeamByteDetector
from gravcomm.sequence_decoder import CFFSK

class BeamDetectorTests(unittest.TestCase):
    def setUp(self):
        self.model=CFFSK(2.,.4,.1,4)
        self.taps=np.array([.8,.2,.1,.05,.03],complex)
        self.det=BeamByteDetector(self.model,self.taps,100000)
        self.words=np.array([[v//49,v//7%7,v%7] for v in range(256)])-3
        self.signals=np.array([lfilter(self.taps,[1],np.r_[np.ones(4),self.model.waveform(np.r_[w,[0,0]])])[4:] for w in self.words])
        rng=np.random.default_rng(33)
        self.obs=self.signals[173]+.6*(rng.normal(size=20)+1j*rng.normal(size=20))
        self.cost=np.sum(abs(self.signals)**2,axis=1)-2*np.real(self.signals@self.obs.conj())
    def test_exhaustive_hard_and_soft(self):
        r=self.det.decode(self.obs,2)
        np.testing.assert_array_equal(r['symbols'],self.words[np.argmin(self.cost)])
        expected=[]
        for bit in range(8):
            ones=(np.arange(256)>>(7-bit))&1
            expected.append(self.cost[ones==1].min()-self.cost[ones==0].min())
        np.testing.assert_allclose(r['llr'],np.clip(expected,-100,100),atol=1e-11)
        self.assertAlmostEqual(r['path_cost'],self.cost.min(),places=10)
    def test_forced_bit_search_cost(self):
        for bit in range(8):
            for value in (0,1):
                r=self.det.decode(self.obs,2,{bit:value},False)
                allowed=((np.arange(256)>>(7-bit))&1)==value
                self.assertAlmostEqual(r['path_cost'],self.cost[allowed].min(),places=10)
    def test_noiseless_long_memory(self):
        taps=np.r_[self.taps,np.zeros(12)]
        det=BeamByteDetector(self.model,taps,512)
        obs=lfilter(taps,[1],np.r_[np.ones(16),self.model.waveform(np.r_[self.words[173],[0,0]])])[16:]
        np.testing.assert_array_equal(det.decode(obs,2)['symbols'],self.words[173])
    def test_nonzero_long_fir_exhaustive(self):
        taps=np.exp(-np.arange(17)/5)*(1+.1j*np.arange(17))
        det=BeamByteDetector(self.model,taps,100000)
        signals=np.array([lfilter(taps,[1],np.r_[np.ones(16),self.model.waveform(np.r_[w,[0,0]])])[16:] for w in self.words])
        rng=np.random.default_rng(36);obs=signals[91]+.8*(rng.normal(size=20)+1j*rng.normal(size=20))
        costs=np.sum(abs(signals)**2,axis=1)-2*np.real(signals@obs.conj())
        r=det.decode(obs,2)
        np.testing.assert_array_equal(r['symbols'],self.words[np.argmin(costs)])
        expected=[]
        for bit in range(8):
            ones=(np.arange(256)>>(7-bit))&1
            expected.append(costs[ones==1].min()-costs[ones==0].min())
        np.testing.assert_allclose(r['llr'],np.clip(expected,-100,100),atol=1e-10)
if __name__=='__main__':unittest.main()

