import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys, unittest
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from gravcomm.coded_sequence import conv_encode,conv_decode,ShortLDPC,MaxLogByteDetector
from gravcomm.whitened_viterbi import WhitenedViterbi
from gravcomm.sequence_decoder import CFFSK
from desktop_cf_fsk import encode

class CodingChecks(unittest.TestCase):
    def test_convolutional_impulse_and_roundtrip(self):
        # Coefficients read MSB-to-LSB in time-reversed generator order.
        np.testing.assert_array_equal(conv_encode([1]),[1,1,0,1,0,0,1,1,1,1,1,0,1,1])
        rng=np.random.default_rng(55)
        for _ in range(20):
            data=rng.integers(0,2,64,dtype=np.uint8);word=conv_encode(data)
            llr=8*(1.-2*word);llr[13]*=-1
            np.testing.assert_array_equal(conv_decode(llr),data)

    def test_ldpc_published_matrices_and_correction(self):
        code=ShortLDPC()
        self.assertEqual(np.linalg.matrix_rank(code.g.astype(float)),64)
        np.testing.assert_array_equal((code.g@code.h.T)%2,np.zeros((64,64)))
        rng=np.random.default_rng(43)
        for _ in range(20):
            data=rng.integers(0,2,64,dtype=np.uint8);word=code.encode(data)
            llr=8*(1.-2*word);llr[19]*=-1
            decoded,valid,_=code.decode(llr)
            self.assertTrue(valid);np.testing.assert_array_equal(decoded,word)

    def test_longer_published_ldpc_blocks(self):
        rng=np.random.default_rng(92)
        for k in (128,256):
            code=ShortLDPC(k)
            np.testing.assert_array_equal((code.g@code.h.T)%2,np.zeros((k,k)))
            for _ in range(5):
                data=rng.integers(0,2,k,dtype=np.uint8);word=code.encode(data)
                llr=8*(1.-2*word);llr[[7,31]]*=-1
                recovered,valid,_=code.decode(llr)
                self.assertTrue(valid);np.testing.assert_array_equal(recovered,word)

    def test_maxlog_against_exhaustive_waveform_likelihoods(self):
        model=CFFSK(2.,.4,.5,8);taps=np.array([1.,-.2])
        hard=WhitenedViterbi(model,taps);soft=MaxLogByteDetector(hard)
        waves=[]
        for byte in range(256):
            wave=np.r_[np.ones(8),model.waveform(np.r_[encode([byte]),[0,0]])]
            waves.append(lfilter(taps,[1.],wave)[8:])
        waves=np.array(waves);rng=np.random.default_rng(8)
        obs=waves[123]+rng.normal(size=40)+1j*rng.normal(size=40)
        costs=np.sum(abs(waves-obs)**2,axis=1)
        expected=[]
        for k in range(8):
            labels=(np.arange(256)>>(7-k))&1
            expected.append(np.min(costs[labels==1])-np.min(costs[labels==0]))
        np.testing.assert_allclose(soft.decode(obs,2),expected,atol=1e-11)
        hardbits=np.unpackbits(np.array([np.argmin(costs)],dtype=np.uint8))
        np.testing.assert_array_equal(soft.decode(obs,2)<0,hardbits)

if __name__=='__main__':unittest.main()