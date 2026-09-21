import unittest
import numpy as np
from numpy.testing import assert_allclose,assert_array_equal
from scipy.linalg import solve
from scipy.integrate import quad
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.sequence_decoder import CFFSK

A=5.180259430854157e-10
class CausalReceiverTests(unittest.TestCase):
    def setUp(self):self.r=CausalReceiver()
    def test_chunk_invariance_and_no_future_dependence(self):
        rng=np.random.default_rng(5);x=rng.normal(size=1003)+1j*rng.normal(size=1003)
        full=self.r.stream(1j).process(x)
        stream=self.r.stream(1j);parts=[];start=0
        for n in (1,7,13,256,3,511,212):
            parts.append(stream.process(x[start:start+n]));start+=n
        assert_array_equal(np.concatenate(parts),full)
        changed=x.copy();changed[501:]+=1e6
        assert_array_equal(self.r.stream(1j).process(changed)[:126],full[:126])
        assert_array_equal(self.r.stream(1j).process(x[:501]),full[:126])
    def test_alias_sum_equals_high_rate_lag_selection(self):
        n=65536;f=np.fft.fftfreq(n*4,1/16.)
        high=np.fft.ifft(self.r.input_psd(f,A)*abs(self.r.response(f))**2)
        low=self.r.covariance_lags(160,A,n)
        assert_allclose(low,high[np.arange(160)*4],atol=1e-12,rtol=1e-11)
        # Asymmetry is real model information, not forcibly symmetrized away.
        self.assertGreater(np.max(abs(low.imag)),1e-5)
    def test_filter_order_meets_alias_budget(self):
        f=np.linspace(-1,1,4001);parts=self.r.output_psd(f,A,True)
        self.assertLess(np.max(parts[1:].sum(axis=0)/parts[0]),1e-4)
        assert_allclose(abs(self.r.response([0,1]))**2,[1,.5],atol=1e-12)
    def test_covariance_convergence_and_innovation_prefix(self):
        c=self.r.finite_record(128,A,32768);fine=self.r.finite_record(128,A,65536)
        assert_allclose(c.covariance,fine.covariance,atol=1e-12,rtol=1e-10)
        short=self.r.finite_record(64,A)
        assert_allclose(short.factor,c.factor[:64,:64],atol=1e-11,rtol=1e-10)
        rng=np.random.default_rng(7);x=rng.normal(size=128)+1j*rng.normal(size=128)
        stream=c.stream();z=np.r_[stream.process(x[:17]),stream.process(x[17:81]),stream.process(x[81:])]
        assert_allclose(z,c.whiten(x),atol=1e-10,rtol=1e-10)
        modified=x.copy();modified[64:]+=1e6
        assert_allclose(c.whiten(modified)[:64],c.whiten(x)[:64],atol=1e-10,rtol=1e-10)
    def test_exhaustive_finite_window_likelihood(self):
        means=[]
        for byte in range(256):
            word=np.array([byte//49,byte//7%7,byte%7])-3
            signal=sampled_cffsk(np.r_[word,np.zeros(6,int)],2.)
            means.append(self.r.stream(1.).process(signal))
        means=np.array(means).T;model=self.r.finite_record(len(means),A)
        rng=np.random.default_rng(18);noise=model.factor@((rng.normal(size=len(means))+1j*rng.normal(size=len(means)))/np.sqrt(2))
        y=means[:,91]+noise;residual=y[:,None]-means
        innovations=model.whiten(residual);cost=np.sum(abs(innovations)**2,axis=0)
        direct=np.real(np.sum(residual.conj()*solve(model.covariance,residual,assume_a='her'),axis=0))
        assert_allclose(cost,direct,atol=1e-9,rtol=1e-10)
        self.assertEqual(np.argmin(cost),np.argmin(direct))
        for bit in range(8):
            ones=(np.arange(256)>>(7-bit))&1
            self.assertAlmostEqual(cost[ones==1].min()-cost[ones==0].min(),direct[ones==1].min()-direct[ones==0].min(),places=8)
    def test_sampler_matches_analytic_symbol_model(self):
        symbols=np.array([-3,2,0,1,-2,3]);m=CFFSK(2.,.4,.1,32)
        assert_allclose(sampled_cffsk(symbols,2.),m.waveform(symbols),atol=1e-13)
        # No reset of receiver sample clock when symbol duration is fractional.
        x=sampled_cffsk(symbols,256/126)
        self.assertEqual(len(x),int(np.ceil(6*(256/126)*16)))
        duration=256/126;transition=.4*duration;spacing=.1/(.6*duration)
        def frequency(t):
            index=min(int(t/duration),len(symbols)-1)
            previous=symbols[index-1] if index else 0
            u=min((t-index*duration)/transition,1.)
            return spacing*(previous+(symbols[index]-previous)*(10*u**3-15*u**4+6*u**5))
        boundaries=sorted([i*duration+d for i in range(len(symbols)) for d in (0,transition)])
        for index in range(0,len(x),7):
            t=index/16
            phase=quad(frequency,0,t,points=[v for v in boundaries if 0<v<t],epsabs=1e-12)[0]
            assert_allclose(x[index],np.exp(2j*np.pi*phase),atol=1e-11)
    def test_continuous_filter_state_across_packets(self):
        # A continuous transmitter phase history split at a packet boundary.
        symbols=np.array([3,-2,1,0,-3,2, -1,3,0,2,-2,1])
        x=sampled_cffsk(symbols,256/126)
        boundary=int(np.ceil(6*(256/126)*16))
        stream=self.r.stream(1.)
        first=stream.process(x[:boundary]);second=stream.process(x[boundary:])
        assert_array_equal(np.r_[first,second],self.r.stream(1.).process(x))
        self.assertEqual(stream.samples_seen,len(x))
        # Resetting a filter at the second packet changes its physical response.
        reset=self.r.stream(1.).process(x[boundary:])
        self.assertGreater(np.max(abs(second[:min(len(second),len(reset))]-reset[:min(len(second),len(reset))])),.01)

if __name__=='__main__':unittest.main()
