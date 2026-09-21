"""Numerical regression checks for the separate range-revision harness."""
import sys, unittest
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from range_revision import CenteredReceiver,duty,make_decoder,FC
from gravcomm import FiniteRotor
from gravcomm.receivers import StructuralOscillator
from gravcomm.whitened_viterbi import noise_psd,receiver_noise
from dataclasses import asdict,replace
from desktop_cf_fsk import encode

class RangeRevisionChecks(unittest.TestCase):
    def test_physical_carrier_translation_preserves_noise(self):
        physical=replace(StructuralOscillator(),resonance_hz=25.15)
        receiver=CenteredReceiver(physical,25.15)
        fs,a,n=2.,1e-10,1024
        frequencies=25.15+np.fft.fftfreq(n,1/fs)
        np.testing.assert_allclose(noise_psd(fs,a,n,receiver=receiver),
            2*fs*physical.acceleration_noise_asd(frequencies)**2/a**2,rtol=1e-12)

    def test_exact_message_duty_reproduces_reference(self):
        result=duty(2.,.4,.1)
        self.assertTrue(result['feasible'])
        self.assertAlmostEqual(result['worst_valid_message_rms_nm'],.4474016719053217,places=7)
        self.assertAlmostEqual(result['peak_torque_nm'],1.2136335596129677,places=5)
        self.assertFalse(duty(1.75,.4,.1)['feasible'])

    def test_similarity_for_finite_source(self):
        p=np.array([[.2,0,.02],[-.2,0,-.02],[.05,.07,0],[-.05,-.07,0]])
        m=np.array([2.,2.,.5,.5]);source=FiniteRotor(p,m,.25)
        s=2.4;scaled=FiniteRotor(s*p,s**3*m,s*.25)
        phases=np.linspace(0,2*np.pi,23)
        np.testing.assert_allclose(scaled.acceleration_component(2.,phases),
            s*source.acceleration_component(2./s,phases),rtol=1e-13)

    def test_custom_receiver_path_matches_original_at_baseline(self):
        a=5.180259430854158e-10
        case=dict(receiver=asdict(StructuralOscillator()),carrier_hz=FC,amplitude_m_s2=a)
        model,decoder,receiver,fit=make_decoder(case,2.,.1)
        symbols=encode(np.array([0,255,42,251,128,7,99,201]))
        original=decoder.observe_and_decode(symbols,np.random.default_rng(61842),a,tail_symbols=6,byte_mapping=True)
        prefix=decoder.memory_symbols*model.samples_per_symbol
        waveform=np.r_[np.ones(prefix,complex),model.waveform(np.r_[symbols,np.zeros(6,int)])]
        waveform+=receiver_noise(np.random.default_rng(61842),len(waveform),8.,a,receiver=receiver)
        adapted=decoder.decode_whitened(lfilter(decoder.taps,[1.],waveform)[prefix:],6,True)[:24]
        np.testing.assert_array_equal(adapted,original)

if __name__=='__main__':unittest.main()
