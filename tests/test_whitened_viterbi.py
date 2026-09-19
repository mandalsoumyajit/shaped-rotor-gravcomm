import itertools
import numpy as np
from scipy.signal import lfilter
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import WhitenedViterbi,fit_whitener,receiver_noise,noise_psd


def test_precomputed_trellis_matches_exhaustive_whitened_metric():
    model=CFFSK(samples_per_symbol=4,alphabet=(-1,0,1))
    taps=np.array([1.,-.3+.2j,.1]);decoder=WhitenedViterbi(model,taps)
    rng=np.random.default_rng(23)
    for _ in range(3):
        y=rng.normal(size=12)+1j*rng.normal(size=12)
        candidates=[]
        for symbols in itertools.product(model.alphabet,repeat=3):
            x=np.r_[np.ones(4),model.waveform(symbols)]
            filtered=lfilter(taps,[1.],x)[4:]
            candidates.append((sum(abs(y-filtered)**2),symbols))
        expected=min(candidates)[1]
        np.testing.assert_array_equal(decoder.decode_whitened(y),expected)


def test_full_alphabet_noiseless_packet():
    model=CFFSK(symbol_s=4,transition_fraction=.5,spacing_dwell_product=.5,samples_per_symbol=8)
    taps,_=fit_whitener(2,5.180259430854158e-10,16)
    decoder=WhitenedViterbi(model,taps)
    rng=np.random.default_rng(12);symbols=rng.choice(model.alphabet,40)
    np.testing.assert_array_equal(decoder.observe_and_decode(symbols,rng,5.180259430854158e-10,noise=False),symbols)


def test_known_trailer_matches_exhaustive_payload_only_search():
    model=CFFSK(2.,.5,.0625,4,(-1,0,1));taps=np.array([1.,-.2j])
    decoder=WhitenedViterbi(model,taps);rng=np.random.default_rng(99)
    y=rng.normal(size=16)+1j*rng.normal(size=16)
    candidates=[]
    for payload in itertools.product(model.alphabet,repeat=2):
        x=np.r_[np.ones(4),model.waveform(payload+(0,0))]
        prediction=lfilter(taps,[1.],x)[4:]
        candidates.append((sum(abs(y-prediction)**2),payload+(0,0)))
    np.testing.assert_array_equal(decoder.decode_whitened(y,tail_symbols=2),min(candidates)[1])


def test_byte_constrained_viterbi_matches_all_256_messages():
    model=CFFSK(2.,.4,.1,4);taps=np.array([1.,-.2j])
    decoder=WhitenedViterbi(model,taps);rng=np.random.default_rng(671)
    y=rng.normal(size=20)+1j*rng.normal(size=20)
    candidates=[]
    for byte in range(256):
        symbols=(byte//49-3,(byte//7)%7-3,byte%7-3,0,0)
        signal=np.r_[np.ones(4),model.waveform(symbols)]
        filtered=lfilter(taps,[1.],signal)[4:]
        candidates.append((sum(abs(y-filtered)**2),symbols))
    np.testing.assert_array_equal(decoder.decode_whitened(y,2,True),min(candidates)[1])


def test_whitener_converges_and_noise_normalization():
    A=5.180259430854158e-10
    short,s=fit_whitener(2,A,8);long,l=fit_whitener(2,A,16)
    assert l['rms_psd_error']<.005
    assert l['rms_psd_error']<s['rms_psd_error']/5
    rng=np.random.default_rng(123)
    noise=receiver_noise(rng,100000,2,A,minimum_fft=262144)
    expected=np.mean(noise_psd(2,A,262144))
    assert abs(np.mean(abs(noise)**2)/expected-1)<.025
    white=lfilter(long,[1.],noise)[100:]
    assert abs(np.mean(abs(white)**2)-1)<.025
    assert abs(np.mean(white[1:]*white[:-1].conj()))<.025


def test_custom_receiver_noise_and_whitener_use_same_psd():
    from gravcomm.receivers import CARTER_OSCILLATOR
    class ScaledReceiver:
        def acceleration_noise_asd(self, f):
            return 2*CARTER_OSCILLATOR.acceleration_noise_asd(f)
    receiver=ScaledReceiver(); A=5.180259430854158e-10
    a=receiver_noise(np.random.default_rng(47),1000,4,A)
    b=receiver_noise(np.random.default_rng(47),1000,4,A,receiver=receiver)
    np.testing.assert_allclose(b,2*a,rtol=1e-13,atol=1e-13)
    h,_=fit_whitener(4,A,24)
    g,_=fit_whitener(4,A,24,receiver=receiver)
    np.testing.assert_allclose(g,h/2,rtol=1e-11,atol=1e-13)
