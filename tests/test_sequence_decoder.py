"""Independent exhaustive likelihood and waveform checks for sequence detection."""
import itertools
import numpy as np
import pytest
from gravcomm.sequence_decoder import CFFSK, viterbi_white, viterbi_fir


def oracle(model,y,K,H,previous,phase):
    candidates=[]
    for symbols in itertools.product(model.alphabet,repeat=len(y)//model.samples_per_symbol):
        error=y-H@model.waveform(symbols,previous,phase)
        metric=float(np.vdot(error,np.linalg.solve(K,error)).real)
        candidates.append((metric,symbols))
    return min(candidates)


def test_boundary_phase_matches_independent_frequency_integration():
    model=CFFSK();assert model.phase_lattice()[0]==8
    t=np.linspace(0,model.symbol_s,100001);tr=model.transition_fraction*model.symbol_s
    u=np.minimum(t/tr,1);shape=10*u**3-15*u**4+6*u**5
    for a,b in itertools.product(model.alphabet,repeat=2):
        numeric=np.trapz(model.spacing_hz*(a+(b-a)*shape),t)
        assert abs(numeric-model.phase_increment_cycles(a,b))<1e-10


@pytest.mark.parametrize('seed',range(4))
def test_viterbi_matches_exhaustive_seven_tone_packets(seed):
    rng=np.random.default_rng(seed);model=CFFSK(samples_per_symbol=8)
    phase=.371;previous=-2;truth=rng.choice(model.alphabet,3)
    signal=model.waveform(truth,previous,phase);N=len(signal)
    y=signal+1.8*(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
    expected=oracle(model,y,np.eye(N),np.eye(N),previous,phase)
    result=viterbi_white(model,y,previous=previous,phase_rad=phase)
    np.testing.assert_array_equal(result['symbols'],expected[1])
    assert result['metric']==pytest.approx(expected[0],rel=1e-12)


def test_long_noiseless_packet_preserves_phase():
    model=CFFSK(samples_per_symbol=12);truth=np.random.default_rng(31).choice(model.alphabet,80)
    y=model.waveform(truth,-3,1.193)
    np.testing.assert_array_equal(viterbi_white(model,y,previous=-3,phase_rad=1.193)['symbols'],truth)


def test_waveform_matches_existing_desktop_generator():
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
    import desktop_cf_fsk as original
    symbols,previous=original.all_transitions()
    _,_,_,_,phase=original.trajectory(symbols,previous,8.)
    model=CFFSK(samples_per_symbol=800)
    np.testing.assert_allclose(model.waveform(symbols,previous),np.exp(1j*phase[:-1]),atol=2e-11)


@pytest.mark.parametrize('seed',range(3))
def test_colored_noise_viterbi_matches_direct_packet_likelihood(seed):
    rng=np.random.default_rng(seed);model=CFFSK(samples_per_symbol=4)
    N=12;taps=np.array([1.,-.7*np.exp(.24j)])
    W=np.eye(N,dtype=complex)+taps[1]*np.eye(N,k=-1)
    inverse=np.linalg.inv(W);K=inverse@inverse.conj().T
    truth=rng.choice(model.alphabet,3);phase=.47;previous=1
    noise=inverse@(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
    y=model.waveform(truth,previous,phase)+noise
    expected=oracle(model,y,K,np.eye(N),previous,phase)
    result=viterbi_fir(model,y,previous=previous,phase_rad=phase,whitening_taps=taps)
    np.testing.assert_array_equal(result['symbols'],expected[1])
    assert result['metric']==pytest.approx(expected[0],rel=1e-12)


def test_multisymbol_filter_memory_matches_exhaustive():
    model=CFFSK(samples_per_symbol=2,alphabet=(-1,0,1));truth=[1,-1,0,1]
    taps=np.array([1.,-.2,.1,0,.05j]);N=8
    W=sum(t*np.eye(N,k=-j) for j,t in enumerate(taps));inv=np.linalg.inv(W)
    K=inv@inv.conj().T
    rng=np.random.default_rng(45)
    y=model.waveform(truth)+inv@(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
    expected=oracle(model,y,K,np.eye(N),0,0.)
    result=viterbi_fir(model,y,whitening_taps=taps)
    np.testing.assert_array_equal(result['symbols'],expected[1])
    assert result['metric']==pytest.approx(expected[0],rel=1e-12)
