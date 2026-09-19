"""Spectral whitening and precomputed standard CF-FSK Viterbi trellis.

The discrete complex-baseband PSD is 2*fs*S_a(fc+f)/A**2. Signals are
normalized to the harmonic amplitude A. Samples span [-fs/2,fs/2).
Finite FIR whitening is an approximation checked against the original PSD.
Numba accelerates only add-compare-select and traceback; it is optional.
"""
import itertools
from math import gcd
import numpy as np
from scipy.signal import lfilter
from .receivers import CARTER_OSCILLATOR
try:
    from numba import njit
except ImportError:
    def njit(*args,**kwargs):
        return lambda f:f


def noise_psd(fs, amplitude, size, receiver=CARTER_OSCILLATOR, carrier=50.3):
    f=np.fft.fftfreq(size,1/fs)
    return 2*fs*receiver.acceleration_noise_asd(carrier+f)**2/amplitude**2


def fit_whitener(fs, amplitude, memory, size=32768, *, receiver=CARTER_OSCILLATOR):
    """Minimum-phase spectral factor of inverse noise PSD, truncated causally."""
    psd=noise_psd(fs,amplitude,size,receiver=receiver)
    cepstrum=np.fft.ifft(-.5*np.log(psd))
    causal=np.zeros_like(cepstrum);causal[0]=cepstrum[0]
    causal[1:size//2]=2*cepstrum[1:size//2];causal[size//2]=cepstrum[size//2]
    impulse=np.fft.ifft(np.exp(np.fft.fft(causal)))
    taps=impulse[:memory+1].copy()
    ratio=abs(np.fft.fft(taps,size))**2*psd
    variance=float(np.mean(ratio));taps/=np.sqrt(variance);ratio/=variance
    f=np.fft.fftfreq(size,1/fs);central=abs(f)<.8*fs/2
    report=dict(sample_rate_hz=fs,taps=len(taps),memory_s=memory/fs,
        rms_psd_error=float(np.sqrt(np.mean((ratio-1)**2))),
        max_psd_error=float(max(abs(ratio-1))),
        central_max_psd_error=float(max(abs(ratio[central]-1))))
    return taps,report


def receiver_noise(rng, length, fs, amplitude, minimum_fft=32768, *, receiver=CARTER_OSCILLATOR):
    """Stationary colored noise from the physical PSD, independent of FIR fit.

    A segment shorter than half the FFT record avoids periodic packet wrapping.
    Standard proper complex Gaussian FFT coefficients use the unitary scaling
    provided by FFT(white); each sample has variance mean(discrete PSD).
    """
    size=1<<int(np.ceil(np.log2(max(minimum_fft,2*length))))
    white=(rng.normal(size=size)+1j*rng.normal(size=size))/np.sqrt(2)
    colored=np.fft.ifft(np.fft.fft(white)*np.sqrt(noise_psd(fs,amplitude,size,receiver=receiver)))
    return colored[size//4:size//4+length]


@njit(cache=True,nogil=True)
def _decode(correlation,energy,phase_factors,increments,histories,alphabet,initial_phase,initial_history,stride,alpha,forced,byte_payload_symbols):
    length=correlation.shape[0];phases=len(phase_factors);states=(phases//stride)*histories
    alphabet_size=len(alphabet)
    costs=np.full(states,np.inf);costs[initial_phase*histories+initial_history]=0.
    parents=np.empty((length,states),np.int32)
    decisions=np.empty((length,states),np.int8)
    for step in range(length):
        new=np.full(states,np.inf)
        for state in range(states):
            if not np.isfinite(costs[state]):continue
            h=state%histories;old=h%alphabet_size
            p=((state//histories)*stride-alpha*alphabet[old])%phases
            for current in range(alphabet_size):
                if forced[step]>=0 and current!=forced[step]:continue
                if step<byte_payload_symbols and step%3==2:
                    first=(h//alphabet_size)%alphabet_size
                    if first*49+old*7+current>255:continue
                edge=h*alphabet_size+current
                metric=costs[state]+energy[edge]-2*(phase_factors[p]*correlation[step,edge]).real
                next_phase=(p+increments[old,current])%phases
                reduced_phase=((next_phase+alpha*alphabet[current])%phases)//stride
                dest=reduced_phase*histories+(h*alphabet_size+current)%histories
                if metric<new[dest]:
                    new[dest]=metric;parents[step,dest]=state;decisions[step,dest]=current
        offset=np.min(new);costs=new-offset
    state=np.argmin(costs);decoded=np.empty(length,np.int64)
    for step in range(length-1,-1,-1):
        decoded[step]=decisions[step,state];state=parents[step,state]
    return decoded


class WhitenedViterbi:
    """Precompute the exact finite-FIR trellis; never prune survivors.

    Packet boundaries use an acquired zero-tone, zero-phase prehistory.
    Observations must already include that prehistory when FIR filtering.
    """
    def __init__(self,model,taps,maximum_states=100000):
        self.model=model;self.taps=np.asarray(taps);n=model.samples_per_symbol
        self.memory_symbols=int(np.ceil((len(taps)-1)/n))
        self.history_length=1+self.memory_symbols
        alphabet=model.alphabet;M=len(alphabet)
        self.phases,inc=model.phase_lattice()
        self.increments=np.array([[inc[a,b] for b in alphabet] for a in alphabet],np.int64)
        # Reachable phase obeys p+alpha*last = 0 modulo gcd(D,alpha+beta).
        # This removes unreachable states only; no survivor pruning occurs.
        self.alpha=inc[1,0] if 1 in alphabet else 0
        self.stride=gcd(self.phases,inc[1,1]) if 1 in alphabet else 1
        self.histories=M**self.history_length
        self.state_count=(self.phases//self.stride)*self.histories
        if self.state_count>maximum_states:raise ValueError('trellis state budget exceeded')
        templates=[]
        for history in itertools.product(alphabet,repeat=self.history_length):
            # End of these historical transitions is current boundary phase 0.
            initial=-2*np.pi*sum(model.phase_increment_cycles(a,b) for a,b in zip(history[:-1],history[1:]))
            prefix=model.waveform(history[1:],history[0],initial)
            for current in alphabet:
                waveform=np.r_[prefix,model.branch(history[-1],current,0.)]
                templates.append(lfilter(taps,[1.],waveform)[-n:])
        self.templates=np.array(templates)
        self.energy=np.sum(abs(self.templates)**2,axis=1)
        self.phase_factors=np.exp(2j*np.pi*np.arange(self.phases)/self.phases)
        zero=alphabet.index(0)
        self.initial_history=sum(zero*M**j for j in range(self.history_length))

    def decode_whitened(self,packet,tail_symbols=0,byte_mapping=False):
        blocks=np.asarray(packet).reshape(-1,self.model.samples_per_symbol)
        # |y-s|^2 drops the same observation energy from all branch metrics.
        correlation=blocks.conj()@self.templates.T
        forced=np.full(len(blocks),-1,np.int64)
        if tail_symbols:forced[-tail_symbols:]=self.model.alphabet.index(0)
        payload=len(blocks)-tail_symbols
        if byte_mapping and (self.model.alphabet!=(-3,-2,-1,0,1,2,3) or self.history_length<2 or payload%3):
            raise ValueError('byte mapping requires seven tones, complete triples and two-tone history')
        decoded=_decode(correlation,self.energy,self.phase_factors,self.increments,
                        self.histories,np.array(self.model.alphabet),0,self.initial_history,self.stride,self.alpha,forced,payload if byte_mapping else 0)
        return np.array(self.model.alphabet)[decoded]

    def observe_and_decode(self,symbols,rng,amplitude,noise=True,tail_symbols=0,byte_mapping=False):
        n=self.model.samples_per_symbol;prefix=self.memory_symbols*n
        transmitted=np.r_[symbols,np.zeros(tail_symbols,dtype=int)]
        waveform=np.r_[np.ones(prefix,complex),self.model.waveform(transmitted)]
        if noise:
            waveform+=receiver_noise(rng,len(waveform),n/self.model.symbol_s,amplitude)
        received=lfilter(self.taps,[1.],waveform)[prefix:]
        return self.decode_whitened(received,tail_symbols,byte_mapping)[:len(symbols)]
