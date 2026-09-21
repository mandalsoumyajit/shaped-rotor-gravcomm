"""Causal sampled receiver and matched finite-record Gaussian innovations.

Boundary: 16 complex samples/s of calibrated acceleration envelope. The
reference input PSD is specified at this sampled boundary; analog acquisition
and its anti-aliasing before this boundary are prerequisites, not simulated.
A stable SOS low-pass and fixed-phase decimation process signal and noise alike.
No FFT projection, future samples, or packet-wise noise normalization is used
in receiver processing. FFTs below only evaluate PSD/covariance or generate
stationary input noise for offline tests.
"""
from dataclasses import dataclass
import numpy as np
from scipy.signal import butter,sosfilt,sosfilt_zi,sosfreqz
from scipy.linalg import toeplitz,cholesky,solve_triangular
from .receivers import StructuralOscillator

@dataclass(frozen=True)
class CausalReceiver:
    input_rate_hz: float = 16.
    output_rate_hz: float = 4.
    cutoff_hz: float = 1.
    order: int = 6
    reference_center_hz: float = 50.3
    carrier_hz: float | None = None
    background: object | None = None
    def __post_init__(self):
        if self.background is not None:
            if self.carrier_hz is None or not np.isfinite(self.carrier_hz):
                raise ValueError('Environmental noise requires a physical carrier')
            # Validate the entire acquisition band, not only the nominal passband.
            self.background.acceleration_psd(np.array([self.carrier_hz-self.input_rate_hz/2,
                                                       self.carrier_hz+self.input_rate_hz/2]))
        ratio=self.input_rate_hz/self.output_rate_hz
        if ratio!=int(ratio) or ratio<1:raise ValueError('Integer decimation required')
        if not 0<self.cutoff_hz<self.output_rate_hz/2:raise ValueError('Cutoff below output Nyquist required')
        if self.order<1 or int(self.order)!=self.order:raise ValueError('Positive integer order required')
    @property
    def decimation(self):return int(self.input_rate_hz/self.output_rate_hz)
    @property
    def sos(self):return butter(self.order,self.cutoff_hz,fs=self.input_rate_hz,output='sos')
    def response(self,offset_hz):
        return sosfreqz(self.sos,worN=np.atleast_1d(offset_hz),fs=self.input_rate_hz)[1]
    def input_psd(self,offset_hz,amplitude_m_s2=1.):
        """Discrete proper-complex PSD, normalized by the chosen signal amplitude.

        Physical receiver ASD is unchanged. amplitude is a unit conversion,
        not a tunable noise floor. PSD average is per-sample variance.
        """
        if amplitude_m_s2<=0:raise ValueError('Positive normalization amplitude required')
        f=np.asarray(offset_hz)
        psd=StructuralOscillator().acceleration_noise_asd(self.reference_center_hz+f)**2
        if self.background is not None:
            psd=psd+self.background.acceleration_psd(self.carrier_hz+f)
        return 2*self.input_rate_hz*psd/amplitude_m_s2**2
    def output_psd(self,offset_hz,amplitude_m_s2=1.,components=False):
        """Aliased discrete output PSD: sum of D high-rate terms divided by D."""
        f=np.atleast_1d(offset_hz);terms=[]
        for k in range(self.decimation):
            shifted=(f+k*self.output_rate_hz+self.input_rate_hz/2)%self.input_rate_hz-self.input_rate_hz/2
            terms.append(abs(self.response(shifted))**2*self.input_psd(shifted,amplitude_m_s2)/self.decimation)
        terms=np.asarray(terms)
        return terms if components else terms.sum(axis=0)
    def covariance_lags(self,count,amplitude_m_s2=1.,nfft=65536):
        if nfft<2*count:raise ValueError('FFT must resolve requested covariance lags')
        f=np.fft.fftfreq(nfft,1/self.output_rate_hz)
        r=np.fft.ifft(self.output_psd(f,amplitude_m_s2))[:count]
        r[0]=r[0].real
        return r
    def finite_record(self,count,amplitude_m_s2=1.,nfft=65536):
        return FiniteRecordInnovations(self.covariance_lags(count,amplitude_m_s2,nfft))
    def stream(self,steady_input=0j):return ReceiverStream(self,steady_input)
    def stationary_noise(self,count,rng,amplitude_m_s2=1.,burn_s=32.,nfft=None):
        """Offline stationary-noise generator; filtering/decimation is causal.

        Burn-in belongs to generator initialization, not per-packet receiver
        latency. Physical use assumes an already operating stationary receiver.
        """
        burn=int(np.ceil(burn_s*self.output_rate_hz))*self.decimation
        needed=burn+count*self.decimation
        size=1<<int(np.ceil(np.log2(max(65536,4*needed)))) if nfft is None else nfft
        if size<2*needed:raise ValueError('Noise FFT record too short')
        f=np.fft.fftfreq(size,1/self.input_rate_hz)
        white=(rng.normal(size=size)+1j*rng.normal(size=size))/np.sqrt(2)
        noise=np.fft.ifft(np.fft.fft(white)*np.sqrt(self.input_psd(f,amplitude_m_s2)))
        segment=noise[size//4:size//4+needed]
        result=self.stream().process(segment)
        return result[burn//self.decimation:burn//self.decimation+count]

class ReceiverStream:
    def __init__(self,model,steady_input=0j):
        self.model=model;self.sos=model.sos
        # Known deterministic mean state only. Random stationary noise state
        # must be supplied by prehistory/continuous operation, never reset.
        self.state=sosfilt_zi(self.sos).astype(complex)*steady_input
        self.samples_seen=0
    def process(self,samples):
        samples=np.asarray(samples,dtype=complex)
        if samples.ndim!=1:raise ValueError('One stream required')
        if not len(samples):return np.empty(0,complex)
        filtered,self.state=sosfilt(self.sos,samples,zi=self.state)
        first=(-self.samples_seen)%self.model.decimation
        result=filtered[first::self.model.decimation]
        self.samples_seen+=len(samples)
        return result

class FiniteRecordInnovations:
    """Exact finite-window covariance metric for the specified stationary model.

    Cholesky is lower triangular: each innovation depends only on current/past
    observations. The coefficient table depends on model and horizon, not data.
    Full candidate mean histories must undergo the same transform. This is NOT
    a claim that an existing finite-memory trellis supplies exact metrics.
    """
    def __init__(self,lags):
        self.covariance=toeplitz(np.asarray(lags,dtype=complex))
        self.factor=cholesky(self.covariance,lower=True)
    def whiten(self,samples):
        x=np.asarray(samples,dtype=complex)
        if x.shape[0]!=len(self.factor):raise ValueError('Wrong observation length')
        return solve_triangular(self.factor,x,lower=True,check_finite=True)
    def metric(self,observed,mean):
        residual=np.asarray(observed)-np.asarray(mean)
        z=self.whiten(residual)
        return np.sum(abs(z)**2,axis=0)
    def stream(self):return InnovationStream(self.factor)

class InnovationStream:
    def __init__(self,factor):
        self.factor=factor;self.past=np.empty(len(factor),complex);self.seen=0
    def process(self,samples):
        x=np.asarray(samples,dtype=complex)
        if x.ndim!=1 or self.seen+len(x)>len(self.factor):raise ValueError('Innovation horizon exceeded')
        result=np.empty(len(x),complex)
        for j,value in enumerate(x):
            i=self.seen;z=(value-self.factor[i,:i]@self.past[:i])/self.factor[i,i]
            self.past[i]=z;self.seen+=1;result[j]=z
        return result

def sampled_cffsk(symbols,symbol_s,input_rate_hz=16.,transition_fraction=.4,q=.1):
    """Continuous-phase analytic waveform sampled on a fixed receiver clock.

    Arbitrary symbol durations need not have an integer number of samples.
    Initial offset tone and phase are zero. End time is excluded.
    """
    symbols=np.asarray(symbols,dtype=int)
    if not len(symbols) or symbol_s<=0:raise ValueError('Nonempty positive-duration waveform required')
    r=transition_fraction;spacing=q/((1-r)*symbol_s)
    before=np.r_[0,symbols[:-1]]
    increments=spacing*symbol_s*(r*before/2+(1-r/2)*symbols)
    phase0=np.r_[0,np.cumsum(increments)[:-1]]
    count=int(np.ceil(len(symbols)*symbol_s*input_rate_hz-1e-12))
    t=np.arange(count)/input_rate_hz;index=np.minimum((t/symbol_s).astype(int),len(symbols)-1)
    local=t-index*symbol_s;tr=r*symbol_s;u=np.minimum(local/tr,1)
    primitive=tr*(2.5*u**4-3*u**5+u**6)+np.maximum(local-tr,0)
    phase=phase0[index]+spacing*(before[index]*local+(symbols[index]-before[index])*primitive)
    return np.exp(2j*np.pi*phase)
