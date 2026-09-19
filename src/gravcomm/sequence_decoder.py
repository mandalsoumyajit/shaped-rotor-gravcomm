"""Full-waveform CF-FSK sequence detection with explicit phase memory.

Standard Viterbi uses accumulated phase, previous tone and any finite
whitening-filter memory. Branches include transition and dwell samples.
Initial tone, phase, symbol timing and frequency reference are acquired inputs.
These routines do not implement acquisition or infer an actuator/channel model.
Noise covariance uses E[n n^H] for proper complex samples, in signal units.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import lcm
import numpy as np


@dataclass(frozen=True)
class CFFSK:
    symbol_s: float = 8.
    transition_fraction: float = .2
    spacing_dwell_product: float = 1.
    samples_per_symbol: int = 32
    alphabet: tuple = (-3, -2, -1, 0, 1, 2, 3)

    def __post_init__(self):
        if not np.isfinite(self.symbol_s) or self.symbol_s <= 0:
            raise ValueError('symbol_s must be positive')
        if not 0 < self.transition_fraction < 1:
            raise ValueError('transition fraction must be between zero and one')
        if not np.isfinite(self.spacing_dwell_product) or self.spacing_dwell_product <= 0:
            raise ValueError('spacing must be positive')
        if not isinstance(self.samples_per_symbol, int) or self.samples_per_symbol < 2:
            raise ValueError('samples_per_symbol must be an integer >= 2')
        if not self.alphabet or len(set(self.alphabet)) != len(self.alphabet) or any(int(x)!=x for x in self.alphabet):
            raise ValueError('alphabet must contain distinct integers')

    @property
    def spacing_hz(self):
        return self.spacing_dwell_product/((1-self.transition_fraction)*self.symbol_s)

    def phase_increment_cycles(self, previous, current):
        r=self.transition_fraction
        return self.spacing_hz*self.symbol_s*(r*previous/2+(1-r/2)*current)

    def branch(self, previous, current, phase_rad=0.):
        """Complex unit-amplitude envelope, including the entire transition.

        Samples use left endpoints, matching desktop_cf_fsk. The next branch
        starts at the analytically integrated boundary phase, without a reset.
        """
        t=np.arange(self.samples_per_symbol)*self.symbol_s/self.samples_per_symbol
        tr=self.transition_fraction*self.symbol_s;u=np.minimum(t/tr,1)
        integral=tr*(2.5*u**4-3*u**5+u**6)+np.maximum(t-tr,0)
        phase=phase_rad+2*np.pi*self.spacing_hz*(previous*t+(current-previous)*integral)
        return np.exp(1j*phase)

    def waveform(self, symbols, previous=0, phase_rad=0.):
        if previous not in self.alphabet or any(x not in self.alphabet for x in symbols):
            raise ValueError('unknown tone')
        blocks=[]
        for current in symbols:
            blocks.append(self.branch(previous,current,phase_rad))
            phase_rad=(phase_rad+2*np.pi*self.phase_increment_cycles(previous,current))%(2*np.pi)
            previous=current
        return np.concatenate(blocks) if blocks else np.empty(0,complex)

    def phase_lattice(self, maximum_states=4096):
        """Exact rational lattice; reject excessive state counts, never round."""
        r=Fraction(str(self.transition_fraction));q=Fraction(str(self.spacing_dwell_product))
        increments={(a,b):q/(1-r)*(r*a/2+(1-r/2)*b) for a in self.alphabet for b in self.alphabet}
        size=lcm(*(x.denominator for x in increments.values()))
        if size*len(self.alphabet)>maximum_states:
            raise ValueError('phase lattice too large; use full-history decoding')
        return size,{key:int(value*size)%size for key,value in increments.items()}


def _observations(model, received, previous, phase_rad):
    y=np.asarray(received,dtype=complex)
    if y.ndim!=1 or not len(y) or len(y)%model.samples_per_symbol or np.any(~np.isfinite(y)):
        raise ValueError('received must contain complete finite symbols')
    if previous not in model.alphabet or not np.isfinite(phase_rad):
        raise ValueError('invalid initial state')
    return y


def viterbi_fir(model, received, *, previous=0, phase_rad=0., amplitude=1.,
                whitening_taps=(1.,), noise_variance=1., maximum_states=100000):
    """Standard add-compare-select Viterbi and predecessor traceback.

    A supplied causal FIR whitens the noise. Its intersymbol memory is included
    in the trellis state together with phase and tone history. Initial filter
    residual state is zero (or known prepacket contributions were removed).
    noise_variance is E[|innovation|^2]. Exact ML requires white innovations;
    fitting a finite whitener to a measured spectrum is a separate operation.
    """
    from scipy.signal import lfilter
    y=_observations(model,received,previous,phase_rad)
    taps=np.asarray(whitening_taps,dtype=complex)
    if taps.ndim!=1 or not len(taps) or np.any(~np.isfinite(taps)) or taps[0]==0:
        raise ValueError('invalid whitening filter')
    if not np.isfinite(noise_variance) or noise_variance<=0 or not np.isfinite(amplitude):
        raise ValueError('invalid amplitude or noise variance')
    size,increments=model.phase_lattice();n=model.samples_per_symbol
    memory=len(taps)-1
    # One preceding tone is needed to construct the earliest retained branch.
    history_length=1+(memory+n-1)//n
    initial=(0,(previous,))
    survivors={initial:(0.,np.empty(0,complex))}
    backpointers=[];white_y=lfilter(taps,[1.],y)
    for start in range(0,len(y),n):
        next_survivors={};parents={}
        for state,(metric,tail) in survivors.items():
            p,history=state;old=history[-1]
            for current in model.alphabet:
                block=amplitude*model.branch(old,current,phase_rad+2*np.pi*p/size)
                joined=np.r_[tail,block]
                predicted=lfilter(taps,[1.],joined)[-n:]
                score=metric+float(np.sum(abs(white_y[start:start+n]-predicted)**2)/noise_variance)
                destination=((p+increments[old,current])%size,(history+(current,))[-history_length:])
                if destination not in next_survivors or score<next_survivors[destination][0]:
                    next_survivors[destination]=(score,joined[-memory:] if memory else np.empty(0,complex))
                    parents[destination]=(state,current)
        if len(next_survivors)>maximum_states:
            raise ValueError('trellis exceeds maximum_states; no survivor pruning performed')
        survivors=next_survivors;backpointers.append(parents)
    state=min(survivors,key=lambda key:survivors[key][0]);score=survivors[state][0]
    decoded=[]
    for parents in reversed(backpointers):
        state,symbol=parents[state];decoded.append(symbol)
    return dict(symbols=np.array(decoded[::-1]),metric=score,phase_states=size,
                final_survivor_states=len(survivors),whitening_memory_samples=memory,exact=True)


def viterbi_white(model, received, **kwargs):
    """White-noise special case of the same standard Viterbi recursion."""
    return viterbi_fir(model,received,**kwargs)
