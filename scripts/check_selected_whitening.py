"""Paired filter checks; these are convergence diagnostics, not BER estimates."""
import json
import numpy as np
from scipy.signal import lfilter
import desktop_cf_fsk as cf
from search_sequence_cf_fsk import A,OUT
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import fit_whitener,WhitenedViterbi,receiver_noise

model=CFFSK(2.,.4,.1,4);rng=np.random.default_rng(728314)
short,rs=fit_whitener(2,A,12);long,rl=fit_whitener(2,A,16)
ds=WhitenedViterbi(model,short,maximum_states=1000000)
dl=WhitenedViterbi(model,long,maximum_states=1000000)
counts=dict(packets=500,short_bit_errors=0,long_bit_errors=0,disagreeing_payload_symbols=0)
prefix=16
for _ in range(counts['packets']):
    data=rng.integers(0,256,8,dtype=np.uint8);symbols=cf.encode(data)
    signal=np.r_[np.ones(prefix),model.waveform(np.r_[symbols,np.zeros(6,int)])]
    received=signal+receiver_noise(rng,len(signal),2,A)
    a=ds.decode_whitened(lfilter(short,[1.],received)[prefix:],6,True)[:24]
    b=dl.decode_whitened(lfilter(long,[1.],received)[prefix:],6,True)[:24]
    counts['short_bit_errors']+=int(np.unpackbits(data^cf.decode(a)).sum())
    counts['long_bit_errors']+=int(np.unpackbits(data^cf.decode(b)).sum())
    counts['disagreeing_payload_symbols']+=int(np.count_nonzero(a!=b))
counts.update(short_filter=rs,long_filter=rl)
(OUT/'whitening_convergence.json').write_text(json.dumps(counts,indent=2)+'\n')
print(json.dumps(counts,indent=2),flush=True)

model=CFFSK(2.,.4,.1,8);taps,fit=fit_whitener(4,A,24)
decoder=WhitenedViterbi(model,taps,maximum_states=1000000)
bit_errors=0
for _ in range(500):
    data=rng.integers(0,256,8,dtype=np.uint8);symbols=cf.encode(data)
    recovered=cf.decode(decoder.observe_and_decode(symbols,rng,A,tail_symbols=6,byte_mapping=True))
    bit_errors+=int(np.unpackbits(data^recovered).sum())
record=dict(packets=500,bits=32000,bit_errors=bit_errors,whitener=fit,
    purpose='Sample-rate convergence diagnostic; low counts do not estimate BER precisely')
(OUT/'sampling_convergence.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
