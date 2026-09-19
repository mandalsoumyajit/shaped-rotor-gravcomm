"""Independent error-count confirmation at finer sampling and added noise.

Run under WSL (fork workers). Retains every packet error count and supports
safe resumption with an identical configuration. No bits within a packet are
assumed independent for confidence intervals.
"""
import argparse
import json
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import numpy as np
from scipy.signal import lfilter
import desktop_cf_fsk as cf
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import fit_whitener, receiver_noise, WhitenedViterbi
from gravcomm.packet_confidence import packet_ber_interval
from gravcomm.receivers import CARTER_OSCILLATOR

A=5.180259430854158e-10
OUT=Path(__file__).resolve().parents[1]/'results/communications_checks'
OUT.mkdir(exist_ok=True)

class AddedNoise:
    def __init__(self, asd): self.asd=asd
    def acceleration_noise_asd(self, f):
        return np.hypot(CARTER_OSCILLATOR.acceleration_noise_asd(f), self.asd)

def trial(index):
    rng=np.random.default_rng(np.random.SeedSequence([246819,index]))
    data=rng.integers(0,256,8,dtype=np.uint8)
    symbols=cf.encode(data)
    signal=np.r_[np.ones(PREFIX), MODEL.waveform(np.r_[symbols,np.zeros(6,int)])]
    received=signal+receiver_noise(rng,len(signal),ARGS.fs,A,receiver=RECEIVER)
    decoded=DECODER.decode_whitened(lfilter(TAPS,[1.],received)[PREFIX:],6,True)[:24]
    return int(np.unpackbits(data^cf.decode(decoded)).sum())

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--fs',type=float,default=4)
    p.add_argument('--memory',type=float,default=6)
    p.add_argument('--technical-asd',type=float,default=0)
    p.add_argument('--max-packets',type=int,default=150000)
    p.add_argument('--workers',type=int,default=4)
    ARGS=p.parse_args()
    assert ARGS.max_packets%500==0
    MODEL=CFFSK(2.,.4,.1,round(2*ARGS.fs))
    assert MODEL.samples_per_symbol/2==ARGS.fs
    RECEIVER=AddedNoise(ARGS.technical_asd)
    TAPS,fit=fit_whitener(ARGS.fs,A,round(ARGS.memory*ARGS.fs),receiver=RECEIVER)
    DECODER=WhitenedViterbi(MODEL,TAPS,maximum_states=1000000)
    PREFIX=DECODER.memory_symbols*MODEL.samples_per_symbol
    config=dict(fs=ARGS.fs,memory=ARGS.memory,technical_asd=ARGS.technical_asd,
                symbol_s=2.,transition_fraction=.4,q=.1,seed=246819)
    name=f'fs{ARGS.fs:g}_mem{ARGS.memory:g}_noise{ARGS.technical_asd:g}'.replace('.','p')
    path=OUT/(name+'.json'); errors=[]; previous=0.
    if path.exists():
        old=json.loads(path.read_text()); assert old['config']==config
        errors=old['bit_errors_per_packet'];previous=old['elapsed_s']
        if old['adequate_error_count']: raise SystemExit('Already complete')
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=ARGS.workers,mp_context=get_context('fork')) as pool:
        for end in range(len(errors)+500,ARGS.max_packets+1,500):
            errors.extend(pool.map(trial,range(len(errors),end)))
            count=sum(errors); failed=int(np.count_nonzero(errors))
            adequate=count>=100 and failed>=30
            record=dict(config=config,whitener=fit,states=DECODER.state_count,
                        packets=end,bits=end*64,bit_errors=count,errored_packets=failed,
                        observed_ber=count/(end*64),adequate_error_count=adequate,
                        interval95=packet_ber_interval(errors),bit_errors_per_packet=errors,
                        complete=adequate or end==ARGS.max_packets,
                        elapsed_s=previous+time.monotonic()-started)
            path.write_text(json.dumps(record,indent=2)+'\n')
            print(json.dumps({k:v for k,v in record.items() if k not in ('bit_errors_per_packet','config','whitener')}),flush=True)
            if adequate: break
