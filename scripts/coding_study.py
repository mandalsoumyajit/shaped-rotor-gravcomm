"""Initial coding screen, separate from uncoded validation; no rare-BER claim."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,argparse,zlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.signal import lfilter
from scipy.stats import beta
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from range_revision import make_decoder,duty,save,OUT as RANGE
from desktop_cf_fsk import encode,decode
from gravcomm.whitened_viterbi import receiver_noise
from gravcomm.coded_sequence import MaxLogByteDetector,ShortLDPC,conv_encode,conv_decode
from gravcomm.packet_confidence import packet_ber_interval
OUT=ROOT/'results/coding_study_2026_09_20'


def run(case,ts,q,count,workers,k=64):
    mechanics=duty(ts,.4,q,inertia=case['inertia_kg_m2'],fc=case['carrier_hz'],**case['drive'])
    header=dict(case=case['name'],symbol_s=ts,q=q,information_bits=k,mechanics=mechanics)
    if not mechanics['feasible']:return dict(header,status='mechanically_excluded')
    model,detector,receiver,fit=make_decoder(case,ts,q)
    soft=MaxLogByteDetector(detector);ldpc=ShortLDPC(k)
    prefix=detector.memory_symbols*model.samples_per_symbol;n=model.samples_per_symbol
    # Fixed within-codeword interleaver disperses sequential detector errors.
    permutations={length:np.random.default_rng(4000+length).permutation(length) for length in (2*k,2*(k+6))}
    def one(index,scheme,noiseless=False):
        rng=np.random.default_rng(np.random.SeedSequence([902104,zlib.crc32(case['name'].encode()),index]))
        data=rng.integers(0,2,k,dtype=np.uint8)
        word=ldpc.encode(data) if scheme=='ldpc' else conv_encode(data) if scheme=='conv171133' else data
        perm=permutations.get(len(word),np.arange(len(word)))
        symbols=encode(np.packbits(word[perm]))
        waveform=np.r_[np.ones(prefix,complex),model.waveform(np.r_[symbols,np.zeros(6,int)])]
        if not noiseless:waveform+=receiver_noise(rng,len(waveform),n/ts,case['amplitude_m_s2'],receiver=receiver)
        observed=lfilter(detector.taps,[1.],waveform)[prefix:]
        valid=True;iterations=0
        if scheme=='uncoded':
            recovered=np.unpackbits(decode(detector.decode_whitened(observed,6,True)[:-6]))
            raw=int(np.count_nonzero(recovered!=data))
        else:
            channel=soft.decode(observed)[:len(word)];llr=np.empty_like(channel);llr[perm]=channel
            raw=int(np.count_nonzero((llr<0)!=word))
            if scheme=='ldpc':
                recovered,valid,iterations=ldpc.decode(llr);recovered=recovered[:k]
            else:recovered=conv_decode(llr)
        return dict(errors=int(np.count_nonzero(recovered!=data)),raw_errors=raw,valid=bool(valid),iterations=int(iterations))
    rows=[]
    for scheme,length in [('uncoded',k),('conv171133',2*(k+6)),('ldpc',2*k)]:
        assert one(0,scheme,True)['errors']==0
        start=time.monotonic()
        with ThreadPoolExecutor(max_workers=workers) as pool:records=list(pool.map(lambda i:one(i,scheme),range(count)))
        errors=[r['errors'] for r in records];failures=sum(r['errors']>0 or not r['valid'] for r in records)
        duration=(3*((length+7)//8)+6)*ts
        per=[float(beta.ppf(.025,failures,count-failures+1)) if failures else 0.,
             float(beta.ppf(.975,failures+1,count-failures)) if failures<count else 1.]
        row=dict(scheme=scheme,packets=count,coded_bits=length,padding_bits=(-length)%8,packet_duration_s=duration,
                 payload_bit_s=k/duration,observed_goodput_bit_s=k/duration*(1-failures/count),
                 observed_ber=sum(errors)/(k*count),packet_failure_rate=failures/count,packet_failure_ci95=per,
                 packet_aware_ber_interval=packet_ber_interval(errors,bits_per_packet=k),raw_detector_ber=sum(r['raw_errors'] for r in records)/(count*length),
                 declared_failures=sum(not r['valid'] for r in records),
                 undetected_wrong_packets=sum(r['valid'] and r['errors']>0 for r in records) if scheme=='ldpc' else None,
                 packet_records=records,elapsed_s=time.monotonic()-start)
        rows.append(row);print(json.dumps(dict(case=case['name'],ts=ts,q=q,**{k:v for k,v in row.items() if k not in ('packet_records',)})),flush=True)
    return dict(header,status='screen_only',whitener=fit,states=detector.state_count,rows=rows,
                nominal_tone_span_hz=6*model.spacing_hz,
                assumptions='Acquired timing/phase; stationary physical colored noise; finite FIR whitening; max-log soft metrics; no iterative detector/decoder; interleaver within codeword only; no CRC or acquisition overhead. Tone span is not occupied bandwidth. Undetected count uses known simulated payload; syndrome pass does not certify correctness.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--packets',type=int,default=100);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--tag',default='initial');p.add_argument('--block-lengths',nargs='+',type=int,default=[64]);p.add_argument('--length-study',action='store_true');args=p.parse_args()
    cases={c['name']:c for c in json.loads((RANGE/'scenarios.json').read_text())['cases']}
    specs=[('desktop_0.5m',2.,.1),('scale2_2m_same_receiver',12.6916,.25),
           ('scale2_2m_same_receiver',6.3458,.25),('scale2_2m_slow_retuned_receiver',7.18008,.25),
           ('scale2_2m_slow_retuned_receiver',3.59004,.25)]
    if args.length_study: specs=[('scale2_2m_same_receiver',9.5,.25),('scale2_2m_slow_retuned_receiver',5.4,.25)]
    path=OUT/(args.tag+'_screen.json')
    if path.exists():raise FileExistsError('Use a fresh tag to preserve results')
    result=dict(packets_per_point=args.packets,seed=902104,status='exploratory_screen_not_validation',results=[])
    for name,ts,q in specs:
        for k in args.block_lengths:
            result['results'].append(run(cases[name],ts,q,args.packets,args.workers,k));save(path,result)

if __name__=='__main__':main()