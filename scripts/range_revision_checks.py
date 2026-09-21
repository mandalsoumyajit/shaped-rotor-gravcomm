"""Independent sampling/filter diagnostics for frozen range candidates.

Diagnostics are not new optimized candidates or accurate rare-event BER claims.
Filter comparisons use the identical received waveform for both filters.
"""
import argparse,json,time,zlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.signal import lfilter
from range_revision import OUT,TAIL,make_decoder,trials,save
from desktop_cf_fsk import encode,decode
from gravcomm.whitened_viterbi import receiver_noise


def main():
    p=argparse.ArgumentParser();p.add_argument('--cases',nargs='+',default=['desktop_2m','scale2_2m_slow_retuned_receiver'])
    p.add_argument('--tag',default='');p.add_argument('--from-validation',action='store_true')
    p.add_argument('--packets',type=int,default=1000);p.add_argument('--workers',type=int,default=4);args=p.parse_args()
    cases=json.loads((OUT/'scenarios.json').read_text())['cases']
    screen=json.loads((OUT/'screen.json').read_text())
    records=[]
    for case in cases:
        if case['name'] not in args.cases:continue
        eligible=[r for r in screen if r['case']==case['name'] and r.get('observed_ber',1)<=.0005]
        eligible.sort(key=lambda r:(r['symbol_s'],r['observed_ber'],r['states']))
        if not eligible:continue
        row=eligible[0]
        if args.from_validation:
            suffix=('_'+args.tag) if args.tag else ''
            row=json.loads((OUT/(case['name']+suffix+'_validation.json')).read_text())
        ts,q=row['symbol_s'],row['q']
        model,short,rec,fitshort=make_decoder(case,ts,q,16,3)
        _,long,_,fitlong=make_decoder(case,ts,q,16,4)
        prefix=4*16;seed=zlib.crc32(case['name'].encode());a=case['amplitude_m_s2']
        def one(index):
            rng=np.random.default_rng(np.random.SeedSequence([910732,seed,index]))
            data=rng.integers(0,256,8,dtype=np.uint8);symbols=encode(data)
            signal=np.r_[np.ones(prefix,complex),model.waveform(np.r_[symbols,np.zeros(TAIL,int)])]
            observed=signal+receiver_noise(rng,len(signal),16/ts,a,receiver=rec)
            recovered=[]
            for dec in (short,long):
                decoded=dec.decode_whitened(lfilter(dec.taps,[1.],observed)[prefix:],TAIL,True)[:24]
                recovered.append(decode(decoded))
            return [int(np.unpackbits(data^v).sum()) for v in recovered]+[int(np.unpackbits(recovered[0]^recovered[1]).sum())]
        started=time.monotonic()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:results=np.array(list(pool.map(one,range(args.packets))))
        highmodel,high,highrec,fithigh=make_decoder(case,ts,q,32,3)
        higher=trials(case,highmodel,high,highrec,0,args.packets,982014,args.workers)
        record=dict(case=case['name'],symbol_s=ts,q=q,packets=args.packets,payload_bits=64*args.packets,
                    paired_filter_bit_errors_3symbols=int(results[:,0].sum()),paired_filter_bit_errors_4symbols=int(results[:,1].sum()),
                    paired_filter_disagreeing_bits=int(results[:,2].sum()),independent_double_sampling_bit_errors=sum(higher),
                    short_filter=fitshort,long_filter=fitlong,double_sample_filter=fithigh,
                    elapsed_s=time.monotonic()-started,
                    scope='Paired memory check and independently seeded sample-rate diagnostic; no rare-event confidence claim.')
        records.append(record);save(OUT/('sampling_filter_checks'+('_'+args.tag if args.tag else '')+'.json'),records);print(json.dumps(record),flush=True)

if __name__=='__main__':main()
