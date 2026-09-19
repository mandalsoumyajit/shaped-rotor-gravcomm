"""Reproducible full-packet Viterbi CF-FSK search with physical colored noise.

Mapped rates exclude acquisition/coding overhead. The screen is followed by
independent seeds and longer whitening/sample-rate convergence checks. No
binomial confidence interval treats correlated symbol errors as independent.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.stats import beta
import desktop_cf_fsk as cf
from optimize_cf_fsk import trajectory,mechanics
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import fit_whitener,WhitenedViterbi

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/sequence_cf_fsk';OUT.mkdir(exist_ok=True)
A=5.180259430854158e-10


def packet_trials(model,decoder,count,seed,bytes_per_packet=8):
    rng=np.random.default_rng(seed);errors=[];byte_errors=[];probes=[]
    tail=max(6,decoder.memory_symbols+2);by_position=np.zeros(3*bytes_per_packet,int)
    for _ in range(count):
        data=rng.integers(0,256,bytes_per_packet)
        symbols=cf.encode(data);decoded=decoder.observe_and_decode(symbols,rng,A,tail_symbols=tail)
        bad=decoded!=symbols
        errors.append(int(sum(bad)))
        by_position+=bad
        byte_errors.append(int(np.sum(np.any(bad.reshape(-1,3),axis=1))))
        # One random symbol per independent packet supplies a valid binomial
        # sample of average symbol-error probability, despite within-packet memory.
        probes.append(int(bad[rng.integers(len(bad))]))
    failed=int(np.count_nonzero(errors));probe_errors=sum(probes)
    upper=lambda k:float(beta.ppf(.95,k+1,count-k)) if k<count else 1.
    return dict(packets=count,bytes_per_packet=bytes_per_packet,tail_symbols=tail,
        packet_payload_bit_s=8*bytes_per_packet/((3*bytes_per_packet+tail)*model.symbol_s),
        errors_by_position=by_position.tolist(),symbol_errors=sum(errors),
        symbols=count*3*bytes_per_packet,ser=sum(errors)/(count*3*bytes_per_packet),
        packet_errors=failed,packet_error_upper95=upper(failed),byte_errors=sum(byte_errors),
        independent_symbol_probe_errors=probe_errors,average_ser_upper95=upper(probe_errors),
        symbol_errors_per_packet=errors)


def evaluate(ts,r,q,packets,seed,filter_seconds=6.,sample_rate=2.):
    n=int(round(ts*sample_rate));model=CFFSK(ts,r,q,n)
    fs=n/ts;memory=int(np.ceil(filter_seconds*fs))
    taps,fit=fit_whitener(fs,A,memory)
    decoder=WhitenedViterbi(model,taps,maximum_states=1000000)
    result=packet_trials(model,decoder,packets,seed)
    return dict(ts=ts,transition_fraction=r,spacing_dwell_product=q,spacing_hz=model.spacing_hz,
        rate=8/(3*ts),whitener=fit,states=decoder.state_count,trial=result)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    parser.add_argument('--refine',action='store_true')
    parser.add_argument('--packets',type=int,default=100);args=parser.parse_args()
    if args.verify:
        rows=json.loads((OUT/'screen.json').read_text())
        viable=[x for x in rows if x['trial']['ser']<=.003]
        viable.sort(key=lambda x:(x['ts'],x['trial']['ser'],x['states']))
        selected=[]
        # Several candidates survive the small preliminary screen. Their next
        # independently seeded stage determines the final fixed candidate.
        for index,row in enumerate(viable[:8]):
            result=evaluate(row['ts'],row['transition_fraction'],row['spacing_dwell_product'],args.packets,87122+index,8.)
            result['mechanics']=row['mechanics'];selected.append(result)
            print(json.dumps(result),flush=True)
            (OUT/'validation.json').write_text(json.dumps(selected,indent=2)+'\n')
        return
    rows=[];started=time.monotonic()
    periods=np.arange(2.,3.51,.25) if args.refine else np.arange(2.,8.01,.5)
    fractions=(.4,.5) if args.refine else (.2,.4,.5)
    spacings=(.0625,.09375,.125,.1875,.25) if args.refine else (.125,.25,.5,1.)
    filename='refinement.json' if args.refine else 'screen.json'
    for ts in periods:
        for r in fractions:
            for q in spacings:
                # Keep tone centres inside the original 1.5625 Hz receiver band.
                if 3*q/((1-r)*ts)>.78125:continue
                t,f,df,_,_=trajectory(ts,r,q);mech=mechanics(ts,r,q,t,f,df)
                if not mech['feasible']:continue
                try:result=evaluate(float(ts),r,q,8,17876)
                except ValueError as exc:
                    if 'state' in str(exc):continue
                    raise
                result['mechanics']=mech;rows.append(result)
                print(json.dumps({k:result[k] for k in ('ts','transition_fraction','spacing_dwell_product','states','trial')}),flush=True)
                (OUT/filename).write_text(json.dumps(rows,indent=2)+'\n')
        print(f'Completed Ts={ts:g}; elapsed {time.monotonic()-started:.1f}s',flush=True)

if __name__=='__main__':main()
