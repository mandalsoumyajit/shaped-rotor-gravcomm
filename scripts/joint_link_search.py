"""Checkpointed joint-link pilot. Conditional models; no raw-BER eligibility gate."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,sys,time,platform
import scipy,numba
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
from scipy.stats import beta
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from range_revision import make_decoder,duty,save
from gravcomm.coded_sequence import MaxLogByteDetector,ShortLDPC,conv_encode,conv_decode
from gravcomm.whitened_viterbi import receiver_noise
from gravcomm import convolutional_candidates as tb
from desktop_cf_fsk import encode,decode
OUT=ROOT/'results/joint_link_study_2026_09_20'


def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()


def source_hashes():
    paths=['scripts/joint_link_search.py','scripts/desktop_cf_fsk.py','scripts/range_revision.py','src/gravcomm/coded_sequence.py',
           'src/gravcomm/whitened_viterbi.py','src/gravcomm/sequence_decoder.py','src/gravcomm/receivers.py']
    extra=ROOT/'src/gravcomm/convolutional_candidates.py'
    if extra.exists():paths.append(str(extra.relative_to(ROOT)))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def manifest(bank,tag,packets,seed,crossed=False):
    configs=[]
    for case in bank['cases']:
        if crossed and case['receiver']['quality_factor']<600000:continue
        cap=case.get('gaussian_benchmark_bit_s',0)
        if cap<=0:raise ValueError('Receiver bank must supply positive benchmark as scale-only starting seed')
        base=max(1.,8/(3*cap))
        # Coarse starting seeds, not a capacity-ratio target or optimized grid.
        for q,factor in (((.1,1.45),(.25,1.7)) if crossed else ((.1,.9),(.25,1.7))):
            ts=max(1.,base*factor)
            for _ in range(40):
                mech=duty(ts,.4,q,inertia=case['inertia_kg_m2'],fc=case['carrier_hz'],**case['drive'])
                if mech['feasible']:break
                ts*=1.2
            else:continue
            for k in ((64,128) if crossed else ((64,) if q==.1 else (128,))):
                for scheme in ('uncoded','ldpc','tb_1/2','tb_2/3','tb_3/4'):
                    c=dict(case=case,symbol_s=float(f'{ts:.8g}'),q=q,information_bits=k,scheme=scheme,
                           samples=16,memory=3,tail=6,seed=seed,packets=packets)
                    c['id']=digest(c)[:16];configs.append(c)
    return dict(tag=tag,receiver_bank_sha256=digest(bank),source_sha256=source_hashes(),configs=configs,
                stage='crossed_screen' if crossed else 'pilot_screen',crossed=crossed,notes='No raw-BER gate. Periods seeded by capacity only for scale, then mechanical pruning; no optimality claims. Receiver statuses distinguish sensitivity designs. No CRC, acquisition, achieved-drive or phase-drift model yet.')


def interval(events,n):
    return [float(beta.ppf(.025,events,n-events+1)) if events else 0.,
            float(beta.ppf(.975,events+1,n-events)) if events<n else 1.]


def run_config(c,workers,context):
    case=c['case'];ts=c['symbol_s'];q=c['q'];k=c['information_bits'];scheme=c['scheme']
    key=(digest(case),ts,q,c['samples'],c['memory'])
    if context.get('key')!=key:
        model,det,rec,fit=make_decoder(case,ts,q,c['samples'],c['memory'])
        context.clear();context.update(key=key,model=model,det=det,rec=rec,fit=fit,soft=MaxLogByteDetector(det))
    model,det,rec,soft=context['model'],context['det'],context['rec'],context['soft']
    ldpc=ShortLDPC(k) if scheme=='ldpc' else None
    n=model.samples_per_symbol;prefix=det.memory_symbols*n
    code_n=tb.metadata(k,scheme[3:])['n'] if scheme.startswith('tb_') else 2*k if scheme=='ldpc' else 2*(k+6) if scheme=='conv_terminated' else k
    perm=np.random.default_rng(4000+code_n).permutation(code_n) if scheme!='uncoded' else np.arange(k)
    def one(index,noiseless=False):
        # Same payload/base noise across schemes at a given case/k/index, for screening only.
        rng=np.random.default_rng(np.random.SeedSequence([c['seed'],int(digest(case)[:8],16),k,index]))
        bits=rng.integers(0,2,k,dtype=np.uint8)
        word=tb.encode(bits,scheme[3:]) if scheme.startswith('tb_') else ldpc.encode(bits) if scheme=='ldpc' else conv_encode(bits) if scheme=='conv_terminated' else bits
        symbols=encode(np.packbits(word[perm]));total=np.r_[symbols,np.zeros(c['tail'],int)]
        wave=np.r_[np.ones(prefix,complex),model.waveform(total)]
        if not noiseless:wave+=receiver_noise(rng,len(wave),n/ts,case['amplitude_m_s2'],receiver=rec)
        observed=lfilter(det.taps,[1.],wave)[prefix:]
        valid=True;iterations=0
        if scheme=='uncoded':
            recovered=np.unpackbits(decode(det.decode_whitened(observed,c['tail'],True)[:len(symbols)]))
            raw=np.count_nonzero(recovered!=bits)
        else:
            channel=soft.decode(observed,c['tail'])[:len(word)];llr=np.empty_like(channel);llr[perm]=channel
            raw=np.count_nonzero((llr<0)!=word)
            if scheme=='ldpc':
                recovered,valid,iterations=ldpc.decode(llr);recovered=recovered[:k]
            else:recovered=tb.decode(llr,scheme[3:],k) if scheme.startswith('tb_') else conv_decode(llr)
        errors=int(np.count_nonzero(recovered!=bits))
        return dict(errors=errors,raw_errors=int(raw),accepted=bool(valid),failure=bool(errors or not valid),iterations=int(iterations))
    assert one(0,True)['errors']==0
    start=time.monotonic()
    with ThreadPoolExecutor(max_workers=workers) as pool:records=list(pool.map(one,range(c['packets'])))
    count=len(records);fail=sum(r['failure'] for r in records);accepted=sum(r['accepted'] for r in records)
    duration=(3*((code_n+7)//8)+c['tail'])*ts
    return dict(config=c,packets=count,records=records,postdecode_ber=sum(r['errors'] for r in records)/(k*count),
                accepted_packet_fraction=accepted/count,
                accepted_ber=sum(r['errors'] for r in records if r['accepted'])/(k*accepted) if accepted else None,
                packet_failure_rate=fail/count,packet_failure_ci95=interval(fail,count),
                wrong_accepted_packets=sum(r['accepted'] and r['errors']>0 for r in records),
                error_detection='LDPC syndrome only' if scheme=='ldpc' else 'none; all outputs accepted',
                raw_ber=sum(r['raw_errors'] for r in records)/(code_n*count),
                payload_bit_s=k/duration,observed_correct_payload_bit_s=k/duration*(1-fail/count),packet_duration_s=duration,
                coded_bits=code_n,padding_bits=(-code_n)%8,crc_bits=0,states=det.state_count,whitener=context['fit'],
                mechanics=duty(ts,.4,q,inertia=case['inertia_kg_m2'],fc=case['carrier_hz'],**case['drive']),
                elapsed_s=time.monotonic()-start,status='pilot_only_not_reliability_validation')


def progress(run,m,done,running,status):
    rows=[]
    for path in sorted(run.glob('case_*.json')):
        r=json.loads(path.read_text());c=r['config']
        rows.append(dict(id=c['id'],distance_m=c['case']['distance_m'],receiver=c['case'].get('receiver_id'),
                         receiver_status=c['case'].get('receiver_status'),scheme=c['scheme'],information_bits=c['information_bits'],
                         symbol_s=c['symbol_s'],q=c['q'],postdecode_ber=r['postdecode_ber'],packet_failure_rate=r['packet_failure_rate'],
                         payload_bit_s=r['payload_bit_s'],packet_duration_s=r['packet_duration_s'],packets=r['packets']))
    save(run/'progress.json',dict(status=status,total=len(m['configs']),completed=done,running_id=running,rows=rows))
    lines=['# Joint-link pilot run progress','',f'Status: {status}; {done}/{len(m["configs"])} configurations complete.',
           '', 'All results are exploratory acquired-packet screens, not validated reliability claims.', '',
           '| Distance (m) | Receiver | Code | Payload bits | Ts (s) | BER after decoding | Packet failure | Payload bit/s |',
           '|---:|---|---|---:|---:|---:|---:|---:|']
    for r in rows:lines.append(f'| {r["distance_m"]} | {r["receiver"]} | {r["scheme"]} | {r["information_bits"]} | {r["symbol_s"]:.4g} | {r["postdecode_ber"]:.3g} | {r["packet_failure_rate"]:.3g} | {r["payload_bit_s"]:.4g} |')
    (run/'PROGRESS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser();p.add_argument('--tag',default='pilot_v1');p.add_argument('--packets',type=int,default=16)
    p.add_argument('--crossed',action='store_true');p.add_argument('--workers',type=int,default=4);p.add_argument('--limit',type=int);p.add_argument('--seed',type=int,default=20370920)
    args=p.parse_args()
    if not args.tag.isidentifier():raise ValueError('Safe identifier required')
    bank=json.loads((OUT/'receiver_bank.json').read_text());run=OUT/args.tag;run.mkdir(exist_ok=True)
    path=run/'manifest.json'
    if path.exists():
        if args.limit:raise ValueError('--limit is creation-only; omit on resume')
        m=json.loads(path.read_text())
        if m['receiver_bank_sha256']!=digest(bank) or m['source_sha256']!=source_hashes():raise ValueError('Inputs/code changed; use new tag')
        if any(c['packets']!=args.packets or c['seed']!=args.seed for c in m['configs']):raise ValueError('Run settings changed')
    else:
        m=manifest(bank,args.tag,args.packets,args.seed,args.crossed)
        if args.limit:m['configs']=m['configs'][:args.limit]
        save(path,m)
    save(run/'runtime.json',dict(python=sys.version,platform=platform.platform(),host=platform.node(),numpy=np.__version__,scipy=scipy.__version__,numba=numba.__version__,workers=args.workers,source_sha256=source_hashes()))
    context={};done=0
    progress(run,m,done,None,'running')
    for c in m['configs']:
        output=run/('case_'+c['id']+'.json')
        if output.exists():done+=1;continue
        progress(run,m,done,c['id'],'running')
        r=run_config(c,args.workers,context);save(output,r);done+=1
        print(json.dumps(dict(done=done,total=len(m['configs']),case=c['case']['name'],scheme=c['scheme'],k=c['information_bits'],ber=r['postdecode_ber'],per=r['packet_failure_rate'])),flush=True)
    progress(run,m,done,None,'pilot_complete')

if __name__=='__main__':main()