"""Pilot coded packets with exact fixed-band observations; beam convergence audit."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json,binascii,time,hashlib,argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from gravcomm.beam_sequence import BeamByteDetector
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import noise_psd
from gravcomm.coded_sequence import ShortLDPC
from gravcomm import convolutional_candidates as tb
from gravcomm.packet_confidence import packet_ber_interval
from desktop_cf_fsk import encode,decode
from range_revision import save,duty
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/band_coded_pilot_v1'

def crc(x):return binascii.crc_hqx(bytes(x),0xffff)
def run(job):
    c,index=job;before=time.monotonic();scheme=c['scheme'];k=c['information_bits'];bits_payload=k-32
    rng=np.random.default_rng(np.random.SeedSequence([20370928,bits_payload,index]))
    payload=rng.integers(0,256,bits_payload//8,dtype=np.uint8);body=b'\x00\x01'+bytes(payload)
    info=np.unpackbits(np.frombuffer(body+crc(body).to_bytes(2,'big'),dtype=np.uint8))
    ldpc=ShortLDPC(k) if scheme=='ldpc' else None
    word=info if scheme=='uncoded' else ldpc.encode(info) if ldpc else tb.encode(info,scheme[3:])
    perm=np.random.default_rng(4000+len(word)).permutation(len(word)) if scheme!='uncoded' else np.arange(len(word))
    symbols=encode(np.packbits(word[perm]));ts=bits_payload/(len(symbols)+18)
    model=CFFSK(ts,.4,.1,16);N=32768;fs=16/ts;f=np.fft.fftfreq(N,1/fs)
    psd=noise_psd(fs,c['source']['amplitude_m_s2'],N);mask=(abs(f)<1).astype(float)
    cep=np.fft.ifft(.5*np.log(np.maximum(mask,1e-6)/psd));ca=np.zeros_like(cep);ca[0]=cep[0];ca[1:N//2]=2*cep[1:N//2];ca[N//2]=cep[N//2]
    H=np.exp(np.fft.fft(ca))*mask;imp=np.fft.ifft(H)
    transmitted=model.waveform(np.r_[symbols,np.zeros(6,int)]);start=8192
    wave=np.ones(N,complex);wave[start:start+len(transmitted)]=transmitted;tones=np.r_[0,symbols,np.zeros(6,int)]
    wave[start+len(transmitted):]=np.exp(2j*np.pi*sum(model.phase_increment_cycles(a,b) for a,b in zip(tones[:-1],tones[1:])))
    white=(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
    observed=np.fft.ifft((np.fft.fft(wave)+np.fft.fft(white)*np.sqrt(psd))*H)[start:start+len(transmitted)]
    def detect(memory,width):
        det=BeamByteDetector(model,imp[:memory*16+1],width);r=det.decode(observed,6)
        llr=np.empty(len(word));llr[perm]=r['llr'][:len(word)]
        syndrome=True
        if scheme=='uncoded':recovered=np.unpackbits(decode(r['symbols']))[:k]
        elif ldpc:recovered,syndrome,it=ldpc.decode(llr);recovered=recovered[:k]
        else:recovered=tb.decode(llr,scheme[3:],k)
        data=np.packbits(recovered);accepted=bool(syndrome and bytes(data[:2])==b'\x00\x01' and crc(data[:-2])==int.from_bytes(bytes(data[-2:]),'big'))
        errors=int(np.unpackbits(data[2:-2]^payload).sum())
        return dict(errors=errors,accepted=accepted,wrong_accepted=bool(errors and accepted),failure=bool(errors or not accepted),raw_errors=int(np.count_nonzero((llr<0)!=word)),missing=r['missing_bit_hypotheses'],repaired=r['repaired_hypotheses'],llr=llr,decoded=recovered)
    result=detect(18,8192)
    checks=[]
    if index==0:
        for memory,width in ((19,8192),(18,32768)):
            other=detect(memory,width)
            checks.append(dict(memory=memory,width=width,relative_llr_change=float(np.linalg.norm(other['llr']-result['llr'])/max(np.linalg.norm(other['llr']),1e-30)),llr_sign_changes=int(np.count_nonzero((other['llr']<0)!=(result['llr']<0))),decoded_bits_changed=int(np.count_nonzero(other['decoded']!=result['decoded'])),errors=other['errors'],accepted=other['accepted']))
    del result['llr'];del result['decoded']
    return dict(config=c,index=index,**result,convergence=checks,symbol_s=ts,payload_bits=bits_payload,coded_bits=len(word),packet_airtime_s=bits_payload,nominal_payload_bit_s=1.,elapsed_s=time.monotonic()-before)

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=12);p.add_argument('--packets',type=int,default=8);args=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    fields=json.loads((ROOT/'results/fixed_rate_sizing_2026_09_20/reconciliation_v1/fields.json').read_text())
    configs=[]
    for scale in (.95,1.):
        source=next(r for r in fields if r['scale']==scale)
        for scheme,k in [('uncoded',288),('tb_3/4',288),('tb_2/3',288),('tb_1/2',256),('ldpc',256)]:configs.append(dict(id=f's{scale}_{scheme.replace("/","_")}',source=source,scheme=scheme,information_bits=k))
    paths=['scripts/band_coded_pilot.py','src/gravcomm/beam_sequence.py','src/gravcomm/coded_sequence.py','src/gravcomm/convolutional_candidates.py']
    manifest=dict(configs=configs,packets=args.packets,seed=20370928,source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},scope='Exact common 2 Hz projection; fixed noise versus offset, tunable carrier. Acquired phase/timing, acquisition budget only. Approximate beam max-log detector; full-packet convergence checks required. Pilot not BER validation.')
    if (OUT/'manifest.json').exists():
        assert json.loads((OUT/'manifest.json').read_text())==manifest,'Inputs changed; use fresh run directory'
    else:save(OUT/'manifest.json',manifest)
    jobs=[(c,i) for c in configs for i in range(args.packets) if not (OUT/f'{c["id"]}_{i}.json').exists()]
    done=len(configs)*args.packets-len(jobs)
    def progress():save(OUT/'progress.json',dict(completed=done,total=len(configs)*args.packets,stage='complete' if done==len(configs)*args.packets else 'running',scope='pilot, not validated code ranking'))
    progress()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();save(OUT/f'{r["config"]["id"]}_{r["index"]}.json',r);done+=1;progress();print(json.dumps(dict(done=done,total=len(configs)*args.packets,id=r['config']['id'],index=r['index'],errors=r['errors'],seconds=r['elapsed_s'])),flush=True)
    summaries=[]
    for c in configs:
        rows=[json.loads((OUT/f'{c["id"]}_{i}.json').read_text()) for i in range(args.packets)];n=rows[0]['payload_bits'];ts=rows[0]['symbol_s'];s=c['source']['scale']
        summaries.append(dict(config=c,packets=len(rows),ber=sum(r['errors'] for r in rows)/(n*len(rows)),ber_ci95=packet_ber_interval([r['errors'] for r in rows],n),packet_failure_rate=sum(r['failure'] for r in rows)/len(rows),wrong_accepted=sum(r['wrong_accepted'] for r in rows),convergence=rows[0]['convergence'],mechanics=duty(ts,.4,.1,inertia=c['source']['inertia_kg_m2'],fc=c['source']['carrier_hz'],drag=.05*s**3,max_rpm=1800/s),symbol_s=ts,airtime_s=n))
    save(OUT/'summary.json',summaries);print('BAND_CODED_PILOT_COMPLETE',flush=True)
if __name__=='__main__':main()
