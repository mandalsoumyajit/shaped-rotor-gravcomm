"""Audit exact 2 Hz Fourier projection against finite-memory sequence metrics.
No normalization to unit per-sample variance: projected noise has unit PSD
in retained bins and zero PSD outside. Its covariance is a projection.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from gravcomm.whitened_viterbi import noise_psd,WhitenedViterbi
from gravcomm.sequence_decoder import CFFSK
from gravcomm.coded_sequence import MaxLogByteDetector
from desktop_cf_fsk import encode,decode
from range_revision import save
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/bandwidth_audit_v1'
A=5.180259430854157e-10

def design(ts,memory,nfft=32768):
    fs=16/ts;f=np.fft.fftfreq(nfft,1/fs);psd=noise_psd(fs,A,nfft)
    band=(abs(f)<1).astype(float)
    ideal=np.sqrt(band/psd)
    impulse=np.fft.ifft(ideal);half=memory*16//2
    taps=impulse[np.arange(-half,half+1)%nfft]
    exact=ideal*np.exp(-2j*np.pi*f*half/fs)
    return taps,exact,psd,band,fs

def main():
    OUT.mkdir(exist_ok=True)
    rows=[];packet_checks=[]
    for ts in (256/126,224/147,224/210):
        model=CFFSK(ts,.4,.1,16)
        rng=np.random.default_rng(20370926)
        symbols=encode(rng.integers(0,256,256,dtype=np.uint8))
        body=model.waveform(np.r_[symbols,np.zeros(6,int)])
        N=32768;start=8192
        wave=np.ones(N,complex);wave[start:start+len(body)]=body
        phase=sum(model.phase_increment_cycles(a,b) for a,b in zip(np.r_[0,symbols,np.zeros(6,int)][:-1],np.r_[0,symbols,np.zeros(6,int)][1:]))
        wave[start+len(body):]=np.exp(2j*np.pi*phase)
        for memory in (3,4,5,8,12,16):
            taps,H,psd,band,fs=design(ts,memory)
            finite=np.fft.fft(taps,N)
            exact_signal=np.fft.ifft(np.fft.fft(wave)*H)
            finite_signal=np.fft.ifft(np.fft.fft(wave)*finite)
            region=slice(start,start+len(body))
            signal_error=float(np.linalg.norm((finite_signal-exact_signal)[region])/np.linalg.norm(exact_signal[region]))
            projected_noise_psd=abs(finite)**2*psd
            rows.append(dict(symbol_s=ts,memory_symbols=memory,states_estimate=12*7**(memory+1),relative_signal_l2_error=signal_error,relative_noise_projection_l2_error=float(np.linalg.norm(projected_noise_psd-band)/np.linalg.norm(band)),out_of_band_noise_fraction=float(projected_noise_psd[band==0].sum()/projected_noise_psd.sum())))
        # Exact Fourier-domain construction checks: no out-of-band output;
        # whitening+projection is applied to BOTH source and physical noise.
        taps,H,psd,band,fs=design(ts,3)
        white=(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
        rawnoise=np.fft.ifft(np.fft.fft(white)*np.sqrt(psd))
        filterednoise=np.fft.ifft(np.fft.fft(rawnoise)*H)
        direct=np.fft.ifft(np.fft.fft(white)*np.sqrt(band)*np.exp(-2j*np.pi*np.fft.fftfreq(N)*24))
        np.testing.assert_allclose(filterednoise,direct,atol=1e-12,rtol=1e-11)
        # Parseval verifies likelihood norm on retained frequency bins.
        test=(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
        lhs=np.vdot(np.fft.ifft(np.fft.fft(test)*H),np.fft.ifft(np.fft.fft(test)*H)).real
        rhs=float((abs(np.fft.fft(test))**2*band/psd).sum()/N)
        np.testing.assert_allclose(lhs,rhs,rtol=1e-12)
    save(OUT/'filter_metrics.json',rows)
    # Small independent, exact-band observation tests. Full source and noise
    # are projected first; only trellis templates approximate channel memory.
    for ts in (256/126,224/210):
        model=CFFSK(ts,.4,.1,16);rng=np.random.default_rng(20370927)
        payload=rng.integers(0,256,8,dtype=np.uint8);symbols=encode(payload)
        body=model.waveform(np.r_[symbols,np.zeros(6,int)]);N=32768;start=8192
        wave=np.ones(N,complex);wave[start:start+len(body)]=body
        total=np.r_[0,symbols,np.zeros(6,int)]
        wave[start+len(body):]=np.exp(2j*np.pi*sum(model.phase_increment_cycles(a,b) for a,b in zip(total[:-1],total[1:])))
        white=(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
        records=[]
        for memory in (3,4):
            taps,H,psd,band,fs=design(ts,memory)
            decoder=WhitenedViterbi(model,taps,maximum_states=2000000)
            soft=MaxLogByteDetector(decoder)
            exact_clean=np.fft.ifft(np.fft.fft(wave)*H)[start:start+len(body)]
            recovered=decode(decoder.decode_whitened(exact_clean,6,True)[:-6]);assert np.array_equal(recovered,payload)
            obs=np.fft.ifft((np.fft.fft(wave)+np.fft.fft(white)*np.sqrt(psd))*H)[start:start+len(body)]
            hard=decode(decoder.decode_whitened(obs,6,True)[:-6]);llr=soft.decode(obs,6)
            records.append(dict(memory=memory,states=decoder.state_count,errors=int(np.unpackbits(hard^payload).sum()),hard=hard.tolist(),llr=llr.tolist()))
        a,b=records
        packet_checks.append(dict(symbol_s=ts,records=records,hard_equal=a['hard']==b['hard'],llr_sign_disagreements=int(np.count_nonzero((np.array(a['llr'])<0)!=(np.array(b['llr'])<0))),relative_llr_change=float(np.linalg.norm(np.array(a['llr'])-b['llr'])/np.linalg.norm(b['llr']))))
        print(json.dumps({k:v for k,v in packet_checks[-1].items() if k!='records'}),flush=True)
    result=dict(filter_metrics=rows,packet_checks=packet_checks,checks=['CRC-independent numerical audit','Exact common projection of physical signal and physical colored noise','Fourier whitening identity','Parseval likelihood identity','Noiseless exact-band observations decode at memories 3 and 4'],status='audit_only_no_code_ranking',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    save(OUT/'audit.json',result)
    lines=['# Fixed 2 Hz receive-band audit','','Signal and colored noise undergo the same exact Fourier projection and noise weighting. The retained baseband is (-1,1) Hz. Output noise has unit spectral eigenvalues in-band and zero outside; it must not be incorrectly normalized to independent unit-variance samples. Parseval and direct noise-projection checks pass. Receiver PSD remains fixed versus carrier offset.','','Finite-memory sequence templates approximate the resulting long impulse response. This audit measures that approximation before any code ranking.','','| Symbol duration (s) | Memory (symbols) | Relative waveform L2 error | Relative noise-projection L2 error |','|---:|---:|---:|---:|']
    for r in rows:lines.append(f"| {r['symbol_s']:.6g} | {r['memory_symbols']} | {r['relative_signal_l2_error']:.4g} | {r['relative_noise_projection_l2_error']:.4g} |")
    lines+=['','| Symbol duration (s) | Hard decisions agree, memory 3/4 | LLR sign changes | Relative LLR change |','|---:|---|---:|---:|']
    for r in packet_checks:lines.append(f"| {r['symbol_s']:.6g} | {r['hard_equal']} | {r['llr_sign_disagreements']} | {r['relative_llr_change']:.4g} |")
    lines+=['','These small diagnostics are not BER validation. Exact hard bandwidth truncation introduces long channel memory; the previous three-symbol whitening trellis is not automatically adequate. Code comparisons require metric convergence or a different sequence-detection algorithm. No coded ranking is generated by this audit.']
    (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('BANDWIDTH_AUDIT_COMPLETE',flush=True)
if __name__=='__main__':main()
