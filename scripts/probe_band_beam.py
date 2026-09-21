import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json,time
from pathlib import Path
import numpy as np
sys.path[:0]=['AIP_Advances_Submission/src','AIP_Advances_Submission/scripts']
from gravcomm.beam_sequence import BeamByteDetector
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import noise_psd
from desktop_cf_fsk import encode,decode
from range_revision import save
out=Path('AIP_Advances_Submission/results/fixed_rate_sizing_2026_09_20/bandwidth_audit_v1')
records=[]
for ts in (256/126,224/210):
 N=32768;fs=16/ts;f=np.fft.fftfreq(N,1/fs);psd=noise_psd(fs,5.180259430854157e-10,N);mask=(abs(f)<1).astype(float)
 cep=np.fft.ifft(.5*np.log(np.maximum(mask,1e-6)/psd));ca=np.zeros_like(cep);ca[0]=cep[0];ca[1:N//2]=2*cep[1:N//2];ca[N//2]=cep[N//2]
 H=np.exp(np.fft.fft(ca))*mask;imp=np.fft.ifft(H)
 model=CFFSK(ts,.4,.1,16);rng=np.random.default_rng(28);data=rng.integers(0,256,8,dtype=np.uint8);symbols=encode(data);body=model.waveform(np.r_[symbols,np.zeros(6,int)])
 start=8192;wave=np.ones(N,complex);wave[start:start+len(body)]=body;total=np.r_[0,symbols,np.zeros(6,int)]
 wave[start+len(body):]=np.exp(2j*np.pi*sum(model.phase_increment_cycles(a,b) for a,b in zip(total[:-1],total[1:])))
 white=(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
 observed=np.fft.ifft((np.fft.fft(wave)+np.fft.fft(white)*np.sqrt(psd))*H)[start:start+len(body)]
 for mem in (8,12):
  for width in (128,512,2048):
   before=time.monotonic();d=BeamByteDetector(model,imp[:mem*16+1],width);r=d.decode(observed,6);hard=decode(r['symbols'])
   row=dict(ts=ts,memory=mem,width=width,errors=int(np.unpackbits(hard^data).sum()),missing=r['missing_bit_hypotheses'],llr=r['llr'].tolist(),hard=hard.tolist(),seconds=time.monotonic()-before)
   records.append(row);save(out/'beam_probe.json',records);print(json.dumps({k:v for k,v in row.items() if k not in ('llr','hard')}),flush=True)
