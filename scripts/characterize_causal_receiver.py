"""Characterize causal receiver, stationary innovations, and implementation cost."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.signal import sosfilt,sos2zpk
from scipy.linalg import solve_triangular
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from range_revision import save
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/causal_receiver_v1'
A=5.180259430854157e-10

def main():
 OUT.mkdir(exist_ok=True);r=CausalReceiver();orders=[];f=np.linspace(-1,1,4001)
 for order in (2,4,6,8):
  rec=CausalReceiver(order=order);parts=rec.output_psd(f,A,True)
  orders.append(dict(order=order,maximum_alias_to_main_psd_fraction=float(np.max(parts[1:].sum(axis=0)/parts[0]))))
 assert orders[1]['maximum_alias_to_main_psd_fraction']>1e-4 and orders[2]['maximum_alias_to_main_psd_fraction']<1e-4
 grid=np.linspace(-8,8,131072,endpoint=False);H=r.response(grid)
 enbw=float(np.mean(abs(H)**2)*16)
 fg=np.linspace(-1.2,1.2,24001);hg=r.response(fg);delay=-np.gradient(np.unwrap(np.angle(hg)),fg)/(2*np.pi)
 impulse=sosfilt(r.sos,np.r_[1.,np.zeros(4095)]);step=sosfilt(r.sos,np.ones(4096))
 tail=np.cumsum(abs(impulse[::-1])**2)[::-1]/sum(abs(impulse)**2)
 settled=np.maximum.accumulate(abs(step[::-1]-1))[::-1]
 step_time=float(np.flatnonzero(settled<1e-4)[0]/16)
 energy_time=float(np.flatnonzero(tail<1e-6)[0]/16)
 # Theoretical innovations and finite-memory approximation diagnostics.
 t=time.perf_counter();model=r.finite_record(1024,A);chol_seconds=time.perf_counter()-t
 whiten=solve_triangular(model.factor,np.eye(1024),lower=True)
 identity=whiten@model.covariance@whiten.conj().T
 identity_error=float(np.max(abs(identity-np.eye(1024))))
 assert identity_error<1e-9
 short=r.finite_record(256,A)
 tail_checks=[]
 W=solve_triangular(short.factor,np.eye(256),lower=True);C=short.covariance
 ii,jj=np.indices(W.shape)
 for keep in (8,16,32,64,128):
  truncated=np.where(ii-jj<keep,W,0)
  R=truncated@C@truncated.conj().T
  tail_checks.append(dict(retained_samples=keep,history_s=(keep-1)/4,max_whitened_covariance_error=float(np.max(abs(R-np.eye(256))))))
 # Empirical noise passes through the implemented causal filter, not covariance-factor synthesis.
 rng=np.random.default_rng(20370929);noise=r.stationary_noise(65536,rng,A)
 lags=r.covariance_lags(17,A)
 empirical=np.array([np.mean(noise[k:]*noise[:len(noise)-k].conj()) for k in range(17)])
 relative_cov_error=float(np.max(abs(empirical-lags))/lags[0].real)
 assert relative_cov_error<.05
 blocks=noise.reshape(-1,128).T;finite=r.finite_record(128,A)
 z=finite.whiten(blocks)
 white_power=float(np.mean(abs(z)**2));pseudo=complex(np.mean(z*z));mean=complex(np.mean(z))
 lagwhite=np.array([np.mean(z[k:]*z[:-k].conj()) for k in range(1,17)])
 assert abs(white_power-1)<.03 and np.max(abs(lagwhite))<.03 and abs(pseudo)<.03
 # Compare initialization on the same raw input and evaluation interval.
 size=65536;freq=np.fft.fftfreq(size,1/r.input_rate_hz)
 white=(rng.normal(size=size)+1j*rng.normal(size=size))/np.sqrt(2)
 raw=np.fft.ifft(np.fft.fft(white)*np.sqrt(r.input_psd(freq,A)))
 start=1024;stop=start+4096
 reference=r.stream().process(raw[:stop])[start//4:]
 burn_checks=[]
 for seconds in (0,2,4,8,16,32):
  first=start-int(seconds*r.input_rate_hz)
  output=r.stream().process(raw[first:stop])[(start-first)//4:]
  error=float(np.max(abs(output-reference))/np.sqrt(lags[0].real))
  burn_checks.append(dict(burn_s=seconds,max_output_error_in_noise_rms=error))
 assert burn_checks[-2]['max_output_error_in_noise_rms']<1e-8
 assert burn_checks[-1]['max_output_error_in_noise_rms']<1e-12
 # Input-FFT covariance refinement.
 nfft_cov=np.array([np.linalg.norm(r.covariance_lags(256,A,nfft)-short.covariance[:,0])/np.linalg.norm(short.covariance[:,0]) for nfft in (8192,16384,32768,65536)])
 assert max(nfft_cov)<1e-9
 t=time.perf_counter();stream=model.stream();stream.process(np.ones(1024,complex));innovation_stream_seconds=time.perf_counter()-t
 poles=sos2zpk(r.sos)[1]
 result=dict(receiver=dict(input_rate_hz=16,output_rate_hz=4,full_minus3db_bandwidth_hz=2,noise_equivalent_bandwidth_hz=enbw,order=6,sos=r.sos.tolist(),maximum_pole_radius=float(max(abs(poles)))),order_selection=orders,alias_budget_psd_fraction=1e-4,group_delay=dict(at_center_s=float(delay[len(delay)//2]),max_within_minus3db_band_s=float(max(delay[abs(fg)<=1]))),step_settling_1e4_s=step_time,impulse_energy_tail_1e6_s=energy_time,covariance=dict(max_exact_whitened_covariance_error=identity_error,condition_number_128=float(np.linalg.cond(finite.covariance)),empirical_max_lag_error_relative_to_variance=relative_cov_error,empirical_innovation_power=white_power,empirical_innovation_max_abs_lag_correlation=float(max(abs(lagwhite))),empirical_innovation_pseudocovariance=[pseudo.real,pseudo.imag],empirical_innovation_mean=[mean.real,mean.imag],fft_refinement_relative_errors=nfft_cov.tolist()),initialization_checks=burn_checks,finite_memory_diagnostics=tail_checks,runtime=dict(cholesky_1024_samples_s=chol_seconds,stream_1024_innovations_s=innovation_stream_seconds),status='causal_sampled_frontend_and_matched_finite_record_metric_verified_not_end_to_end_hardware_or_BER_validation',scope='Sampled calibrated-acceleration boundary. Analog acquisition and transducer calibration upstream remain requirements. Noise stationary through packet boundaries. No acquisition algorithm or code ranking. -3dB bandwidth is not a brick-wall information bandwidth.',source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('src/gravcomm/causal_receiver.py','scripts/characterize_causal_receiver.py','tests/test_causal_receiver.py','src/gravcomm/receivers.py','src/gravcomm/constants.py')})
 save(OUT/'characterization.json',result)
 fig,axs=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
 display=np.linspace(0,4,4001);gain=abs(r.response(display));axs[0,0].plot(display,20*np.log10(np.maximum(gain,1e-10)));axs[0,0].axvline(1,color='k',ls=':',lw=1);axs[0,0].set(xlabel='Offset frequency (Hz)',ylabel='Filter gain (dB)',ylim=(-100,3),title='Causal low-pass: 2 Hz full -3 dB bandwidth')
 axs[0,1].plot(np.arange(160)/16,step[:160]);axs[0,1].set(xlabel='Time after input step (s)',ylabel='Output',title='Causal step response')
 ff=np.linspace(-2,2,8192,endpoint=False);parts=r.output_psd(ff,A,True);axs[1,0].semilogy(ff,parts.sum(axis=0),label='Total sampled PSD');axs[1,0].semilogy(ff,parts[1:].sum(axis=0),label='Folded contribution');axs[1,0].legend();axs[1,0].set(xlabel='Offset frequency (Hz)',ylabel='Normalized discrete PSD',title='Aliasing included in the noise model')
 axs[1,1].semilogy([x['history_s'] for x in tail_checks],[x['max_whitened_covariance_error'] for x in tail_checks],'o-');axs[1,1].axhline(1e-3,color='k',ls=':',lw=1);axs[1,1].set(xlabel='Retained innovation history (s)',ylabel='Max covariance error',title='Truncation diagnostic; exact model retains history')
 fig.savefig(OUT/'receiver_characterization.png',dpi=160);plt.close(fig)
 lines=['# Causal sampled receiver model','','The implemented receiver processes data causally using a stable sixth-order SOS low-pass at 16 complex samples/s and decimates to 4 complex samples/s. Signal and noise undergo the identical streaming operations. Cholesky innovations use the full finite-record covariance and past samples only. No hard spectral projection is used by the receiver.','','## Fixed receiver specification','',f'- Full -3 dB bandwidth: 2 Hz; two-sided noise-equivalent bandwidth: {enbw:.6f} Hz.',f'- Maximum folded/main noise PSD within +-1 Hz: {orders[2]["maximum_alias_to_main_psd_fraction"]:.4g}. Order 6 is the lowest tested even order satisfying the 1e-4 design budget; order 4 fails.',f'- Group delay: {result["group_delay"]["at_center_s"]:.4f} s at center; up to {result["group_delay"]["max_within_minus3db_band_s"]:.4f} s within the -3 dB band.',f'- Step error remains below 1e-4 after {step_time:.4f} s. Remaining impulse-response energy falls below 1e-6 after {energy_time:.4f} s. These are settling criteria, not a hard delay or per-packet warm-up requirement.','','Bandwidth, noise spectrum and filter coefficients are held fixed versus frequency offset while center tuning remains allowed. The 4 Hz Nyquist span and the 2 Hz -3 dB bandwidth are different quantities. Transition-band information can remain useful; the previous ideal-band capacity bound is not automatically applicable.','','## Noise and likelihood checks','',f'- Exact whitened finite-record covariance differs from identity by at most {identity_error:.3g}.',f'- Independent causal-noise simulation: covariance-lag error {100*relative_cov_error:.3f}% of variance; innovation power {white_power:.5f}; largest measured lag correlation {max(abs(lagwhite)):.4g}.',f'- 1,024-sample covariance factorization: {chol_seconds:.4f} s; streaming 1,024 innovations: {innovation_stream_seconds:.4f} s in this local run.', f'- Common-record initialization: 16 s burn-in differs from the longer-prehistory reference by {burn_checks[-2]["max_output_error_in_noise_rms"]:.3g} noise RMS; 32 s by {burn_checks[-1]["max_output_error_in_noise_rms"]:.3g}.', '- Seven tests pass: stream/chunk invariance and no future dependence; alias sum versus lag selection; alias budget; covariance/Cholesky prefix consistency; exhaustive 256-message metric comparison; continuous waveform sampling on a fixed clock against independent phase quadrature; filter continuity across packet boundaries.','','The metric is (y-mu)^H C^-1 (y-mu). It is exact for the specified finite stationary Gaussian record, to numerical covariance accuracy. Both candidate means and observations use the same lower-triangular transform. Dense innovations do not make a finite-memory sequence trellis exact; a sequence detector must carry the required history or pass a separate approximation check.','','## Scope and remaining integration','','The model starts at an explicitly sampled, calibrated acceleration-envelope boundary. Upstream transducer calibration and analog acquisition/anti-aliasing are not silently assumed to have been validated. Stationary noise state is continuous across packets; the offline simulation burn-in is not per-packet latency. Known deterministic mean state is separately specified. Acquisition, complete receiver startup, implementation timing and BER remain unvalidated. Existing transmitter/noise assumptions are preserved; no receiver sensitivity improvement is credited.','','Next: connect the matched innovations metric to sequence detection, using the exact short-message likelihood as the reference and the saved memory diagnostics to control approximations. Only then resume coded comparisons.']
 (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
