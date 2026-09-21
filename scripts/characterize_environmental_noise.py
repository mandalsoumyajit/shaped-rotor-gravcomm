"""Reproducible environmental scenarios; no optimized link or BER claim."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.receivers import StructuralOscillator
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/environmental_background_v1';OUT.mkdir(exist_ok=True)
f=np.geomspace(.1,100,2001);floor=float(StructuralOscillator().acceleration_noise_asd(50.3));rows=[];scenarios=[]
for seismic in ('NLNM','NHNM'):
 for resonance in (.15,.3,.6):
  model=SeismicBackground(seismic_model=seismic,isolation_frequency_hz=resonance)
  c=model.components(f);background=np.sqrt(c['background_upper_psd'])
  # Last crossing: all higher grid frequencies also meet this center-floor comparison.
  bad=np.flatnonzero(background>floor)
  index=bad[-1]+1 if len(bad) else 0
  crossover=float(f[index]) if index<len(f) else None
  scenarios.append(dict(seismic=seismic,isolation_frequency_hz=resonance,background_below_reference_center_floor_above_hz=crossover))
  if resonance==.3:
   for frequency in (.1,.2,1.,2.,5.,10.,20.,50.3):
    v=model.components(frequency)
    rows.append(dict(seismic=seismic,frequency_hz=frequency,newtonian_nGal_rtHz=float(v['newtonian_acceleration_asd']/1e-11),support_nGal_rtHz=float(v['support_acceleration_asd']/1e-11),total_with_instrument_center_nGal_rtHz=float(np.sqrt(floor**2+v['background_upper_psd'])/1e-11),extrapolated=bool(v['extrapolated'])))
fig,axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
for ax,seismic in zip(axes,('NLNM','NHNM')):
 c=SeismicBackground(seismic_model=seismic).components(f)
 ax.loglog(f,c['newtonian_acceleration_asd']/1e-11,label='Newtonian gravity')
 ax.loglog(f,c['support_acceleration_asd']/1e-11,label='Residual support motion')
 ax.loglog(f,np.sqrt(floor**2+c['background_upper_psd'])/1e-11,label='Total at tuned center')
 ax.axhline(floor/1e-11,color='k',ls=':',label='Fixed instrumental floor')
 ax.axvspan(10,100,alpha=.1,color='gray',label='Assumed extension >10 Hz')
 ax.set(xlabel='Absolute frequency (Hz)',ylabel='Acceleration ASD (nGal / sqrt(Hz))',title=seismic+' seismic scenario',ylim=(1e-7,1e9));ax.legend(fontsize=8)
fig.savefig(OUT/'background_spectra.png',dpi=160);plt.close(fig)
result=dict(reference_instrument_center_asd=floor,nominal=dict(seismic_model='NHNM',density_kg_m3=2000,rayleigh_gamma=.8,height_m=0,isolation_stages=2,isolation_frequency_hz=.3,isolation_damping=.05,above_10_hz='constant ground acceleration ASD at 10 Hz; assumed, not Peterson data'),frequency_hz=f.tolist(),table=rows,sensitivity=scenarios,sources=['https://doi.org/10.3133/ofr93322','https://link.springer.com/article/10.1007/lrr-2015-3','https://search.r-project.org/CRAN/refmans/IRISSeismic/html/noiseModels.html'],source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('src/gravcomm/environmental_noise.py','src/gravcomm/causal_receiver.py','src/gravcomm/noise.py','src/gravcomm/receivers.py','src/gravcomm/constants.py','scripts/characterize_environmental_noise.py','tests/test_environmental_noise.py')})
(OUT/'characterization.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# Assumed environmental background for carrier selection','','Adopt NHNM plus an explicit isolation model as a provisional seismic-background baseline; NLNM is the quiet sensitivity case. This is a reproducible design scenario, not a measured laboratory spectrum or a bound on all disturbances. Instrumental sensitivity and 2 Hz receive-filter bandwidth remain fixed; background noise is evaluated at absolute physical frequency.','','## Model and provenance','','Peterson ground-acceleration PSD: 10^((A+B log10(1/f))/10) in (m/s^2)^2/Hz, using the existing tabulated coefficients in src/gravcomm/noise.py. Model data end at 10 Hz. From 10 to 100 Hz we explicitly hold ground-acceleration ASD at its 10 Hz value; this is an assumed engineering continuation, not published Peterson data or a universal upper bound. The scenario domain is 0.1-100 Hz. Sources: https://doi.org/10.3133/ofr93322 and https://search.r-project.org/CRAN/refmans/IRISSeismic/html/noiseModels.html.','','Convert vertical ground acceleration to displacement by dividing ASD by (2*pi*f)^2. For a single horizontal test mass in an isotropic Rayleigh field, gravitational ASD is (2*pi*G*rho*gamma/sqrt(2))*exp(-2*pi*f*h/c_R)*displacement_ASD, following Harms, Eq. (109) of the 2015 review: https://link.springer.com/article/10.1007/lrr-2015-3. We assume rho=2000 kg/m^3, gamma=0.8 and h=0; these are declared ground/site assumptions, not fitted measurements. Rayleigh speed 200 m/s is immaterial at h=0. This does not treat all seismic wave types as measured Rayleigh waves.','','For support motion, use an explicit hypothetical two-stage horizontal isolation chain: each stage has resonance 0.3 Hz and damping ratio 0.05. Its base-to-mass transfer is (1+2*i*zeta*r)/(1-r^2+2*i*zeta*r), r=f/f_iso; multiply the two transfers. This cascade is an idealized feed-forward stage model, not a validated coupled suspension. Assume horizontal ground ASD equals the vertical seismic envelope for this scenario. Tilt, actuator/control noise and the coupled six-degree-of-freedom isolation dynamics require a later hardware budget. These assumptions must accompany any numerical operating-point result.','','Because support and Newtonian disturbances can share a seismic origin, use the conservative coherent ASD sum before squaring: S_background=(N_Newtonian+N_support)^2. Add independent instrumental PSD afterward. This is a worst-phase envelope, not a claim of statistically independent seismic paths.','','## Calculated levels at the tuned carrier center','','| Seismic case | Frequency (Hz) | Newtonian ASD (nGal/sqrt Hz) | Support ASD | Total with instrumental floor |','|---|---:|---:|---:|---:|']
for row in rows:lines.append(f"| {row['seismic']} | {row['frequency_hz']:g} | {row['newtonian_nGal_rtHz']:.4g} | {row['support_nGal_rtHz']:.4g} | {row['total_with_instrument_center_nGal_rtHz']:.4g} |")
lines+=['','All entries above 10 Hz use the stated extrapolation. The fixed instrumental center floor is '+f'{floor/1e-11:.4f} nGal/sqrt(Hz).','', '## Implications and next calculation','','The Rayleigh gravitational component alone is small relative to the current instrumental floor at frequencies of a few hertz and above in these scenarios. Residual support motion is the dominant frequency-dependent penalty of the assumed terrestrial receiver. Do not assert that Newtonian noise necessarily sets an optimum in the tens-of-hertz band. A frequency optimum must follow the actual noise budget, rather than an imposed low-frequency penalty.','','Varying isolation resonance to 0.15, 0.3 and 0.6 Hz is a hardware sensitivity study. The crossings below compare background with the instrumental center floor only; they are not operating-frequency optima or BER thresholds:','']
for row in scenarios:lines.append(f"- {row['seismic']}, isolation resonance {row['isolation_frequency_hz']:g} Hz: background remains below center floor above {row['background_below_reference_center_floor_above_hz']} Hz on the evaluated grid.")
lines+=['','The optional CausalReceiver(background=SeismicBackground(), carrier_hz=fc) interface now uses S_total(nu;fc)=S_instrument_reference(50.3+nu)+S_background(fc+nu) through the same causal filter and alias sum, preserving the fixed instrumental profile. It validates the entire acquisition-frequency interval and rejects out-of-domain carriers. Four environmental/integration tests pass. Evaluate full waveform likelihoods and power at candidate carriers; do not use a center-only SNR to choose a carrier. Equal-amplitude simulation reuse is now conditional on carrier and environmental scenario.','','The present 16 complex samples/s envelope interface spans offsets +-8 Hz. At fc<=8.1 Hz it reaches below this scenario domain, and for fc<=8 Hz crosses zero physical frequency. Do not substitute abs(f) or clip silently: low-carrier packet simulations require a consistent real-passband/acquisition model or a separately validated narrower sampled-envelope interface. A 2 Hz modulation band also requires a positive lower band edge. The scalar spectrum calculation here is valid independently of that unresolved interface.','','Atmospheric gravity, nearby moving objects, coherent lines and technical/control disturbances are not quantified by this seismic scenario. It is a defensible assumed baseline, not complete environmental qualification. Add explicit spectra/coherent disturbance cases when data are available; do not credit chopping with rejecting genuine gravity noise. No code ranking, chosen carrier or post-decoding BER claim follows from these levels.']
(OUT/'MODEL.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(table=rows,sensitivity=scenarios),indent=2))
