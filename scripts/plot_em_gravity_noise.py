from pathlib import Path
import sys,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.macnae_noise import macnae_asd
from gravcomm.environmental_noise import SeismicBackground
OUT=ROOT/'results/noise_comparison_2026_09_21';OUT.mkdir(parents=True,exist_ok=True)
def fit(name,func,lo,hi):
 f=np.geomspace(lo,hi,300);y=func(f);s,b=np.polyfit(np.log10(f),np.log10(y),1)
 residual=np.log10(y)-(b+s*np.log10(f))
 return dict(name=name,band_Hz=[lo,hi],alpha=float(-s),ASD_at_1Hz=float(10**b),rms_residual_dB=float(np.sqrt(np.mean((20*residual)**2))))
bg=SeismicBackground()
fits=[fit('Magnetic',macnae_asd,.1,1),fit('Residual support',lambda f:bg.components(f)['support_acceleration_asd'],1,10),fit('Rayleigh gravity',lambda f:bg.components(f)['newtonian_acceleration_asd'],1,10)]
(OUT/'fits.json').write_text(json.dumps(fits,indent=2)+'\n')
fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
f=np.geomspace(.01,100,1200);axs[0].loglog(f,macnae_asd(f)*1e12,label='Macnae nominal background',color='C0')
f=np.geomspace(.1,100,1200);c=bg.components(f)
for key,label,col in [('support_acceleration_asd','Residual support motion','C1'),('newtonian_acceleration_asd','Rayleigh gravitational noise','C2')]:axs[1].loglog(f,c[key]/1e-11,label=label,color=col)
axs[1].axvspan(10,100,color='gray',alpha=.12,label='Ground ASD continuation')
for j,r in enumerate(fits):
 ax=axs[0 if j==0 else 1];x=np.geomspace(*r['band_Hz'],200);unit=1e12 if j==0 else 1e11
 ax.loglog(x,r['ASD_at_1Hz']*x**(-r['alpha'])*unit,'--',color='black',lw=1,label=fr'Fit: $\alpha={r["alpha"]:.2f}$ ({r["band_Hz"][0]:g}–{r["band_Hz"][1]:g} Hz)')
for ax in axs:
 for fc in (18,24):ax.axvline(fc,color='gray',ls=':',lw=.8)
 ax.set_xlabel('Frequency (Hz)');ax.grid(which='both',alpha=.15);ax.legend(fontsize=7,loc='best')
axs[0].set_ylabel(r'Magnetic ASD (pT/$\sqrt{\mathrm{Hz}}$)');axs[1].set_ylabel(r'Acceleration ASD (nGal/$\sqrt{\mathrm{Hz}}$)')
axs[0].set_title('(a) Environmental magnetic field');axs[1].set_title('(b) Gravitational-receiver environment')
for ext in ('pdf','png'):fig.savefig(OUT/f'noise_comparison.{ext}',dpi=180)
for d in (ROOT/'figures',):
 d.mkdir(exist_ok=True);shutil.copy2(OUT/'noise_comparison.pdf',d/'noise_comparison.pdf')
print(json.dumps(fits,indent=2))
