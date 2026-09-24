"""Conductivity-range screening at gravitational-design power ceilings."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.magnetic_link import *
OUT=ROOT/'results/magnetic_comparison_2026_09_21'
DESIGNS=[dict(distance=.5,mass=7.027334750750038,radius=.2402963533160248,power=120.49236454609739),dict(distance=1.,mass=37.69173275922457,radius=.4206287956499422,power=1980.2436500888432),dict(distance=2.,mass=227.3153385558173,radius=.7656434061793719,power=29680.143)]
RINFO=1-(-.001*np.log2(.001)-.999*np.log2(.999))

def grid(n=2048):
 edges=np.geomspace(.01,1000,n+1);return np.sqrt(edges[:-1]*edges[1:]),np.diff(edges)

def calc(sigma,d,design,ambient=30e-15,sensor='pcb',knee=0,n=2048):
 f,df=grid(n);kw=dict(cavity_radius=design['radius'],coil_fraction=.5,copper_mass=.1*design['mass'],ambient=ambient,sensor=sensor,ambient_knee=knee)
 g=log_gain(f,sigma,d,**kw)
 return minimum_capacity_logpower(g.max(axis=0),df,RINFO)

def crossing(design,d,method='capacity',ambient=30e-15,sensor='pcb',knee=0,powerfactor=1):
 def fun(x):
  if method=='capacity':v=calc(10**x,d,design,ambient,sensor,knee)
  else:v=bpsk_logpower(10**x,d,rate=256/224,cavity_radius=design['radius'],coil_fraction=.5,copper_mass=.1*design['mass'],ambient=ambient,sensor=sensor,ambient_knee=knee)[0]
  return v-np.log(design['power']*powerfactor)
 lo,hi=fun(-4),fun(9)
 if lo>=0:return None
 if hi<=0:return 1e9
 return 10**brentq(fun,-4,9,xtol=1e-5)

def main():
 OUT.mkdir(parents=True,exist_ok=True);results=[]
 for design in DESIGNS:
  d=design['distance'];kw=dict(cavity_radius=design['radius'],coil_fraction=.5,copper_mass=.1*design['mass'])
  row=dict(design=design,noise_scenarios=[])
  for ambient,sensor,knee,label in [(30e-15,'pcb',0,'quiet PCB extrapolation'),(1e-12,'pcb',0,'1 pT background'),(10e-12,'pcb',0,'10 pT stress'),(30e-15,'pcb',10,'quiet with 10 Hz ambient knee'),(1e-15,'ideal_flat',0,'optimistic 1 fT flat total')]:
   rec=dict(label=label,ambient_T_rtHz=ambient,sensor=sensor,knee_Hz=knee,capacity_failure_sigma=crossing(design,d,ambient=ambient,sensor=sensor,knee=knee),bpsk_budget_sigma=crossing(design,d,'bpsk',ambient,sensor,knee))
   row['noise_scenarios'].append(rec)
  row['points']=[]
  for sigma in (0.,.01,4.,1e4,1e6,1e7,5.8e7):
   power,fc,axis=bpsk_logpower(sigma,d,rate=256/224,**kw)
   row['points'].append(dict(sigma=sigma,required_average_W=float(np.exp(min(power,700))),carrier_Hz=fc,orientation=['axial','transverse'][axis]))
  base=row['noise_scenarios'][0]['capacity_failure_sigma'];row['half_power_capacity_sigma']=crossing(design,d,powerfactor=.5)
  row['capacity_grid_refinement_logpower_difference']=calc(base,d,design,n=8192)-calc(base,d,design,n=2048)
  print(json.dumps(row),flush=True);results.append(row)
 (OUT/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
 sigmas=np.geomspace(1e-4,1e8,121);fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained');maps=[]
 for ax,design,row in zip(axs,DESIGNS,results):
  ds=np.geomspace(max(.3,design['radius']*1.1),100,91);Z=np.empty((len(sigmas),len(ds)))
  for i,sigma in enumerate(sigmas):
   for j,d in enumerate(ds):Z[i,j]=(calc(sigma,d,design)-np.log(design['power']))/np.log(10)*10
  maps.append(dict(distance=ds.tolist(),sigma=sigmas.tolist(),required_power_margin_dB=Z.tolist()))
  im=ax.pcolormesh(ds,sigmas,np.clip(Z,-100,100),cmap='RdBu_r',vmin=-100,vmax=100,shading='auto');ax.contour(ds,sigmas,Z,levels=[0],colors='black',linewidths=1.5)
  ax.scatter([design['distance']],[row['noise_scenarios'][0]['capacity_failure_sigma']],c='black',s=22)
  ax.axhline(4,color='gray',ls=':',lw=.8);ax.set(xscale='log',yscale='log',xlabel='Center distance (m)',title=f"{design['power']:.3g} W ceiling; cavity {design['radius']:.3g} m")
 axs[0].set_ylabel('Conductivity (S/m)');fig.colorbar(im,ax=axs,label='Minimum ideal-channel power / power ceiling (dB)',shrink=.9)
 fig.savefig(OUT/'conductivity_range.png',dpi=150);fig.savefig(OUT/'conductivity_range.pdf');plt.close(fig)
 (OUT/'maps.json').write_text(json.dumps(maps))
 # Two frequency diagnostics distinguish propagation attenuation from optimal operating frequency.
 f=np.geomspace(.01,1000,1000);fig,axs=plt.subplots(1,2,figsize=(10,3.7),layout='constrained');design=DESIGNS[2]
 for sigma in (4,1e4,1e6,1e7):
  gain=log_gain(f,sigma,2,cavity_radius=design['radius'],copper_mass=.1*design['mass']).max(axis=0)
  axs[0].semilogx(f,10*gain/np.log(10),label=f'{sigma:g} S/m')
 for ambient,label in [(30e-15,'30 fT ambient'),(1e-12,'1 pT ambient'),(10e-12,'10 pT ambient')]:axs[1].loglog(f,noise_asd(f,ambient)*1e15,label=label)
 axs[0].set(xlabel='Frequency (Hz)',ylabel='Channel gain / noise PSD (dB Hz/W)',ylim=(-50,200));axs[0].legend()
 axs[1].set(xlabel='Frequency (Hz)',ylabel='Total ASD (fT / sqrt(Hz))');axs[1].legend()
 fig.savefig(OUT/'frequency_noise.png',dpi=150);fig.savefig(OUT/'frequency_noise.pdf');plt.close(fig)
 paths=[Path(__file__),ROOT/'src/gravcomm/magnetic_link.py',ROOT/'tests/test_magnetic_link.py']
 (OUT/'manifest.json').write_text(json.dumps(dict(source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},designs=DESIGNS,capacity_information_rate=RINFO,frequency_domain_Hz=[.01,1000],conductivity_map_domain=[1e-4,1e8],power_note='Magnetic average-power relaxation equals gravitational peak mechanical demand; an optimistic magnetic bound, not an equal-average-electrical claim.'),indent=2)+'\n')
 print('COMPLETE',flush=True)
if __name__=='__main__':main()
