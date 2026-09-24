"""Conducting barrier fills the link: finite coil close to each wall face."""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json
import numpy as np
from scipy.optimize import brentq
from scipy.special import erfcinv,logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.magnetic_wall import wall_channel
from gravcomm.magnetic_link import noise_asd,minimum_capacity_logpower
from magnetic_conductivity_study import DESIGNS,RINFO
OUT=ROOT/'results/magnetic_wall_finite_hardware_2026_09_21'
F_EDGES=np.geomspace(1e-4,1000,2049);F=np.sqrt(F_EDGES[:-1]*F_EDGES[1:]);DF=np.diff(F_EDGES)

class Channel:
 def __init__(self,design,d,gap=.01,n=640):
  self.design=design;self.d=d;self.gap=gap
  self.thickness=d-design['radius']-.0127-2*gap
  if self.thickness<=0:raise ValueError('No space for wall after gravitational hardware and clearances')
  # Square copper winding cross-section with 60% packing; add 1 cm housing clearance.
  volume=.1*design['mass']/8960
  self.tx_half=.5*np.sqrt(volume/(2*np.pi*design['radius']*.6))
  for _ in range(12):self.tx_half=.5*np.sqrt(volume/(2*np.pi*(design['radius']-self.tx_half)*.6))
  self.coil_radius=design['radius']-self.tx_half
  self.z1=gap+self.tx_half;self.z2=gap+.005
  self.magnetic_distance=self.z1+self.thickness+self.z2
  self.rx_radius=np.sqrt(36.388/(100*np.pi))
  self.product=np.geomspace(1e-9,1e14,n)
  H,L,E=wall_channel(self.product,1.,self.magnetic_distance,thickness=self.thickness,coil_radius=self.coil_radius,copper_mass=.1*design['mass'],stand_off=self.z1,receive_radius=self.rx_radius)
  _,_,Er=wall_channel(self.product,1.,self.magnetic_distance,thickness=self.thickness,coil_radius=self.rx_radius,copper_mass=1.,stand_off=self.z2)
  self.logEr=np.log(np.maximum(Er/self.product,1e-300))
  self.logH=np.log(np.maximum(abs(H),1e-300));self.logE=np.log(np.maximum(E/self.product,1e-300));self.copper=float(L[0]-E[0])
 def gain(self,f,sigma,ambient=30e-15,sensor='pcb',knee=0):
  x=np.log(np.maximum(sigma*f,self.product[0]));h=np.interp(x,np.log(self.product),self.logH);e=np.exp(np.interp(x,np.log(self.product),self.logE))*f
  wallnoise=4*1.380649e-23*300*np.exp(np.interp(x,np.log(self.product),self.logEr))/(4*np.pi**2*f)
  return 2*h-np.log(self.copper+e)-np.log(noise_asd(f,ambient,sensor,knee)**2+wallnoise)
 def capacity(self,sigma,ambient=30e-15,sensor='pcb',knee=0):return minimum_capacity_logpower(self.gain(F,sigma,ambient,sensor,knee),DF,RINFO)
 def bpsk(self,sigma,ambient=30e-15,sensor='pcb',knee=0):
  rate=256/224;alpha=.25;half=rate*(1+alpha)/2;carriers=np.geomspace(1.,999.,601)
  z=(np.arange(96)+.5)/96*2*half-half;pulse=np.ones_like(z)/rate;edge=abs(z)>(1-alpha)*rate/2;pulse[edge]=.5/rate*(1+np.cos(np.pi/(alpha*rate)*(abs(z[edge])-(1-alpha)*rate/2)));weights=pulse/pulse.sum()
  g=self.gain(carriers[:,None]+z,sigma,ambient,sensor,knee)
  powers=np.log(erfcinv(.002)**2*rate)+logsumexp(np.log(weights)-g,axis=-1);i=int(np.argmin(powers))
  return float(powers[i]),float(carriers[i])
 def crossing(self,method='capacity',**kw):
  def fun(x):
   result=getattr(self,method)(10**x,**kw);return (result[0] if isinstance(result,tuple) else result)-np.log(self.design['power'])
  if fun(-4)>0:return None
  if fun(12)<0:return 1e12
  return float(10**brentq(fun,-4,12,xtol=1e-5))

def main():
 OUT.mkdir(parents=True,exist_ok=True);summary=[]
 for design in DESIGNS:
  c=Channel(design,design['distance']);r=dict(design=design,gap_m=.01,wall_thickness_m=c.thickness,magnetic_center_distance_m=c.magnetic_distance,tx_coil_half_thickness_m=c.tx_half,tx_coil_mean_radius_m=c.coil_radius,rx_radius_m=c.rx_radius,noise_scenarios=[],points=[])
  for ambient,sensor,knee,label in [(30e-15,'pcb',0,'quiet PCB extrapolation'),(1e-12,'pcb',0,'1 pT ambient'),(10e-12,'pcb',0,'10 pT stress'),(30e-15,'pcb',10,'10 Hz ambient knee'),(1e-15,'ideal_flat',0,'optimistic flat 1 fT total')]:
   r['noise_scenarios'].append(dict(label=label,capacity_failure_sigma=c.crossing('capacity',ambient=ambient,sensor=sensor,knee=knee),bpsk_budget_sigma=c.crossing('bpsk',ambient=ambient,sensor=sensor,knee=knee)))
  for sigma in (.01,4.,1e4,1e6,1e7,5.8e7):
   power,f=c.bpsk(sigma);r['points'].append(dict(conductivity=sigma,bpsk_required_average_W=float(np.exp(min(power,700))),carrier_Hz=f,capacity_minimum_average_W=float(np.exp(min(c.capacity(sigma),700)))))
  # Validate actual channel at frequencies relevant to the conductivity crossing.
  sigma=r['noise_scenarios'][0]['capacity_failure_sigma'];freq=np.geomspace(.001,10,40)
  H,L,E=wall_channel(freq,sigma,c.magnetic_distance,thickness=c.thickness,coil_radius=c.coil_radius,copper_mass=.1*design['mass'],stand_off=c.z1,order=256,loss_refinement=2,receive_radius=c.rx_radius)
  _,_,Er=wall_channel(freq,sigma,c.magnetic_distance,thickness=c.thickness,coil_radius=c.rx_radius,copper_mass=1.,stand_off=c.z2,order=256,loss_refinement=2)
  gd=2*np.log(np.maximum(abs(H),1e-300))-np.log(L)-np.log(noise_asd(freq)**2+4*1.380649e-23*300*Er/(2*np.pi*freq)**2)
  mask=gd>max(gd)-40
  r['gain_validation_max_abs_dB']=float(max(abs(gd[mask]-c.gain(freq,sigma)[mask]))/np.log(10)*10)
  for gap in (.005,.02):
   alt=Channel(design,design['distance'],gap=gap);r[f'gap_{gap}_capacity_sigma']=alt.crossing()
  summary.append(r);print(json.dumps(r),flush=True)
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 # Boundaries for each source size and numerical power ceiling, with wall filling separation.
 sigmas=np.geomspace(1e-4,1e10,101);maps=[];fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
 for ax,design in zip(axs,DESIGNS):
  ds=np.geomspace(max(.3,design['radius']+.10),30,41);Z=np.empty((len(sigmas),len(ds)))
  for j,d in enumerate(ds):
   c=Channel(design,d,n=384)
   for i,sigma in enumerate(sigmas):Z[i,j]=(c.capacity(sigma)-np.log(design['power']))/np.log(10)*10
  im=ax.pcolormesh(ds,sigmas,np.clip(Z,-100,100),cmap='RdBu_r',vmin=-100,vmax=100,shading='auto');ax.contour(ds,sigmas,Z,levels=[0],colors='black',linewidths=1.5)
  ax.axvline(design['distance'],color='black',ls=':',lw=.8);ax.axhline(4,color='gray',ls=':',lw=.8)
  ax.set(xscale='log',yscale='log',xlabel='Gravitational axis-to-sensor distance (m)',title=f"{design['power']:.3g} W; coil radius {design['radius']:.3g} m")
  maps.append(dict(distance_m=ds.tolist(),conductivity_S_m=sigmas.tolist(),power_margin_dB=Z.tolist()))
 axs[0].set_ylabel('Conductivity (S/m)');fig.colorbar(im,ax=axs,label='Ideal minimum average power / stated ceiling (dB)',shrink=.9)
 fig.savefig(OUT/'wall_boundary.png',dpi=160);fig.savefig(OUT/'wall_boundary.pdf');plt.close(fig)
 (OUT/'maps.json').write_text(json.dumps(maps));print('COMPLETE',flush=True)
if __name__=='__main__':main()
