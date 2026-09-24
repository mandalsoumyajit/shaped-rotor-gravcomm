import sys,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'src'))
from magnetic_wall_study import Channel,DESIGNS,F,DF,RINFO
from gravcomm.magnetic_wall import wall_channel
out=[]
for design in DESIGNS:
 c=Channel(design,design['distance']);sigma=5.8e7
 fs=np.array([.1,1.,18.,24.]);delta=1/np.sqrt(np.pi*fs*4*np.pi*1e-7*sigma)
 kw=dict(thickness=c.thickness,coil_radius=c.coil_radius,copper_mass=.1*design['mass'],stand_off=c.z1,receive_radius=c.rx_radius,order=256,loss_refinement=2)
 h=wall_channel(fs,sigma,c.magnetic_distance,**kw)[0];h0=wall_channel(fs,0,c.magnetic_distance,**kw)[0]
 g=c.gain(F,sigma,ambient=1e-12)
 def cap(L):return np.sum(DF*np.maximum(0,L+g))/np.log(2)
 L=brentq(lambda L:cap(L)-RINFO,-1000,1000)
 ps=np.maximum(0,np.exp(L)-np.exp(np.minimum(700,-g)))
 watts=ps*DF;bits=DF*np.maximum(0,L+g)/np.log(2)
 active=ps>0
 quant=lambda v:[float(F[np.searchsorted(np.cumsum(v)/np.sum(v),q)]) for q in (.05,.5,.95)]
 out.append(dict(distance=design['distance'],thickness=c.thickness,skin_checks=[dict(f_Hz=float(f),skin_depth_m=float(de),bulk_amplitude=float(np.exp(-c.thickness/de)),finite_loop_ratio=float(abs(x/y))) for f,de,x,y in zip(fs,delta,h,h0)],capacity_W=float(sum(watts)),capacity_bits_s=float(sum(bits)),active_band_Hz=[float(min(F[active])),float(max(F[active]))],power_percentiles_Hz=quant(watts),information_percentiles_Hz=quant(bits)))
p=(ROOT/'results/magnetic_wall_elf_background_2026_09_21/skin_depth_audit.json');p.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
