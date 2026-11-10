"""Rerun the aligned-carrier magnetic-wall comparison for the NB-LDPC GF(512) joint design (hop-limited CF-FSK, 0.58 / 0.697 A0).

Reuses the released model code unchanged (Channel, AtomicChannel, bpsk_power,
capacity_power, crossing). Only the design dictionaries (rotor radius, rotor mass,
power ceiling) change; all other assumptions are those of the 2026-09-21 studies:
coil copper mass 10% of rotor mass, coil outer radius = rotor radius, wall thickness
t = d - R - 0.0327 m, PCB and atomic receivers, 0.38/1/3.14/10 pT backgrounds,
BPSK benchmark, capacity window fc +/- 1 Hz.

Usage (from the release root):
  PYTHONPATH=src python scripts/magnetic_wall_nb_2026_11_05.py
"""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
import json
import numpy as np
from scipy.optimize import brentq
from magnetic_wall_study import Channel,DESIGNS as OLD_DESIGNS,ROOT
from gravcomm.magnetic_wall import wall_channel
from magnetic_wall_aligned_carriers import bpsk_power,capacity_power,crossing
from magnetic_wall_atomic_receiver import AtomicChannel

TAG='2026_11_05_nb'
OUT_FH=ROOT/f'results/magnetic_wall_finite_hardware_{TAG}'
OUT_AL=ROOT/f'results/magnetic_wall_aligned_carriers_{TAG}'
OUT_AT=ROOT/f'results/magnetic_wall_atomic_receiver_{TAG}'
# Hop-limited CF-FSK designs: rotor radius = 0.25 m * scale (diameter 0.5 m * scale).
NEW_DESIGNS=[
 dict(distance=.5,scale=0.9010437566657981,radius=.25*0.9010437566657981,mass=5.791239,inertia=0.187497,carrier=24.,power=38.9548),
 dict(distance=1.,scale=1.5761135310793883,radius=.25*1.5761135310793883,mass=30.982996,inertia=3.074508,carrier=24.,power=638.768),
 dict(distance=2.,scale=2.854447700523453,radius=.25*2.854447700523453,mass=184.048840,inertia=59.892758,carrier=18.,power=9367.954)]
for d in OLD_DESIGNS:d.setdefault('carrier',18. if d['distance']==2 else 24.)
COPPER=5.8e7;MU0=4*np.pi*1e-7;BACKGROUNDS=(.38,1.,3.14,10.)

def at_distance(design,d):
 x=dict(design);x['distance']=d;return x

def pcb(design,d=None):return Channel(design,design['distance'] if d is None else d)
def atomic(design,floor=35e-15,d=None):return AtomicChannel(design if d is None else at_distance(design,d),floor)

def critical_copper_thickness(design,make,ambient=1e-12):
 """Copper thickness at which the capacity-bound minimum power equals the ceiling,
 holding source, receiver, clearances, carrier and band fixed."""
 fc=design['carrier'];off=design['radius']+.0127+.02
 f=lambda t:capacity_power(make(design,t+off),COPPER,fc,ambient)-np.log(design['power'])
 lo,hi=.005,3.
 if f(lo)>0:return None
 return float(brentq(f,lo,hi,xtol=1e-5))

def geometry(design):
 c=pcb(design)
 return dict(distance_m=design['distance'],rotor_radius_m=design['radius'],rotor_diameter_m=2*design['radius'],rotor_mass_kg=design['mass'],copper_mass_kg=.1*design['mass'],power_ceiling_W=design['power'],carrier_Hz=design['carrier'],wall_thickness_m=c.thickness,tx_coil_half_thickness_m=float(c.tx_half),tx_coil_mean_radius_m=float(c.coil_radius),magnetic_center_distance_m=float(c.magnetic_distance),pcb_rx_radius_m=float(c.rx_radius))

def free_band(design):
 """Released finite-hardware (freely optimized band) crossings; supplementary."""
 c=pcb(design);out=[]
 for ambient,sensor,knee,label in [(30e-15,'pcb',0,'quiet PCB extrapolation'),(1e-12,'pcb',0,'1 pT ambient'),(10e-12,'pcb',0,'10 pT stress'),(30e-15,'pcb',10,'10 Hz ambient knee'),(1e-15,'ideal_flat',0,'optimistic flat 1 fT total')]:
  out.append(dict(label=label,capacity_failure_sigma=c.crossing('capacity',ambient=ambient,sensor=sensor,knee=knee),bpsk_budget_sigma=c.crossing('bpsk',ambient=ambient,sensor=sensor,knee=knee)))
 return out

def aligned(design):
 c=pcb(design);fc=design['carrier'];rows=[]
 for pt in BACKGROUNDS:
  a=pt*1e-12;bf=lambda s:bpsk_power(c,s,fc,a);cf=lambda s:capacity_power(c,s,fc,a)
  bs=crossing(bf,design['power']);cs=crossing(cf,design['power'])
  rows.append(dict(distance_m=design['distance'],wall_thickness_m=c.thickness,carrier_Hz=fc,capacity_band_Hz=[fc-1,fc+1],background_pT_sqrtHz=pt,power_ceiling_W=design['power'],bpsk_boundary_S_m=bs,capacity_boundary_S_m=cs,copper_bpsk_W=float(np.exp(bf(COPPER))),copper_capacity_W=float(np.exp(cf(COPPER))),quadrature_change_dB=dict(bpsk=float((bpsk_power(c,bs,fc,a,n=1024)-bf(bs))*10/np.log(10)),capacity=float((capacity_power(c,cs,fc,a,n=4096)-cf(cs))*10/np.log(10)))))
 return rows

def atomic_rows(design):
 fc=design['carrier'];rows=[]
 for floor,label in [(35e-15,'optimistic compensated SERF'),(5e-12,'Mx benchmark')]:
  c=atomic(design,floor)
  for pt in BACKGROUNDS:
   a=pt*1e-12;cf=lambda s:capacity_power(c,s,fc,a)
   rows.append(dict(distance_m=design['distance'],carrier_Hz=fc,sensor=label,instrument_ASD_T_sqrtHz=floor,background_pT_sqrtHz=pt,bpsk_boundary_S_m=crossing(lambda s:bpsk_power(c,s,fc,a),design['power']),capacity_boundary_S_m=crossing(cf,design['power']),copper_capacity_W=float(np.exp(cf(COPPER)))))
 return rows

def skin(design):
 c=pcb(design);fs=np.array([18.,24.]);delta=1/np.sqrt(np.pi*fs*MU0*COPPER)
 kw=dict(thickness=c.thickness,coil_radius=c.coil_radius,copper_mass=.1*design['mass'],stand_off=c.z1,receive_radius=c.rx_radius,order=256,loss_refinement=2)
 h=wall_channel(fs,COPPER,c.magnetic_distance,**kw)[0];h0=wall_channel(fs,0,c.magnetic_distance,**kw)[0]
 return [dict(f_Hz=float(f),skin_depth_m=float(de),bulk_amplitude=float(np.exp(-c.thickness/de)),finite_loop_ratio=float(abs(x/y))) for f,de,x,y in zip(fs,delta,h,h0)]

def evaluate(designs,label):
 res=dict(label=label,designs=[])
 for d in designs:
  r=dict(geometry=geometry(d),skin_depth=skin(d),aligned_pcb=aligned(d),atomic=atomic_rows(d),
   critical_copper_thickness_m=dict(pcb_1pT=critical_copper_thickness(d,lambda x,dd:pcb(x,dd)),serf_1pT=critical_copper_thickness(d,lambda x,dd:atomic(x,35e-15,dd))))
  res['designs'].append(r);print(label,d['distance'],'done',flush=True)
 return res

def main():
 for o in (OUT_FH,OUT_AL,OUT_AT):o.mkdir(parents=True,exist_ok=True)
 old=evaluate(OLD_DESIGNS,'released 2026-09-21 designs (reproduction)')
 new=evaluate(NEW_DESIGNS,'hop-limited CF-FSK designs 2026-11-05')
 for r,d in zip(new['designs'],NEW_DESIGNS):r['free_band_supplementary']=free_band(d)
 OUT_FH.joinpath('summary.json').write_text(json.dumps(dict(new=new,old_reproduction=old,new_designs=NEW_DESIGNS,old_designs=OLD_DESIGNS),indent=2)+'\n')
 OUT_AL.joinpath('summary.json').write_text(json.dumps([x for r in new['designs'] for x in r['aligned_pcb']],indent=2)+'\n')
 OUT_AT.joinpath('summary.json').write_text(json.dumps([x for r in new['designs'] for x in r['atomic']],indent=2)+'\n')
 L=['# Magnetic-wall comparison for hop-limited CF-FSK designs (2026-11-05)','',
 'Released model code unchanged; only rotor radius (0.25 m x scale), rotor mass and power ceiling (new peak mechanical power) change. Coil copper mass 10% of rotor mass, coil outer radius = rotor radius, t = d - R - 0.0327 m, receivers, noise levels, BPSK benchmark and the fc +/- 1 Hz capacity window are as in the 2026-09-21 studies. The old columns are recomputed here and match the archived values.','',
 '| d (m) | R old/new (m) | t old/new (m) | ceiling old/new (W) | fc (Hz) | PCB 1 pT old/new (1e6 S/m) | SERF 1 pT old/new (1e6 S/m) |','|---:|---:|---:|---:|---:|---:|---:|']
 pick=lambda rows,pt,s=None:[x for x in rows if x['background_pT_sqrtHz']==pt and (s is None or x['sensor']==s)][0]
 for o,n in zip(old['designs'],new['designs']):
  go,gn=o['geometry'],n['geometry']
  L.append(f"| {gn['distance_m']:g} | {go['rotor_radius_m']:.4f}/{gn['rotor_radius_m']:.4f} | {go['wall_thickness_m']:.4f}/{gn['wall_thickness_m']:.4f} | {go['power_ceiling_W']:.2f}/{gn['power_ceiling_W']:.2f} | {gn['carrier_Hz']:g} | {pick(o['aligned_pcb'],1.)['capacity_boundary_S_m']/1e6:.3g}/{pick(n['aligned_pcb'],1.)['capacity_boundary_S_m']/1e6:.3g} | {pick(o['atomic'],1.,'optimistic compensated SERF')['capacity_boundary_S_m']/1e6:.3g}/{pick(n['atomic'],1.,'optimistic compensated SERF')['capacity_boundary_S_m']/1e6:.3g} |")
 L+=['','## New designs: aligned-carrier PCB receiver','','| d (m) | fc | background (pT) | BPSK crossing (S/m) | capacity crossing (S/m) | copper BPSK W | copper capacity W |','|---:|---:|---:|---:|---:|---:|---:|']
 for n in new['designs']:
  for x in n['aligned_pcb']:L.append(f"| {x['distance_m']:g} | {x['carrier_Hz']:g} | {x['background_pT_sqrtHz']:g} | {x['bpsk_boundary_S_m']:.4g} | {x['capacity_boundary_S_m']:.4g} | {x['copper_bpsk_W']:.3g} | {x['copper_capacity_W']:.3g} |")
 L+=['','## New designs: atomic receiver','','| d (m) | sensor | background (pT) | BPSK (S/m) | capacity (S/m) | copper capacity W |','|---:|---|---:|---:|---:|---:|']
 for n in new['designs']:
  for x in n['atomic']:L.append(f"| {x['distance_m']:g} | {x['sensor']} | {x['background_pT_sqrtHz']:g} | {x['bpsk_boundary_S_m']:.4g} | {x['capacity_boundary_S_m']:.4g} | {x['copper_capacity_W']:.3g} |")
 L+=['','## Copper (5.8e7 S/m) checks','','| d (m) | t (m) old/new | bulk amp 24 Hz old/new | bulk amp 18 Hz old/new | finite-loop ratio at fc old/new | critical Cu thickness PCB 1 pT old/new (m) | critical Cu thickness SERF 1 pT old/new (m) |','|---:|---:|---:|---:|---:|---:|---:|']
 for o,n in zip(old['designs'],new['designs']):
  fc=n['geometry']['carrier_Hz'];s=lambda r,f:[x for x in r['skin_depth'] if x['f_Hz']==f][0]
  L.append(f"| {n['geometry']['distance_m']:g} | {o['geometry']['wall_thickness_m']:.4f}/{n['geometry']['wall_thickness_m']:.4f} | {s(o,24.)['bulk_amplitude']:.3g}/{s(n,24.)['bulk_amplitude']:.3g} | {s(o,18.)['bulk_amplitude']:.3g}/{s(n,18.)['bulk_amplitude']:.3g} | {s(o,fc)['finite_loop_ratio']:.3g}/{s(n,fc)['finite_loop_ratio']:.3g} | {o['critical_copper_thickness_m']['pcb_1pT']:.4f}/{n['critical_copper_thickness_m']['pcb_1pT']:.4f} | {o['critical_copper_thickness_m']['serf_1pT']:.4f}/{n['critical_copper_thickness_m']['serf_1pT']:.4f} |")
 L+=['','Critical copper thickness: wall thickness at which the capacity-bound minimum source power at 5.8e7 S/m equals the design ceiling, holding source coil, receiver, clearances, carrier and 2 Hz window fixed (new diagnostic; the manuscript quotes the evaluated wall thicknesses).','','Reproduce: PYTHONPATH=src python scripts/magnetic_wall_nb_2026_11_05.py. Raw data: summary.json here and in the aligned_carriers / atomic_receiver folders with the same tag.']
 OUT_FH.joinpath('RESULTS.md').write_text('\n'.join(L)+'\n',encoding='utf-8')
 for o in (OUT_AL,OUT_AT):o.joinpath('RESULTS.md').write_text(f'See ../magnetic_wall_finite_hardware_{TAG}/RESULTS.md. Raw data: summary.json. Reproduce: PYTHONPATH=src python scripts/magnetic_wall_nb_2026_11_05.py\n',encoding='utf-8')
 print('\n'.join(L))
if __name__=='__main__':main()
