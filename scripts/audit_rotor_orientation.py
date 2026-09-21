"""Finite-volume orientation audit, maximizing a real single-axis harmonic projection."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json
from pathlib import Path
import numpy as np
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm import FiniteRotor
from gravcomm.constants import G
from streaming_dynamics_runtime import OUT,save
rotor=FiniteRotor.from_npz(ROOT/'fea/results/shaped_waveforms/medium_rounded.npz')
def one(spec):
 distance,elevation=spec;theta=np.deg2rad(elevation);R=distance*np.array([np.cos(theta),0,np.sin(theta)]);phase=np.arange(64)*2*np.pi/64;p=rotor.positions_m;m=rotor.masses_kg
 h=np.zeros(3,complex)
 for phi in phase:
  c,s=np.cos(phi),np.sin(phi);delta=np.column_stack((c*p[:,0]-s*p[:,1]-R[0],s*p[:,0]+c*p[:,1],p[:,2]-R[2]));w=G*m/np.sum(delta*delta,axis=1)**1.5
  h+=2/len(phase)*np.sum(delta*w[:,None],axis=0)*np.exp(-2j*phi)
 vals,axes=np.linalg.eigh(np.outer(h.real,h.real)+np.outer(h.imag,h.imag))
 return dict(equivalent_distance_m=distance,elevation_deg=elevation,max_single_axis_amplitude=float(np.sqrt(max(vals[-1],0))),radial_amplitude=float(abs(h@ (R/distance))),axis=axes[:,-1].tolist(),phasor_real=h.real.tolist(),phasor_imag=h.imag.tolist())
if __name__=='__main__':
 specs=[(d,a) for d in (.5,.571,.625) for a in (0,5,15,30,45,60,75,90)];rows=[]
 with ThreadPoolExecutor(max_workers=3) as pool:
  for f in as_completed([pool.submit(one,s) for s in specs]):
   x=f.result();rows.append(x);save(OUT/'orientation_audit.json',rows);print(json.dumps(x),flush=True)
