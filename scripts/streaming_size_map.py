import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys,json,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from scipy.optimize import brentq
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from gravcomm import FiniteRotor
from range_revision import MASS,INERTIA
from streaming_link_search import A0,save
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/streaming_link_search_v1';rotor=FiniteRotor.from_npz(ROOT/'fea/results/shaped_waveforms/medium_rounded.npz')
def one(spec):
 distance,multiplier=spec;target=A0*multiplier;maximum=(distance-.0627)/.25
 def field(s,n=64):return s*rotor.harmonic_amplitude(distance/s,samples=n)
 lo=.1;hi=min(6.,maximum)
 if field(hi)<target:return dict(distance_m=distance,amplitude_scale=multiplier,excluded='no field bracket within geometry search')
 s=brentq(lambda x:field(x)-target,lo,hi,xtol=1e-8);refined=field(s,128)
 assert abs(refined/target-1)<1e-6
 return dict(distance_m=distance,amplitude_scale=multiplier,amplitude_m_s2=refined,scale=s,diameter_m=.5*s,mass_kg=MASS*s**3,inertia_kg_m2=INERTIA*s**5,radial_gap_m=distance-.25*s-.0127,field_refinement_relative=abs(refined/target-1))
if __name__=='__main__':
 rows=[]
 with ThreadPoolExecutor(max_workers=3) as pool:
  futures=[pool.submit(one,(d,a)) for d in (.5,1.,2.) for a in (.5,.75,1.,1.25,2.,4.)]
  for f in as_completed(futures):
   row=f.result();rows.append(row);save(OUT/'amplitude_size_map.json',sorted(rows,key=lambda x:(x['distance_m'],x['amplitude_scale'])));print(row,flush=True)
