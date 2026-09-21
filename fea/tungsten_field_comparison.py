"""Matched finite-volume substitution: shorten insert z, preserve mass weights."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT.parent/'src'))
from gravcomm import FiniteRotor
from tungsten_link_fea import OUT,save
source=ROOT/'results/shaped_waveforms/medium_rounded.npz';z=np.load(source);p=z['positions_m'].copy();m=z['masses_kg'];mask=z['is_steel'].astype(bool)
p[mask,2]*=7850/18000
steel=FiniteRotor.from_npz(source);tungsten=FiniteRotor(p,m,float(z['bounding_radius_m']))
assert abs(np.sum(m*np.sum(p[:,:2]**2,axis=1))/np.sum(m*np.sum(z['positions_m'][:,:2]**2,axis=1))-1)<1e-14
rows=[]
for ch in json.loads((OUT/'frozen_designs.json').read_text())['channels']:
 s=ch['source_designs'][0]['scale'];d=ch['distance_m'][0];a={}
 for label,rotor in [('steel',steel),('tungsten',tungsten)]:
  v=[s*rotor.harmonic_amplitude(d/s,samples=n) for n in (64,128)]
  assert abs(v[0]/v[1]-1)<1e-7
  a[label+'_a2_m_s2']=v[1];a[label+'_phase_relative_error']=abs(v[0]/v[1]-1)
 a.update(distance_m=d,relative_signal_change_percent=100*(a['tungsten_a2_m_s2']/a['steel_a2_m_s2']-1));rows.append(a);print(json.dumps(a),flush=True)
save(OUT/'field_comparison.json',dict(rows=rows,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),method='Exact axial coordinate transform of existing insert quadrature, preserving every mass weight. Same faceted radial geometry and carrier as steel. Fourier phase refinement checked; no new BER qualification.'))
