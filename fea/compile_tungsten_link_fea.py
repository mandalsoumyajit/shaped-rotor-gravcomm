"""Audit linear load bases, mesh convergence and scaled rotor response."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import json,hashlib
from pathlib import Path
import numpy as np
from tungsten_link_fea import ROOT,OUT,save
from postprocess_carrier import read_input,read_displacements
from scipy.signal import lsim,TransferFunction

def vm(s):return np.sqrt(((s[...,0]-s[...,1])**2+(s[...,1]-s[...,2])**2+(s[...,2]-s[...,0])**2)/2+3*np.sum(s[...,3:]**2,axis=-1))
def fields(folder):
 z=np.load(folder/'fields.npz');u=z['u'];s=z['stress'];u=u[np.argsort(u[:,0])];s=s[np.lexsort((s[:,1],s[:,0]))];return u,s,z['rf']
def main(meshes=(.018,.013)):
 suffix='' if len(meshes)>1 else f'_h{meshes[0]:g}'
 manifest=json.loads((OUT/'manifest.json').read_text())
 for rel,digest in manifest['source_sha256'].items():assert hashlib.sha256((ROOT.parent/rel).read_bytes()).hexdigest()==digest,rel
 results=[];checks=[];modal=[]
 for h in meshes:
  prefix=OUT/f'h{h:g}';geo=np.load(OUT/f'geometry_h{h:g}.npz');ids=geo['ids'];xyz=geo['xyz'];stats=json.loads((OUT/f'mesh_h{h:g}.json').read_text())
  src=OUT/f'reference_h{h:g}/model.inp';nodes,elements,materials=read_input(src);eid=np.array(sorted(elements));conn=np.array([elements[int(e)] for e in eid]);lookup={int(n):i for i,n in enumerate(ids)};ix=np.array([[lookup[int(n)] for n in c] for c in conn]);X=xyz[ix];vol=np.linalg.det(np.stack((X[:,1]-X[:,0],X[:,2]-X[:,0],X[:,3]-X[:,0]),axis=2))/6
  bary=np.full((4,4),.138196601125011);np.fill_diagonal(bary,1-3*.138196601125011);points=np.einsum('ij,ejk->eik',bary,X[:,:4]).reshape(-1,3);weights=np.repeat(vol/4,4);alu=np.repeat([materials[int(e)]=='ALUMINUM' for e in eid],4)
  bases={axis:fields(prefix/f'basis_{axis}') for axis in ('alpha','gx','gy','gz')}
  for axis,(u,stress,rf) in bases.items():
   assert np.array_equal(u[:,0],ids);assert np.array_equal(stress[:,0],np.repeat(eid,4))
   resultant=np.sum(rf[:,1:4],axis=0);positions=np.array([nodes[int(n)] for n in rf[:,0]]);moment=np.sum(np.cross(positions,rf[:,1:4]),axis=0)
   expected=np.zeros(3)
   if axis!='alpha':expected[['gx','gy','gz'].index(axis)]=-stats['mesh_mass_kg']
   err=np.linalg.norm(resultant-expected)/max(stats['mesh_mass_kg'],1)
   terr=abs(moment[2]/stats['mesh_inertia_kg_m2']-1) if axis=='alpha' else None
   assert err<.002 and (terr is None or terr<.002),(axis,err,terr)
   checks.append(dict(mesh=h,basis=axis,force_relative_error=err,torque_relative_error=terr))
  # Identical dimensionless physical windows across scaled meshes.
  width=.03;hub=np.sqrt(.03**2-(width/2)**2);boss=.2-np.sqrt(.05**2-(width/2)**2);rim=np.sqrt(.23**2-(width/2)**2)
  centers=[(sx*a,sy*b) for a,b in [(hub,.015),(boss,.015),(.015,hub),(.015,rim)] for sx in (-1,1) for sy in (-1,1)]
  windows=[alu & (np.sum((points[:,:2]-c)**2,axis=1)<=r*r) for r in (.005,.01) for c in centers]
  assert all(np.any(mask) for mask in windows)
  for i,ch in enumerate(manifest['frozen']['channels'],1):
   source=ch['source_designs'][0];scale=source['scale'];folder=prefix/f'design{i}_high';u,stress,rf=fields(folder);center=fields(prefix/f'design{i}_center');ss=stress[:,2:8];uu=u[:,1:4]
   fc=ch['config']['carrier_hz'];omega_max=np.pi*(fc+3*.08/(.4*(224/162)));alpha=source['peak_torque_nm']/source['inertia_kg_m2']
   # Triangle-inequality envelopes include speed/acceleration maxima independently: deliberately conservative.
   astress=bases['alpha'][1][:,2:8]*scale**2*alpha;au=bases['alpha'][0][:,1:4]*scale**3*alpha
   gs=[bases[a][1][:,2:8]*scale*9.80665 for a in ('gx','gy','gz')];gu=[bases[a][0][:,1:4]*scale**2*9.80665 for a in ('gx','gy','gz')]
   gvbound=np.sqrt(sum(vm(g)**2 for g in gs));gubound=np.sqrt(sum(np.sum(g*g,axis=1) for g in gu))
   bound=vm(ss)+vm(astress)+gvbound;ubound=np.linalg.norm(uu,axis=1)+np.linalg.norm(au,axis=1)+gubound
   # Actual phase combinations: full extreme-to-extreme hop in both directions, independent axial and sampled in-plane attitudes.
   phase_peaks=[];transition=.6*224/162;spacing=.08/(.4*(224/162))
   for sign in (-1,1):
    for t in np.linspace(0,1,21):
     smooth=10*t**3-15*t**4+6*t**5;f=fc+sign*3*spacing*(2*smooth-1);a=sign*np.pi*6*spacing*(30*t*t*(1-t)**2)/transition
     coeff=(np.pi*f/omega_max)**2;S=ss*coeff+bases['alpha'][1][:,2:8]*scale**2*a;U=uu*coeff+bases['alpha'][0][:,1:4]*scale**3*a
     # Norm bounds over gravity attitude at this physical speed/acceleration pair.
     phase_peaks.append((float(max(vm(S)+gvbound)),float(max(np.linalg.norm(U,axis=1)+gubound))))
   row=dict(design=i,mesh=h,mass_kg=source['mass_kg'],scale=scale,mesh_mass_kg=stats['mesh_mass_kg']*scale**3,mesh_inertia_kg_m2=stats['mesh_inertia_kg_m2']*scale**5,centrifugal_raw_alu_vm_mpa=float(max(vm(ss)[alu])/1e6),centrifugal_raw_tungsten_vm_mpa=float(max(vm(ss)[~alu])/1e6),centrifugal_max_u_um=float(max(np.linalg.norm(uu,axis=1))*1e6),combined_conservative_alu_vm_mpa=float(max(bound[alu])/1e6),combined_conservative_tungsten_vm_mpa=float(max(bound[~alu])/1e6),combined_conservative_u_um=float(max(ubound)*1e6),sampled_transition_gravity_bound_vm_mpa=max(v[0] for v in phase_peaks)/1e6,sampled_transition_gravity_bound_u_um=max(v[1] for v in phase_peaks)*1e6,angular_only_max_u_um=float(max(np.linalg.norm(au,axis=1))*1e6),gravity_any_attitude_u_bound_um=float(max(np.sqrt(sum(np.sum(g*g,axis=1) for g in gu))))*1e6,regional_10mm_scaled_vm_mpa=float(max(np.average(bound[mask],weights=weights[mask]) for mask in windows[16:])/1e6))
   horizontal_stress=vm(ss)+vm(astress)+np.sqrt(vm(gs[0])**2+vm(gs[1])**2)
   horizontal_u=np.linalg.norm(uu,axis=1)+np.linalg.norm(au,axis=1)+np.sqrt(np.sum(gu[0]**2+gu[1]**2,axis=1))
   row.update(horizontal_conservative_alu_vm_mpa=float(max(horizontal_stress[alu])/1e6),horizontal_conservative_u_um=float(max(horizontal_u)*1e6),angular_only_alu_vm_mpa=float(max(vm(astress)[alu])/1e6),horizontal_gravity_alu_vm_bound_mpa=float(max(np.sqrt(vm(gs[0])**2+vm(gs[1])**2)[alu])/1e6))
   results.append(row)
   for speed in ('low','center','high'):
    summary=json.loads((prefix/f'design{i}_{speed}'/'summary.json').read_text());freq=np.array(summary['frequencies_hz']);modal.append(dict(design=i,mesh=h,speed=speed,rpm=summary['rpm'],frequencies_hz=freq.tolist()))
   print(json.dumps(row),flush=True)
 save(OUT/f'static_compiled{suffix}.json',results);save(OUT/f'equilibrium_checks{suffix}.json',checks);save(OUT/f'modal_frequencies{suffix}.json',modal)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--mesh',type=float);a=p.parse_args();main((a.mesh,) if a.mesh else (.018,.013))
