"""Consistent-mass modal coupling and bounded diagnostic forced responses."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,math
import numpy as np
from scipy.signal import lsim,TransferFunction
from scipy.optimize import linear_sum_assignment
from tungsten_link_fea import ROOT,OUT,save
from postprocess_carrier import read_input,read_displacements,EDGES

def mass_shape():
 terms=[]
 for i in range(4):
  a=np.zeros(4,int);a[i]=2;b=np.zeros(4,int);b[i]=1;terms.append([(2,a),(-1,b)])
 for i,j in EDGES:
  a=np.zeros(4,int);a[i]=a[j]=1;terms.append([(4,a)])
 return np.array([[sum(c*d*6*math.prod(math.factorial(int(x)) for x in a+b)/math.factorial(int(sum(a+b))+3) for c,a in row for d,b in col) for col in terms] for row in terms])
def main(meshes=(.018,.013)):
 suffix='' if len(meshes)>1 else f'_h{meshes[0]:g}'
 frozen=json.loads((OUT/'manifest.json').read_text())['frozen'];out=[];reference={};C=mass_shape();assert abs(C.sum()-1)<1e-14
 for h in meshes:
  geo=np.load(OUT/f'geometry_h{h:g}.npz');ids=geo['ids'];xyz=geo['xyz'];lookup={int(n):i for i,n in enumerate(ids)}
  _,els,mats=read_input(OUT/f'reference_h{h:g}/model.inp');ix=np.array([[lookup[n] for n in c] for c in els.values()]);X=xyz[ix];vol=np.linalg.det(np.stack((X[:,1]-X[:,0],X[:,2]-X[:,0],X[:,3]-X[:,0]),axis=2))/6;dm=vol*np.array([2700 if mats[e]=='ALUMINUM' else 18000 for e in els]);I=json.loads((OUT/f'mesh_h{h:g}.json').read_text())['mesh_inertia_kg_m2']
  for i,ch in enumerate(frozen['channels'],1):
   s=ch['source_designs'][0]['scale']
   for speed in (('center',) if h==.018 else ('center','low','high')):
    folder=OUT/f'h{h:g}'/f'design{i}_{speed}';summary=json.loads((folder/'summary.json').read_text());modes=read_displacements(folder/'model.frd');keys=sorted(k for k in modes if k>0);assert len(keys)==12
    phi=np.stack([np.array([modes[k][int(n)] for n in ids]) for k in keys],axis=2);del modes
    V=phi[ix];W=np.einsum('ij,ejkn->eikn',C,V,optimize=True)*(dm*s**3)[:,None,None,None];del V
    Mphi=np.zeros_like(phi)
    for j in range(10):np.add.at(Mphi,ix[:,j],W[:,j])
    del W
    gram=np.einsum('nik,nil->kl',phi,Mphi);mass=np.diag(gram);assert max(abs(mass-1))<.01,mass.tolist()
    qa=np.einsum('nik,ni->k',phi,geo['angular_force']*s**4);fractions=qa**2/(mass*I*s**5);assert sum(fractions)<1.01
    components=np.einsum('nik,nik->ik',phi,Mphi)/mass[None,:]
    translation_fraction=np.sum(Mphi,axis=0)**2/(mass[None,:]*sum(dm)*s**3)
    insert_means={}
    for sign in (-1,1):
     selected=np.array([mats[e]=='TUNGSTEN' for e in els]) & (np.mean(X[:,:,0],axis=1)*sign>0)
     meanN=np.r_[np.full(4,-.05),np.full(6,.2)]
     insert_means[str(sign)]=(np.einsum('j,ejkn,e->kn',meanN,phi[ix[selected]],dm[selected],optimize=True)/sum(dm[selected])).tolist()
    key=(h,i)
    if speed=='center':reference[key]=(phi.copy(),Mphi.copy(),mass.copy())
    ref,refM,refmass=reference[key];cross=np.einsum('nik,nil->kl',ref,Mphi);mac=cross**2/(refmass[:,None]*mass[None,:]);ri,ci=linear_sum_assignment(-mac)
    row=dict(translation_effective_mass_fraction=translation_fraction.tolist(),component_mass_fractions=components.tolist(),insert_mean_displacements=insert_means,mesh=h,design=i,speed=speed,frequencies_hz=summary['frequencies_hz'],modal_mass=mass.tolist(),angular_effective_inertia_fraction=fractions.tolist(),matched_center_modes=[int(x+1) for x in ci],matched_mac=[float(mac[a,b]) for a,b in zip(ri,ci)])
    if speed=='center':
     T=.6*224/162;dt=.0005;t=np.arange(0,2*T+dt/2,dt);u=np.minimum(t/T,1);pulse=16*u*u*(1-u)**2;diagnostics=[]
     for damping in (.005,.02,.05):
      gains=[]
      for f in summary['frequencies_hz']:
       w=2*np.pi*f;_,y,_=lsim(TransferFunction([w*w],[1,2*damping*w,w*w]),U=pulse,T=t);gains.append(float(max(abs(y))))
      diagnostics.append(dict(damping_ratio=damping,isolated_transition_modal_peak_over_static_peak=gains,gravity_rotation_modal_amplification=[float(1/np.sqrt((1-(ch['config']['carrier_hz']/2/f)**2)**2+(2*damping*ch['config']['carrier_hz']/2/f)**2)) for f in summary['frequencies_hz']]))
     row['forced_response_diagnostics']=diagnostics
    out.append(row);save(OUT/f'modal_response{suffix}.json',out);print(json.dumps(dict(mesh=h,design=i,speed=speed,min_mac=min(row['matched_mac']),dominant_torsional_mode=int(np.argmax(fractions)+1),inertia_fraction=float(max(fractions)))),flush=True)
  reference.clear()
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--mesh',type=float);a=p.parse_args();main((a.mesh,) if a.mesh else (.018,.013))
