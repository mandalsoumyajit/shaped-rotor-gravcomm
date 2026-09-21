"""Reproducible fixed-bore FEA of same-mass tungsten-insert rotors."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,hashlib,time,shutil,tempfile,subprocess,argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from postprocess_carrier import read_input,shape
from run_calculix import parse_dat
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/tungsten_comparison_2026_09_20'
FROZEN=OUT/'frozen_designs.json'
def save(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def build_reference(h):
 from design_sweep import mesh_builder
 from shaped_rotor import nset
 original=json.loads((ROOT/'results/design_sweep/medium_rounded/h0.013_rpm1800/metadata.json').read_text())['spec']
 spec=dict(original);spec['insert_size']=original['insert_size']*7850/18000
 folder=OUT/f'reference_h{h:g}';folder.mkdir(parents=True,exist_ok=True)
 nodes,elements,stats=mesh_builder(spec)(h,folder);fixed=stats.pop('_fixed_nodes')
 stats['tungsten_mass_kg']=stats.pop('steel_mass_kg')*18000/7850
 stats['tungsten_polar_inertia_kg_m2']=stats.pop('steel_polar_inertia_kg_m2')*18000/7850
 assert abs(stats['aluminum_mass_kg']+stats['tungsten_mass_kg']-7.913509491619676)<1e-7
 # Geometry builder reports its original steel-weighted quadrupole; correct that metadata.
 stats['quadrupole_x2_minus_y2_kg_m2'] += stats['tungsten_mass_kg']*(1-7850/18000)*spec['center']**2
 lines=['*HEADING','Same-mass 95 percent tungsten heavy alloy inserts, SI','*NODE,NSET=NALL']
 lines += [f'{n},'+','.join(f'{v:.12g}' for v in p) for n,p in sorted(nodes.items())]
 for name,group in elements.items():
  name='TUNGSTEN' if name=='STEEL' else name
  lines += [f'*ELEMENT,TYPE=C3D10,ELSET={name}']+[f'{e},'+','.join(map(str,c)) for e,c in group]
 lines+=['*ELSET,ELSET=EALL','ALUMINUM,TUNGSTEN']+nset('FIXED',fixed)
 for name,E,nu,rho in [('ALUMINUM',70e9,.33,2700.),('TUNGSTEN',330e9,.28,18000.)]:
  lines += [f'*MATERIAL,NAME={name}','*ELASTIC',f'{E},{nu}','*DENSITY',str(rho),f'*SOLID SECTION,ELSET={name},MATERIAL={name}']
 lines+=['*BOUNDARY','FIXED,1,3']
 (folder/'model.inp').write_text('\n'.join(lines)+'\n')
 save(folder/'metadata.json',dict(spec=spec,cad=stats,insert_material=dict(density_kg_m3=18000,E_pa=330e9,poisson_ratio=.28,poisson_ratio_status='assumed')))
 print(json.dumps(dict(reference_mesh=h,nodes=len(nodes),elements=sum(map(len,elements.values())),cad_mass_kg=stats['aluminum_mass_kg']+stats['tungsten_mass_kg'])),flush=True)

def geometry(h):
 source=OUT/f'reference_h{h:g}/model.inp'
 if not source.exists():build_reference(h)
 nodes,els,mats=read_input(source);ids=np.array(sorted(nodes));xyz=np.array([nodes[int(n)] for n in ids]);lookup={int(n):i for i,n in enumerate(ids)}
 conn=np.array([[lookup[n] for n in c] for c in els.values()]);rho=np.array([2700. if mats[e]=='ALUMINUM' else 18000. for e in els]);X=xyz[conn]
 det=np.linalg.det(np.stack((X[:,1]-X[:,0],X[:,2]-X[:,0],X[:,3]-X[:,0]),axis=2));assert min(det)>0
 # All archived tetrahedra have straight midsides: degree-three integration is exact for N*r.
 for k,(i,j) in enumerate(((0,1),(1,2),(2,0),(0,3),(1,3),(2,3)),4):assert np.max(abs(X[:,k]-(X[:,i]+X[:,j])/2))<1e-10
 force=np.zeros_like(xyz);mass=0.;inertia=0.
 points=[(np.full(4,.25),-.8)]+[(np.array([.5 if j==i else 1/6 for j in range(4)]),.45) for i in range(4)]
 for L,w in points:
  N,_=shape(L);q=np.einsum('j,ijk->ik',N,X);dm=rho*det/6*w;mass+=sum(dm);inertia+=sum(dm*np.sum(q[:,:2]**2,axis=1))
  accel=np.column_stack((q[:,1],-q[:,0],np.zeros(len(q)))) # Euler inertia for +alpha
  for j in range(10):np.add.at(force,conn[:,j],N[j]*dm[:,None]*accel)
 torque=sum(np.cross(xyz,force)[:,2]);assert abs(torque/inertia+1)<1e-10
 return source,ids,xyz,force,dict(mesh_mass_kg=mass,mesh_inertia_kg_m2=inertia,min_jacobian= float(min(det)),unit_acceleration_torque_nm=float(torque),nodes=len(ids),elements=len(els))
def prepare():
 OUT.mkdir(parents=True,exist_ok=True);frozen=json.loads(FROZEN.read_text());jobs=[]
 for h in (.018,.013):
  source,ids,xyz,force,stats=geometry(h);text=source.read_text().split('*STEP')[0];prefix=[];node=False
  # Preserve all material, connectivity and CAD-selected bore constraints.
  for line in text.splitlines():
   if line.startswith('*'):node=line.startswith('*NODE,')
   if not node or line.startswith('*'):prefix.append(line)
   else:prefix.append(None)
  def deck(scale):
   it=iter([f'{n},'+','.join(f'{v:.12g}' for v in point*scale) for n,point in zip(ids,xyz)])
   return '\n'.join(next(it) if line is None else line for line in prefix)+'\n'
  save(OUT/f'mesh_h{h:g}.json',dict(stats,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
  np.savez_compressed(OUT/f'geometry_h{h:g}.npz',ids=ids,xyz=xyz,angular_force=force)
  specs=[dict(name=f'basis_{axis}',scale=1.,kind='basis',axis=axis) for axis in ('alpha','gx','gy','gz')]
  spacing=.08/(.4*(224/162))
  for i,ch in enumerate(frozen['channels']):
   for label,fc in [('low',ch['config']['carrier_hz']-3*spacing),('center',ch['config']['carrier_hz']),('high',ch['config']['carrier_hz']+3*spacing)]:specs.append(dict(name=f'design{i+1}_{label}',scale=ch['source_designs'][0]['scale'],kind='centrifugal',rpm=30*fc,design=i+1))
  for spec in specs:
   folder=OUT/f'h{h:g}'/spec['name'];folder.mkdir(parents=True,exist_ok=True);lines=[deck(spec['scale']),'*STEP','*STATIC']
   if spec['kind']=='centrifugal':lines+=['*DLOAD',f"EALL,CENTRIF,{(spec['rpm']*np.pi/30)**2:.12g},0,0,0,0,0,1"]
   elif spec['axis']=='alpha':
    lines+=['*CLOAD'];lines += [f'{n},{j+1},{v:.12g}' for n,f in zip(ids,force) for j,v in enumerate(f) if abs(v)>1e-18]
   else:
    axis=['gx','gy','gz'].index(spec['axis']);direction=[int(i==axis) for i in range(3)];lines+=['*DLOAD','EALL,GRAV,1,'+','.join(map(str,direction))]
   lines+=['*NODE PRINT,NSET=NALL','U','*NODE PRINT,NSET=FIXED','RF','*EL PRINT,ELSET=EALL','S','*NODE FILE','U','*END STEP']
   if spec['kind']=='centrifugal':lines+=['*STEP,PERTURBATION','*FREQUENCY','12','*NODE FILE','U','*END STEP']
   inp='\n'.join(lines)+'\n';(folder/'model.inp').write_text(inp);job=dict(spec,mesh=h,folder=str(folder),input_sha256=hashlib.sha256(inp.encode()).hexdigest());jobs.append(job)
 manifest=dict(frozen=frozen,jobs=jobs,source_sha256={str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'run_calculix.py',ROOT/'postprocess_carrier.py',ROOT/'design_sweep.py',ROOT/'shaped_rotor.py',FROZEN]},assumptions='linear elastic, fixed bore, bonded inserts; gravity/tangential static bases omit geometric stiffness; modal step includes centrifugal prestress only')
 save(OUT/'manifest.json',manifest)
 return jobs
def solve(job):
 folder=Path(job['folder']);assert hashlib.sha256((folder/'model.inp').read_bytes()).hexdigest()==job['input_sha256']
 if (folder/'summary.json').exists():return json.loads((folder/'summary.json').read_text())
 t=time.time()
 with tempfile.TemporaryDirectory(prefix='gravcomm-scaled-') as d:
  work=Path(d);shutil.copyfile(folder/'model.inp',work/'model.inp')
  try:
   p=subprocess.run(['ccx','-i','model'],cwd=work,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=1800,env={**os.environ,'OMP_NUM_THREADS':'2'})
   (work/'solver.log').write_text(p.stdout)
  finally:
   for path in work.iterdir():
    if path.is_file():shutil.copyfile(path,folder/path.name)
 assert p.returncode==0 and '*ERROR' not in p.stdout and 'Job finished' in p.stdout,folder
 blocks,freq=parse_dat(folder/'model.dat');u=np.array(blocks['displacements']);stress=np.array(blocks['stresses']);rf=np.array(blocks['forces']);assert len(u)>0 and len(stress)>0
 assert job['kind']!='centrifugal' or (len(freq)==12 and min(freq)>0)
 np.savez_compressed(folder/'fields.npz',u=u,stress=stress,rf=rf)
 s=stress[:,2:8];vm=np.sqrt(((s[:,0]-s[:,1])**2+(s[:,1]-s[:,2])**2+(s[:,2]-s[:,0])**2)/2+3*np.sum(s[:,3:]**2,axis=1))
 row=dict(job,frequencies_hz=freq,max_u_m=float(max(np.linalg.norm(u[:,1:4],axis=1))),raw_max_vm_pa=float(max(vm)),elapsed_s=time.time()-t);save(folder/'summary.json',row);return row
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--workers',type=int,default=3);a=p.parse_args()
 jobs=prepare() if a.prepare else json.loads((OUT/'manifest.json').read_text())['jobs'];completed=[]
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  fs={pool.submit(solve,j):j for j in jobs}
  for f in as_completed(fs):
   row=f.result();completed.append(row);save(OUT/'progress.json',dict(completed=len(completed),total=len(jobs),workers_running=len(completed)<len(jobs),latest=row));print(json.dumps(dict(done=len(completed),total=len(jobs),name=row['name'],mesh=row['mesh'],u=row['max_u_m'],frequencies=row['frequencies_hz'][:3])),flush=True)
 save(OUT/'solver_results.json',completed)
