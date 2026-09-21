"""Third-mesh centrifugal convergence check for the tungsten source."""
import json,hashlib
import numpy as np
from tungsten_link_fea import OUT,geometry,solve,save
if __name__=='__main__':
 h=.010;src,ids,xyz,force,stats=geometry(h)
 ch=json.loads((OUT/'frozen_designs.json').read_text())['channels'][0];s=ch['source_designs'][0]['scale'];rpm=30*(ch['config']['carrier_hz']+3*.08/(.4*(224/162)))
 lines=[];node=False;it=iter([f'{n},'+','.join(f'{v:.12g}' for v in p*s) for n,p in zip(ids,xyz)])
 for line in src.read_text().splitlines():
  if line.startswith('*'):node=line.startswith('*NODE,')
  lines.append(next(it) if node and not line.startswith('*') else line)
 lines+=['*STEP','*STATIC','*DLOAD',f'EALL,CENTRIF,{(rpm*np.pi/30)**2:.12g},0,0,0,0,0,1','*NODE PRINT,NSET=NALL','U','*NODE PRINT,NSET=FIXED','RF','*EL PRINT,ELSET=EALL','S','*NODE FILE','U','*END STEP','*STEP,PERTURBATION','*FREQUENCY','12','*NODE FILE','U','*END STEP']
 folder=OUT/'refinement_h0.01_design1_high';folder.mkdir(exist_ok=True);inp='\n'.join(lines)+'\n';(folder/'model.inp').write_text(inp)
 job=dict(name='design1_high_refinement',scale=s,kind='centrifugal',rpm=rpm,design=1,mesh=h,folder=str(folder),input_sha256=hashlib.sha256(inp.encode()).hexdigest());save(OUT/'refinement_manifest.json',dict(job=job,mesh=stats));row=solve(job);print(json.dumps(row),flush=True)
