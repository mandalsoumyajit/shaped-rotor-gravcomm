"""Locate static stress maxima by physical region, including fillet studies."""
import json
import numpy as np
from design_sweep import ROOT,OUT
from postprocess_carrier import read_input,quadrature,tensor,write_csv
from run_calculix import parse_dat


def stress_averages(tensors,volumes):
    """Return <VM(sigma)> and VM(<sigma>) with physical volume weights."""
    tensors=np.asarray(tensors);volumes=np.asarray(volumes)
    if len(volumes)==0 or np.any(volumes<=0):raise ValueError('Positive sampled volumes required')
    dev=tensors-np.trace(tensors,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
    mean_vm=np.average(np.sqrt(1.5*np.sum(dev*dev,axis=(1,2))),weights=volumes)
    mean_tensor=np.average(tensors,axis=0,weights=volumes)
    dev_mean=mean_tensor-np.trace(mean_tensor)*np.eye(3)/3
    return float(mean_vm),float(np.sqrt(1.5*np.sum(dev_mean*dev_mean)))


def regions(folder):
    metadata=json.loads((folder/'metadata.json').read_text());spec=metadata['spec']
    c=spec.get('center',.2);insert_r=spec.get('insert_radius_m',.03)
    nodes,elements,materials=read_input(folder/'model.inp')
    blocks,_=parse_dat(folder/'model.dat')
    stresses={(int(v[0]),int(v[1])):v[2:8] for v in blocks['stresses']}
    rows=[];samples=[]
    for eid,ip,N,p,v,rho in quadrature(nodes,elements,materials):
        if materials[eid]!='ALUMINUM':continue
        S=tensor(stresses[eid,ip]);dev=S-np.trace(S)*np.eye(3)/3
        vm=np.sqrt(1.5*np.sum(dev*dev));radius=np.hypot(p[0],p[1])
        pocket=np.hypot(abs(p[0])-c,p[1])
        region='rim_or_spokes'
        if radius<.045:region='hub_junction'
        if pocket<insert_r+.018:region='insert_boss'
        if radius<.01:region='fixed_bore'
        rows.append((region,vm,*p,v))
        samples.append((p,v,vm,S))
    summary=[]
    for region in ('all','hub_junction','insert_boss','fixed_bore','rim_or_spokes'):
        subset=[r for r in rows if region=='all' or r[0]==region]
        peak=max(subset,key=lambda r:r[1])
        summary.append(dict(case=folder.parent.name,mesh=metadata['summary']['mesh_size_m'],region=region,
            vm_max_mpa=peak[1]/1e6,x_mm=peak[2]*1000,y_mm=peak[3]*1000,z_mm=peak[4]*1000))
    # Fixed physical windows, centered on the nominal sharp-geometry junction.
    # This keeps the location identical between sharp and rounded variants.
    width=spec.get('spoke',.025);rim=spec.get('rim',.035);boss=insert_r+.015
    hub_x=np.sqrt(.03**2-(width/2)**2)
    boss_x=c-np.sqrt(boss**2-(width/2)**2)
    rim_y=np.sqrt((.25-rim)**2-(width/2)**2)
    centers={
        'weighted_hub':[(sx*hub_x,sy*width/2) for sx in (-1,1) for sy in (-1,1)],
        'inner_insert_boss':[(sx*boss_x,sy*width/2) for sx in (-1,1) for sy in (-1,1)],
        'unweighted_hub':[(sx*width/2,sy*hub_x) for sx in (-1,1) for sy in (-1,1)],
        'unweighted_rim':[(sx*width/2,sy*rim_y) for sx in (-1,1) for sy in (-1,1)],
        'midarc_rim':[((.25-rim/2)*np.cos(a),(.25-rim/2)*np.sin(a)) for a in np.pi/4+np.arange(4)*np.pi/2]}
    p=np.array([s[0] for s in samples]);v=np.array([s[1] for s in samples]);vm=np.array([s[2] for s in samples]);S=np.array([s[3] for s in samples])
    averages=[]
    for region,locations in centers.items():
        for radius in (.005,.010):
            for index,center in enumerate(locations,1):
                mask=np.sum((p[:,:2]-center)**2,axis=1)<=radius**2
                if not any(mask):raise ValueError('Empty physical averaging window')
                mean_vm,vm_mean_tensor=stress_averages(S[mask],v[mask])
                averages.append(dict(case=folder.parent.name,mesh=metadata['summary']['mesh_size_m'],
                    region=region,window_index=index,window_radius_mm=radius*1000,
                    center_x_mm=center[0]*1000,center_y_mm=center[1]*1000,
                    sampled_volume_mm3=float(sum(v[mask])*1e9),sample_count=int(sum(mask)),
                    mean_vm_mpa=mean_vm/1e6,
                    vm_of_mean_tensor_mpa=vm_mean_tensor/1e6))
    return summary,averages


def main():
    rows=[];averages=[]
    for folder in [ROOT/'results'/'shaped_rotor'/'h0.013_rpm1800']+sorted(OUT.glob('*rounded/h*_rpm1800')):
        if not (folder/'metadata.json').exists():continue
        result,average=regions(folder);rows+=result;averages+=average
        print(json.dumps(result),flush=True)
    write_csv(OUT/'stress_regions.csv',rows)
    write_csv(OUT/'averaged_stress_regions.csv',averages)
    paired=[]
    lookup={(r['case'],r['mesh'],r['region'],r['window_index'],r['window_radius_mm']):r for r in averages}
    for r in averages:
        if r['mesh']!=.013:continue
        coarse=lookup.get((r['case'],.018,r['region'],r['window_index'],r['window_radius_mm']))
        if coarse is None:continue
        paired.append(dict(case=r['case'],region=r['region'],window_index=r['window_index'],
            window_radius_mm=r['window_radius_mm'],fine_mean_vm_mpa=r['mean_vm_mpa'],
            coarse_mean_vm_mpa=coarse['mean_vm_mpa'],
            mean_vm_change_pct=100*(r['mean_vm_mpa']/coarse['mean_vm_mpa']-1),
            sampled_volume_change_pct=100*(r['sampled_volume_mm3']/coarse['sampled_volume_mm3']-1)))
    if paired:write_csv(OUT/'averaged_stress_convergence.csv',paired)


if __name__=='__main__':main()
