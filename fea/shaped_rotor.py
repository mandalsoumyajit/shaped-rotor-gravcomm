"""Fixed-bore, bonded-insert approximation of cad/bench_source_disk.scad.

No shaft/bearing model, fasteners, pocket clearance or frictional contact.
Units m, kg, s, Pa. Rotation about z. Geometry and materials are explicit.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import shutil

import gmsh
import numpy as np
from run_calculix import nset, parse_dat

ROOT = Path(__file__).resolve().parent
SPEC = dict(outer_radius_m=.25, thickness_m=.02, rim_width_m=.035,
            spoke_width_m=.025, hub_radius_m=.03, bore_radius_m=.006,
            insert_radius_m=.03, insert_height_m=.09, insert_center_radius_m=.2,
            boss_radius_m=.045, aluminum_density=2700., steel_density=7850.,
            aluminum_E_pa=70e9, aluminum_nu=.33, steel_E_pa=200e9, steel_nu=.3,
            rpm=1800., insert_interface='bonded conformal; clearance removed',
            support='all bore-surface translations fixed; no shaft or bearings',
            omitted='setscrew, fasteners, contact, fillets, manufacturing tolerances')


def create_mesh(h, folder):
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0)
        gmsh.option.setNumber('General.NumThreads',1)
        gmsh.model.add('shaped_rotor')
        occ = gmsh.model.occ
        def cyl(x,r,z=-.01,height=.02):
            return (3,occ.addCylinder(x,0,z,0,0,height,r))
        outer, inner = cyl(0,.25), cyl(0,.215)
        rim,_ = occ.cut([outer],[inner])
        parts = [cyl(0,.03),cyl(.2,.045),cyl(-.2,.045),
                 (3,occ.addBox(-.216,-.0125,-.01,.432,.025,.02)),
                 (3,occ.addBox(-.0125,-.216,-.01,.025,.432,.02))]
        carrier,_ = occ.fuse(rim,parts)
        holes = [cyl(0,.006,z=-.011,height=.022),cyl(.2,.03,z=-.011,height=.022),
                 cyl(-.2,.03,z=-.011,height=.022)]
        carrier,_ = occ.cut(carrier,holes)
        inserts = [cyl(.2,.03,z=-.045,height=.09),cyl(-.2,.03,z=-.045,height=.09)]
        _, mapping = occ.fragment(carrier,inserts)
        aluminum = {tag for entries in mapping[:len(carrier)] for dim,tag in entries if dim==3}
        steel = {tag for entries in mapping[len(carrier):] for dim,tag in entries if dim==3}
        assert aluminum.isdisjoint(steel)
        occ.synchronize()
        stats = {}
        for name, volumes, rho in [('aluminum',aluminum,2700.),('steel',steel,7850.)]:
            mass = inertia = 0.
            for tag in volumes:
                v = occ.getMass(3,tag)
                com = occ.getCenterOfMass(3,tag)
                tensor = np.asarray(occ.getMatrixOfInertia(3,tag)).reshape(3,3)
                mass += rho*v
                inertia += rho*(tensor[2,2]+v*(com[0]**2+com[1]**2))
            stats[name+'_mass_kg'] = mass
            stats[name+'_polar_inertia_kg_m2'] = inertia
        gmsh.option.setNumber('Mesh.MeshSizeMin',h)
        gmsh.option.setNumber('Mesh.MeshSizeMax',h)
        # Smaller bore mesh resolves the actual 12 mm attachment.
        points = gmsh.model.getEntities(0)
        for dim,tag in points:
            p = gmsh.model.getValue(dim,tag,[])
            if np.hypot(p[0],p[1]) < .007:
                gmsh.model.mesh.setSize([(dim,tag)],min(h,.003))
        gmsh.option.setNumber('Mesh.MeshSizeMin',min(h,.003))
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)
        gmsh.write(str(folder/'mesh.msh'))
        gmsh.write(str(folder/'geometry.brep'))
        tags, xyz,_ = gmsh.model.mesh.getNodes()
        nodes = dict(zip(map(int,tags),xyz.reshape(-1,3)))
        elements = {'ALUMINUM':[],'STEEL':[]}
        for name,volumes in [('ALUMINUM',aluminum),('STEEL',steel)]:
            for vol in volumes:
                types,ids,conns = gmsh.model.mesh.getElements(3,vol)
                for typ,eids,conn in zip(types,ids,conns):
                    _,dim,_,n,local,primary = gmsh.model.mesh.getElementProperties(typ)
                    assert n==10 and primary==4
                    local = np.asarray(local).reshape(n,dim)
                    order = list(range(4))
                    for i,j in ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3)):
                        delta = np.linalg.norm(local-(local[i]+local[j])/2,axis=1)
                        assert min(delta)<1e-12
                        order.append(int(np.argmin(delta)))
                    elements[name] += [(int(e),list(map(int,c[order]))) for e,c in zip(eids,conn.reshape(-1,n))]
        return nodes,elements,stats
    finally:
        gmsh.finalize()


def solve(h,rpm,output,mesh_builder=None,spec=None):
    folder = output/f'h{h:g}_rpm{rpm:g}'
    folder.mkdir(parents=True,exist_ok=True)
    nodes,elements,stats = (mesh_builder or create_mesh)(h,folder)
    fixed = stats.pop('_fixed_nodes',None)
    if fixed is None:fixed = [n for n,p in nodes.items() if abs(np.hypot(p[0],p[1])-.006)<1e-7]
    assert len(fixed)>12
    lines = ['*HEADING','Shaped carrier with bonded steel inserts, SI units','*NODE,NSET=NALL']
    lines += [f'{n},'+','.join(f'{x:.12g}' for x in p) for n,p in nodes.items()]
    for name,group in elements.items():
        lines += [f'*ELEMENT,TYPE=C3D10,ELSET={name}']
        lines += [f'{e},'+','.join(map(str,c)) for e,c in group]
    lines += ['*ELSET,ELSET=EALL','ALUMINUM,STEEL']+nset('FIXED',fixed)
    for name,E,nu,rho in [('ALUMINUM',70e9,.33,2700.),('STEEL',200e9,.3,7850.)]:
        lines += [f'*MATERIAL,NAME={name}','*ELASTIC',f'{E},{nu}','*DENSITY',str(rho),
                  f'*SOLID SECTION,ELSET={name},MATERIAL={name}']
    lines += ['*BOUNDARY','FIXED,1,3']
    omega = rpm*2*np.pi/60
    if rpm:
        lines += ['*STEP','*STATIC','*DLOAD',f'EALL,CENTRIF,{omega**2:.12g},0,0,0,0,0,1',
                  '*NODE PRINT,NSET=NALL','U','*EL PRINT,ELSET=ALUMINUM','S',
                  '*NODE FILE','U','*EL FILE','S','*END STEP']
    lines += ['*STEP,PERTURBATION' if rpm else '*STEP','*FREQUENCY','8','*NODE FILE','U','*END STEP']
    (folder/'model.inp').write_text('\n'.join(lines)+'\n')
    # Linux scratch avoids thousands of small writes through /mnt/c. Preserve
    # all solver output, including failure evidence, in the workspace afterward.
    with tempfile.TemporaryDirectory(prefix='gravcomm-ccx-') as temporary:
        work=Path(temporary)
        shutil.copyfile(folder/'model.inp',work/'model.inp')
        try:
            result = subprocess.run(['ccx','-i','model'],cwd=work,text=True,stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'2'},timeout=600)
        finally:
            for path in work.iterdir():
                if path.is_file():shutil.copyfile(path,folder/path.name)
    (folder/'solver.log').write_text(result.stdout)
    if result.returncode or '*ERROR' in result.stdout or 'Job finished' not in result.stdout:
        raise RuntimeError(f'Solver failed: {folder}')
    blocks,freqs = parse_dat(folder/'model.dat')
    assert len(freqs)==8 and all(f>0 for f in freqs)
    stress = np.array(blocks['stresses'])
    vm = np.zeros(1)
    if stress.size:
        s = stress[:,2:8]
        vm = np.sqrt(((s[:,0]-s[:,1])**2+(s[:,1]-s[:,2])**2+(s[:,2]-s[:,0])**2)/2+3*np.sum(s[:,3:]**2,axis=1))
    disp = np.array(blocks['displacements'])
    row = dict(mesh_size_m=h,rpm=rpm,nodes=len(nodes),elements=sum(map(len,elements.values())),
               first_frequency_hz=freqs[0],max_aluminum_integration_point_vm_pa=float(max(vm)),
               max_displacement_m=float(max(np.linalg.norm(disp[:,1:4],axis=1))) if disp.size else 0.,
               **stats)
    (folder/'metadata.json').write_text(json.dumps(dict(spec=spec or SPEC,summary=row,frequencies_hz=freqs),indent=2))
    print(json.dumps(row),flush=True)
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--quick',action='store_true')
    args = parser.parse_args()
    output = ROOT/'results'/'shaped_rotor'
    output.mkdir(parents=True,exist_ok=True)
    (output/'spec.json').write_text(json.dumps(SPEC,indent=2))
    rows=[]
    for h in ([.025] if args.quick else [.025,.018,.013]):
        for rpm in (0.,1800.):
            rows.append(solve(h,rpm,output))
            with (output/'summary.csv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)


if __name__=='__main__':
    main()
