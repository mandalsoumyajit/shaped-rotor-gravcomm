"""Gmsh/CalculiX verification study. Run with Ubuntu system Python under WSL.

SI units throughout. A fixed-root half rotor represents one of two identical
arms; it does not model the shaft, bearings or coupling between the arms.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess

import gmsh
import numpy as np

ROOT = Path(__file__).resolve().parent
E, NU, RHO = 200e9, .3, 7850.


def geometry(kind):
    if kind == 'bar':
        return dict(root=.03, end=.3, width=.01, mass=0., center=.3)
    # End block is exactly 25 kg; its center matches the Bench point mass.
    side = (25/RHO)**(1/3)
    return dict(root=.03, end=.3-side/2, width=.12, mass=25., center=.3, side=side)


def mesh(kind, size, folder):
    g = geometry(kind)
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.option.setNumber('General.NumThreads', 1)
        gmsh.model.add(kind)
        w = g['width']
        arm = gmsh.model.occ.addBox(g['root'], -w/2, -w/2, g['end']-g['root'], w, w)
        if g['mass']:
            s = g['side']
            block = gmsh.model.occ.addBox(g['end'], -s/2, -s/2, s, s, s)
            gmsh.model.occ.fuse([(3, arm)], [(3, block)])
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber('Mesh.MeshSizeMin', size)
        gmsh.option.setNumber('Mesh.MeshSizeMax', size)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)
        gmsh.write(str(folder / 'mesh.msh'))
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        nodes = dict(zip(map(int, tags), xyz.reshape(-1, 3)))
        types, etags, conns = gmsh.model.mesh.getElements(3)
        elements = []
        for typ, ids, connectivity in zip(types, etags, conns):
            name, dim, order, n, local, primary = gmsh.model.mesh.getElementProperties(typ)
            if n != 10 or primary != 4:
                raise RuntimeError(f'Expected quadratic tetrahedra, got {name}')
            # Derive edge-node permutation from reference coordinates rather
            # than assuming Gmsh and CalculiX use the same tetra10 ordering.
            loc = np.asarray(local).reshape(n, dim)
            order_ccx = list(range(4))
            for i, j in ((0,1), (1,2), (2,0), (0,3), (1,3), (2,3)):
                distances = np.linalg.norm(loc-(loc[i]+loc[j])/2, axis=1)
                index = int(np.argmin(distances))
                assert distances[index] < 1e-12
                order_ccx.append(index)
            for eid, conn in zip(ids, connectivity.reshape(-1, n)):
                elements.append((int(eid), list(map(int, conn[order_ccx]))))
        return g, nodes, elements
    finally:
        gmsh.finalize()


def nset(name, ids):
    ids = list(ids)
    return [f'*NSET,NSET={name}'] + [','.join(map(str, ids[i:i+16])) for i in range(0,len(ids),16)]


def deck(g, nodes, elements, omega, modal):
    fixed = [n for n, p in nodes.items() if abs(p[0]-g['root']) < 1e-8]
    xmax = max(p[0] for p in nodes.values())
    tip = [n for n,p in nodes.items() if abs(p[0]-xmax) < 1e-8]
    lines = ['*HEADING', 'Rigid-root radial arm; SI units; quadratic tetrahedra', '*NODE,NSET=NALL']
    lines += [f'{n},'+','.join(f'{x:.12g}' for x in p) for n,p in nodes.items()]
    lines += ['*ELEMENT,TYPE=C3D10,ELSET=EALL']
    lines += [f'{eid},'+','.join(map(str, conn)) for eid,conn in elements]
    lines += nset('FIXED',fixed) + nset('TIP',tip)
    lines += ['*MATERIAL,NAME=STEEL','*ELASTIC',f'{E},{NU}','*DENSITY',str(RHO),
              '*SOLID SECTION,ELSET=EALL,MATERIAL=STEEL','*BOUNDARY','FIXED,1,3']
    # Linear static benchmark: body force evaluated on undeformed geometry.
    if omega:
        lines += ['*STEP','*STATIC','*DLOAD',f'EALL,CENTRIF,{omega**2:.12g},0,0,0,0,0,1',
                  '*NODE PRINT,NSET=FIXED','RF','*NODE PRINT,NSET=TIP','U',
                  '*EL PRINT,ELSET=EALL','S','*NODE FILE','U','*EL FILE','S','*END STEP']
    if modal:
        lines += ['*STEP,PERTURBATION' if omega else '*STEP','*FREQUENCY','8',
                  '*NODE FILE','U','*END STEP']
    return '\n'.join(lines)+'\n'


def parse_dat(path):
    blocks = {'forces': [], 'displacements': [], 'stresses': []}
    mode = None
    frequencies = []
    for line in path.read_text().splitlines():
        low = line.lower().strip()
        # Later modal stress/reaction tables are mass-normalized mode shapes,
        # not physical static stresses. Stop before those tables begin.
        if 'participationfactors' in ''.join(low.split()):
            break
        if 'eigenvalue' in low and 'frequency' in low:
            mode = 'frequency'
            continue
        for key in blocks:
            if low.startswith(key+' '):
                mode = key
        fields = line.split()
        if not fields or not fields[0].isdigit():
            continue
        try:
            values = [float(x.replace('D','E')) for x in fields]
        except ValueError:
            continue
        if mode in blocks:
            blocks[mode].append(values)
        elif mode == 'frequency' and len(values) >= 4:
            if len(values)>=5 and values[4]!=0:
                raise ValueError('Nonzero imaginary eigenfrequency: inspect stability before reporting modes')
            frequencies.append(values[3])
    return blocks, frequencies


def run(kind, size, omega, modal, output):
    name = f'{kind}_h{size:g}_w{omega:.6g}'
    folder = output / name
    folder.mkdir(parents=True, exist_ok=True)
    g, nodes, elements = mesh(kind,size,folder)
    (folder/'model.inp').write_text(deck(g,nodes,elements,omega,modal))
    process = subprocess.run(['ccx','-i','model'], cwd=folder, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             env={**os.environ,'OMP_NUM_THREADS':'2'}, timeout=300)
    (folder/'solver.log').write_text(process.stdout)
    if process.returncode or '*ERROR' in process.stdout or 'Job finished' not in process.stdout:
        raise RuntimeError(f'CalculiX failed: {folder}/solver.log')
    blocks, freqs = parse_dat(folder/'model.dat')
    if modal and (len(freqs) != 8 or not all(f > 0 for f in freqs)):
        raise RuntimeError(f'Expected eight positive real frequencies: {freqs}')
    length = g['end']-g['root']
    area = g['width']**2
    force = omega**2*(g['mass']*g['center']+RHO*area*(g['end']**2-g['root']**2)/2)
    row = dict(case=name, kind=kind, mesh_size_m=size, nodes=len(nodes), elements=len(elements),
               omega_rad_s=omega, analytical_root_force_n=force,
               fea_root_force_n=-sum(v[1] for v in blocks['forces']),
               nominal_root_stress_pa=force/area,
               mean_tip_displacement_m=float(np.mean([v[1] for v in blocks['displacements']])) if omega else 0.,
               first_frequency_hz=freqs[0] if freqs else None)
    if kind == 'bar':
        row['analytical_tip_displacement_m'] = RHO*omega**2/(2*E)*(g['end']**2*length-(g['end']**3-g['root']**3)/3)
        row['beam_first_frequency_hz'] = 1.8751040687**2/(2*math.pi*length**2)*math.sqrt(E*g['width']**2/(12*RHO))
    else:
        row['analytical_tip_displacement_m'] = None
        row['beam_first_frequency_hz'] = None
    stress = np.asarray(blocks['stresses'])
    if stress.size:
        s = stress[:,2:8]
        vm = np.sqrt(((s[:,0]-s[:,1])**2+(s[:,1]-s[:,2])**2+(s[:,2]-s[:,0])**2)/2
                     +3*np.sum(s[:,3:]**2,axis=1))
        row['max_integration_point_von_mises_pa'] = float(max(vm))
    else:
        row['max_integration_point_von_mises_pa'] = None
    (folder/'metadata.json').write_text(json.dumps(dict(geometry=g,youngs_modulus_pa=E,
        poisson_ratio=NU,density_kg_m3=RHO,frequencies_hz=freqs,summary=row),indent=2))
    print(json.dumps(row),flush=True)
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--quick',action='store_true',help='one benchmark mesh')
    args = parser.parse_args()
    output = ROOT / 'results' / 'benchmark'
    output.mkdir(parents=True,exist_ok=True)
    versions = subprocess.run(['dpkg-query','-W','calculix-ccx','python3-gmsh'],capture_output=True,text=True)
    (output/'versions.txt').write_text(versions.stdout+f'Gmsh Python: {gmsh.__version__}\n')
    rows = []
    for h in ([.008] if args.quick else [.01,.006,.004]):
        for omega in (0.,250/.3):
            rows.append(run('bar',h,omega,True,output))
    with (output/'summary.csv').open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
