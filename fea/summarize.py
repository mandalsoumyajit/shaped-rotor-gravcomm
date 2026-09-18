"""Assemble analytical checks and plots from completed CalculiX runs."""
import csv
import json
from pathlib import Path
import shutil

import gmsh
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'


def main():
    benchmark=OUT/'benchmark'
    benchmark.mkdir(exist_ok=True)
    # Relocate only this task's completed benchmark outputs from the initial
    # mixed exploratory run, keeping all input decks and solver evidence.
    for h in (.01,.006,.004):
        for omega in (0.,250/.3):
            name=f'bar_h{h:g}_w{omega:.6g}'
            src=(OUT/name).resolve()
            dst=(benchmark/name).resolve()
            assert src.parent==OUT.resolve() and dst.parent==benchmark.resolve()
            if src.exists() and not dst.exists():
                src.rename(dst)
    if (OUT/'versions.txt').exists():
        shutil.copyfile(OUT/'versions.txt',benchmark/'versions.txt')
    bars=[json.loads(p.read_text())['summary'] for p in sorted(benchmark.glob('bar_*/metadata.json'))]
    assert len(bars)==6
    with (benchmark/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(bars[0]))
        writer.writeheader();writer.writerows(bars)
    spinning=sorted([r for r in bars if r['omega_rad_s']],key=lambda r:r['mesh_size_m'])
    b=spinning[0]
    stationary=min((r for r in bars if not r['omega_rad_s']),key=lambda r:r['mesh_size_m'])
    force_error=abs(b['fea_root_force_n']/b['analytical_root_force_n']-1)
    displacement_error=abs(b['mean_tip_displacement_m']/b['analytical_tip_displacement_m']-1)
    frequency_error=abs(stationary['first_frequency_hz']/stationary['beam_first_frequency_hz']-1)
    assert force_error<.001 and displacement_error<.005 and frequency_error<.005
    shaped=list(csv.DictReader((OUT/'shaped_rotor'/'summary.csv').open()))
    assert len(shaped)==6, 'Wait for all shaped-rotor simulations to finish'
    rotating=sorted([r for r in shaped if float(r['rpm'])>0],key=lambda r:float(r['mesh_size_m']))
    finest,previous=rotating[:2]
    fig,axes=plt.subplots(1,3,figsize=(12,3.7))
    for ax,key,label,scale in zip(axes,['first_frequency_hz','max_displacement_m','max_aluminum_integration_point_vm_pa'],
                                 ['First prestressed mode (Hz)','Maximum displacement (micrometers)',
                                  'Raw peak aluminium stress (MPa)'],[1,1e6,1e-6]):
        ax.plot([float(r['mesh_size_m'])*1000 for r in rotating],
                [float(r[key])*scale for r in rotating],'o-')
        ax.set_xlabel('Nominal mesh size (mm)');ax.set_ylabel(label);ax.grid(alpha=.3)
    fig.suptitle('Shaped carrier at 1800 rpm: fixed bore, ideal bonded inserts')
    fig.tight_layout();fig.savefig(OUT/'shaped_rotor'/'convergence.png',dpi=160);plt.close(fig)
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal',0)
    gmsh.open(str(OUT/'shaped_rotor'/'h0.013_rpm1800'/'mesh.msh'))
    tags,xyz,_=gmsh.model.mesh.getNodes()
    points=dict(zip(map(int,tags),xyz.reshape(-1,3)))
    faces=[];colors=[]
    types,ids,conns=gmsh.model.mesh.getElements(2)
    for typ,conn in zip(types,conns):
        _,_,_,n,_,primary=gmsh.model.mesh.getElementProperties(typ)
        if primary!=3: continue
        for c in conn.reshape(-1,n):
            pts=np.array([points[int(tag)] for tag in c[:3]])
            if np.max(abs(pts[:,2]-.01))<1e-7 or np.max(abs(pts[:,2]-.045))<1e-7:
                faces.append(pts[:,:2]*1000)
                colors.append('#aabacb' if abs(pts[0,2]-.01)<1e-7 else '#6d737c')
    gmsh.finalize()
    fig,ax=plt.subplots(figsize=(7,7))
    ax.add_collection(PolyCollection(faces,facecolors=colors,edgecolors='#34404c',linewidths=.2))
    ax.autoscale();ax.set_aspect('equal');ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)')
    ax.set_title('Shaped rotor: top surface mesh\nSteel inserts project above the aluminium carrier')
    fig.tight_layout();fig.savefig(OUT/'shaped_rotor'/'geometry_mesh.png',dpi=160);plt.close(fig)
    changes={key:abs(float(finest[key])/float(previous[key])-1)*100 for key in
             ['first_frequency_hz','max_displacement_m','max_aluminum_integration_point_vm_pa']}
    report=f'''# Initial CalculiX results

Follow-up qualification: local refinement around explicit fillets changed the
peak stress substantially. The original 38.9 MPa value is a coarse-local-mesh
diagnostic, not a converged strength result. See STRESS_MODES_AND_DESIGN.md.

The active design is the shaped spoked carrier from `cad/bench_source_disk.scad`.
The earlier arm-and-block calculations are superseded and excluded.

## Solver verification

Three quadratic-tetrahedral meshes were run for the rotating-bar benchmark.
At the finest mesh, centrifugal reaction differs from the analytical resultant
by {100*force_error:.4f}%, extension by {100*displacement_error:.4f}%, and the
stationary first bending frequency from Euler-Bernoulli theory by
{100*frequency_error:.4f}%. All verification thresholds pass
(0.1% force, 0.5% extension, 0.5% frequency).

## Shaped rotor

Three meshes, each at zero speed and 1800 rpm, completed successfully. The model
uses a fixed 12 mm bore, a 20 mm aluminium carrier and ideal bonded steel inserts.
The total mass is {float(finest['aluminum_mass_kg'])+float(finest['steel_mass_kg']):.5f} kg;
polar inertia is {float(finest['aluminum_polar_inertia_kg_m2'])+float(finest['steel_polar_inertia_kg_m2']):.6f} kg m^2.

Finest-mesh preliminary results at 1800 rpm:

- First prestressed structural mode: {float(finest['first_frequency_hz']):.3f} Hz.
- Maximum displacement: {float(finest['max_displacement_m'])*1e6:.3f} micrometers.
- Raw maximum aluminium integration-point von Mises stress:
  {float(finest['max_aluminum_integration_point_vm_pa'])*1e-6:.3f} MPa.

Changes from the intermediate to the finest mesh:

- First mode: {changes['first_frequency_hz']:.3f}%.
- Maximum displacement: {changes['max_displacement_m']:.3f}%.
- Raw peak stress: {changes['max_aluminum_integration_point_vm_pa']:.3f}%.

The sharp-junction peak stress is a diagnostic, not a certified strength margin.
Further local refinement and physical fillets are needed before interpreting
peak stress. The fixed-bore modes do not validate the earlier shaft-whirl estimate.
The CAD clearance, attachment hardware, shaft, bearings, gyroscopic mode analysis,
modulation loading, fatigue and contact are not modeled in this initial study.

The intended shaft bore is cut after the carrier union; unlike the original
OpenSCAD ordering, it is not refilled by the spokes. Pocket clearance is removed
only for the explicitly bonded-interface approximation. No source CAD is edited.

See `README.md` for reproducible commands, solver references and assumptions;
`results/shaped_rotor/summary.csv` for all six runs; `geometry_mesh.png` and
`convergence.png` for the geometry and refinement plots. Meshes, input decks,
raw results and solver logs are retained for inspection.
'''
    (ROOT/'INITIAL_RESULTS.md').write_text(report)
    print(report)


if __name__=='__main__': main()
