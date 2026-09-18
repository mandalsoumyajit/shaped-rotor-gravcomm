"""Explicit 3 mm junction fillets with insert mass adjusted to fixed total mass."""
import argparse
import json
import shutil
from design_sweep import OUT,design,mesh_builder,gravity
from shaped_rotor import solve
from postprocess_carrier import write_csv


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--geometry-only',action='store_true')
    args=parser.parse_args()
    cases={'baseline_rounded':dict(shape='cylinder',radius=.03,center=.2,spoke=.025,rim=.035,fillet=.003,faceted_geometry=True),
           'medium_rounded':dict(shape='cylinder',radius=.035,center=.2,spoke=.03,rim=.02,fillet=.003,faceted_geometry=True),
           'wide_rounded':dict(shape='cylinder',radius=.035,center=.2,spoke=.04,rim=.02,fillet=.003,faceted_geometry=True)}
    chosen={name:design(spec) for name,spec in cases.items()}
    assert all(chosen.values()),'Filleted geometry violates the fixed-mass envelope'
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'filleted_assumptions.json').write_text(json.dumps(chosen,indent=2))
    print(json.dumps(chosen,indent=2),flush=True)
    if args.geometry_only:return
    rows=[]
    for name,candidate in chosen.items():
        for h in (.018,.013):
            folder=OUT/name/f'h{h:g}_rpm1800'
            if all((folder/file).exists() for file in ('metadata.json','model.inp','model.dat','model.frd')):
                saved=json.loads((folder/'metadata.json').read_text())
                assert saved['spec']==candidate,'Refusing to reuse a different geometry'
                row=saved['summary']
            else:
                if (folder/'solver.log').exists():shutil.copyfile(folder/'solver.log',folder/'solver_failed_initial.log')
                row=solve(h,1800.,OUT/name,mesh_builder(candidate),spec=candidate)
            row.update(gravity(folder))
            row.update(case=name,shape=candidate['shape'],center_m=candidate['center'],
                spoke_width_m=candidate['spoke'],rim_width_m=candidate['rim'],fillet_radius_m=.003)
            rows.append(row);write_csv(OUT/'filleted_screen.csv',rows)


if __name__=='__main__':main()
