"""Refine promising fixed-mass candidates and assess a sphere near its radial limit."""
import csv
import json
from scipy.optimize import brentq
from design_sweep import ROOT,OUT,design,mesh_builder,gravity
from shaped_rotor import solve
from postprocess_carrier import write_csv


def main():
    assumptions=json.loads((OUT/'assumptions.json').read_text())
    chosen=assumptions['selected']
    chosen['medium_spokes']=design(dict(shape='cylinder',radius=.035,center=.2,spoke=.03,rim=.02))
    def extent(c):
        r=design(dict(shape='sphere',radius=.03,center=c,spoke=.02,rim=.02),enforce_envelope=False)
        # Outer extent is max(.25,...), so solve against actual boss edge.
        return c+r['insert_radius_m']+.015-.25
    center=brentq(extent,.18,.2,xtol=1e-8)
    chosen['sphere_outer_limit']=design(dict(shape='sphere',radius=.03,center=center,spoke=.02,rim=.02))
    assert chosen['sphere_outer_limit'] is not None
    assumptions['selected']=chosen
    (OUT/'assumptions.json').write_text(json.dumps(assumptions,indent=2))
    rows=[]
    for name,h in [('medium_spokes',.018),('sphere_outer_limit',.018),('best_quadrupole',.013),
                   ('medium_spokes',.013),('wide_spokes',.013),('sphere_outer_limit',.013)]:
        candidate=chosen[name]
        row=solve(h,1800.,OUT/name,mesh_builder(candidate),spec=candidate)
        row.update(gravity(OUT/name/f'h{h:g}_rpm1800'))
        row.update(case=name,shape=candidate['shape'],center_m=candidate['center'],
                   spoke_width_m=candidate['spoke'],rim_width_m=candidate['rim'])
        rows.append(row);write_csv(OUT/'refined_screen.csv',rows)
    # Fine original baseline is an independent mesh of the same geometry.
    folder=ROOT/'results'/'shaped_rotor'/'h0.013_rpm1800'
    base=json.loads((folder/'metadata.json').read_text())['summary']
    base.update(gravity(folder));base['case']='baseline'
    (OUT/'fine_baseline.json').write_text(json.dumps(base,indent=2))
    # Temporal quadrature check on a representative finite-volume waveform.
    folder=OUT/'medium_spokes'/'h0.013_rpm1800'
    check=gravity(folder,nphase=256)
    reference=next(r for r in rows if r['case']=='medium_spokes' and r['mesh_size_m']==.013)
    errors={k:abs(check[k]/reference[k]-1) for k in check if k.startswith('a2')}
    assert max(errors.values())<1e-10
    (OUT/'gravity_phase_convergence.json').write_text(json.dumps(errors,indent=2))


if __name__=='__main__':main()
