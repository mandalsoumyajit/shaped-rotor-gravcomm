"""Compare sphere and cylinder at the same center, spoke and rim dimensions."""
import json
from design_sweep import OUT,design,mesh_builder,gravity
from shaped_rotor import solve
from postprocess_carrier import write_csv


def main():
    rows=[]
    candidate=design(dict(shape='cylinder',radius=.035,center=.18,spoke=.02,rim=.02))
    name='cylinder_shape_control'
    row=solve(.018,1800.,OUT/name,mesh_builder(candidate),spec=candidate)
    row.update(gravity(OUT/name/'h0.018_rpm1800'));row['case']=name;rows.append(row)
    folder=OUT/'best_sphere'/'h0.018_rpm1800'
    row=json.loads((folder/'metadata.json').read_text())['summary']
    row.update(gravity(folder));row['case']='sphere_shape_control';rows.append(row)
    write_csv(OUT/'matched_shape_pair.csv',rows)


if __name__=='__main__':main()
