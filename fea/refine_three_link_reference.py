import json
from pathlib import Path
from three_link_fea import ROOT,OUT
from design_sweep import mesh_builder
from shaped_rotor import solve
spec=json.loads((ROOT/'results/design_sweep/medium_rounded/h0.013_rpm1800/metadata.json').read_text())['spec']
if __name__=='__main__':
 row=solve(.010,1800.,OUT/'refinement_reference',mesh_builder(spec),spec=spec)
 print(json.dumps(row),flush=True)
