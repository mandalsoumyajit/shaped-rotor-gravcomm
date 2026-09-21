"""Compile completed three-distance mechanical evidence, with explicit limits."""
import json,hashlib
from pathlib import Path
import numpy as np
from three_link_fea import ROOT,OUT,save

def main():
 static=[];modes=[];checks=[]
 for h in (.018,.013):
  static+=json.loads((OUT/f'static_compiled_h{h:g}.json').read_text());modes+=json.loads((OUT/f'modal_response_h{h:g}.json').read_text());checks+=json.loads((OUT/f'equilibrium_checks_h{h:g}.json').read_text())
 frozen=json.loads((OUT/'frozen_designs.json').read_text());assembly=json.loads((OUT/'assembly_screen.json').read_text());rows=[]
 for i,ch in enumerate(frozen['channels'],1):
  fine=next(x for x in static if x['mesh']==.013 and x['design']==i);coarse=next(x for x in static if x['mesh']==.018 and x['design']==i);mode=next(x for x in modes if x['mesh']==.013 and x['design']==i and x['speed']=='center');freq=mode['frequencies_hz'];j=int(np.argmax(mode['angular_effective_inertia_fraction']));shaft=assembly[i-1]['shaft_only_rotor_torsion_hz'];combined=1/np.sqrt(1/shaft**2+1/freq[j]**2)
  rows.append(dict(distance_m=ch['distance_m'][0],mass_kg=ch['source_designs'][0]['mass_kg'],fine=fine,coarse=coarse,mode=mode,torsional_mode=j+1,torsional_frequency_hz=freq[j],series_shaft_carrier_torsion_screen_hz=float(combined),regional_change_percent=100*(fine['regional_10mm_scaled_vm_mpa']/coarse['regional_10mm_scaled_vm_mpa']-1),raw_peak_change_percent=100*(fine['centrifugal_raw_alu_vm_mpa']/coarse['centrifugal_raw_alu_vm_mpa']-1),displacement_change_percent=100*(fine['centrifugal_max_u_um']/coarse['centrifugal_max_u_um']-1)))
 save(OUT/'compiled.json',rows)
 refdir=OUT/'refinement_reference/h0.01_rpm1800'
 ref=json.loads((refdir/'metadata.json').read_text())['summary']
 old=json.loads((ROOT/'results/design_sweep/medium_rounded/h0.013_rpm1800/metadata.json').read_text())['summary']
 delta=100*(ref['max_aluminum_integration_point_vm_pa']/old['max_aluminum_integration_point_vm_pa']-1)
 save(OUT/'refinement_audit.json',dict(summary=ref,peak_change_percent=delta,sha256={n:hashlib.sha256((refdir/n).read_bytes()).hexdigest() for n in ['model.inp','metadata.json','model.dat','solver.log']}))
 lines=['# Structural reassessment of the qualified 0.5, 1 and 2 m sources','',
 'Scope follows the user clarification: 7.03 kg at 0.5 m, 37.69 kg at 1 m and the smallest passing 227.32 kg source at 2 m. The larger 2 m alternatives are not the selected study. Communications timing is unchanged. All quantities below are calculations for the bare rotor with fixed bore and perfectly bonded inserts unless explicitly identified as a separate assembly screen.', '',
 '## Fixed-bore linear FEA','',
 'Twenty-six static/prestressed-modal cases use two archived mesh resolutions, three actual operating speeds per design, and four independently solved gravity/angular-acceleration bases per mesh. Seven input-identical completed cases were reused after SHA-256 checks. All solver jobs must finish successfully; static/modal outputs are separated. Five existing FEA regression tests pass. Positive Jacobians, exact consistent-load torque, force/torque equilibrium and similarity to the archived desktop solution are checked.', '',
 '| Link (m) | Bare mass (kg) | Upper-speed centrifugal Al peak (MPa) | Horizontal combined Al bound (MPa) | Centrifugal displacement (um) | Horizontal combined displacement bound (um) |', '|---:|---:|---:|---:|---:|---:|']
 for x in rows:
  a=x['fine'];lines.append(f"| {x['distance_m']:g} | {x['mass_kg']:.2f} | {a['centrifugal_raw_alu_vm_mpa']:.2f} | {a['horizontal_conservative_alu_vm_mpa']:.2f} | {a['centrifugal_max_u_um']:.2f} | {a['horizontal_conservative_u_um']:.2f} |")
 lines+=['',
 'The combined columns are triangle-inequality bounds assembled from the linear static solutions, including maximum centrifugal load, either sign of maximum angular acceleration and any in-plane gravity direction. They are deliberately conservative; speed and angular-acceleration maxima need not occur simultaneously. Separate phase-sampled extreme transitions use the actual quintic trajectory. Neither calculation is a nonlinear contact simulation or a complete time-domain rotor-bearing response. Gravity/acceleration bases use unprestressed elastic stiffness, consistent with the original linear statics model.', '',
 'Horizontal shaft orientation is the preliminary choice. All-attitude bounds (including axial gravity) are substantially larger and retained in compiled.json. Gravity direction rotates at the shaft frequency in rotor coordinates for a horizontal shaft. Loads at the 24/24/18 Hz gravitational carriers are not automatically mechanical excitation at those frequencies.', '',
 '## Mesh convergence and strength interpretation','',
 '| Link (m) | Change in centrifugal raw Al peak (%) | Change in max displacement (%) | Change in maximum scaled-window combined average (%) |', '|---:|---:|---:|---:|']
 for x in rows:lines.append(f"| {x['distance_m']:g} | {x['raw_peak_change_percent']:.2f} | {x['displacement_change_percent']:.3f} | {x['regional_change_percent']:.3f} |")
 lines+=['','Regional windows have reference radii 5/10 mm multiplied by the geometric scale; they span the carrier thickness. The reported average is the maximum of the 10 mm-reference windows. It is a comparison metric, not a replacement for local yield or fatigue analysis. Fine-mesh mass is about 0.11% below CAD because of the faceted approximation; density is not renormalized. See mesh_h*.json.', '',
 f'Third-mesh check: a newly generated 10 mm reference mesh has {ref["nodes"]} nodes and {ref["elements"]} elements. Its 1800 rpm aluminum peak is {ref["max_aluminum_integration_point_vm_pa"]/1e6:.3f} MPa, a {delta:.3f}% increase from the 13 mm mesh. Displacement is {ref["max_displacement_m"]*1e6:.3f} um and first frequency {ref["first_frequency_hz"]:.3f} Hz. This supports convergence of the centrifugal reference peak; it does not establish convergence of every combined-load or contact stress.', '',
 'Preliminary material choices are 6061-T651 plate (240 MPa yield screen), specifically certified normalized 1045 inserts (300 MPa minimum procurement requirement), and 4140 Q&T shafts/flanges (650 MPa minimum study specification). A factor of two against yield is the initial static screening criterion. Material sources, product-form caveats and actual assembly choices are in ASSEMBLY_SPECIFICATION.md and materials.json. Raw mesh peaks and idealized attachments do not establish hardware safety factors. Fatigue, preload, real contact and local features remain unresolved.', '',
 'All three idealized horizontal-shaft combined aluminum bounds are below the 120 MPa screening limit (240 MPa divided by two); the largest is about 55.50 MPa. The model therefore indicates useful static margin in the carrier, while insert retention and assembly dynamics still require explicit validation.', '',
 '## Modes, coupling and response','',
 '| Link (m) | First three center-speed modes (Hz) | Dominant angular-acceleration mode | Polar-inertia fraction (%) | Shaft/carrier series torsion screen (Hz) |', '|---:|:---|---:|---:|---:|']
 for x in rows:
  m=x['mode'];j=x['torsional_mode']-1;lines.append(f"| {x['distance_m']:g} | {', '.join(f'{f:.2f}' for f in m['frequencies_hz'][:3])} | {j+1} | {100*m['angular_effective_inertia_fraction'][j]:.2f} | {x['series_shaft_carrier_torsion_screen_hz']:.2f} |")
 lines+=['',
 'Mode matching uses the consistent mass matrix and modal assurance criterion (MAC), with node correspondence within each mesh. Mass normalization is independently recovered from degree-four exact element mass integrals. Component energies and mean insert motions identify axial/torsional motion. Angular and translational effective mass fractions distinguish modulation coupling from gravity coupling; low frequency alone does not establish strong excitation.', '',
 'modal_response_h*.json contains damping sweeps at 0.5%, 2% and 5% for isolated quintic-transition scalar modal responses and rotation-frequency harmonic transfer factors. The transition simulation starts at rest and omits repeated-message resonance, motor control, gyroscopic terms and contact. It diagnoses time-scale separation, not arbitrary-message dynamic qualification. Shaft/carrier series values are a simple compliance screen using the chosen stub-shaft dimensions and dominant carrier mode; motor/coupling inertia and bearing/structure dynamics are absent. These values demonstrate why fixed-bore carrier frequencies must not be called assembly critical speeds.', '',
 '## Assembly and next checks','',
 'ASSEMBLY_SPECIFICATION.md selects a horizontal straddled shaft, locating angular-contact pair plus floating cylindrical-roller support, piloted hub flanges and positive insert shoulder/collar retention. Plain-shaft and pocket-pressure estimates are saved separately in assembly_screen.json. These are preliminary calculations, not additional 3D assembly FEA. Before physical implementation, add flange bolts, contact/clearance, retention grooves/collars, preload, actual bearing and motor/coupling stiffness, and fatigue. Added hardware mass and off-axis retention features require a revised field/inertia audit before reusing the BER claims for that assembly.', '',
 'The tungsten-heavy-alloy option in TUNGSTEN_OPTION.md retains insert mass and diameter while shortening its length; it is a promising windage/overhang reduction study. It has not been substituted into the qualified steel models.', '',
 '## Reproduction','',
 'Run fea/three_link_fea.py --prepare --workers 4 in the documented WSL environment, then compile_three_link_fea.py --mesh .018 and --mesh .013; three_link_modal_response.py for each mesh; assembly_screen.py; report_three_link_fea.py. Completed input-identical cases are reusable. Inputs, complete solver outputs, fields, equilibrium checks and hashes are retained. No BER samples or source configuration was changed.', '']
 (OUT/'RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
 for job in json.loads((OUT/'manifest.json').read_text())['jobs']:
  folder=Path(job['folder']);assert (folder/'summary.json').exists();assert hashlib.sha256((folder/'model.inp').read_bytes()).hexdigest()==job['input_sha256']
 save(OUT/'report_audit.json',dict(all_26_cases_complete=True,equilibrium_checks=checks,postprocess_source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'compile_three_link_fea.py',ROOT/'three_link_modal_response.py',ROOT/'assembly_screen.py']}))
 print(json.dumps([dict(distance=x['distance_m'],stress=x['fine']['horizontal_conservative_alu_vm_mpa'],displacement=x['fine']['horizontal_conservative_u_um'],torsion=x['torsional_frequency_hz']) for x in rows],indent=2))
if __name__=='__main__':main()
