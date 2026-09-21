"""Compile audited steel / tungsten comparison without changing steel results."""
import json,hashlib
from pathlib import Path
import numpy as np
from tungsten_link_fea import ROOT,OUT,save
BASE=ROOT/'results/three_link_designs_2026_09_20'
def read(p):return json.loads(p.read_text())
def main():
 rows=[];audit=[]
 for h in (.018,.013):
  steel=read(BASE/f'static_compiled_h{h:g}.json');tungsten=read(OUT/f'static_compiled_h{h:g}.json');sm=read(BASE/f'modal_response_h{h:g}.json');tm=read(OUT/f'modal_response_h{h:g}.json')
  for i,(a,b) in enumerate(zip(steel,tungsten),1):
   assert a['design']==b['design']==i
   sa=next(x for x in sm if x['design']==i and x['speed']=='center');tb=next(x for x in tm if x['design']==i and x['speed']=='center')
   rows.append(dict(design=i,distance_m=(.5,1.,2.)[i-1],mesh=h,steel=a,tungsten=b,steel_modes=sa,tungsten_modes=tb))
  # Independent cross-design scale check (identical mesh, different size and speed).
  jobs=read(OUT/'manifest.json')['jobs'];r=[]
  for b in tungsten:
   job=next(j for j in jobs if j['mesh']==h and j['name']==f"design{b['design']}_high")
   r.append([b['centrifugal_raw_alu_vm_mpa']/b['scale']**2/job['rpm']**2,b['centrifugal_max_u_um']/b['scale']**3/job['rpm']**2])
  err=np.max(abs(np.array(r)/r[0]-1));assert err<.002
  audit.append(dict(mesh=h,similarity_relative_error=float(err),equilibrium=read(OUT/f'equilibrium_checks_h{h:g}.json')))
 manifest=read(OUT/'manifest.json')
 for job in manifest['jobs']:
  folder=Path(job['folder']);assert (folder/'summary.json').exists();assert hashlib.sha256((folder/'model.inp').read_bytes()).hexdigest()==job['input_sha256']
 for rel,digest in manifest['source_sha256'].items():assert hashlib.sha256((ROOT.parent/rel).read_bytes()).hexdigest()==digest
 save(OUT/'comparison.json',rows);save(OUT/'audit.json',dict(solver_cases=len(manifest['jobs']),checks=audit,postprocess_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'compile_tungsten_link_fea.py',ROOT/'tungsten_modal_response.py',ROOT/'tungsten_field_comparison.py',ROOT/'report_tungsten_comparison.py']}))
 fine=[x for x in rows if x['mesh']==.013]
 ref=read(OUT/'refinement_h0.01_design1_high/summary.json')
 from postprocess_carrier import read_input
 from compile_tungsten_link_fea import vm
 _,_,materials=read_input(OUT/'refinement_h0.01_design1_high/model.inp')
 z=np.load(OUT/'refinement_h0.01_design1_high/fields.npz');S=z['stress'];mask=np.array([materials[int(e)]=='ALUMINUM' for e in S[:,0]])
 peak=float(max(vm(S[mask,2:8]))/1e6);change=100*(peak/fine[0]['tungsten']['centrifugal_raw_alu_vm_mpa']-1)
 save(OUT/'refinement_comparison.json',dict(fine_alu_peak_mpa=fine[0]['tungsten']['centrifugal_raw_alu_vm_mpa'],refined_alu_peak_mpa=peak,change_percent=change,refined_displacement_um=ref['max_u_m']*1e6,refined_frequencies_hz=ref['frequencies_hz'],input_sha256=ref['input_sha256']))
 assert hashlib.sha256((OUT/'refinement_h0.01_design1_high/model.inp').read_bytes()).hexdigest()==ref['input_sha256']
 lines=['# Same-mass tungsten insert comparison','','The principal advantage is reduced insert protrusion. Carrier stresses are essentially unchanged (fine-mesh differences about 0.03%, below mesh sensitivity), combined horizontal displacement decreases about 0.4-0.5%, and insert stress decreases about 23%. Low bending frequencies rise modestly while the dominant torsional mode stays nearly constant. The carrier and insert bounds remain below their provisional factor-of-two yield screening limits. These are idealized rotor results, not a validated support/retention assembly.','','## Scope and material','',
 'The 0.5, 1 and smallest qualified 2 m rotor carriers, insert radii/centers, bare masses, spin speeds and modulation timing are unchanged. Each insert is shortened to 43.61% of its steel length. CAD mass and polar inertia are preserved; conformal meshes are regenerated. This comparison repeats all 26 static/prestressed-modal cases of the steel baseline, with fixed bore and perfectly bonded inserts. It is not an assembly contact model or a new BER qualification.','',
 'Use illustrative 95% tungsten heavy alloy K1800/ET95: density 18000 kg/m3, E=330 GPa from [Elmet manufacturer data](https://www.elmettechnologies.com/wp-content/uploads/2025/09/Tungsten-Heavy-Alloy-Product-Datasheet.pdf). Poisson ratio 0.28 is an explicit modeling assumption. The table lists yield >=517 MPa; a conservative 500 MPa study specification requires procurement certification. Aluminum remains 6061-T651 with a 240 MPa screening value. Static yield screening uses a factor of two.','',
 '## Matched fine-mesh comparison','',
 '| Link (m) | Aluminum combined bound, steel / W (MPa) | Combined displacement, steel / W (um) | Tungsten combined bound (MPa), all attitudes |','|---:|---:|---:|---:|']
 for x in fine:
  a,b=x['steel'],x['tungsten'];lines.append(f"| {x['distance_m']:g} | {a['horizontal_conservative_alu_vm_mpa']:.2f} / {b['horizontal_conservative_alu_vm_mpa']:.2f} | {a['horizontal_conservative_u_um']:.2f} / {b['horizontal_conservative_u_um']:.2f} | {b['combined_conservative_tungsten_vm_mpa']:.2f} |")
 lines+=['','Combined aluminum/displacement columns assume a horizontal shaft (gravity in the rotor plane). They add conservative norms of centrifugal, maximum angular-acceleration and gravity solutions; simultaneous maxima are not implied. Insert stress column additionally bounds axial gravity. All-attitude and phase-sampled results are retained in comparison.json.','',
 '## Prestressed modes at center speed','',
 '| Link (m) | Steel first three (Hz) | Tungsten first three (Hz) | W dominant angular-coupled mode / effective inertia (%) |','|---:|:---|:---|:---|']
 for x in fine:
  a,b=x['steel_modes'],x['tungsten_modes'];j=int(np.argmax(b['angular_effective_inertia_fraction']));lines.append(f"| {x['distance_m']:g} | {', '.join(f'{v:.2f}' for v in a['frequencies_hz'][:3])} | {', '.join(f'{v:.2f}' for v in b['frequencies_hz'][:3])} | {j+1} / {100*b['angular_effective_inertia_fraction'][j]:.2f} |")
 lines+=['','The frequency ordering is a comparison, not proof of matching shapes across different meshes/materials. Component energies, insert motions and angular coupling are saved in modal_response_h*.json. Within each tungsten mesh, center/low/high modes are matched by consistent-mass MAC. Isolated-transition damping sweeps are diagnostic; they omit repeated-message resonance, shaft/bearing compliance, gyroscopic terms and motor dynamics.','',
 '## Mesh and load checks','',
 '| Link (m) | W centrifugal Al peak, coarse / fine (MPa) | W displacement change (%) | W regional combined stress change (%) |','|---:|---:|---:|---:|']
 for x in fine:
  c=next(r['tungsten'] for r in rows if r['design']==x['design'] and r['mesh']==.018);b=x['tungsten'];lines.append(f"| {x['distance_m']:g} | {c['centrifugal_raw_alu_vm_mpa']:.2f} / {b['centrifugal_raw_alu_vm_mpa']:.2f} | {100*(b['centrifugal_max_u_um']/c['centrifugal_max_u_um']-1):.3f} | {100*(b['regional_10mm_scaled_vm_mpa']/c['regional_10mm_scaled_vm_mpa']-1):.3f} |")
 lines+=['','Positive element Jacobians, exact consistent-load torque, force/moment equilibrium, cross-design static scaling and independent modal mass normalization are checked. Fine/coarse raw stress peaks can remain sensitive to local discretization; regional averages do not replace local strength assessment. The steel reference third mesh is not evidence of convergence for tungsten.','',
 f'Additional 10 mm reference mesh, scaled to the 0.5 m rotor at its upper operating speed: aluminum centrifugal peak {peak:.4f} MPa, change {change:.3f}% from the 13 mm mesh. Displacement {ref["max_u_m"]*1e6:.4f} um. This checks tungsten centrifugal convergence only; combined-load and contact convergence are not established by it.', '',
 '## Signal and windage implications','',
 '| Link (m) | Steel second harmonic (m/s2) | W second harmonic (m/s2) | Change (%) |','|---:|---:|---:|---:|']
 for f in read(OUT/'field_comparison.json')['rows']:lines.append(f"| {f['distance_m']:g} | {f['steel_a2_m_s2']:.6e} | {f['tungsten_a2_m_s2']:.6e} | {f['relative_signal_change_percent']:.3f} |")
 lines+=['','Signal uses the established steel field quadrature with only insert axial coordinates compressed, keeping every mass weight and carrier point unchanged. This isolates the finite-height effect at the same equatorial receiver position; Fourier phase refinement is checked. Leading quadrupole and polar inertia stay constant, but finite-range signal rises slightly as mass moves closer to the receiver plane. No packet simulations are repeated.','',
 'Exposed insert height and projected side area outside the carrier fall 75.4%. This is a geometric windage proxy, not a CFD prediction or a 75.4% whole-rotor power saving. Rotor rim, spokes, surface friction, flow interactions and the enclosure also matter. Acceleration torque and radial centrifugal forces do not fall at fixed mass/radius/speed. Real insert retention and material-specific fabrication remain to be modeled.','',
 '## Horizontal support drawing','',
 '![Horizontal-shaft concept](horizontal_shaft_support.png)','',
 'The rotor plane is vertical, like a wheel. Two stationary pedestals support horizontal stub shafts on opposite faces of its hub; the motor is outside the drive-side support. Bearings and pedestals lie outside the axial swept envelope. One bearing arrangement locates the rotor axially, while the opposite permits expansion. The front and edge views show rotor angles 90 degrees apart. The 0.5 m rotor dimensions illustrate the concept; supports/drive are schematic and insert retention/enclosure are omitted. SVG and PDF versions are included.','',
 '## Reproduction','',
 'Run tungsten_link_fea.py --prepare --workers 4; compile_tungsten_link_fea.py --mesh .018 and --mesh .013; tungsten_modal_response.py for both meshes; tungsten_field_comparison.py; draw_horizontal_support.py; refine_tungsten_reference.py; report_tungsten_comparison.py. Inputs, solver outputs, material metadata and hashes are retained in this folder.','']
 (OUT/'RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
 print(json.dumps([dict(distance=x['distance_m'],stress=x['tungsten']['horizontal_conservative_alu_vm_mpa'],u=x['tungsten']['horizontal_conservative_u_um']) for x in fine]),flush=True)
if __name__=='__main__':main()
