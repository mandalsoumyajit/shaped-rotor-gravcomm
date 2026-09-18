"""Build the final comparison report and plots from completed result files."""
import csv
import json
from pathlib import Path
import numpy as np
import gmsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from design_sweep import ROOT,OUT


def curves(path):
    gmsh.initialize();gmsh.option.setNumber('General.Terminal',0)
    try:
        gmsh.open(str(path));result=[]
        for dim,tag in gmsh.model.getEntities(1):
            lo,hi=gmsh.model.getParametrizationBounds(dim,tag)
            points=np.asarray(gmsh.model.getValue(dim,tag,np.linspace(lo[0],hi[0],80))).reshape(-1,3)
            if np.max(abs(points[:,2]-.01))<1e-7:result.append(points[:,:2]*1000)
        return result
    finally:gmsh.finalize()


def main():
    rounded=list(csv.DictReader((OUT/'filleted_screen.csv').open()))
    assert len(rounded)==6
    fine={r['case']:r for r in rounded if float(r['mesh_size_m'])==.013}
    wave=list(csv.DictReader((ROOT/'results'/'shaped_waveforms'/'summary.csv').open()))
    signal={r['case']:float(r['a2_m_s2']) for r in wave if float(r['distance_m'])==1.}
    # Independent waveform implementation and 256 phases must reproduce the
    # screening integrator's 128-phase second harmonic on the same geometry.
    for name,r in fine.items():
        np.testing.assert_allclose(signal[name],float(r['a2_at_1m_m_s2']),rtol=1e-9,atol=0)
    base=signal['baseline_rounded']
    mass_errors=[100*(float(r['quadrature_mass_kg'])/7.913509491619677-1) for r in fine.values()]
    averages=list(csv.DictReader((OUT/'averaged_stress_regions.csv').open()))
    convergence=list(csv.DictReader((OUT/'averaged_stress_convergence.csv').open()))
    worst10=max(abs(float(r['mean_vm_change_pct'])) for r in convergence if float(r['window_radius_mm'])==10)
    worst5=max(abs(float(r['mean_vm_change_pct'])) for r in convergence if float(r['window_radius_mm'])==5)
    def average(name,h,radius):
        selected=[v for v in averages if v['case']==name and float(v['mesh'])==h and float(v['window_radius_mm'])==radius]
        return max(float(v['mean_vm_mpa']) for v in selected)
    mode=list(csv.DictReader((ROOT/'results'/'stress_modes'/'mode_summary.csv').open()))
    stress=json.loads((ROOT/'results'/'stress_modes'/'h0.013_rpm1800.json').read_text())
    shape=list(csv.DictReader((OUT/'matched_shape_pair.csv').open()))
    shape_ratio=float(shape[1]['a2_at_1m_m_s2'])/float(shape[0]['a2_at_1m_m_s2'])-1
    point_comparison=json.loads((ROOT/'results'/'shaped_waveforms'/'point_mass_comparison.json').read_text())
    point_change=100*(point_comparison['full_shaped_a2_at_1m']/point_comparison['insert_only_point_a2_at_1m']-1)
    fig,axes=plt.subplots(1,2,figsize=(11,4.7))
    for name,r in fine.items():
        gain=100*(signal[name]/base-1)
        label={'baseline_rounded':'Baseline','medium_rounded':'30 mm spokes','wide_rounded':'40 mm spokes'}[name]
        axes[0].scatter(gain,float(r['first_frequency_hz']),s=60)
        axes[0].annotate(label,(gain,float(r['first_frequency_hz'])),xytext=(-4,5) if gain>20 else (4,5),ha='right' if gain>20 else 'left',textcoords='offset points',fontsize=8)
        axes[1].scatter(gain,average(name,.013,10),s=60)
        axes[1].annotate(label,(gain,average(name,.013,10)),xytext=(-4,5) if gain>20 else (4,5),ha='right' if gain>20 else 'left',textcoords='offset points',fontsize=8)
    for ax in axes:ax.set_xlabel('Signal amplitude gain at 1 m (%)');ax.grid(alpha=.3);ax.margins(.3)
    axes[0].set_ylabel('First prestressed mode (Hz)');axes[1].set_ylabel('Largest region-average von Mises (MPa; 10 mm radius)')
    fig.suptitle('Fixed 7.91 kg mass: three carrier designs with 3 mm junction fillets')
    fig.tight_layout();fig.savefig(OUT/'design_tradeoffs.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    rounded_curves=curves(OUT/'medium_rounded'/'h0.013_rpm1800'/'geometry.brep')
    sharp_curves=curves(OUT/'medium_spokes'/'h0.013_rpm1800'/'geometry.brep')
    for p in rounded_curves:
        for ax in axes:ax.plot(p[:,0],p[:,1],color='steelblue',lw=1)
    for p in sharp_curves:axes[1].plot(p[:,0],p[:,1],'r--',lw=.7,alpha=.65)
    axes[0].set_title('30 mm spokes, 20 mm rim, 3 mm fillets')
    axes[1].set_title('Hub junction detail: rounded / original dashed')
    axes[1].set_xlim(17,43);axes[1].set_ylim(6,30)
    for ax in axes:ax.set_aspect('equal');ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)')
    fig.tight_layout();fig.savefig(OUT/'rounded_geometry.png',dpi=170);plt.close(fig)
    table=[]
    for name,r in fine.items():
        coarse=next(v for v in rounded if v['case']==name and float(v['mesh_size_m'])==.018)
        change=100*(average(name,.013,10)/average(name,.018,10)-1)
        table.append(f"| {name} | {1000*float(r['spoke_width_m']):.0f} | {1000*float(r['rim_width_m']):.0f} | "
                     f"{float(r['steel_mass_kg'])/2:.3f} | {signal[name]/1e-11:.4f} | {100*(signal[name]/base-1):+.2f}% | "
                     f"{float(r['first_frequency_hz']):.2f} | {average(name,.013,5):.2f} | {average(name,.013,10):.2f} | {change:+.2f}% |")
    report='''# Matched stresses, carrier modes and fixed-mass design study

## Scope and assumptions

Objective chosen by the user: maximize transmitted gravitational signal at fixed
total rotor mass. The study holds CAD mass at 7.91350949 kg, outer radius at
at most 250 mm, axial envelope at at most 90 mm, carrier thickness at 20 mm,
speed at 1800 rpm, bore diameter at 12 mm and nominal boss wall at 15 mm.
The signal objective is radial second-harmonic amplitude at 1 m from the axis;
0.5 m and 2 m are also recomputed. These are source fields, not receiver outputs.
The hub is fixed and inserts perfectly bonded. Contact/fasteners, shaft/bearings,
fatigue, damping, gyroscopic response and practical motor drive limits still
require separate work. This is a sampled trade study, not a global optimum.

## Location-matched analytical checks on the original shaped model

- Volume-weighted hoop stress in four mid-arc rim sectors (radius 220-245 mm;
  angle within 15 degrees of 45, 135, 225, 315 degrees) is 5.974-5.988 MPa.
  The existing thin-ring estimate using outer radius is 5.996 MPa. This validates
  the scale of nominal mid-arc hoop stress, not every point on the rim: sampled
  local stresses span about 3.2-8.7 MPa. A free-ring estimate evaluated at mean
  radius would differ; this loaded spoked ring is not an isolated thin ring.
- Each insert needs 14.1950 kN from m*omega^2*r. Consistent steel-element nodal
  force recovery gives 14.1945 and 14.1946 kN, within 0.004%. It converges from
  about 0.037% error on the coarse mesh. The recovery integrates B-transpose
  times elastic stress and subtracts consistent centrifugal nodal loading.
- A geometry-matched projected area is insert diameter times engagement
  thickness: 60*20=1200 mm^2. The associated NOMINAL load/area is 11.829 MPa.
  The older 20.279 MPa estimate used rim width*thickness=700 mm^2, a different
  gross-area proxy. Neither number is a local bonded/contact peak stress.
- Extrapolating elemental stresses onto the curved interface gave force errors
  of several percent and lacked satisfactory convergence. Those traction
  estimates are retained in JSON as rejected diagnostics, not evidence of
  validated local contact stress. The consistent-force method avoids that
  extrapolation, but still does not turn a bonded model into a contact model.

## Carrier modes

Modes below are for the original fixed-bore carrier at 1800 rpm, including
centrifugal prestress. They are not shaft-whirl critical speeds. Shapes are
mass-normalized and their amplitudes are arbitrary, not predicted vibration.

| Mode | Hz | Identification |
|---|---:|---|
| 1 | 78.58 | Predominantly axial bending: opposite inserts move in opposite axial directions |
| 2 | 93.73 | In-plane torsional motion of carrier relative to fixed hub |
| 3 | 102.98 | Predominantly axial bending: inserts move together |
| 4 | 123.94 | Axial bending dominated by the unweighted transverse spokes/rim |
| 5 | 232.66 | Higher axial bending with opposing weighted/unweighted sectors |
| 6 | 267.84 | Higher axial twisting/saddle-like carrier deformation |
| 7 | 467.79 | Higher in-plane bending with inserts moving mainly together tangentially |
| 8 | 569.40 | Higher in-plane rim/spoke deformation |

The second mode accounts for about 99.86% of effective polar inertia for a
distributed angular-acceleration excitation in this fixed-hub model. It is
therefore more directly relevant to speed modulation than the lowest axial
mode. This is a coupling calculation, not a forced-response or BER prediction.
Quadrature reconstruction of all eight generalized modal masses agrees with
the solver normalization to about 1e-7. See `results/stress_modes/carrier_modes.png`.

## Corners and local mesh limitations

Explicit 3 mm in-plane junction fillets were added. Their small mass change is
compensated by adjusting insert height, retaining fixed total CAD mass. The
through-bore and bonded pocket axial edges are not fully detailed hardware.
An initial curved quadratic mesh contained inverted elements and was rejected.
Final rounded cases use quadratic displacement elements with straight midside
geometry, local junction refinement and resolved circular boundaries. Every
accepted mesh passes positive-Jacobian checks. All bore nodes are constrained
by CAD surface membership, not by testing sampled node radius.

The earlier 38.9 MPa peak is NOT a converged stress validation: introducing local
refinement changed the peak substantially. Comparison with it must not be used
to claim that rounding increases stress. Global mesh-size agreement had missed
local stress sensitivity. Raw peaks are retained only as diagnostics in
`results/design_sweep/stress_regions.csv`.

The primary comparison uses volume-weighted average von Mises stress in fixed
5 mm and 10 mm radius circular neighbourhoods, spanning the full plate thickness
and including aluminium only. Windows are centered on nominal hub/spoke,
inner insert-boss/spoke and unweighted rim/spoke junctions, plus mid-arc rim
locations. Four symmetry-related windows per group are evaluated separately;
the table reports the largest regional average, not a pointwise maximum.
Centers remain fixed in physical coordinates between meshes. The CSV also
reports actual sampled volumes, sample counts and von Mises stress of the
averaged tensor. The latter can be lower because differently oriented stresses
cancel; it must not be substituted silently for the average of von Mises stress.
See `results/design_sweep/averaged_stress_regions.csv` for all windows.
Regional averages improve robustness for design comparisons but can conceal
local yielding or fatigue initiation; they are not alone a strength criterion.

## Fixed-mass rounded-design comparisons

All dimensions below are mm, insert mass is per insert. Signal includes the
complete finite aluminium and steel geometry. CAD mass is exact by construction;
quadrature/mesh mass error is separately recorded in the raw CSVs.

| Design | Spoke | Rim | Insert kg | A2 at 1 m, nGal | Gain | First mode Hz | Mean VM, 5 mm MPa | Mean VM, 10 mm MPa | 10 mm mean change, coarse to fine |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(table)+f'''

The controlling windows in this table are at the weighted hub/spoke junctions;
their 10 mm averages change by less than 0.2%. This does NOT imply that every
window has converged equally well. Across all individual windows, the largest
change is {worst10:.2f}% for 10 mm windows and {worst5:.2f}% for 5 mm windows.
Some small windows outside the locally refined junctions have sparse quadrature
sampling and appreciable sampled-volume changes. The 10 mm averages are the
primary comparison; `averaged_stress_convergence.csv` records each paired window
and its sampled-volume change so these limitations remain visible.

The final faceted meshes have mass errors from {min(mass_errors):+.3f}% to
{max(mass_errors):+.3f}% relative to exact CAD mass. Gravity weights are not
renormalized. This discretization error should be kept in mind when comparing
small signal changes; the larger design gains are well above that scale.

The 30 and 40 mm spoke alternatives use 70 mm diameter cylinders at 200 mm
radius, versus 60 mm cylinders in the original design. Reducing rim width
frees material budget for the signal-bearing steel; widening spokes can regain
stiffness. The 20 mm spoke, 20 mm rim sharp candidate gives about 32% more signal
but its first mode falls to about 68 Hz and its raw peak rises substantially.
It is not automatically the preferred design.

Sphere/cylinder isolation check: with the same 180 mm center radius, 20 mm spokes,
20 mm rim and total mass, spherical inserts change A2 at 1 m by {100*shape_ratio:+.3f}%
relative to the cylinder control. This matched comparison shows essentially
the same signal, with finite-size corrections and carrier/pocket geometry
accounted for. Their larger radial footprint also limits how far outward
they fit with a 15 mm boss wall. A separately tested sphere near its radial
limit has center radius about 192.35 mm, not 200 mm.

For further optimization, vary continuous placement, cylinder aspect ratio,
spoke taper/thickness, rim width and fillet radius together. Preserve the full
finite-volume signal objective and include constraints on locally converged
stress, deformation, modal response and attachment design. Maximizing nominal
quadrupole alone pushes toward a light flexible carrier and is insufficient.

## Redone gravitational waveforms

`gravcomm.FiniteRotor` integrates the rotating full-body Newtonian field using
positive element quadrature masses. The carrier is not assumed silent. Its
fourfold structure contributes higher harmonics; bosses/pockets also affect
the second harmonic. The rotor's polar inertia and signal quadrupole anisotropy
are distinct quantities. Complex material harmonic coefficients must be summed
before taking magnitude.

For the original shaped geometry at 1 m, the full-body second-harmonic amplitude
is {point_comparison['full_shaped_a2_at_1m']/1e-11:.4f} nGal, compared with
{point_comparison['insert_only_point_a2_at_1m']/1e-11:.4f} nGal for the two-point
insert approximation ({point_change:+.2f}%). The aluminium contribution at the
second harmonic partly opposes the steel contribution; its fourth harmonic
adds to the steel fourth harmonic. The correction is small in this example,
but it is now calculated rather than assumed away.

`results/shaped_waveforms/` contains full phase traces and spectra at 0.5, 1 and
2 m, 256 mechanical phases, reusable NPZ source models, the separate material
harmonics and a comparison plot. Rotation is about z and the receiver lies on x.
These waveforms are for undeformed geometry at constant 1800 rpm; the same
FiniteRotor can be passed to the actuator trace to evaluate its achieved phase.
They do not include a receiver transfer function, noise or a demonstrated bit rate.

Validation includes the point-mass limit, phase/sign convention, axisymmetric
ring cancellation, actuator compatibility, C3D10 affine-strain reconstruction,
stress-component ordering, modal parsing, mass checks, mechanical mesh comparisons,
and a 128/256 phase comparison for the volume-integrated signal. The existing
model suite passes 51 tests and the FEA parsing/interpolation/averaging suite passes five.

Sources for implementation conventions: [CalculiX manual](https://www.dhondt.de/ccx_2.21.pdf),
[CalculiX Gauss order](https://github.com/Dhondtguido/CalculiX/blob/master/src/gauss.f),
[Gmsh API](https://gmsh.info/doc/texinfo/).
'''
    (ROOT/'STRESS_MODES_AND_DESIGN.md').write_text(report)
    print('\n'.join(table))


if __name__=='__main__':main()
