"""Regenerate complete finite-body waveforms/spectra for source design variants."""
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from postprocess_carrier import read_input,quadrature,write_csv
from design_sweep import ROOT,OUT

sys.path.insert(0,str(ROOT.parent/'src'))
from gravcomm import FiniteRotor,TwoMassRotor


def model(folder,destination,name):
    nodes,elements,materials=read_input(folder/'model.inp')
    q=quadrature(nodes,elements,materials)
    p=np.array([v[3] for v in q]);m=np.array([v[4]*v[5] for v in q])
    steel=np.array([materials[v[0]]=='STEEL' for v in q])
    np.savez_compressed(destination/(name+'.npz'),positions_m=p,masses_kg=m,
                        bounding_radius_m=.250000001,is_steel=steel)
    return FiniteRotor(p,m,.250000001),steel


def main():
    destination=ROOT/'results'/'shaped_waveforms';destination.mkdir(exist_ok=True)
    sources={'original_shaped':ROOT/'results'/'shaped_rotor'/'h0.013_rpm1800'}
    sources.update({name:OUT/name/'h0.013_rpm1800' for name in
                    ('baseline_rounded','medium_rounded','wide_rounded','sphere_outer_limit','best_quadrupole')})
    fig,axes=plt.subplots(1,2,figsize=(12,4.5));summaries=[]
    for name,folder in sources.items():
        if not (folder/'metadata.json').exists():raise RuntimeError(f'Missing completed model: {folder}')
        rotor,steel=model(folder,destination,name)
        for d in (.5,1.,2.):
            phase,ac=rotor.ac_signal(d,samples=256)
            coeff=2*np.fft.rfft(ac)/len(ac)
            harmonics=[dict(harmonic=n,frequency_hz=30*n,peak_amplitude_m_s2=abs(coeff[n]),
                            real_m_s2=coeff[n].real,imag_m_s2=coeff[n].imag) for n in range(1,17)]
            write_csv(destination/f'{name}_{d:g}m_harmonics.csv',harmonics)
            a2=np.real(coeff[2]*np.exp(2j*phase))
            write_csv(destination/f'{name}_{d:g}m_waveform.csv',[
                dict(phase_rad=p,time_s=p/(60*np.pi),ac_acceleration_m_s2=a,second_harmonic_m_s2=b)
                for p,a,b in zip(phase,ac,a2)])
            summaries.append(dict(case=name,distance_m=d,mass_kg=rotor.total_mass,
                polar_inertia_kg_m2=rotor.polar_inertia,
                quadrupole_anisotropy_kg_m2=abs(rotor.quadrupole_anisotropy),
                a2_m_s2=abs(coeff[2]),a4_m_s2=abs(coeff[4]),total_rms_m_s2=np.sqrt(np.mean(ac*ac)),
                higher_even_to_a2_rms_ratio=float(np.sqrt(sum(abs(coeff[n])**2 for n in range(4,17,2)))/abs(coeff[2]))))
            if d==1.:
                axes[0].plot(phase*180/np.pi,ac/1e-11,label=name)
                axes[1].semilogy(range(2,13,2),[abs(coeff[n])/1e-11 for n in range(2,13,2)],'o-',label=name)
        if name=='original_shaped':
            # Separate aluminium/steel fields to show carrier harmonics explicitly.
            separate=[]
            phi=np.arange(256)*2*np.pi/256
            for mask,label in ((steel,'steel'),(~steel,'aluminum')):
                part=FiniteRotor(rotor.positions_m[mask],rotor.masses_kg[mask],.250000001)
                _,ac=part.ac_signal(1,samples=256)
                for n in (2,4,6):
                    coefficient=2*np.mean(ac*np.exp(-1j*n*phi))
                    separate.append(dict(component=label,harmonic=n,coefficient_real_m_s2=coefficient.real,
                                         coefficient_imag_m_s2=coefficient.imag,amplitude_m_s2=abs(coefficient)))
            write_csv(destination/'original_material_harmonics_1m.csv',separate)
            m=float(sum(rotor.masses_kg[steel]))
            point=TwoMassRotor.balanced_quadrupole(m,.2)
            (destination/'point_mass_comparison.json').write_text(json.dumps(dict(
                insert_only_point_a2_at_1m=point.harmonic_amplitude(1.),
                full_shaped_a2_at_1m=next(r['a2_m_s2'] for r in summaries if r['case']==name and r['distance_m']==1.)),indent=2))
    write_csv(destination/'summary.csv',summaries)
    axes[0].set_xlabel('Mechanical phase (degrees)');axes[0].set_ylabel('Radial AC acceleration (nGal)')
    axes[1].set_xlabel('Mechanical harmonic number');axes[1].set_ylabel('Peak amplitude (nGal)')
    for ax in axes:ax.grid(alpha=.25)
    axes[0].legend(fontsize=7);fig.suptitle('Full shaped-rotor gravity at 1 m: carrier and finite inserts included')
    fig.tight_layout();fig.savefig(destination/'waveforms_and_spectra.png',dpi=170);plt.close(fig)
    (destination/'README.md').write_text('''# Finite-volume gravitational waveforms

NPZ files contain SI quadrature positions and positive mass weights for the
complete aluminium carrier and steel inserts. Load with
`gravcomm.FiniteRotor.from_npz(path)`. Rotation is about z, receiver at (d,0,0),
and x acceleration is radial. AC traces subtract their spatial-phase mean;
the field method itself returns the full field including DC. No receiver noise,
filtering or deformation correction is included. Rotor speed is 1800 rpm.

The selected second-harmonic signal is at 60 Hz. Higher even harmonics are
reported separately; they must not be folded into its signal power. Tiny odd
harmonics can reflect mesh asymmetry. The material CSV gives COMPLEX component
coefficients: add these before taking magnitudes. The carrier is not generally
gravitationally silent. These files replace the point-insert approximation for
this shaped demonstrator; they do not redefine the separate deployment classes.

Quadrature: 4 points per quadratic tetrahedron, including curved Jacobians.
Waveforms: 256 mechanical phases per turn, distances 0.5, 1 and 2 m from axis.
The mesh-convergence tables and 128/256 phase check are in ../design_sweep/.
''')


if __name__=='__main__':main()
