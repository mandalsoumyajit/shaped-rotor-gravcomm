"""Preliminary equal-moving-mass coherent array comparison at 2 m.

Finite-volume fields at actual rotor positions, not co-located approximations.
This is a source/noise benchmark, NOT a validated array communication result.
All rotors within a configuration send the same waveform with chosen phases.
Independent parallel subchannels and synchronization errors are not simulated.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import json, sys
from dataclasses import replace
from pathlib import Path
import numpy as np
from range_revision import ROOT,OUT,MASS,INERTIA,FC,save
from gravcomm import FiniteRotor
from gravcomm.receivers import StructuralOscillator
from gravcomm.constants import G
from paper_results import gaussian_waterfill


def coefficient(rotor,receiver,samples):
    phase=np.arange(samples)*2*np.pi/samples
    p=rotor.positions_m; mass=rotor.masses_kg; acc=[]
    for phi in phase:
        dx=np.cos(phi)*p[:,0]-np.sin(phi)*p[:,1]-receiver[0]
        dy=np.sin(phi)*p[:,0]+np.cos(phi)*p[:,1]-receiver[1]
        dz=p[:,2]-receiver[2]
        acc.append(G*np.sum(mass*dx/(dx*dx+dy*dy+dz*dz)**1.5))
    return 2*np.mean(np.asarray(acc)*np.exp(-2j*phase))


def main():
    rotor=FiniteRotor.from_npz(ROOT/'fea/results/shaped_waveforms/medium_rounded.npz')
    layouts={1:[(0.,0.,0.)],2:[(0.,y,0.) for y in (-.42,.42)],
             4:[(0.,y,z) for y in (-.34,.34) for z in (-.10,.10)],
             8:[(0.,y,z) for y in (-.30,.30) for z in (-.36,-.12,.12,.36)]}
    receiver=np.array([2.,0.,0.]);records=[]
    for count,centers in layouts.items():
        scale=(8/count)**(1/3);carrier=FC/scale
        radius=.25*scale; halfheight=.0792404*scale/2
        # Nonintersection of swept cylinders: sufficient conservative check.
        for i,a in enumerate(centers):
            for b in centers[i+1:]:
                delta=np.abs(np.array(a)-b)
                assert np.hypot(delta[0],delta[1])>2*radius or delta[2]>2*halfheight
        coefficients=[];refinement=[]
        for center in centers:
            relative=(receiver-center)/scale
            c128=scale*coefficient(rotor,relative,128)
            c256=scale*coefficient(rotor,relative,256)
            refinement.append(abs(c128-c256)/abs(c256));coefficients.append(c256)
        assert max(refinement)<1e-7
        if count==1:
            expected=2*rotor.harmonic_coefficient(1.,samples=256)
            np.testing.assert_allclose(coefficients[0],expected,rtol=1e-10,atol=0)
        aligned=float(sum(abs(c) for c in coefficients));same=float(abs(sum(coefficients)))
        # Each retuned case is a separate receiver assumption; fixed receiver
        # column is the like-for-like noise comparison across configurations.
        capacities={}
        for name,rec in [('fixed_50p3Hz_receiver',StructuralOscillator()),
                         ('retuned_receiver',replace(StructuralOscillator(),resonance_hz=carrier))]:
            band=2.;f=carrier+band*((np.arange(65536)+.5)/65536-.5)
            capacities[name]=gaussian_waterfill(rec.acceleration_noise_asd(f)**2,band,aligned**2/2)
        c=np.array(centers)
        extent=np.ptp(c,axis=0)+[2*radius,2*radius,2*halfheight]
        rec=dict(rotors=count,scale_per_rotor=scale,total_moving_mass_kg=count*MASS*scale**3,
                 individual_mass_kg=MASS*scale**3,diameter_per_rotor_m=2*radius,
                 total_inertia_kg_m2=count*INERTIA*scale**5,carrier_hz=carrier,mechanical_rpm=30*carrier,
                 rotor_centers_m=centers,receiver_m=receiver.tolist(),swept_source_extent_xyz_m=extent.tolist(),
                 individual_complex_coefficients=[[float(v.real),float(v.imag)] for v in coefficients],
                 ideal_receiver_aligned_amplitude_m_s2=aligned,equal_mechanical_phase_amplitude_m_s2=same,
                 ideal_mechanical_phase_offsets_rad=[float(-np.angle(v)/2) for v in coefficients],
                 gaussian_benchmarks_bit_s=capacities,max_phase_refinement_relative=max(refinement),
                 total_rotational_energy_j=.5*count*INERTIA*scale**5*(np.pi*carrier)**2,
                 scope='Equal moving mass and centrifugal stress scale; unequal footprint. Drives/support mass, power limits, phase errors, packet decoding and environmental noise not included. Retuned column changes receiver resonance only, retaining assumed mass/Q/readout.')
        records.append(rec);print(json.dumps(rec),flush=True)
        save(OUT/'coherent_array_benchmarks.json',records)
    np.testing.assert_allclose([r['total_moving_mass_kg'] for r in records],8*MASS,rtol=1e-12)
    np.testing.assert_allclose([r['total_rotational_energy_j'] for r in records],records[0]['total_rotational_energy_j'],rtol=1e-12)

if __name__=='__main__':main()
