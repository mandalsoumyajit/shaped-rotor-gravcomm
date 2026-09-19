"""Independent cross-check of manuscript numbers against retained result files.

Does not modify the simulation or regenerate its Monte Carlo outcomes. Writes
an audit record with closed-form mechanical/noise checks and data provenance.
"""
import csv
import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with (ROOT/path).open() as stream:
        return list(csv.DictReader(stream))


def main():
    result = {}
    design = {r['case']: r for r in rows('fea/results/design_sweep/filleted_screen.csv')
              if float(r['mesh_size_m']) == .013}
    selected = design['medium_rounded']
    I = float(selected['aluminum_polar_inertia_kg_m2']) + float(selected['steel_polar_inertia_kg_m2'])
    mass = float(selected['aluminum_mass_kg']) + float(selected['steel_mass_kg'])
    result['cad'] = dict(mass_kg=mass, inertia_kg_m2=I,
                         quadrupole_over_inertia=float(selected['quadrupole_x2_minus_y2_kg_m2'])/I)
    result['signal_gains_percent'] = {name: 100*(float(r['a2_at_1m_m_s2']) /
        float(design['baseline_rounded']['a2_at_1m_m_s2'])-1) for name, r in design.items()}
    field = [r for r in rows('fea/results/shaped_waveforms/summary.csv') if r['case']=='medium_rounded']
    result['field_checks'] = []
    for r in field:
        d = float(r['distance_m'])
        trace = np.genfromtxt(ROOT/f'fea/results/shaped_waveforms/medium_rounded_{d:g}m_waveform.csv', delimiter=',', names=True)
        ac = trace['ac_acceleration_m_s2']; phi = trace['phase_rad']
        coeff = 2*np.mean(ac*np.exp(-2j*phi))
        rms = np.sqrt(np.mean(ac**2))
        np.testing.assert_allclose([abs(coeff), rms], [float(r['a2_m_s2']), float(r['total_rms_m_s2'])], rtol=1e-11)
        np.testing.assert_allclose(abs(coeff), float(selected[f'a2_at_{d:g}m_m_s2']), rtol=1e-9)
        result['field_checks'].append(dict(distance_m=d,a2_nGal=abs(coeff)/1e-11,rms_nGal=rms/1e-11))
    averages = rows('fea/results/design_sweep/averaged_stress_regions.csv')
    result['stress_checks'] = {}
    for name in design:
        means = {h: max(float(r['mean_vm_mpa']) for r in averages if r['case']==name
                       and float(r['mesh'])==h and float(r['window_radius_mm'])==10)
                 for h in (.018,.013)}
        result['stress_checks'][name] = dict(fine_mean_mpa=means[.013],
             change_pct=100*(means[.013]/means[.018]-1))
    cf = json.loads((ROOT/'results/desktop_cf_fsk/summary.json').read_text())
    ts = cf['symbol_s']; td = .8*ts; tr = .2*ts
    peak = (15/8)*math.pi*I*6/td/tr
    np.testing.assert_allclose(peak,cf['peak_inertial_torque_nm'],rtol=1e-10)
    result['operating_point'] = dict(payload_bit_s=8/(3*ts),peak_inertial_torque_nm=peak,
        low_rpm=30*(50.3-3/td),high_rpm=30*(50.3+3/td),
        centre_energy_j=I*(math.pi*50.3)**2/2,
        max_energy_j=I*(math.pi*(50.3+3/td))**2/2,
        study_energy_j=I*(60*math.pi)**2/2,
        scaled_centrifugal_mean_mpa=result['stress_checks']['medium_rounded']['fine_mean_mpa']*((50.3+3/td)/60)**2)
    result['rate_torque_table'] = [dict(symbol_s=t,rate_bit_s=8/(3*t),spacing_hz=1/(.8*t),
        peak_inertial_torque_nm=(15/8)*math.pi*I*6/(.16*t*t)) for t in (16,8,4)]
    f = np.array(cf['tones_hz']); m=.0031; Q=637000; T=300; f0=50.3; sx=1e-13
    thermal2=8*np.pi*1.380649e-23*T*f0*f0/(m*f*Q)
    readout2=sx*sx*(2*np.pi)**4*((f0*f0-f*f)**2+f0**4/Q**2)
    noise=np.sqrt(thermal2+readout2)
    np.testing.assert_allclose(noise,cf['receiver']['noise_asd_at_tones'],rtol=1e-12)
    result['receiver'] = dict(tone_asd_si=noise.tolist(),amplitude_decay_s=Q/(np.pi*f0),
                              broadband_0p1_ng_in_nGal=.1*9.80665e-9/1e-11)
    # Independently check the retained continuous-window covariance and its
    # frequency-grid convergence. Quantify the sampled-window approximation.
    def covariance(points, sampled=False):
        nu = np.linspace(-5/td, 5/td, points)
        df = nu[None, :] - np.arange(-3, 4)[:, None]/td
        if sampled:
            count = round(100*td)
            W = np.sinc(df*count/100)/np.sinc(df/100)*np.exp(1j*np.pi*df*(count-1)/100)
        else:
            W = np.sinc(df*td)*np.exp(1j*np.pi*df*td)
        freq = f0+nu
        psd = (8*np.pi*1.380649e-23*T*f0*f0/(m*freq*Q)
               +sx*sx*(2*np.pi)**4*((f0*f0-freq*freq)**2+f0**4/Q**2))
        weights = np.full(points, nu[1]-nu[0]); weights[[0,-1]] *= .5
        return (W*(2*psd*weights))@W.conj().T
    saved_cov = np.load(ROOT/'results/desktop_cf_fsk/detector_inputs.npz')['covariance']
    cov = covariance(12001); refined = covariance(24001)
    np.testing.assert_allclose(cov, saved_cov, rtol=1e-10, atol=1e-34)
    scale = np.linalg.norm(refined)
    grid_error = np.linalg.norm(cov-refined)/scale
    sampled_error = np.linalg.norm(covariance(24001, True)-refined)/scale
    assert grid_error < 1e-6
    assert sampled_error < .01
    result['covariance'] = dict(relative_grid_change=float(grid_error),
        relative_sampled_window_difference=float(sampled_error),
        minimum_eigenvalue=float(np.linalg.eigvalsh(cov)[0]))
    errors=np.array(cf['monte_carlo']['errors_per_context']); n=100000
    upper=float(np.max(beta.ppf(1-.05/49,errors+1,n-errors)))
    np.testing.assert_allclose(upper,cf['monte_carlo']['worst_simultaneous_95pct_ser_upper'],rtol=1e-12)
    diag=cf['diagnostics']
    residual=diag['motoring_work_j']-diag['absorbed_work_j']-diag['drag_heat_j']-diag['kinetic_energy_change_j']
    assert abs(residual)<1e-6
    result['statistics_and_work'] = dict(total_errors=int(errors.sum()),decisions=49*n,
        mean_ser=float(errors.sum()/(49*n)),simultaneous_upper=upper,work_residual_j=residual)
    out=ROOT/'results/manuscript_audit';out.mkdir(exist_ok=True)
    (out/'numerical_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
