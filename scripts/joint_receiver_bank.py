"""Build coupled receiver sensitivity bank; no hardware qualification is implied."""
import json
import sys
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gravcomm.receivers import StructuralOscillator
from paper_results import gaussian_waterfill
OUT = ROOT / 'results/joint_link_study_2026_09_20'


def build():
    original = StructuralOscillator()
    previous = json.loads((ROOT/'results/range_revision_2026_09_20/scenarios.json').read_text())
    receivers = []
    for resonance in (50.3, 40., 25.15):
        scale = original.resonance_hz / resonance
        for loss_multiplier in ((1.,) if scale == 1 else (1., 4.)):
            parameters = asdict(original)
            parameters.update(mass_kg=original.mass_kg*scale**3, resonance_hz=resonance,
                              quality_factor=original.quality_factor/loss_multiplier,
                              maximum_frequency_hz=original.maximum_frequency_hz/scale)
            model = StructuralOscillator(**parameters)
            identity = f'similar_f{resonance:g}_loss{loss_multiplier:g}'
            baseline = scale == 1 and loss_multiplier == 1
            receivers.append(dict(id=identity, status='published_parameter_model' if baseline else 'conditional_design_sensitivity',
                parameters=parameters, mechanical_design=dict(linear_scale=scale,
                stiffness_n_m=model.mass_kg*(2*np.pi*resonance)**2,
                mechanical_diameter_m=.0254*scale, first_higher_mode_estimate_hz=209/scale,
                structural_loss_multiplier=loss_multiplier,
                assumptions='Uniform geometric scaling, fixed material/density. Q and unchanged optical readout need verification; loss multiplier 4 is a stress test, not a confidence bound.'),
                diagnostics=dict(near_resonance_amplitude_decay_s=model.quality_factor/(np.pi*resonance),
                approximate_mechanical_linewidth_hz=resonance/model.quality_factor,
                resonance_acceleration_asd_m_s2_rt_hz=float(model.acceleration_noise_asd(resonance)),
                resonance_displacement_per_acceleration_s2=float(model.displacement_per_acceleration(resonance)))))
    cases=[]
    for distance in (.5, 1., 2.):
        base=next(c for c in previous['cases'] if c['scale']==1 and c['distance_m']==distance)
        for receiver in receivers:
            p=receiver['parameters']; physical=StructuralOscillator(**p)
            carriers=[p['resonance_hz']]

            for carrier in carriers:
                case=deepcopy(base)
                case.update(name=f'desktop_{distance:g}m__{receiver["id"]}__fc{carrier:g}',
                    carrier_hz=carrier,receiver=p.copy(),receiver_id=receiver['id'],receiver_status=receiver['status'],
                    receiver_redesigned=receiver['mechanical_design']['linear_scale']!=1,
                    carrier_model_band_hz=[p['minimum_frequency_hz'],p['maximum_frequency_hz']],
                    maximum_symmetric_baseband_sample_rate_hz=2*min(carrier-p['minimum_frequency_hz'],p['maximum_frequency_hz']-carrier))
                frequencies=carrier+2*((np.arange(65536)+.5)/65536-.5)
                noise=physical.acceleration_noise_asd(frequencies)**2
                case['gaussian_benchmark_bit_s']=gaussian_waterfill(noise,2,case['amplitude_m_s2']**2/2)
                case['wideband_energy_upper_bit_s']=float(case['amplitude_m_s2']**2/(2*np.min(noise)*np.log(2)))
                factor=(carrier/60)**2
                case.update(centrifugal_stress_factor_vs_base1800=factor,
                    illustrative_regional_mean_mpa=25.48*factor,illustrative_unconverged_peak_mpa=58.13*factor,
                    rotor_edge_to_receiver_envelope_gap_m=distance-.25-receiver['mechanical_design']['mechanical_diameter_m']/2,
                    carrier_steady_state_displacement_m=float(case['amplitude_m_s2']*physical.displacement_per_acceleration(carrier)),
                    assumptions='Stationary acquired-packet model, no environmental/control noise or hardware range qualification. Geometrically coupled receiver scaling is conditional; optical readout and loss require verification. Drive ratings fixed across links; drag coefficient is the baseline carrier-relative assumption.')
                assert case['rotor_edge_to_receiver_envelope_gap_m']>0
                cases.append(case)
    # Coupling and reproduction checks, independent of field recomputation.
    for receiver in receivers:
        p=receiver['parameters']; s=receiver['mechanical_design']['linear_scale']
        k=receiver['mechanical_design']['stiffness_n_m']
        np.testing.assert_allclose(k/(original.mass_kg*(2*np.pi*original.resonance_hz)**2),s)
        np.testing.assert_allclose(np.sqrt(k/p['mass_kg'])/(2*np.pi),p['resonance_hz'])
        f=p['resonance_hz']; physical=StructuralOscillator(**p)
        np.testing.assert_allclose(physical.acceleration_noise_asd(f)**2,
            physical.thermal_acceleration_asd(f)**2+physical.readout_acceleration_asd(f)**2)
    for d in (.5,1.,2.):
        ref=next(c for c in previous['cases'] if c['distance_m']==d and c['scale']==1)
        current=next(c for c in cases if c['distance_m']==d and c['receiver_status']=='published_parameter_model' and c['carrier_hz']==50.3)
        np.testing.assert_allclose(current['gaussian_benchmark_bit_s'],ref['gaussian_benchmark_bit_s'],rtol=1e-12)
    return dict(schema_version=1, status='deterministic_candidate_bank_not_optimized_links',
        source_sha256=previous['source_sha256'],source_field_provenance='Accepted finite-volume fields copied without renormalization from range_revision_2026_09_20/scenarios.json.',
        common_resource_envelope=dict(source_mass_kg=cases[0]['mass_kg'], source_diameter_m=.5,
            receiver_maximum_proof_mass_kg=.0248, receiver_maximum_mechanical_diameter_m=.0508,
            receiver_temperature_k=300,assumed_displacement_readout_asd_m_rt_hz=1e-13),
        evidence_url='https://www.nature.com/articles/s41598-024-68623-0',
        checks='stiffness/geometry similarity, resonance identity, thermal+readout PSD sum, all three original benchmark reproductions; positive envelope clearance',
        receivers=receivers,cases=cases)


def main():
    bank=build(); OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'receiver_bank.json').write_text(json.dumps(bank,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    rows=['| Receiver | Proof mass (g) | f0 (Hz) | Q | Decay time (s) | ASD at f0 (nGal/rtHz) |', '|---|---:|---:|---:|---:|---:|']
    for r in bank['receivers']:
        p=r['parameters']; d=r['diagnostics']
        rows.append(f'| {r["id"]} | {1000*p["mass_kg"]:.3f} | {p["resonance_hz"]:.2f} | {p["quality_factor"]:.0f} | {d["near_resonance_amplitude_decay_s"]:.1f} | {1e11*d["resonance_acceleration_asd_m_s2_rt_hz"]:.3f} |')
    note='''# Coupled receiver candidate bank

20 September 2026. This is a deterministic starting bank, not an optimized receiver or qualified hardware design.

The fixed published-parameter oscillator plus the manuscript's idealized flat readout is the evidence-anchored control. Carter et al. describe flexure geometry, proof mass and competing loss mechanisms; changing resonance alone cannot establish a new receiver. Their demonstrated complete sensor noise must not be equated with this thermal-plus-flat-readout model. [Primary source](https://www.nature.com/articles/s41598-024-68623-0).

## Coupled design hypothesis

Uniform geometric similarity at fixed material/density gives m=s^3 m0, stiffness k=s k0, f0=f00/s, and approximate higher-mode frequencies divided by s. The mechanical envelope scales with s as well. This is an elastic-similarity hypothesis, not a flexure redesign/FEA result. Gravity loading, stress, loss dilution, surface/thermoelastic losses and optical alignment require redesign checks. The bank uses s=1, 50.3/40 and 2 under one common maximum receiver envelope at every distance: 24.8 g proof mass and 50.8 mm mechanical diameter. Housing, optical cavity and isolation resources are not accounted for by that mechanical diameter.

Q cannot be derived from similarity alone. For each scaled size we evaluate the original Q and Q/4, with the latter a loss sensitivity, not a probabilistic uncertainty bound. The same 300 K and 100 fm/rtHz optical displacement readout are assumed; no cooling or readout improvement is silently introduced. No scaled receiver is labeled experimentally supported. The original receiver supplies the 50.3 Hz control. Scaled designs start at their resonance; carrier detuning and neighboring carrier refinement remain required.

## Recomputed quantities

All receiver variants recompute response, thermal/readout noise, linewidth, decay time, benchmark and carrier displacement together. The stored 2 Hz Gaussian result is only a screening reference. Linewidth is not imposed as a brick-wall information bandwidth, and decay time is not inserted as per-symbol dead time. Actual startup, saturation, causal tracking and receive filtering are still unmodeled; static displacement diagnostics are not dynamic-range qualification. Source amplitude is unchanged with carrier in this quasistatic model. Source centrifugal stress summaries are updated for carrier, while drive torque limits remain common.

'''+ '\n'.join(rows)+'''

## Interface and checks

`receiver_bank.json` has `receivers` (identity, status, StructuralOscillator parameters, mechanical design, diagnostics) and `cases`. Each case retains the range harness fields, so it can feed the existing carrier-aware decoder directly. Use receiver identity AND full parameter values in hashes. Cases cover the same five receiver/carrier combinations at each of 0.5, 1 and 2 m (15 cases). Expanded combinations should be pruned/balanced before packet Monte Carlo; this is not a request to multiply all 15 by the entire code/waveform grid.

The generator asserts stiffness/resonance consistency, PSD sum identity, positive radial envelope clearance, and reproduction of all three original 2 Hz benchmarks. Fields are copied from the accepted finite-volume cache and retain its source hash. Rerun with `python scripts/joint_receiver_bank.py`. No new field quadrature, receiver FEA, experimental validation or end-to-end coding performance is claimed.
'''
    (OUT/'RECEIVER_DESIGN.md').write_text(note,encoding='utf-8')
    print(json.dumps(dict(receivers=len(bank['receivers']),cases=len(bank['cases']),checks=bank['checks'])))

if __name__=='__main__': main()
