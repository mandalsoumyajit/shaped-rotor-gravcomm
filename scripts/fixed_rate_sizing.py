"""Deterministic fixed-noise, tunable-center Tx sizing; no decoded-link claims."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import sys,json,hashlib,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from range_revision import MASS,INERTIA,FC,duty,save
from gravcomm import FiniteRotor
from gravcomm.receivers import StructuralOscillator
from gravcomm.convolutional_candidates import metadata
from paper_results import gaussian_waterfill
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20'
SOURCE=ROOT/'fea/results/shaped_waveforms/medium_rounded.npz'

def main():
    OUT.mkdir(exist_ok=True)
    rotor=FiniteRotor.from_npz(SOURCE); rec=StructuralOscillator()
    offset=2*((np.arange(32768)+.5)/32768-.5)
    psd=rec.acceleration_noise_asd(FC+offset)**2
    h=-.001*np.log2(.001)-.999*np.log2(.999)
    specs=[]
    for d in (.5,1.,2.):
        for s in (.5,.75,1.,1.25,1.5,1.75,2.,2.5,3.,3.5,4.):
            specs.append(dict(distance_m=d,radial_scale=s,axial_scale=s,family='uniform'))
        for z in (2.,4.,8.):specs.append(dict(distance_m=d,radial_scale=1.,axial_scale=z,family='axial'))
    eligible=[c for c in specs if c['distance_m']-.25*c['radial_scale']-.0127>=.05]
    excluded=[c for c in specs if c not in eligible]
    def status(done,stage):
        save(OUT/'progress.json',dict(stage=stage,completed_fields=done,total_fields=len(eligible),geometry_excluded=len(excluded),packet_simulations=0))
        (OUT/'PROGRESS.md').write_text(f'# Fixed-rate sizing progress\n\nPlan saved with fixed receiver noise/bandwidth and tunable center.\n\n{stage}: {done}/{len(eligible)} finite-volume cases; {len(excluded)} geometry exclusions.\n\nPacket simulations: not started. Deterministic feasibility only.\n',encoding='utf-8')
    status(0,'running_fields')
    def field(c):
        s,z,d=c['radial_scale'],c['axial_scale'],c['distance_m']
        name=f"{c['family']}_d{d:g}_r{s:g}_z{z:g}"
        path=OUT/(name+'.json')
        if path.exists():return json.loads(path.read_text())
        if s==z:
            a64=s*rotor.harmonic_amplitude(d/s,samples=64)
            a128=s*rotor.harmonic_amplitude(d/s,samples=128)
        else:
            scaled=FiniteRotor(rotor.positions_m*np.array([s,s,z]),rotor.masses_kg*s*s*z,rotor.bounding_radius_m*s)
            a64=scaled.harmonic_amplitude(d,samples=64);a128=scaled.harmonic_amplitude(d,samples=128)
        rel=abs(a64/a128-1)
        if rel>1e-7:
            a256=(s*rotor.harmonic_amplitude(d/s,samples=256) if s==z else scaled.harmonic_amplitude(d,samples=256))
            rel=abs(a128/a256-1);a128=a256
        assert rel<1e-7
        cap=gaussian_waterfill(psd,2.,a128*a128/2)
        capfine=gaussian_waterfill(rec.acceleration_noise_asd(FC+2*((np.arange(65536)+.5)/65536-.5))**2,2.,a128*a128/2)
        assert abs(cap/capfine-1)<1e-4
        row=dict(c,id=name,mass_kg=MASS*s*s*z,inertia_kg_m2=INERTIA*s**4*z,diameter_m=.5*s,axial_extent_m=float(np.ptp(rotor.positions_m[:,2]))*z,carrier_hz=FC/s if s==z else FC,amplitude_m_s2=a128,phase_refinement_relative=rel,gaussian_band_benchmark_bit_s=capfine,ber_distortion_information_threshold_bit_s=1-h,passes_relaxed_information_screen=bool(capfine>=1-h),model='Fixed PSD versus offset, tunable center, 2 Hz band; no physical tuning/noise verification.')
        save(path,row);return row
    rows=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(field,c) for c in eligible]
        for future in as_completed(futures):
            row=future.result();rows.append(row);status(len(rows),'running_fields')
            print(json.dumps(dict(done=len(rows),total=len(eligible),id=row['id'],capacity=row['gaussian_band_benchmark_bit_s'])),flush=True)
    rows.sort(key=lambda c:(c['distance_m'],c['mass_kg'],c['family']))
    status(len(rows),'mechanical_screen')
    templates=[]
    for k in (64,128,256):
        payload=k-32
        for scheme in ('uncoded','ldpc','tb_1/2','tb_2/3','tb_3/4'):
            n=k if scheme=='uncoded' else 2*k if scheme=='ldpc' else metadata(k,scheme[3:])['n']
            symbols=3*((n+7)//8)+18
            ts=payload/symbols
            for q in (.1,.25):
                # Unit inertia, zero drag gives scalable exact valid-byte RMS.
                unit=duty(ts,.4,q,inertia=1.,fc=FC,drag=0.)
                templates.append(dict(scheme=scheme,information_bits=k,payload_bits=payload,coded_bits=n,padding_bits=(-n)%8,symbol_s=ts,packet_airtime_s=float(payload),q=q,unit_peak=unit['peak_torque_nm'],unit_rms=unit['worst_valid_message_rms_nm']))
    mechanics=[]
    for row in rows:
        if not row['passes_relaxed_information_screen']:continue
        for t in templates:
            spacing=t['q']/(.6*t['symbol_s']);fmax=row['carrier_hz']+3*spacing
            peak=row['inertia_kg_m2']*t['unit_peak'];rms=row['inertia_kg_m2']*t['unit_rms']
            mechanics.append(dict(source_id=row['id'],**t,tone_span_hz=6*spacing,tones_within_band=bool(3*spacing<=1),peak_inertial_torque_nm=peak,rms_inertial_torque_nm=rms,conservative_peak_inertial_power_w=peak*np.pi*fmax,maximum_rpm=30*fmax,stress_scale_vs_reference1800=(row['radial_scale']*fmax/60)**2,within_reference_stress_scale=bool(row['radial_scale']*fmax<=60),nominal_stored_energy_j=.5*row['inertia_kg_m2']*(np.pi*row['carrier_hz'])**2))
    manifest=dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),receiver=dict(reference_center_hz=FC,bandwidth_hz=2,noise_profile='StructuralOscillator default evaluated at 50.3 Hz + offset; translated unchanged with carrier'),target_payload_bit_s=1,postdecode_ber_target=.001,latency_limit_s=300,fields=rows,geometry_excluded=excluded,mechanics=mechanics,limitations='No packet simulations, occupied-bandwidth filter, CRC/acquisition implementation, drag or new FEA. Power is conservative torque-times-maximum-speed bound. RMS is sustained valid-byte maximum and does not include startup. Capacity relaxation and source-distortion bound are not finite-packet achievability.')
    save(OUT/'deterministic_screen.json',manifest)
    lines=['# First deterministic sizing results','','Receiver noise versus offset and 2 Hz bandwidth fixed; center tracks carrier. No coded-link or hardware feasibility claims.','', '| Distance (m) | Family | Diameter (m) | Axial scale | Mass (kg) | Carrier (Hz) | Gaussian benchmark (bit/s) | Pass relaxed information screen |','|---:|---|---:|---:|---:|---:|---:|---|']
    for r in rows:lines.append(f"| {r['distance_m']} | {r['family']} | {r['diameter_m']:.3g} | {r['axial_scale']:g} | {r['mass_kg']:.2f} | {r['carrier_hz']:.3g} | {r['gaussian_band_benchmark_bit_s']:.4g} | {r['passes_relaxed_information_screen']} |")
    lines+=['','The information screen uses 1-H2(0.001), about 0.9886 bit/s, for a uniform binary source with allowed reconstruction errors. Passing only retains a source for further study. Failing excludes that idealized band-limited asymptotic model at the target distortion, not arbitrary source geometries.','',f'{len(mechanics)} framed waveform/mechanics combinations were evaluated. See deterministic_screen.json. Tone containment does not establish occupied-bandwidth compliance. Inertial torque excludes drag; no drive was selected.','', 'Next: implement fixed-band receive filtering and common framing/error accounting, prune mechanically unsuitable candidates, then run decoded-packet screens near the mass boundaries.']
    (OUT/'DETERMINISTIC_RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    status(len(rows),'deterministic_screen_complete')
    print('DETERMINISTIC_SCREEN_COMPLETE',flush=True)
if __name__=='__main__':main()
