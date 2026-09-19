"""Read-only cross-check of current manuscript results and their input records."""
import csv,json,math
from pathlib import Path
import numpy as np
from gravcomm.packet_confidence import packet_ber_interval
from gravcomm.receivers import CARTER_OSCILLATOR

ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads((ROOT/path).read_text())
checks=ROOT/'results/communications_checks'
for path in sorted(checks.glob('fs*.json')):
    r=json.loads(path.read_text());errors=np.array(r['bit_errors_per_packet'])
    assert r['complete'] and len(errors)==r['packets']
    assert r['bits']==64*len(errors) and errors.sum()==r['bit_errors']
    assert np.count_nonzero(errors)==r['errored_packets']
    assert r['adequate_error_count']==bool(errors.sum()>=100 and np.count_nonzero(errors)>=30)
    np.testing.assert_allclose(packet_ber_interval(errors),r['interval95'],rtol=1e-10)
    np.testing.assert_allclose(errors.mean()/64,r['observed_ber'],rtol=1e-13)
selected=read('results/communications_checks/fs8_mem6_noise0.json')
assert selected['bit_errors']==144 and selected['bits']==2912000
assert selected['interval95'][1]<1e-3
spacing=list((ROOT/'results/spacing_sweep').glob('*_fs8_m3.json'));assert len(spacing)==9
for path in spacing:
    r=json.loads(path.read_text());errors=np.array(r['bit_errors_per_packet'])
    assert r['complete'] and errors.sum()==r['bit_errors'] and len(errors)==r['packets']
    np.testing.assert_allclose(packet_ber_interval(errors),r['packet_aware_anytime95_ber_interval'],rtol=1e-10)
designs=list(csv.DictReader((ROOT/'fea/results/design_sweep/filleted_screen.csv').open()))
d=next(r for r in designs if r['case']=='medium_rounded' and float(r['mesh_size_m'])==.013)
inertia=float(d['aluminum_polar_inertia_kg_m2'])+float(d['steel_polar_inertia_kg_m2'])
mass=float(d['aluminum_mass_kg'])+float(d['steel_mass_kg'])
np.testing.assert_allclose([mass,inertia],[7.9135094916,.316071309585],rtol=1e-9)
field=list(csv.DictReader((ROOT/'fea/results/shaped_waveforms/summary.csv').open()))
for r in field:
    if r['case']!='medium_rounded':continue
    distance=float(r['distance_m'])
    trace=np.genfromtxt(ROOT/f'fea/results/shaped_waveforms/medium_rounded_{distance:g}m_waveform.csv',delimiter=',',names=True)
    ac=trace['ac_acceleration_m_s2'];phase=trace['phase_rad']
    np.testing.assert_allclose(abs(2*np.mean(ac*np.exp(-2j*phase))),float(r['a2_m_s2']),rtol=1e-11)
    np.testing.assert_allclose(np.sqrt(np.mean(ac**2)),float(r['total_rms_m_s2']),rtol=1e-11)
    if distance==.5:
        np.testing.assert_allclose(float(r['a2_m_s2']),5.180259430854158e-10,rtol=1e-12)
energy=read('results/communications_checks/revised_energy.json');e=energy['diagnostics']
residual=e['motoring_work_j']-e['absorbed_work_j']-e['drag_heat_j']-e['kinetic_energy_change_j']
assert abs(residual)<2e-9
np.testing.assert_allclose(.5*inertia*(50.55*np.pi)**2,energy['stored_energy_at_max_speed_j'],rtol=1e-10)
peak=15/8*np.pi*inertia*6*(.1/1.2)/.8
assert round(peak,3)==1.164
np.testing.assert_allclose(CARTER_OSCILLATOR.acceleration_noise_asd(50.3),5.1493883736e-11,rtol=1e-9)
print('Verified current BER counts/intervals, nine-point sweep, CAD moments, waveform harmonics, energy, torque and receiver normalization.')
