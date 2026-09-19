"""Trace the original capacity claims and the scaled CF-FSK examples separately."""
from dataclasses import replace
import json
from pathlib import Path
import sys
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from gravcomm.designs import DEPLOYMENTS
from gravcomm.link_budget import capacity_bps
from gravcomm.receivers import StructuralOscillator
from paper_results import gaussian_waterfill
fixed=next(x for x in DEPLOYMENTS if x.key=='fixed')
rotor=fixed.rotor();B=fixed.default_bandwidth
rows=[]
for d,old_asd in ((5.,1e-9),(14.,1e-11)):
    peak=rotor.peak_ac_amplitude(d);a2=rotor.harmonic_amplitude(d)
    f=fixed.signal_frequency+B*((np.arange(32768)+.5)/32768-.5)
    rec=replace(StructuralOscillator(),resonance_hz=fixed.signal_frequency)
    rows.append(dict(distance_m=d,mass_kg=fixed.total_mass,carrier_hz=fixed.signal_frequency,
        bandwidth_hz=B,old_flat_asd=old_asd,old_peak_capacity=float(capacity_bps(peak,old_asd,B)),
        harmonic_corrected_same_noise_capacity=float(capacity_bps(a2,old_asd,B)),
        harmonic_corrected_Carter_flat_reference_capacity=float(capacity_bps(a2,9.80665e-10,B)),
        retuned_structural_receiver_capacity=gaussian_waterfill(rec.acceleration_noise_asd(f)**2,B,a2*a2/2)))
scaled=json.loads((root/'results/scaled_cf_fsk/summary.json').read_text())
for item in scaled:
    B=2*item['receiver_band']['receiver_halfband_hz']
    f=item['carrier_hz']+B*((np.arange(32768)+.5)/32768-.5)
    rec=StructuralOscillator(**item['receiver'])
    item['gaussian_capacity_bit_s']=gaussian_waterfill(rec.acceleration_noise_asd(f)**2,B,item['amplitude_m_s2']**2/2)
record=dict(original_fixed=rows,scaled=[{key:item[key] for key in
    ('scale','mass_kg','distance_m','carrier_hz','payload_bit_s','gaussian_capacity_bit_s')} for item in scaled])
(root/'results/scaled_cf_fsk/rate_comparison.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
