"""Separate alphabet and spectral-allocation contributions to the desktop gap."""
import json
from pathlib import Path
import sys
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from gravcomm.receivers import CARTER_OSCILLATOR
from paper_results import gaussian_waterfill
b=json.loads((root/'results/paper_revision/gaussian_benchmark.json').read_text())
B=b['bandwidth_hz'];P=b['rows'][0]['a2_m_s2']**2/2
f=50.3+B*((np.arange(32768)+.5)/32768-.5)
noise=CARTER_OSCILLATOR.acceleration_noise_asd(f)**2
result=dict(bandwidth_hz=B,gaussian_uniform_psd_bit_s=float(B*np.log2(1+P/B/noise).mean()),gaussian_waterfilled_bit_s=gaussian_waterfill(noise,B,P),seven_tone_alphabet_bit_s=float(np.log2(7)/8),byte_mapped_bit_s=1/3,byte_mapping_efficiency=float(8/(3*np.log2(7))),dwell_fraction=.8,filter_bandwidth_per_spacing=10,cf_spectral_efficiency=(1/3)/B,peak_torque_only_rate_ceiling_bit_s=(1/3)*np.sqrt((2-.05)/1.090906))
(root/'results/scaled_cf_fsk/desktop_capacity_gap.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
