"""Independent similarity, drive-demand and confidence-bound checks."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import beta
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from gravcomm.finite_rotor import FiniteRotor
r=FiniteRotor.from_npz(root/'fea/results/shaped_waveforms/medium_rounded.npz')
phase=np.arange(8)*np.pi/4
baseline=r.acceleration_component(1,phase)
results=json.loads((root/'results/scaled_cf_fsk/summary.json').read_text())
for item in results:
    s=item['scale'];ts=item['symbol_s'];I=item['inertia_kg_m2']
    scaled=FiniteRotor(r.positions_m*s,r.masses_kg*s**3,r.bounding_radius_m*s)
    np.testing.assert_allclose(scaled.acceleration_component(s,phase),s*baseline,rtol=1e-12,atol=0)
    peak=15*np.pi*I*6/(8*.2*.8*ts**2)
    # Difference from inertial peak is consistent with the small added drag.
    assert abs(item['diagnostics']['peak_torque_nm']-peak)<.05*s**3*1.3
    assert item['diagnostics']['peak_torque_nm']<item['actuator']['peak_torque_nm']
    assert item['diagnostics']['peak_motoring_power_w']<item['actuator']['motoring_power_w']
    assert item['diagnostics']['peak_regenerative_power_w']<item['actuator']['regenerative_power_w']
    assert abs(item['diagnostics']['energy_balance_error_j'])/item['maximum_energy_j']<1e-9
    assert item['monte_carlo']['worst_simultaneous_95pct_ser_upper']<1e-3
    np.testing.assert_allclose(item['monte_carlo']['worst_simultaneous_95pct_ser_upper'],max(beta.ppf(1-.05/49,np.array(item['monte_carlo']['errors_per_context'])+1,100000-np.array(item['monte_carlo']['errors_per_context']))))
    fmax=item['carrier_hz']+3/(.8*ts)
    assert fmax<=60/s
    np.testing.assert_allclose(item['centrifugal_regional_mean_at_max_speed_mpa'],25.475687623668655*(s*fmax/60)**2)
print('PASS: direct field similarity, analytic torque, drive limits, work balance, centrifugal scaling, and confidence bounds')

