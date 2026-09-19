"""Geometrically scaled CF-FSK scenarios with centrifugal similitude.

Retunes the modeled receiver resonance, retaining its other noise parameters.
Uses exact Newtonian similarity a_s(d)=s*a_1(d/s), not a far-field estimate.
Outputs separate results; the desktop results and model are unchanged.
"""
from dataclasses import replace, asdict
import json
from pathlib import Path
import numpy as np
import desktop_cf_fsk as cf
from gravcomm.receivers import StructuralOscillator
from gravcomm.actuator import Actuator, simulate_actuator

root=Path(__file__).resolve().parents[1]
out=root/'results/scaled_cf_fsk';out.mkdir(exist_ok=True)
desktop=json.loads((root/'results/desktop_cf_fsk/summary.json').read_text())
field=json.loads((root/'results/paper_revision/gaussian_benchmark.json').read_text())
a1=next(r['a2_m_s2'] for r in field['rows'] if r['distance_m']==1)
symbols,previous=cf.all_transitions()
results=[]
for scale in (4.,10.):
    fc=50.3/scale; inertia=desktop['inertia_kg_m2']*scale**5
    cf.FC=fc;cf.CARTER_OSCILLATOR=replace(StructuralOscillator(),resonance_hz=fc)
    amplitude=scale*a1;screen=[];candidate=None
    for ts in (2.,4.,6.,8.,12.,16.,24.):
        if fc+3/(.8*ts)>60/scale:continue
        t,f,df,ddf,phase=cf.trajectory(symbols,previous,ts)
        means,cov,rec=cf.receiver_statistics(ts,phase,amplitude)
        mc=cf.monte_carlo(means,cov,symbols,10000,6814)
        row=dict(scale=scale,symbol_s=ts,rate=8/(3*ts),errors=mc['total_errors'],worst_ser=mc['worst_observed_ser'])
        screen.append(row);print(json.dumps(row),flush=True)
        if mc['worst_observed_ser']<1e-3:candidate=(ts,t,f,df,ddf,phase,rec);break
    if candidate is None:raise RuntimeError('No candidate')
    ts,t,f,df,ddf,phase,rec=candidate
    # Explicit assumed drive specification, chosen above trajectory demand.
    # Drag torque scales as s^3 for the scenario; no aerodynamic claim is made.
    response=.02;drag_torque=.05*scale**3;drag=drag_torque/(np.pi*fc)
    inertial_peak=1.875*np.pi*inertia*6/(.8*ts)/(.2*ts)
    peak_limit=float(np.ceil((inertial_peak+drag_torque)*1.2/100)*100)
    power_limit=float(np.ceil(peak_limit*np.pi*(fc+3/(.8*ts))/1000)*1000)
    drive=Actuator(inertia,peak_limit,.5*peak_limit,response,power_limit,power_limit,60*np.pi/scale,drag)
    desired=inertia*np.pi*df+drag*np.pi*f
    derivative=inertia*np.pi*ddf+drag*np.pi*df
    trace=simulate_actuator(drive,t,desired+response*derivative,np.pi*f[0],
        initial_demand_nm=desired[0],target_carrier_hz=f,max_step_s=.005)
    actual_phase=2*trace.phase_rad-2*np.pi*fc*t
    means,cov,rec=cf.receiver_statistics(ts,actual_phase,amplitude)
    mc=cf.monte_carlo(means,cov,symbols,100000,39047)
    assert mc['worst_simultaneous_95pct_ser_upper']<1e-3
    assert not trace.diagnostics['speed_limit_exceeded']
    assert not trace.diagnostics['continuous_torque_exceeded']
    result=dict(scale=scale,mass_kg=7.913509491619676*scale**3,diameter_m=.5*scale,
        distance_m=scale,amplitude_m_s2=amplitude,carrier_hz=fc,symbol_s=ts,
        payload_bit_s=8/(3*ts),inertia_kg_m2=inertia,receiver=asdict(cf.CARTER_OSCILLATOR),
        receiver_band=rec,actuator=asdict(drive),diagnostics=trace.diagnostics,
        centrifugal_regional_mean_at_max_speed_mpa=25.475687623668655*(scale*np.max(f)/60)**2,
        maximum_energy_j=inertia*(np.pi*np.max(f))**2/2,
        monte_carlo=mc,screen=screen)
    np.savez_compressed(out/f'scale_{scale:g}_detector.npz',means=means,covariance=cov,symbols=symbols)
    results.append(result)
    (out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(result),flush=True)

