"""Drive integration refinement and paired nominal/achieved noisy decisions."""
import json
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
import desktop_cf_fsk as cf
import optimize_cf_fsk as dyn
from gravcomm.actuator import Actuator,simulate_actuator
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import fit_whitener,receiver_noise,WhitenedViterbi

OUT=Path(__file__).resolve().parents[1]/'results/communications_checks'
OUT.mkdir(exist_ok=True)
A=5.180259430854158e-10
drive=Actuator(dyn.I,2,.5,.02,400,400,60*np.pi,.05/(np.pi*cf.FC))
records=[]
for ts,q in [(2.,.1),(1.75,.1)]:
    fs=8.; model=CFFSK(ts,.4,q,round(ts*fs))
    taps,fit=fit_whitener(fs,A,3*model.samples_per_symbol)
    decoder=WhitenedViterbi(model,taps,maximum_states=1000000)
    prefix=decoder.memory_symbols*model.samples_per_symbol
    for kind in ['worst_duty','random']:
        data=(np.array([251,42]*4,dtype=np.uint8) if kind=='worst_duty' else
              np.random.default_rng(832919).integers(0,256,8,dtype=np.uint8))
        symbols=cf.encode(data); dyn.SYMBOLS=np.r_[symbols,np.zeros(6,int)];dyn.PREVIOUS=0
        nominal=model.waveform(dyn.SYMBOLS); actual=[];diagnostics=[]
        for grid in [200,400]:
            dyn.FS=grid
            t,f,df,ddf,phase=dyn.trajectory(ts,.4,q)
            tau=dyn.I*np.pi*df+drive.viscous_drag_nm_s*np.pi*f
            cmd=tau+.02*(dyn.I*np.pi*ddf+drive.viscous_drag_nm_s*np.pi*df)
            trace=simulate_actuator(drive,t,cmd,np.pi*f[0],initial_demand_nm=tau[0],
                target_carrier_hz=f,max_step_s=1/grid)
            actual.append(np.exp(1j*(2*trace.phase_rad-2*np.pi*cf.FC*t))[:-1:int(grid/fs)])
            diagnostics.append(trace.diagnostics)
            assert abs(trace.diagnostics['energy_balance_error_j'])<1e-4
        delta=float(max(abs(actual[1]-actual[0])))
        signal=np.r_[np.ones(prefix),actual[1]]
        ideal=np.r_[np.ones(prefix),nominal]
        out=decoder.decode_whitened(lfilter(taps,[1.],signal)[prefix:],6,True)[:24]
        np.testing.assert_array_equal(out,symbols)
        bit_errors=[0,0];disagreements=0
        for index in range(500):
            rng=np.random.default_rng(np.random.SeedSequence([786133,index]))
            noise=receiver_noise(rng,len(signal),fs,A)
            decoded=[]
            for j,wave in enumerate([ideal,signal]):
                dec=decoder.decode_whitened(lfilter(taps,[1.],wave+noise)[prefix:],6,True)[:24]
                decoded.append(dec)
                bit_errors[j]+=int(np.unpackbits(data^cf.decode(dec)).sum())
            disagreements+=int(np.any(decoded[0]!=decoded[1]))
        record=dict(symbol_s=ts,q=q,payload_kind=kind,payload_hex=data.tobytes().hex(),
            coarse_grid_hz=200,fine_grid_hz=400,drive_diagnostics=diagnostics,
            coarse_fine_max_envelope_difference=delta,
            achieved_nominal_max_envelope_difference=float(max(abs(actual[1]-nominal))),
            paired_noise_trials=500,nominal_bit_errors=bit_errors[0],achieved_bit_errors=bit_errors[1],
            disagreeing_decoded_packets=disagreements,
            scope='Fixed-message paired sensitivity diagnostic; no random-message BER estimate')
        records.append(record)
        (OUT/'drive_refinement.json').write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps(record),flush=True)

# Closed 49-transition sequence supplies the revised paper energy example.
dyn.SYMBOLS,dyn.PREVIOUS=cf.all_transitions();dyn.FS=400
t,f,df,ddf,phase=dyn.trajectory(2.,.4,.1)
tau=dyn.I*np.pi*df+drive.viscous_drag_nm_s*np.pi*f
trace=simulate_actuator(drive,t,tau+.02*(dyn.I*np.pi*ddf+drive.viscous_drag_nm_s*np.pi*df),
    np.pi*f[0],initial_demand_nm=tau[0],target_carrier_hz=f,max_step_s=.0025)
record=dict(duration_s=float(t[-1]),diagnostics=trace.diagnostics,
            centrifugal_regional_mean_mpa=25.4756876*(1516.5/1800)**2,
            stored_energy_at_max_speed_j=.5*dyn.I*(np.pi*50.55)**2)
assert abs(trace.diagnostics['energy_balance_error_j'])<1e-4
(OUT/'revised_energy.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)
