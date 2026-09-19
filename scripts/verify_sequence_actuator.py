"""Check selected sequence candidate with achieved phase and constrained drive."""
import argparse,json
import numpy as np
from scipy.signal import lfilter
import optimize_cf_fsk as dynamics
import desktop_cf_fsk as cf
from search_sequence_cf_fsk import A,OUT
from gravcomm.actuator import Actuator,simulate_actuator
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import WhitenedViterbi,fit_whitener

parser=argparse.ArgumentParser();parser.add_argument('--selected',action='store_true');args=parser.parse_args()
ts=2.;r,q=(.4,.1) if args.selected else (.5,.0625);I=dynamics.I
t,f,df,ddf,phase=dynamics.trajectory(ts,r,q)
regional=dynamics.mechanics(ts,r,q,t,f,df)
drag=.05/(np.pi*cf.FC)
drive=Actuator(I,2,.5,.02,400,400,60*np.pi,drag)

def run(t,f,df,ddf):
    tau=I*np.pi*df+drag*np.pi*f
    return simulate_actuator(drive,t,tau+.02*(I*np.pi*ddf+drag*np.pi*df),np.pi*f[0],
        initial_demand_nm=tau[0],target_carrier_hz=f,max_step_s=.005)

trace=run(t,f,df,ddf)
result=dict(all_transition_mechanics=regional,all_transition_drive=trace.diagnostics)
assert regional['feasible']
assert not trace.diagnostics['speed_limit_exceeded']
assert not trace.diagnostics['continuous_torque_exceeded']
assert abs(trace.diagnostics['energy_balance_error_j'])<1e-4

rng=np.random.default_rng(89271);data=rng.integers(0,256,8)
symbols=cf.encode(data);tail=6
dynamics.SYMBOLS=np.r_[symbols,np.zeros(tail,int)];dynamics.PREVIOUS=0
t,f,df,ddf,phase=dynamics.trajectory(ts,r,q);packet=run(t,f,df,ddf)
model=CFFSK(ts,r,q,4);taps,fit=fit_whitener(2,A,12)
decoder=WhitenedViterbi(model,taps)
actual_phase=2*packet.phase_rad-2*np.pi*cf.FC*t
actual=np.exp(1j*actual_phase[:-1:50])
ideal=model.waveform(dynamics.SYMBOLS)
assert len(actual)==len(ideal)
prefix=decoder.memory_symbols*model.samples_per_symbol
observed=lfilter(taps,[1.],np.r_[np.ones(prefix),actual])[prefix:]
decoded=decoder.decode_whitened(observed,tail)[:len(symbols)]
np.testing.assert_array_equal(decoded,symbols)
result['packet_drive']=packet.diagnostics
result['maximum_complex_waveform_error']=float(max(abs(actual-ideal)))
result['actual_drive_packet_decodes']=True
result['payload_hex']=data.astype(np.uint8).tobytes().hex()
filename='selected_actuator.json' if args.selected else 'actuator_check.json'
(OUT/filename).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
