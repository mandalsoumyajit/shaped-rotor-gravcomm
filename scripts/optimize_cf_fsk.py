"""Finite-grid desktop CF-FSK search; conditional stationary receiver model.

Keeps seven tones, byte mapping, quintic transitions and a fixed receiver band.
Independent final trials prevent the preliminary zero-error selection rule.
The dwell likelihood uses unknown-phase pure-tone templates, never true labels
or the known preceding symbol. Full packet acquisition remains outside scope.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import i0e
from scipy.stats import beta
import desktop_cf_fsk as cf
from gravcomm.actuator import Actuator, simulate_actuator

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/cf_fsk_search'
OUT.mkdir(exist_ok=True)
BASE=json.loads((ROOT/'results/desktop_cf_fsk/summary.json').read_text())
I=BASE['inertia_kg_m2']
A=5.180259430854158e-10
B=1.5625
FS=100
SYMBOLS,PREVIOUS=cf.all_transitions()

def trajectory(ts,r,q):
    tr=r*ts;td=ts-tr;spacing=q/td
    n=round(ts*FS);local=np.arange(n)/FS;u=np.minimum(local/tr,1)
    s=10*u**3-15*u**4+6*u**5
    ds=(30*u*u-60*u**3+30*u**4)/tr
    dds=(60*u-180*u*u+120*u**3)/tr**2
    primitive=tr*(2.5*u**4-3*u**5+u**6)+np.maximum(local-tr,0)
    f=[];df=[];ddf=[];phase=[];phi=0.;previous=PREVIOUS
    for current in SYMBOLS:
        before=previous*spacing;delta=(current-previous)*spacing
        f.extend(cf.FC+before+delta*s);df.extend(delta*ds);ddf.extend(delta*dds)
        phase.extend(phi+2*np.pi*(before*local+delta*primitive))
        phi+=2*np.pi*(before*ts+delta*(ts-tr/2));previous=current
    t=np.arange(n*len(SYMBOLS)+1)/FS
    return t,np.r_[f,cf.FC+previous*spacing],np.r_[df,0],np.r_[ddf,0],np.r_[phase,phi]

def statistics(ts,r,q,phase):
    td=(1-r)*ts;spacing=q/td;n=round(ts*FS);start=round(r*ts*FS)
    # Grid choices guarantee exact transition and dwell sample boundaries.
    assert abs(start-r*ts*FS)<1e-8
    local=np.arange(n-start)/FS;offsets=np.arange(-3,4)*spacing
    envelope=A*np.exp(1j*phase[:-1]);freq=np.fft.fftfreq(len(envelope),1/FS)
    filtered=np.fft.ifft(np.fft.fft(envelope)*(abs(freq)<=B/2))
    projections=np.exp(-2j*np.pi*offsets[:,None]*local)/len(local)
    means=filtered.reshape(-1,n)[:,start:]@projections.T
    nu=np.linspace(-B/2,B/2,12001)
    W=np.sinc((nu[None,:]-offsets[:,None])*td)*np.exp(1j*np.pi*(nu[None,:]-offsets[:,None])*td)
    psd=cf.CARTER_OSCILLATOR.acceleration_noise_asd(cf.FC+nu)**2
    weights=np.ones(len(nu))*(nu[1]-nu[0]);weights[[0,-1]]*=.5
    cov=(W*(2*psd*weights))@W.conj().T
    # Each row is the seven projections of one possible pure dwell tone.
    templates=A*np.exp(2j*np.pi*offsets[:,None]*local)@projections.T
    return means/A,cov/A**2,templates/A

def trial(means,cov,templates,count,seed):
    rng=np.random.default_rng(seed);L=np.linalg.cholesky(cov)
    inv=np.linalg.inv(cov)
    linear=inv@templates.T.conj()
    energy=np.real(np.sum(templates*(inv@templates.T.conj()).T,axis=1))
    errors={'energy':[],'likelihood':[]}
    for mean,target in zip(means,SYMBOLS+3):
        totals={key:0 for key in errors}
        for start in range(0,count,2000):
            size=min(2000,count-start)
            z=(rng.normal(size=(size,7))+1j*rng.normal(size=(size,7)))/np.sqrt(2)
            y=mean+z@L.T
            x=2*abs(y@linear)
            scores=np.log(i0e(x))+x-energy
            totals['energy']+=int(np.sum(np.argmax(abs(y)**2,axis=1)!=target))
            totals['likelihood']+=int(np.sum(np.argmax(scores,axis=1)!=target))
        for key in errors:errors[key].append(totals[key])
    return {key:dict(errors_per_context=values,total_errors=sum(values),
        worst_observed=max(values)/count,
        upper=float(max(beta.ppf(1-.05/49,np.array(values)+1,count-np.array(values)))),
        trials_per_context=count) for key,values in errors.items()}

def mechanics(ts,r,q,t,f,df):
    torque=I*np.pi*df+.05*f/cf.FC
    n=round(ts*FS)
    # Upper bound on any repeated symbol history: largest per-symbol RMS.
    rms=float(np.sqrt(np.max(np.mean(torque[:-1].reshape(-1,n)**2,axis=1))))
    peak=float(max(abs(torque)));power=float(max(abs(torque*np.pi*f)))
    return dict(peak_nm=peak,worst_symbol_rms_nm=rms,peak_power_w=power,
        feasible=bool(peak<=2 and rms<=.5 and power<=400 and max(f)<=60))

def main():
    rows=[]
    for ts in np.arange(3.,8.01,.5):
        for r in (.1,.2,.3,.4,.5):
            for q in (.5,.65,.8,1.,1.2):
                if 3*q/((1-r)*ts)>B/2:continue
                t,f,df,ddf,phase=trajectory(ts,r,q)
                mech=mechanics(ts,r,q,t,f,df)
                if not mech['feasible']:continue
                means,cov,templates=statistics(ts,r,q,phase)
                result=trial(means,cov,templates,2000,6814)
                row=dict(ts=float(ts),transition_fraction=r,spacing_dwell_product=q,
                    rate=8/(3*ts),mechanics=mech,screen=result)
                rows.append(row)
        print(f'Completed Ts={ts:g}: {len(rows)} feasible candidates',flush=True)
        (OUT/'screen.json').write_text(json.dumps(rows,indent=2)+'\n')
    # Freeze shortlist before independently seeded validation. Include both
    # receiver methods; do not require zero preliminary errors.
    selected=[]
    for method in ('energy','likelihood'):
        eligible=[x for x in rows if x['screen'][method]['worst_observed']<=.002]
        eligible.sort(key=lambda x:(x['ts'],x['screen'][method]['worst_observed']))
        for row in eligible[:5]:
            ts=row['ts'];r=row['transition_fraction'];q=row['spacing_dwell_product']
            t,f,df,ddf,phase=trajectory(ts,r,q)
            drag=.05/(np.pi*cf.FC);tau=I*np.pi*df+drag*np.pi*f
            command=tau+.02*(I*np.pi*ddf+drag*np.pi*df)
            drive=Actuator(I,2,.5,.02,400,400,60*np.pi,drag)
            trace=simulate_actuator(drive,t,command,np.pi*f[0],initial_demand_nm=tau[0],target_carrier_hz=f,max_step_s=.005)
            actual=2*trace.phase_rad-2*np.pi*cf.FC*t
            means,cov,templates=statistics(ts,r,q,actual)
            result=trial(means,cov,templates,100000,39047)
            record=dict(candidate=row,method=method,verification=result,drive=trace.diagnostics)
            selected.append(record)
            print(json.dumps(record),flush=True)
            (OUT/'verified.json').write_text(json.dumps(selected,indent=2)+'\n')
            if result[method]['upper']<1e-3:break

if __name__=='__main__':
    raise SystemExit('Dwell-only search retired. Validate the full-waveform sequence receiver before restarting rate optimization.')
