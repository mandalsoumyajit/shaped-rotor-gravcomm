"""Reproducible seven-tone desktop operating point; conditional model, not hardware.

Uses the full shaped source's computed second harmonic, quintic frequency
transitions, a constrained rigid actuator and calibrated band-limited receiver.
No Shannon-capacity substitution for the symbol detector is made.
"""
import csv
from dataclasses import asdict
import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.stats import beta as beta_distribution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from gravcomm.receivers import CARTER_OSCILLATOR
from gravcomm.actuator import Actuator,simulate_actuator

OUT=ROOT/'results'/'desktop_cf_fsk'
CASE='medium_rounded'
FC=50.3
FS=100
M=7
TR_FRACTION=.2


def all_transitions():
    """Eulerian circuit of the complete directed graph, including self edges."""
    remaining=[list(range(M)) for _ in range(M)];stack=[0];path=[]
    while stack:
        if remaining[stack[-1]]:stack.append(remaining[stack[-1]].pop())
        else:path.append(stack.pop())
    path=path[::-1]
    assert len(set(zip(path[:-1],path[1:])))==M*M
    return np.array(path[1:])-3,path[0]-3


def encode(data):
    values=np.asarray(data,dtype=int)
    if np.any((values<0)|(values>255)):raise ValueError('Bytes required')
    return np.stack([values//49,(values//7)%7,values%7],axis=1).ravel()-3


def decode(symbols):
    digits=np.asarray(symbols,dtype=int).reshape(-1,3)+3
    values=digits@np.array([49,7,1])
    if np.any((digits<0)|(digits>6)) or np.any(values>255):raise ValueError('Invalid codeword')
    return values.astype(np.uint8)


def trajectory(symbols,previous,ts,fs=FS):
    tr=TR_FRACTION*ts;td=ts-tr;spacing=1/td
    n=round(ts*fs);local=np.arange(n)/fs;u=np.minimum(local/tr,1)
    s=10*u**3-15*u**4+6*u**5
    ds=(30*u**2-60*u**3+30*u**4)/tr
    dds=(60*u-180*u**2+120*u**3)/tr**2
    primitive=tr*(2.5*u**4-3*u**5+u**6)+np.maximum(local-tr,0)
    freq=[];df=[];ddf=[];phase=[];phi=0.
    for current in symbols:
        before=previous*spacing;delta=(current-previous)*spacing
        freq.append(FC+before+delta*s);df.append(delta*ds);ddf.append(delta*dds)
        phase.append(phi+2*np.pi*(before*local+delta*primitive))
        phi+=2*np.pi*(before*ts+delta*(ts-tr/2));previous=current
    t=np.arange(n*len(symbols)+1)/fs
    return t,np.r_[np.concatenate(freq),FC+previous*spacing],np.r_[np.concatenate(df),0.],np.r_[np.concatenate(ddf),0.],np.r_[np.concatenate(phase),phi]


def receiver_statistics(ts,phase,amplitude,fs=FS,technical_asd=0.):
    """Complex dwell-bin means and proper complex Gaussian noise covariance.

    Ideal calibrated analytic acceleration, ideal symmetric bandpass; oscillator
    initial-state/calibration uncertainty and technical noise are excluded.
    One-sided real acceleration PSD S_a becomes complex baseband PSD 2*S_a.
    """
    tr=TR_FRACTION*ts;td=ts-tr;spacing=1/td;halfband=5/td
    envelope=amplitude*np.exp(1j*phase[:-1])
    bins=np.fft.fftfreq(len(envelope),1/fs)
    filtered=np.fft.ifft(np.fft.fft(envelope)*(abs(bins)<=halfband))
    n=round(ts*fs);start=round(tr*fs);nd=round(td*fs)
    dwell=filtered.reshape(-1,n)[:,start:]
    offsets=np.arange(-3,4)*spacing
    templates=np.exp(-2j*np.pi*offsets[:,None]*np.arange(nd)[None,:]/fs)/nd
    means=dwell@templates.T
    nu=np.linspace(-halfband,halfband,12001)
    W=np.sinc((nu[None,:]-offsets[:,None])*td)*np.exp(1j*np.pi*(nu[None,:]-offsets[:,None])*td)
    psd=CARTER_OSCILLATOR.acceleration_noise_asd(FC+nu)**2+technical_asd**2
    weights=np.ones(len(nu))*(nu[1]-nu[0]);weights[[0,-1]]*=.5
    cov=(W*(2*psd*weights))@W.conj().T
    np.testing.assert_allclose(cov,cov.conj().T,atol=1e-35)
    assert np.min(np.linalg.eigvalsh(cov))>0
    return means,cov,dict(spacing_hz=spacing,transition_s=tr,dwell_s=td,
        receiver_halfband_hz=halfband,noise_asd_at_tones=CARTER_OSCILLATOR.acceleration_noise_asd(FC+offsets).tolist())


def monte_carlo(means,cov,symbols,trials,seed):
    """Independent noise trials per fixed signal context; not independent symbols
    in a physical packet. Temporal packet-noise correlations are not simulated.
    """
    scale=np.max(abs(means));L=np.linalg.cholesky(cov/scale**2)
    rng=np.random.default_rng(seed);counts=[]
    for mean,target in zip(means/scale,symbols+3):
        errors=0
        for start in range(0,trials,5000):
            count=min(5000,trials-start)
            z=(rng.normal(size=(count,7))+1j*rng.normal(size=(count,7)))/np.sqrt(2)
            noisy=mean+z@L.T
            errors+=np.count_nonzero(np.argmax(abs(noisy)**2,axis=1)!=target)
        counts.append(errors)
    # Bonferroni makes the one-sided 95% confidence limits simultaneous over
    # the tested fixed contexts. It does not establish bounds for all packets.
    upper=beta_distribution.ppf(1-.05/len(counts),np.array(counts)+1,trials-np.array(counts))
    return dict(trials_per_context=trials,total_trials=trials*len(counts),
        errors_per_context=counts,total_errors=sum(counts),
        mean_ser=sum(counts)/(trials*len(counts)),worst_observed_ser=max(counts)/trials,
        worst_simultaneous_95pct_ser_upper=float(max(upper)))


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def plot_saved():
    result=json.loads((OUT/'summary.json').read_text())
    rows=list(csv.DictReader((OUT/'trace.csv').open()))
    values={key:np.array([float(r[key]) for r in rows]) for key in rows[0]}
    t=values['time_s'];show=t<=8*result['symbol_s']
    with plt.rc_context({'font.family':'serif','font.serif':['STIXGeneral'],
                         'mathtext.fontset':'stix','font.size':7,'axes.labelsize':7,
                         'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.5,
                         'xtick.major.size':2.5,'ytick.major.size':2.5,'pdf.fonttype':42}):
        fig,ax=plt.subplots(3,1,figsize=(3.4,2.5),sharex=True)
        ax[0].plot(t[show],values['target_hz'][show],label='Target',lw=1.2)
        ax[0].plot(t[show],values['actual_hz'][show],'--',label='Achieved',lw=.8)
        ax[0].set_ylabel('Signal (Hz)');ax[0].legend(fontsize=6,loc='lower right',frameon=False,ncol=2)
        ax[1].plot(t[show],values['torque_nm'][show],lw=.8);ax[1].set_ylabel('Torque (N m)')
        ax[2].plot(t[show],values['power_w'][show],lw=.8);ax[2].set_ylabel('Power (W)');ax[2].set_xlabel('Time (s)')
        for a in ax:a.grid(alpha=.16,lw=.4)
        fig.tight_layout(pad=.4,h_pad=.35)
        fig.savefig(OUT/'operating_point.png',dpi=300,bbox_inches='tight',pad_inches=.02)
        fig.savefig(OUT/'operating_point.pdf',bbox_inches='tight',pad_inches=.02);plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    source=json.loads((ROOT/'fea'/'results'/'design_sweep'/CASE/'h0.013_rpm1800'/'metadata.json').read_text())
    inertia=source['summary']['aluminum_polar_inertia_kg_m2']+source['summary']['steel_polar_inertia_kg_m2']
    wave=list(csv.DictReader((ROOT/'fea'/'results'/'shaped_waveforms'/'summary.csv').open()))
    amplitude=float(next(r['a2_m_s2'] for r in wave if r['case']==CASE and float(r['distance_m'])==.5))
    symbols,previous=all_transitions();screen=[];selected=None
    for ts in (4.,6.,8.,10.,12.,16.,20.):
        t,f,df,ddf,phase=trajectory(symbols,previous,ts)
        means,cov,receiver=receiver_statistics(ts,phase,amplitude)
        result=monte_carlo(means,cov,symbols,5000,876)
        row=dict(symbol_s=ts,payload_bit_s=8/(3*ts),peak_inertial_torque_nm=inertia*np.pi*6*receiver['spacing_hz']/receiver['transition_s']*1.875,
            worst_observed_ser=result['worst_observed_ser'],mean_ser=result['mean_ser'])
        screen.append(row);print(json.dumps(row),flush=True)
        # Coarse search selects a candidate with zero observed errors in each
        # tested context; larger verification follows with an explicit CI.
        if result['total_errors']==0:selected=(ts,t,f,df,ddf,phase,receiver);break
    if selected is None:raise RuntimeError('No candidate passed screen')
    ts,t,f,df,ddf,phase,receiver=selected
    response=.02;drag=.05/(np.pi*FC)
    drive=Actuator(inertia,2.,.5,response,400.,400.,60*np.pi,drag)
    desired=inertia*np.pi*df+drag*np.pi*f
    derivative=inertia*np.pi*ddf+drag*np.pi*df
    command=desired+response*derivative
    trace=simulate_actuator(drive,t,command,np.pi*f[0],initial_demand_nm=desired[0],target_carrier_hz=f,max_step_s=.005)
    actual_phase=2*trace.phase_rad-2*np.pi*FC*t
    means,cov,receiver=receiver_statistics(ts,actual_phase,amplitude)
    np.savez_compressed(OUT/'detector_inputs.npz',means=means,covariance=cov,symbols=symbols,
                        actual_phase=actual_phase,symbol_s=ts,amplitude=amplitude)
    verification=monte_carlo(means,cov,symbols,100000,17876)
    if verification['worst_simultaneous_95pct_ser_upper']>1e-3:raise RuntimeError('Candidate failed verification')
    sensitivity=[]
    for technical_asd in (0.,5e-11,1e-10,2e-10,5e-10):
        _,extra_cov,_=receiver_statistics(ts,actual_phase,amplitude,technical_asd=technical_asd)
        mc=monte_carlo(means,extra_cov,symbols,20000,714)
        sensitivity.append(dict(additional_technical_asd_m_s2_sqrt_hz=technical_asd,
            mean_ser=mc['mean_ser'],worst_observed_ser=mc['worst_observed_ser'],trials=mc['total_trials']))
    write_csv(OUT/'technical_noise_sensitivity.csv',sensitivity)
    assert not trace.diagnostics['speed_limit_exceeded']
    assert not trace.diagnostics['continuous_torque_exceeded']
    # Independently check sample refinement in trajectory/bin integration.
    _,_,_,_,finephase=trajectory(symbols,previous,ts,fs=2*FS)
    fine_means,_,_=receiver_statistics(ts,finephase,amplitude,fs=2*FS)
    ideal_means,_,_=receiver_statistics(ts,phase,amplitude)
    bin_refinement=float(np.max(abs(fine_means-ideal_means))/amplitude)
    assert bin_refinement<.001
    np.testing.assert_array_equal(decode(encode(np.arange(256))),np.arange(256))
    # Random, noiseless byte packet: padding separates data from FFT wraparound.
    rng=np.random.default_rng(421);data=rng.integers(0,256,256)
    encoded=encode(data);padded=np.r_[np.repeat(encoded[0],8),encoded,np.repeat(encoded[-1],8)]
    _,_,_,_,packet_phase=trajectory(padded,int(padded[0]),ts)
    packet_means,_,_=receiver_statistics(ts,packet_phase,amplitude)
    decisions=np.argmax(abs(packet_means[8:-8])**2,axis=1)-3
    np.testing.assert_array_equal(decode(decisions),data)
    result=dict(case=CASE,inertia_kg_m2=inertia,source_peak_acceleration_m_s2=amplitude,distance_from_axis_m=.5,
        symbol_s=ts,payload_bit_s=8/(3*ts),alphabet_information_bit_s=math.log2(7)/ts,
        bit_mapping='8 bits -> 3 base-7 digits; 256 of 343 words, no FEC; excludes preamble/pilots',
        carrier_hz=FC,tones_hz=(FC+np.arange(-3,4)*receiver['spacing_hz']).tolist(),
        receiver=receiver,oscillator=asdict(CARTER_OSCILLATOR),actuator=asdict(drive),
        diagnostics=trace.diagnostics,monte_carlo=verification,
        peak_inertial_torque_nm=float(np.max(abs(inertia*np.pi*df))),
        worst_alternating_rms_torque_nm=float(np.sqrt((inertia*np.pi*6*receiver['spacing_hz']/receiver['transition_s'])**2*TR_FRACTION*10/7+.05**2)),
        max_phase_error_rad=float(np.max(abs(actual_phase-phase))),
        receiver_bin_refinement_relative=bin_refinement,noiseless_packet_bytes=len(data),
        assumptions=['ideal calibrated acceleration and ideal bandpass',
          'structural thermal noise plus flat displacement-readout model; not measured carrier-band performance',
          'no environmental/crosstalk noise, calibration or initial oscillator-state uncertainty',
          'symbol clock and carrier frequency acquired; noncoherent dwell-bin energy detection',
          'Monte Carlo confidence applies to 49 tested cyclic contexts, not universal packet BER',
          'assumed actuator ratings and drag; rigid rotor; no structural safety certification'])
    write_csv(OUT/'screen.csv',screen)
    write_csv(OUT/'transition_errors.csv',[dict(previous=int(p),current=int(c),errors=e,trials=100000)
        for p,c,e in zip(np.r_[previous,symbols[:-1]],symbols,verification['errors_per_context'])])
    write_csv(OUT/'trace.csv',[dict(time_s=t[i],target_hz=f[i],actual_hz=trace.carrier_hz()[i],
        torque_nm=trace.actual_torque_nm[i],power_w=trace.mechanical_power_w[i]) for i in range(0,len(t),10)])
    (OUT/'summary.json').write_text(json.dumps(result,indent=2))
    plot_saved()
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    if '--plot-only' in sys.argv:plot_saved()
    else:main()
