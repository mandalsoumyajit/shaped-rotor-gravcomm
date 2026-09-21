"""New range study; preserves baseline data. Run with Numba and one BLAS thread.

Screening is NOT validation. Uses the current acquired-packet, byte-constrained
whitened Viterbi receiver. Optional receiver redesign and drive scaling are
explicit scenario assumptions. No new FEA or acquisition simulation is claimed.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name] = '1'
import argparse, csv, hashlib, itertools, json, math, sys, time, zlib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from gravcomm import FiniteRotor
from gravcomm.receivers import StructuralOscillator
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import WhitenedViterbi, fit_whitener, receiver_noise, noise_psd
from gravcomm.packet_confidence import packet_ber_interval
from paper_results import gaussian_waterfill
from desktop_cf_fsk import encode, decode
OUT = ROOT/'results/range_revision_2026_09_20'
OUT.mkdir(parents=True, exist_ok=True)
MASS, INERTIA, FC = 7.913509491619676, .3160713095850955, 50.3
TARGET, TAIL = 1e-3, 6


def save(path, record):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


class CenteredReceiver:
    """Adapt the legacy PSD's fixed 50.3 Hz origin to a physical carrier.

    Only translates its frequency argument. Does not change the physical
    receiver's resonance or any noise parameter.
    """
    def __init__(self, receiver, carrier):
        self.receiver, self.carrier = receiver, carrier
    def acceleration_noise_asd(self, f):
        return self.receiver.acceleration_noise_asd(np.asarray(f)-FC+self.carrier)


def duty(ts, r, q, inertia=INERTIA, fc=FC, drag=.05, peak_limit=2., rms_limit=.5, power_limit=400., max_rpm=1800.):
    t=np.linspace(0,ts,2001); u=np.minimum(t/(r*ts),1)
    smooth=10*u**3-15*u**4+6*u**5
    slope=(30*u**2-60*u**3+30*u**4)/(r*ts)
    spacing=q/((1-r)*ts); costs=np.zeros((7,7)); peak=power=0.
    for a,b in itertools.product(range(7),repeat=2):
        f=fc+spacing*((a-3)+(b-a)*smooth)
        tau=inertia*np.pi*spacing*(b-a)*slope+drag*f/fc
        costs[a,b]=np.trapz(tau*tau,t)/ts
        peak=max(peak,float(np.max(abs(tau))))
        power=max(power,float(np.max(abs(tau*np.pi*f))))
    words=np.array([[b//49,(b//7)%7,b%7] for b in range(256)])
    edges=np.full((7,7),-np.inf)
    for old in range(7):
        for a,b,c in words:
            edges[old,c]=max(edges[old,c],costs[old,a]+costs[a,b]+costs[b,c])
    # Enumerate all simple cycles; maximum cycle mean gives exact sustained
    # worst valid-byte duty (transients require separate peak checks).
    worst=0.
    for length in range(1,8):
        for nodes in itertools.combinations(range(7),length):
            for rest in itertools.permutations(nodes[1:]):
                cycle=(nodes[0],)+rest
                worst=max(worst,sum(edges[a,b] for a,b in zip(cycle,cycle[1:]+cycle[:1]))/(3*length))
    rms=math.sqrt(worst); rpm=30*(fc+3*spacing)
    return dict(peak_torque_nm=peak,worst_valid_message_rms_nm=rms,peak_mechanical_power_w=power,
                maximum_rpm=rpm,feasible=bool(peak<=peak_limit and rms<=rms_limit and power<=power_limit and rpm<=max_rpm),
                peak_limit_nm=peak_limit,continuous_limit_nm=rms_limit,power_limit_w=power_limit,
                maximum_stored_energy_j=.5*inertia*(np.pi*(fc+3*spacing))**2)


def prepare():
    source=ROOT/'fea/results/shaped_waveforms/medium_rounded.npz'
    rotor=FiniteRotor.from_npz(source); cache={}
    reference=json.loads((ROOT/'results/paper_revision/gaussian_benchmark.json').read_text())['rows']
    def field(d):
        key=f'{d:.12g}'
        if key not in cache:
            a128=rotor.harmonic_amplitude(d,samples=128)
            a256=rotor.harmonic_amplitude(d,samples=256)
            cache[key]=dict(distance_m=d,amplitude_m_s2=a256,phase_refinement_relative=abs(a128/a256-1))
            assert cache[key]['phase_refinement_relative']<1e-7
            for row in reference:
                if d==row['distance_m']:
                    np.testing.assert_allclose(a256,row['a2_m_s2'],rtol=1e-10,atol=0)
            print('FIELD',key,a256,flush=True)
        return cache[key]['amplitude_m_s2']
    cases=[]
    specs=[(f'desktop_{d:g}m',d,1.,FC,False) for d in (.5,.75,1.,1.5,2.)]
    specs += [('scale2_2m_same_receiver',2.,2.,FC,False),
              ('scale2_2m_slow_same_receiver',2.,2.,FC/2,False),
              ('scale2_2m_slow_retuned_receiver',2.,2.,FC/2,True)]
    for name,d,s,fc,retune in specs:
        a=s*field(d/s)
        rec=replace(StructuralOscillator(),resonance_hz=fc) if retune else StructuralOscillator()
        band=2.; frequencies=fc+band*((np.arange(65536)+.5)/65536-.5)
        cap=gaussian_waterfill(rec.acceleration_noise_asd(frequencies)**2,band,a*a/2)
        # PSD-min energy bound resolves ultra-small capacity without waterfill grid dependence.
        energy_bound=a*a/(2*float(np.min(rec.acceleration_noise_asd(frequencies)**2))*np.log(2))
        stress_factor=(s*fc/60)**2
        case=dict(name=name,distance_m=d,scale=s,mass_kg=MASS*s**3,diameter_m=.5*s,
                  inertia_kg_m2=INERTIA*s**5,carrier_hz=fc,receiver=asdict(rec),receiver_redesigned=retune,
                  amplitude_m_s2=a,gaussian_benchmark_bit_s=cap,wideband_energy_upper_bit_s=energy_bound,
                  benchmark_bandwidth_hz=band,centrifugal_stress_factor_vs_base1800=stress_factor,
                  illustrative_regional_mean_mpa=25.48*stress_factor,illustrative_unconverged_peak_mpa=58.13*stress_factor,
                  drive=dict(drag=.05*s**3,peak_limit=2*s**3,rms_limit=.5*s**3,power_limit=400*s**2,max_rpm=1800 if s==1 else 1800/s),
                  assumptions='No environmental noise; ideal acquired timing/phase; similarity transfers centrifugal statics only; no new FEA. Enlarged drive ratings scale explicitly, not fixed-hardware comparisons.')
        if s>1 and fc==FC:
            case['drive']['max_rpm']=1800
            case['mechanical_caution']='Same-speed scaling increases centrifugal stress; not a qualified operating design.'
        cases.append(case)
    np.testing.assert_allclose(duty(2.,.4,.1)['worst_valid_message_rms_nm'],.4474016719053217,rtol=2e-6)
    shifted=CenteredReceiver(replace(StructuralOscillator(),resonance_hz=FC/2),FC/2)
    offsets=np.fft.fftfreq(1024,1/2.)
    expected=4*shifted.receiver.acceleration_noise_asd(FC/2+offsets)**2/cases[0]['amplitude_m_s2']**2
    np.testing.assert_allclose(noise_psd(2.,cases[0]['amplitude_m_s2'],1024,receiver=shifted),expected,rtol=1e-12)
    save(OUT/'scenarios.json',dict(cases=cases,field_cache=cache,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
         provenance='Finite-volume similarity a_s(d)=s*a_1(d/s); CAD inertia used for mechanics, accepted quadrature for field.',
         checks='128/256 phase agreement; original 0.5/1/2 m field reproduction; baseline exact-message torque reproduction; physical carrier PSD translation.'))
    print(json.dumps([{k:c[k] for k in ('name','mass_kg','amplitude_m_s2','gaussian_benchmark_bit_s')} for c in cases]),flush=True)


def make_decoder(case, ts, q, samples=16, memory_symbols=3):
    model=CFFSK(float(ts),.4,float(q),samples)
    fs=samples/ts
    if 3*model.spacing_hz>=.8*fs/2: raise ValueError('outer tone outside guard band')
    physical=StructuralOscillator(**case['receiver']); rec=CenteredReceiver(physical,case['carrier_hz'])
    taps,fit=fit_whitener(fs,case['amplitude_m_s2'],memory_symbols*samples,receiver=rec)
    decoder=WhitenedViterbi(model,taps,maximum_states=1000000)
    return model,decoder,rec,fit


def trials(case,model,decoder,rec,start,count,seed,workers):
    a=case['amplitude_m_s2']; n=model.samples_per_symbol; prefix=decoder.memory_symbols*n
    case_seed=zlib.crc32(case['name'].encode())
    def one(index):
        rng=np.random.default_rng(np.random.SeedSequence([seed,case_seed,index]))
        data=rng.integers(0,256,8,dtype=np.uint8); symbols=encode(data)
        waveform=np.r_[np.ones(prefix,complex),model.waveform(np.r_[symbols,np.zeros(TAIL,int)])]
        waveform+=receiver_noise(rng,len(waveform),n/model.symbol_s,a,receiver=rec)
        observed=lfilter(decoder.taps,[1.],waveform)[prefix:]
        recovered=decode(decoder.decode_whitened(observed,TAIL,True)[:len(symbols)])
        return int(np.unpackbits(data^recovered).sum())
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(one,range(start,start+count)))


def candidate(case,ts,q,count,seed,workers,samples=16,memory_symbols=3):
    started=time.monotonic()
    mech=duty(ts,.4,q,inertia=case['inertia_kg_m2'],fc=case['carrier_hz'],**case['drive'])
    record=dict(case=case['name'],symbol_s=ts,q=q,transition_fraction=.4,mapped_bit_s=8/(3*ts),
                packet_payload_bit_s=64/(30*ts),mechanics=mech,status='mechanically_excluded')
    if not mech['feasible']: return record
    try: model,decoder,rec,fit=make_decoder(case,ts,q,samples,memory_symbols)
    except ValueError as exc:
        record.update(status='computationally_excluded',reason=str(exc));return record
    # Exact-byte, no-noise round-trip regression for every candidate.
    symbols=encode(np.array([0,255,42,251,128,7,99,201],dtype=np.uint8))
    np.testing.assert_array_equal(decoder.observe_and_decode(symbols,np.random.default_rng(1),case['amplitude_m_s2'],noise=False,tail_symbols=TAIL,byte_mapping=True),symbols)
    errors=trials(case,model,decoder,rec,0,count,seed,workers)
    record.update(status='screen_only',packets=count,bit_errors=sum(errors),errored_packets=int(np.count_nonzero(errors)),
                  observed_ber=sum(errors)/(64*count),observed_packet_error=float(np.mean(np.array(errors)>0)),
                  bit_errors_per_packet=errors,whitener=fit,states=decoder.state_count,samples_per_symbol=samples,
                  memory_symbols=memory_symbols,seed=seed,elapsed_s=time.monotonic()-started)
    return record


def screen(args,cases):
    path=OUT/'screen.json'; rows=json.loads(path.read_text()) if path.exists() else []
    done={(r['case'],r['symbol_s'],r['q']) for r in rows}
    for case in cases:
        cap=case['gaussian_benchmark_bit_s']
        if case['name']=='scale2_2m_slow_same_receiver': continue  # benchmark diagnoses receiver mismatch first
        if case['name']=='desktop_0.5m': periods=[1.75,2.,2.5,3.]
        else:
            base=max(2.,8/(3*cap)); periods=[float(f'{base*x:.6g}') for x in (1.,2.,4.,8.)]
        for ts in periods:
            for q in (.1,.25,.5,1.):
                if (case['name'],ts,q) in done:continue
                row=candidate(case,ts,q,args.packets,761093,args.workers)
                rows.append(row);save(path,rows)
                print(json.dumps({k:row[k] for k in ('case','symbol_s','q','status','observed_ber','elapsed_s') if k in row}),flush=True)


def validate(args,cases):
    rows=json.loads((OUT/'screen.json').read_text())
    for case in cases:
        eligible=[r for r in rows if r['case']==case['name'] and r.get('observed_ber',1)<=TARGET/2]
        if not eligible: continue
        eligible.sort(key=lambda r:(r['symbol_s'],r['observed_ber'],r['states']))
        if args.symbol_s is not None: eligible=[r for r in eligible if abs(r['symbol_s']-args.symbol_s)<1e-6]
        if args.q is not None: eligible=[r for r in eligible if abs(r['q']-args.q)<1e-10]
        if not eligible: raise ValueError('Requested frozen candidate was not eligible in the screen')
        row=eligible[0]; suffix=('_'+args.tag) if args.tag else ''
        path=OUT/(case['name']+suffix+'_validation.json')
        validation_seed=492871+zlib.crc32(args.tag.encode()) if args.tag else 492871
        prior=json.loads(path.read_text()) if path.exists() else None
        if prior and (prior['symbol_s']!=row['symbol_s'] or prior['q']!=row['q'] or prior['samples_per_symbol']!=args.samples or prior['memory_symbols']!=args.memory or prior.get('confidence_alpha',.05)!=args.alpha or prior['seed']!=validation_seed):
            raise ValueError('Validation configuration changed; use a distinct tag')
        if prior and prior.get('status')=='conditional_ber_target_pass': continue
        ts,q=row['symbol_s'],row['q'];model,decoder,rec,fit=make_decoder(case,ts,q,args.samples,args.memory)
        if prior and (prior['symbol_s']!=ts or prior['q']!=q or prior['samples_per_symbol']!=args.samples or prior['memory_symbols']!=args.memory):
            raise ValueError('Validation configuration changed; use a distinct output name')
        errors=prior['bit_errors_per_packet'] if prior else []
        started=time.monotonic();previous=prior['elapsed_s'] if prior else 0.
        for start in range(len(errors),args.max_packets,250):
            errors.extend(trials(case,model,decoder,rec,start,min(250,args.max_packets-start),validation_seed,args.workers))
            interval=packet_ber_interval(errors,alpha=args.alpha);passed=interval[1]<=TARGET;failed=interval[0]>TARGET
            record=dict(row,status='conditional_ber_target_pass' if passed else 'conditional_ber_target_fail' if failed else 'validation_inconclusive',
                        packets=len(errors),bit_errors=sum(errors),errored_packets=int(np.count_nonzero(errors)),
                        observed_ber=sum(errors)/(64*len(errors)),observed_packet_error=float(np.mean(np.array(errors)>0)),
                        anytime95_packet_aware_ber_interval=interval,bit_errors_per_packet=errors,
                        seed=validation_seed,confidence_alpha=args.alpha,samples_per_symbol=args.samples,memory_symbols=args.memory,whitener=fit,states=decoder.state_count,
                        elapsed_s=previous+time.monotonic()-started,
                        scope='Frozen candidate, independent validation stream; conditional acquired packets. Sample/filter convergence and achieved-drive check still required.')
            save(path,record)
            print(json.dumps({k:record[k] for k in ('case','packets','observed_ber','anytime95_packet_aware_ber_interval','status','elapsed_s')}),flush=True)
            if passed or failed:break


def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','screen','validate'])
    p.add_argument('--cases',nargs='*');p.add_argument('--packets',type=int,default=32);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--symbol-s',type=float);p.add_argument('--q',type=float);p.add_argument('--tag',default='');p.add_argument('--alpha',type=float,default=.05)
    p.add_argument('--max-packets',type=int,default=10000);p.add_argument('--samples',type=int,default=16);p.add_argument('--memory',type=int,default=3)
    args=p.parse_args()
    if args.tag and not args.tag.isidentifier():raise ValueError('tag must be an identifier')
    if not 0<args.alpha<1:raise ValueError('alpha must be in (0,1)')
    if args.mode=='prepare':prepare();return
    cases=json.loads((OUT/'scenarios.json').read_text())['cases']
    if args.cases: cases=[c for c in cases if c['name'] in args.cases]
    if not cases:raise ValueError('No matching cases')
    (screen if args.mode=='screen' else validate)(args,cases)

if __name__=='__main__':main()
