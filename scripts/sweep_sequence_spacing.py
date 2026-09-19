"""Receiver-only spacing sweep with separately evaluated transmitter demands.

Fixed 1.75 s symbols, 40% quintic transitions, eight payload bytes, six trailer
symbols, acquired phase/timing and no added ECC. All points use the same 8 Hz
complex sample rate. No point is rejected because of motor demand.
"""
import argparse,itertools,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import numpy as np
import desktop_cf_fsk as cf
from search_sequence_cf_fsk import A
from optimize_cf_fsk import I
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import fit_whitener,WhitenedViterbi
from gravcomm.packet_confidence import packet_ber_interval
from gravcomm.receivers import CARTER_OSCILLATOR

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/spacing_sweep';OUT.mkdir(exist_ok=True)
TS=1.75;R=.4;TAIL=6

def mechanics(q):
    t=np.linspace(0,TS,10001);u=np.minimum(t/(R*TS),1)
    smooth=10*u**3-15*u**4+6*u**5
    slope=(30*u**2-60*u**3+30*u**4)/(R*TS)
    spacing=q/((1-R)*TS);cost=np.zeros((7,7));peak=power=0.
    for a,b in itertools.product(range(7),repeat=2):
        f=cf.FC+spacing*((a-3)+(b-a)*smooth)
        tau=I*np.pi*spacing*(b-a)*slope+.05*f/cf.FC
        cost[a,b]=np.trapz(tau*tau,t)/TS
        peak=max(peak,float(max(abs(tau))))
        power=max(power,float(max(abs(tau*np.pi*f))))
    words=np.array([[b//49,(b//7)%7,b%7] for b in range(256)])
    edges=np.full((7,7),-np.inf)
    for old in range(7):
        for a,b,c in words:
            edges[old,c]=max(edges[old,c],cost[old,a]+cost[a,b]+cost[b,c])
    worst=0.
    for length in range(1,8):
        for nodes in itertools.combinations(range(7),length):
            for rest in itertools.permutations(nodes[1:]):
                cycle=(nodes[0],)+rest
                worst=max(worst,sum(edges[a,b] for a,b in zip(cycle,cycle[1:]+cycle[:1]))/(3*length))
    first=np.bincount(words[:,0],minlength=7)/256
    last=np.bincount(words[:,2],minlength=7)/256
    random=(last@cost@first+np.mean(cost[words[:,0],words[:,1]]+cost[words[:,1],words[:,2]]))/3
    return dict(peak_torque_nm=peak,peak_power_w=power,worst_message_rms_nm=float(np.sqrt(worst)),
        random_message_rms_nm=float(np.sqrt(random)),max_rpm=30*(cf.FC+3*spacing))

def simulate(index):
    rng=np.random.default_rng(np.random.SeedSequence([673129,int(round(Q*1000000)),index]))
    data=rng.integers(0,256,8,dtype=np.uint8);symbols=cf.encode(data)
    decoded=DECODER.observe_and_decode(symbols,rng,A,tail_symbols=TAIL,byte_mapping=True)
    return int(np.unpackbits(data^cf.decode(decoded)).sum())

parser=argparse.ArgumentParser();parser.add_argument('--sample-rate',type=float,default=8.)
parser.add_argument('--max-packets',type=int,default=10000)
parser.add_argument('--q',type=float,nargs='+',default=[.05,.075,.1,.15,.2,.3,.5,.75,1.])
parser.add_argument('--workers',type=int,default=4)
parser.add_argument('--memory-symbols',type=int,default=3)
args=parser.parse_args()

for Q in args.q:
    name=f'q{Q:g}_fs{args.sample_rate:g}_m{args.memory_symbols}'.replace('.','p')
    path=OUT/(name+'.json')
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('complete'):continue
    n=round(TS*args.sample_rate);fs=n/TS
    model=CFFSK(TS,R,Q,n)
    if 3*model.spacing_hz>=.8*fs/2:raise ValueError('Increase sample rate to retain outer-tone guard band')
    taps,fit=fit_whitener(fs,A,args.memory_symbols*n)
    DECODER=WhitenedViterbi(model,taps,maximum_states=1000000)
    duty=mechanics(Q);errors=[];started=time.monotonic()
    tones=cf.FC+np.arange(-3,4)*model.spacing_hz
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=get_context('fork')) as pool:
        for end in range(250,args.max_packets+1,250):
            errors.extend(pool.map(simulate,range(len(errors),end)))
            count=sum(errors);failed=int(np.count_nonzero(errors));adequate=count>=100 and failed>=30
            interval=packet_ber_interval(errors)
            complete=adequate or end==args.max_packets
            record=dict(symbol_s=TS,transition_fraction=R,spacing_dwell_product=Q,
                tone_spacing_hz=model.spacing_hz,sample_rate_hz=fs,states=DECODER.state_count,
                source_amplitude_m_s2=A,whitener=fit,mechanics=duty,
                tone_noise_asd=CARTER_OSCILLATOR.acceleration_noise_asd(tones).tolist(),
                packets=end,payload_bits=end*64,bit_errors=count,errored_packets=failed,
                observed_ber=count/(end*64),packet_aware_anytime95_ber_interval=interval,
                adequate_error_count=adequate,complete=complete,bit_errors_per_packet=errors,
                elapsed_s=time.monotonic()-started)
            path.write_text(json.dumps(record,indent=2)+'\n')
            print(json.dumps({key:record[key] for key in ('tone_spacing_hz','packets','bit_errors','errored_packets','observed_ber','adequate_error_count','complete','elapsed_s')}),flush=True)
            if complete:break
