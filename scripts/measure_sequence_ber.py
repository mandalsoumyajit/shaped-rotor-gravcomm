"""Error-count-driven BER trials with byte-valid Viterbi paths.

Whole packets are independent. Bit errors within a packet are correlated.
The reported interval is an anytime-valid empirical-Bernstein bound on the
mean packet bit-error fraction, using a summable allocation across checkpoints.
Low-count runs are explicitly bounds, never accurate BER estimates.
"""
import argparse,json,time,sys
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
from multiprocessing import get_context
import numpy as np
from scipy.stats import beta
import desktop_cf_fsk as cf
from search_sequence_cf_fsk import A,OUT
from gravcomm.sequence_decoder import CFFSK
from gravcomm.whitened_viterbi import fit_whitener,WhitenedViterbi

parser=argparse.ArgumentParser()
parser.add_argument('--boundary',action='store_true')
parser.add_argument('--max-packets',type=int,default=10000)
parser.add_argument('--resume',action='store_true')
parser.add_argument('--workers',type=int,default=4)
parser.add_argument('--symbol-s',type=float)
parser.add_argument('--transition-fraction',type=float)
parser.add_argument('--spacing-product',type=float)
parser.add_argument('--output-name')
parser.add_argument('--legacy-prefix',type=int,default=0,
    help='Reproduce an initial sequential-RNG block from the first validation run; multiple of 250')
args=parser.parse_args()
if args.legacy_prefix%250:raise ValueError('legacy prefix must be a multiple of 250')
ts,r,q=(1.75,.4,.075) if args.boundary else (2.,.4,.1)
ts=args.symbol_s if args.symbol_s is not None else ts
r=args.transition_fraction if args.transition_fraction is not None else r
q=args.spacing_product if args.spacing_product is not None else q
n=round(2*ts);model=CFFSK(ts,r,q,n);taps,fit=fit_whitener(n/ts,A,int(np.ceil(6*n/ts)))
decoder=WhitenedViterbi(model,taps,maximum_states=1000000)
rng=np.random.default_rng(904317 if args.boundary else 937119)
bit_errors=[];symbol_errors=[];tail=6;started=time.monotonic();previous_elapsed=0.
filename='boundary_ber.json' if args.boundary else 'selected_ber.json'
if args.output_name:
    if not args.output_name.isidentifier():raise ValueError('output name must be an identifier')
    filename=args.output_name+'.json'
if args.resume:
    prior=json.loads((OUT/filename).read_text())
    bit_errors=prior['bit_errors_per_packet'];symbol_errors=prior['symbol_errors_per_packet']
    previous_elapsed=prior['elapsed_s']

def simulate(index,local_rng=None):
    # Independent deterministic streams make batching reproducible.
    if local_rng is None:
        local_rng=np.random.default_rng(np.random.SeedSequence([319477 if args.boundary else 178117,index]))
    data=local_rng.integers(0,256,8,dtype=np.uint8);symbols=cf.encode(data)
    decoded=decoder.observe_and_decode(symbols,local_rng,A,tail_symbols=tail,byte_mapping=True)
    recovered=cf.decode(decoded)
    return int(np.unpackbits(np.bitwise_xor(data,recovered)).sum()),int(np.count_nonzero(decoded!=symbols))

pool=(ProcessPoolExecutor(max_workers=args.workers,mp_context=get_context('fork'))
      if sys.platform.startswith('linux') else ThreadPoolExecutor(max_workers=args.workers))
for end in range(len(bit_errors)+250,args.max_packets+1,250):
    if end<=args.legacy_prefix:
        results=[simulate(index,rng) for index in range(len(bit_errors),end)]
    else:
        results=list(pool.map(simulate,range(len(bit_errors),end)))
    bit_errors.extend(x[0] for x in results);symbol_errors.extend(x[1] for x in results)
    packet=len(bit_errors)
    count=sum(bit_errors);failed=int(np.count_nonzero(bit_errors));j=packet//250
    fractions=np.array(bit_errors)/64
    # Empirical Bernstein for independent bounded packet fractions; union over
    # the scheduled looks makes coverage valid under the error-count stopping rule.
    delta=.05/(j*(j+1));logterm=np.log(4/delta)
    mean=float(fractions.mean());variance=float(fractions.var(ddof=1))
    margin=np.sqrt(2*variance*logterm/packet)+7*logterm/(3*(packet-1))
    # An independent uniformly chosen bit from each packet is Bernoulli with
    # probability BER. Conditional on the count, its error probability is count/64.
    # A separate fixed stream samples that bit without assuming within-packet
    # independence. Summable confidence spending also permits adaptive stopping.
    probe_rng=np.random.default_rng(719381)
    probes=int(np.sum(probe_rng.random(packet)<fractions))
    probe_upper=float(beta.ppf(1-delta,probes+1,packet-probes)) if probes<packet else 1.
    adequate=count>=100 and failed>=30
    record=dict(ts=ts,transition_fraction=r,spacing_dwell_product=q,spacing_hz=model.spacing_hz,
        nominal_bit_s=8/(3*ts),packet_payload_bit_s=64/(30*ts),packet_bytes=8,tail_symbols=tail,
        packets=packet,bits=64*packet,bit_errors=count,errored_packets=failed,
        observed_ber=mean,ber_estimate_has_required_error_count=adequate,
        anytime95_ber_interval=[max(0.,mean-margin),min(1.,mean+margin)],
        independent_bit_probe_errors=probes,anytime95_bit_probe_upper=probe_upper,
        symbol_errors=sum(symbol_errors),observed_ser=sum(symbol_errors)/(24*packet),
        states=decoder.state_count,whitener=fit,elapsed_s=previous_elapsed+time.monotonic()-started,
        legacy_prefix_packets=args.legacy_prefix,
        bit_errors_per_packet=bit_errors,symbol_errors_per_packet=symbol_errors)
    (OUT/filename).write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({key:record[key] for key in ('packets','bits','bit_errors','errored_packets','observed_ber','ber_estimate_has_required_error_count','elapsed_s')}),flush=True)
    if adequate:break
pool.shutdown()
