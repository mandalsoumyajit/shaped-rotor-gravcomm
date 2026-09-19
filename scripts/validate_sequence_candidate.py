"""Independent fixed-candidate validation, after the parameter screen."""
import argparse,json,time
from search_sequence_cf_fsk import evaluate,OUT

parser=argparse.ArgumentParser()
parser.add_argument('--selected',action='store_true')
args=parser.parse_args()
r,q=(.4,.1) if args.selected else (.5,.0625)
filename='selected_candidate.json' if args.selected else 'final_candidate.json'
cases=[('final',5000,719351 if args.selected else 604192,6.,2.),('longer_filter',200,918271,8.,2.),
       ('double_sample_rate',200,918271,6.,4.)]
records=[]
for name,packets,seed,duration,fs in cases:
    start=time.monotonic()
    print('Starting '+name,flush=True)
    result=evaluate(2.,r,q,packets,seed,duration,fs)
    result.update(name=name,seed=seed,elapsed_s=time.monotonic()-start)
    records.append(result)
    (OUT/filename).write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('name','elapsed_s','trial','whitener')}),flush=True)
