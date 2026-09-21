"""Integration checks, not BER qualification or code ranking."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import sys,time,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.causal_receiver import CausalReceiver,sampled_cffsk
from gravcomm.causal_sequence import CausalBeamByteDetector
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/causal_sequence_v1';OUT.mkdir(exist_ok=True)
A=5.180259430854157e-10
r=CausalReceiver();rng=np.random.default_rng(20260921);rows=[]
def run(nbytes,duration,widths):
 payload=rng.integers(0,256,nbytes,dtype=np.uint8);v=payload.astype(int)
 symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3
 signal=sampled_cffsk(np.r_[symbols,np.zeros(6,int)],duration)
 mean=r.stream(1).process(signal);y=mean+r.stationary_noise(len(mean),rng,A)
 local=[]
 for width in widths:
  t=time.perf_counter();got=CausalBeamByteDetector(r,duration,A,width).decode(y,nbytes)
  elapsed=time.perf_counter()-t
  residual=y-r.stream(1).process(sampled_cffsk(np.r_[got['symbols'],np.zeros(6,int)],duration))
  direct=r.finite_record(len(y),A).metric(residual,np.zeros(len(y)))
  error=abs(float(direct)-got['path_cost']);assert error<1e-7
  row=dict(nbytes=nbytes,symbol_s=duration,observation_s=len(y)/4,width=width,runtime_s=elapsed,exact_search=got['exact_search'],path_cost=got['path_cost'],direct_metric_error=error,bit_errors=int(np.count_nonzero(np.unpackbits(payload)^np.unpackbits(got['bytes']))),missing_bit_hypotheses=got['missing_bit_hypotheses'],llr=got['llr'].tolist(),decoded_bytes=got['bytes'].tolist())
  local.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ('llr','decoded_bytes')}),flush=True)
 reference=local[-1]
 for row in local:
  a=np.asarray(row['llr']);b=np.asarray(reference['llr']);valid=np.isfinite(a)&np.isfinite(b)
  row['max_llr_difference_to_largest_beam']=float(max(abs(a[valid]-b[valid]))) if valid.any() else None
  row['cost_gap_to_largest_beam']=row['path_cost']-reference['path_cost']
 rows.extend(local)
run(2,.5,(64,512,65536))
# 40 bytes = 120 mapped symbols plus 6 tail, totaling 256 seconds.
# This is a runtime/metric diagnostic, not the finalized framed payload protocol.
run(40,256/126,(32,128,512))
result=dict(status='integrated_exact_history_approximate_beam_search_not_BER_qualification',runs=rows,source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('src/gravcomm/causal_sequence.py','src/gravcomm/causal_receiver.py','scripts/validate_causal_sequence.py','tests/test_causal_sequence.py','src/gravcomm/receivers.py','src/gravcomm/constants.py')})
# Missing soft alternatives are null in JSON, explicitly counted above.
for row in rows:row['llr']=[x if np.isfinite(x) else None for x in row['llr']]
(OUT/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
lines=['# Causal sequence detector integration','','The detector preserves full phase, causal SOS filter state, and dense innovation history per candidate. There is no finite-memory state merging. Beam pruning is explicitly approximate. Positive max-log LLR favors bit zero; missing alternatives are NaN in the API and null in saved JSON, never fabricated confidence.','','Timing, phase, carrier, signal amplitude and deterministic prehistory are supplied inputs. Random receiver noise remains stationary; its record covariance is marginal, not conditioned on previous observations. Acquisition and actuator tracking are not implemented.','','Four integration tests pass: forced-opposite soft reference consistency; exhaustive 256-message hard/soft/cost comparison at integer and fractional symbol durations; explicit pruning and forced-bit constraints; supplied prepacket filter state and nonzero decimation-clock phase. The seven front-end tests also pass.','','## Computational diagnostics','','| Bytes | Beam | Exact search | Runtime (s) | Bit errors | Missing bit alternatives | Cost gap to largest beam |','|---:|---:|:---:|---:|---:|---:|---:|']
for row in rows:lines.append(f"| {row['nbytes']} | {row['width']} | {row['exact_search']} | {row['runtime_s']:.3f} | {row['bit_errors']} | {row['missing_bit_hypotheses']} | {row['cost_gap_to_largest_beam']:.6g} |")
lines+=['','Two-byte diagnostics include all 65,536 complete messages at the largest beam. The 40-byte case spans 256 seconds and measures packet-length runtime and metric consistency; it is not the finalized frame or a BER estimate. For every decoded sequence, the incremental path cost is checked independently against a whole-record covariance likelihood. Only one noisy record per length is used.','','The supplied mean state and absolute input sample offset permit continuation of the receiver filter without an artificial packet reset. Unknown initial state, frequency and timing must be acquired or jointly inferred. Soft outputs after pruning are restricted-list approximations, even when both bit alternatives survive. Wider-beam agreement on this diagnostic alone does not certify reliability.','','An optional decode_soft reference performs one forced-opposite search per bit and pools the retained candidates. It supplies finite, internally consistent approximate max-log outputs and matches exhaustive outputs when unpruned. This costs 1 + 8 times the byte count in searches; it is a validation reference, not yet an efficient production soft detector.\n\nNext: establish beam/soft-output convergence over representative noise records and amplitudes, then integrate code decoders with explicit handling of absent alternatives. Do not use this pilot to rank codes or claim the post-decoding BER target.']
(OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
