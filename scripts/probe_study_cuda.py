"""Bounded, reproducible CUDA probe for the current coded sequence detector.

Measures real finite-whitener template correlations, retaining complex128.
GPU kernel timing and host-transfer-inclusive timing are separately reported.
It does not implement or claim GPU acceleration of the max-log trellis.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name]='1'
os.environ['MKL_THREADING_LAYER']='SEQUENTIAL'
import sys,json,time,statistics,subprocess,platform
from pathlib import Path
import torch # Load its OpenMP runtime before NumPy to avoid duplicate initialization.
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from range_revision import make_decoder
from gravcomm.coded_sequence import MaxLogByteDetector,ShortLDPC,_maxlog
OUT=ROOT/'results/joint_link_study_2026_09_20'

def timed(fn,repeats=7):
    times=[]
    for _ in range(repeats):
        t=time.perf_counter();fn();times.append(time.perf_counter()-t)
    return dict(median_s=statistics.median(times),min_s=min(times),repeats=repeats)

def main():
    import torch
    from numba import cuda
    torch.set_num_threads(1)
    report=dict(timestamp=time.strftime('%Y-%m-%dT%H:%M:%S%z'),python=sys.version,
        platform=platform.platform(),torch_version=torch.__version__,torch_cuda_build=torch.version.cuda,
        torch_cuda_available=torch.cuda.is_available(),numba_cuda_available=cuda.is_available(),
        cuda_path=os.environ.get('CUDA_PATH'),precision='complex128; no TF32',cpu_blas_threads=1)
    report['nvidia_smi']=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.used','--format=csv'],text=True).strip()
    if not torch.cuda.is_available():
        OUT.mkdir(parents=True,exist_ok=True);(OUT/'cuda_probe.json').write_text(json.dumps(report,indent=2));return
    report['gpu']=torch.cuda.get_device_name(0)
    report['compute_capability']=torch.cuda.get_device_capability(0)
    cases=json.loads((ROOT/'results/range_revision_2026_09_20/scenarios.json').read_text())['cases']
    case=next(c for c in cases if c['name']=='scale2_2m_slow_retuned_receiver')
    model,d,receiver,fit=make_decoder(case,5.4,.25)
    soft=MaxLogByteDetector(d)
    rng=np.random.default_rng(8821);n=model.samples_per_symbol;length=54
    # Same shape as 64-information-bit, rate-1/2 LDPC and six trailer symbols.
    blocks=(rng.normal(size=(length,n))+1j*rng.normal(size=(length,n))).astype(np.complex128)
    corr=blocks.conj()@d.templates.T
    def core(c):
        return _maxlog(c,d.energy,d.phase_factors,soft.destinations,soft.edges,
                      soft.phase_indices,soft.byte_values,d.initial_history,length-6,3)
    core(corr) # compile before timing
    ldpc=ShortLDPC(64);llr=np.full(128,2.);ldpc.decode(llr)
    report['workload']=dict(scenario=case['name'],symbol_s=5.4,q=.25,packet_symbols=length,
        samples_per_symbol=n,templates=list(d.templates.shape),trellis_states=d.state_count,
        note='Detector kernel profiling uses deterministic random input, not BER estimation.')
    report['cpu_profile']={
        'correlation':timed(lambda: blocks.conj()@d.templates.T,20),
        'maxlog_after_correlation':timed(lambda:core(corr),5),
        'complete_soft_detector':timed(lambda:soft.decode(blocks.ravel()),5),
        'ldpc_easy_one_iteration':timed(lambda:ldpc.decode(llr),20)}
    templates_gpu=torch.as_tensor(d.templates.T.copy(),device='cuda')
    results=[]
    for batch in (1,16,64):
        a=np.tile(blocks,(batch,1,1))
        ag=torch.as_tensor(a,device='cuda')
        cpu=a.conj()@d.templates.T
        def gpu_resident():
            result=ag.conj()@templates_gpu;torch.cuda.synchronize();return result
        def gpu_transfers():
            current=torch.as_tensor(a,device='cuda')
            result=(current.conj()@templates_gpu).cpu().numpy()
            torch.cuda.synchronize();return result
        got=gpu_transfers();gpu_resident()
        results.append(dict(batch=batch,cpu=timed(lambda:a.conj()@d.templates.T,10),
            gpu_resident=timed(gpu_resident,10),gpu_input_output_transfers=timed(gpu_transfers,10),
            max_abs_error=float(np.max(np.abs(cpu-got))),
            relative_l2_error=float(np.linalg.norm(cpu-got)/np.linalg.norm(cpu))))
    report['correlations']=results
    got=(torch.as_tensor(blocks,device='cuda').conj()@templates_gpu).cpu().numpy()
    reference=core(corr);alternative=core(got)
    report['detector_metric_comparison']=dict(max_abs_llr_error=float(np.max(abs(reference-alternative))),
        identical_llr_signs=bool(np.array_equal(reference<0,alternative<0)))
    # Actual hybrid detector for one packet, including transfers and CPU trellis.
    def hybrid():
        c=(torch.as_tensor(blocks,device='cuda').conj()@templates_gpu).cpu().numpy()
        return core(c)
    hybrid();report['hybrid_gpu_correlation_cpu_trellis']=timed(hybrid,5)
    p=report['cpu_profile'];fraction=p['correlation']['median_s']/p['complete_soft_detector']['median_s']
    report['correlation_fraction_of_cpu_detector']=fraction
    report['ideal_speedup_if_correlation_free']=1/(1-fraction)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'cuda_probe.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
