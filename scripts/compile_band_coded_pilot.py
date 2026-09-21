"""Compile exact-band pilot with explicit convergence and reliability limits."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'src')]
from range_revision import duty,save
from gravcomm.packet_confidence import packet_ber_interval
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/band_coded_pilot_v1'

def main():
    manifest=json.loads((OUT/'manifest.json').read_text());summaries=[]
    for c in manifest['configs']:
        files=[OUT/f'{c["id"]}_{i}.json' for i in range(manifest['packets'])]
        rows=[json.loads(p.read_text()) for p in files if p.exists()]
        if not rows:continue
        assert all(r['config']==c for r in rows)
        n=rows[0]['payload_bits'];ts=rows[0]['symbol_s'];s=c['source']['scale']
        tests=[t for r in rows for t in r['convergence']]
        mechanics=duty(ts,.4,.1,inertia=c['source']['inertia_kg_m2'],fc=c['source']['carrier_hz'],drag=.05*s**3,max_rpm=1800/s)
        summaries.append(dict(id=c['id'],diameter_m=c['source']['diameter_m'],scheme=c['scheme'],packets=len(rows),payload_bits=n,ber=sum(r['errors'] for r in rows)/(n*len(rows)),ber_ci95=packet_ber_interval([r['errors'] for r in rows],n),packet_failures=sum(r['failure'] for r in rows),wrong_accepted=sum(r['wrong_accepted'] for r in rows),mechanics=mechanics,convergence=tests,passes_initial_convergence_check=len(tests)==2 and all(t['relative_llr_change']<.01 and t['llr_sign_changes']==0 and t['decoded_bits_changed']==0 and t['accepted']==next(r for r in rows if r['index']==0)['accepted'] for t in tests),max_memory_llr_change=max((t['relative_llr_change'] for t in tests if t['memory']==19),default=None),max_width_llr_change=max((t['relative_llr_change'] for t in tests if t['width']==32768),default=None),compute_seconds=sum(r['elapsed_s'] for r in rows)))
    complete=sum(r['packets'] for r in summaries)==len(manifest['configs'])*manifest['packets']
    save(OUT/'compiled.json',dict(complete=complete,rows=summaries,scope=manifest['scope']))
    closure=' This diagnostic run was intentionally closed after the numerical gate failed; see termination.json. No jobs remain running.' if (OUT/'termination.json').exists() else ''
    lines=['# Exact-band coded pilot','','The same ideal 2 Hz projection and fixed receiver noise spectrum are applied to every waveform. Receiver center tracks rotor speed. This pilot tests numerical integration and detector convergence; eight packets per point cannot establish BER <= 1e-3. Initial timing/phase are assumed known. Acquisition time is reserved but not demonstrated.','',f'Run complete: {complete}.'+closure,'', '| Diameter (m) | Code | Packets | Payload BER | Failed packets | Peak torque (Nm) | RMS torque (Nm) | Within reference drive |','|---:|---|---:|---:|---:|---:|---:|---|']
    for r in summaries:
        m=r['mechanics'];lines.append(f"| {r['diameter_m']:.3f} | {r['scheme']} | {r['packets']} | {r['ber']:.5g} | {r['packet_failures']} | {m['peak_torque_nm']:.3f} | {m['worst_valid_message_rms_nm']:.3f} | {m['feasible']} |")
    lines+=['','All configurations carry 1 nominal useful bit/s including header, CRC, acquisition allowance and trailer. Airtime is 256 s for uncoded/rate-3/4/rate-2/3 and 224 s for the two rate-1/2 examples. This is one length seed per code, not a length-optimized comparison. Reference drive is 2 Nm peak, 0.5 Nm sustained RMS, 400 W, with the reference centrifugal speed scale. BER includes payload errors in rejected packets; packet failure includes rejection and any payload error. CRC wrong-acceptance counts are retained in compiled.json.','','## Full-packet numerical convergence','','The first packet at each point is repeated with memory 19 instead of 18, and width 32768 instead of 8192. The initial screening check is <1% relative LLR change with unchanged LLR signs, decoded frame and acceptance. This is a numerical screening tolerance, not a BER confidence guarantee.','','| Diameter (m) | Code | Memory change | Width change | Initial check passed |','|---:|---|---:|---:|---|']
    for r in summaries:
        a=r['max_memory_llr_change'];b=r['max_width_llr_change'];sa='pending' if a is None else f'{100*a:.3g}%';sb='pending' if b is None else f'{100*b:.3g}%'
        lines.append(f"| {r['diameter_m']:.3f} | {r['scheme']} | {sa} | {sb} | {r['passes_initial_convergence_check']} |")
    lines+=['','Independent review confirmed common-band normalization and beam bookkeeping, but the cropped projected-noise covariance is not identity. The Euclidean sequence metric is mismatched rather than exact finite-window maximum likelihood. Hard band truncation also has an infinite, noncausal impulse response; a realizable finite-delay front end must be specified and checked before latency claims. These are empirical results for an approximate idealized receiver.', '', 'No best code or minimum reliable transmitter is selected from this pilot. If numerical checks fail, address search/metric convergence before interpreting the corresponding code comparison. Even passing points require more packets, multiple waveform and length choices, acquisition and drive validation. Research decoding runtime is not yet an implemented real-time receiver latency.']
    (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(complete=complete,completed_packets=sum(r['packets'] for r in summaries),points=len(summaries),initial_convergence_passes=sum(r['passes_initial_convergence_check'] for r in summaries))))
if __name__=='__main__':main()

