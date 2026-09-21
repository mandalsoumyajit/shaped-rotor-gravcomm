"""Reconcile paper rates; refine 0.5 m sizing; conditional legacy-front-end screen."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json,hashlib,binascii
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.signal import lfilter,welch
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'src')]
from range_revision import MASS,INERTIA,FC,duty,save,make_decoder
from gravcomm import FiniteRotor
from gravcomm.receivers import StructuralOscillator
from gravcomm.whitened_viterbi import receiver_noise
from gravcomm.packet_confidence import packet_ber_interval
from desktop_cf_fsk import encode,decode
from paper_results import gaussian_waterfill
from dataclasses import asdict
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/reconciliation_v1'
SOURCE=ROOT/'fea/results/shaped_waveforms/medium_rounded.npz'
PACKETS=200

def crc(data):return binascii.crc_hqx(bytes(data),0xffff)
def frame(payload):
    # CRC-16/CCITT-FALSE over header + payload; big endian, no final xor.
    body=b'\x00\x01'+bytes(payload)
    return np.frombuffer(body+crc(body).to_bytes(2,'big'),dtype=np.uint8)
def valid(data):return bytes(data[:2])==b'\x00\x01' and crc(data[:-2])==int.from_bytes(bytes(data[-2:]),'big')
def main():
    OUT.mkdir(exist_ok=True)
    assert crc(b'123456789')==0x29b1
    test=frame(bytes(range(32)));assert valid(test)
    broken=test.copy();broken[6]^=1;assert not valid(broken)
    rotor=FiniteRotor.from_npz(SOURCE);rec=StructuralOscillator()
    psd=rec.acceleration_noise_asd(FC+2*((np.arange(65536)+.5)/65536-.5))**2
    def progress(stage,done,total):
        record=dict(stage=stage,completed=done,total=total,packet_screen_packets=PACKETS,scope='Paper-compatible acquired phase/timing; not new fixed-band frontend validation')
        save(OUT/'progress.json',record)
        (OUT.parent/'PROGRESS.md').write_text('# Fixed-rate sizing progress\n\nOriginal deterministic pass: complete (35 fields, 510 mechanics cases).\n\nReconciliation: '+stage+f' {done}/{total}.\n\nFixed receiver noise profile and tunable center retained. New packet screen uses the paper-compatible sampled front end; a strict 2 Hz receive filter and acquisition remain unvalidated.\n',encoding='utf-8')
    scales=np.round(np.arange(.75,1.001,.025),3)
    def field(s):
        a64=s*rotor.harmonic_amplitude(.5/s,samples=64)
        a128=s*rotor.harmonic_amplitude(.5/s,samples=128)
        assert abs(a64/a128-1)<1e-7
        return dict(scale=float(s),diameter_m=.5*s,mass_kg=MASS*s**3,inertia_kg_m2=INERTIA*s**5,carrier_hz=FC/s,amplitude_m_s2=a128,gaussian_band_benchmark_bit_s=gaussian_waterfill(psd,2,a128*a128/2),phase_refinement_relative=abs(a64/a128-1))
    rows=[];progress('refined_fields',0,len(scales))
    with ThreadPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(field,float(s)) for s in scales]):
            row=f.result();rows.append(row);progress('refined_fields',len(rows),len(scales))
    rows.sort(key=lambda r:r['scale']);save(OUT/'fields.json',rows)
    a=next(r for r in rows if r['scale']==1)['amplitude_m_s2']
    np.testing.assert_allclose(a,5.180259430854157e-10,rtol=1e-10)
    paper=duty(2,.4,.1)
    np.testing.assert_allclose([paper['peak_torque_nm'],paper['worst_valid_message_rms_nm']],[1.2136335596129388,.4474016719053107],rtol=1e-10)
    rates=[]
    for payload in (64,128,224,256,280):
        total=3*((payload+32)//8)+18
        for ts in (2.,payload/total):
            rates.append(dict(payload_bits=payload,header_crc_bits=32,acquisition_symbols=12,tail_symbols=6,symbol_s=ts,airtime_s=ts*total,payload_bit_s=payload/(ts*total),mechanics=duty(ts,.4,.1)))
    save(OUT/'rate_accounting.json',dict(paper=dict(mapped_bit_s=8/6,packet_bit_s=64/60,packet_s=60,acquisition_symbols=0,header_crc_bits=0,mechanics=paper),framed=rates))
    # Each noise profile is reference PSD versus offset, independent of tuned carrier.
    # Use original 50.3 Hz numerical reference with unchanged StructuralOscillator.
    # Physical mechanical carrier is used only by duty, never to change receiver PSD.
    configs=[dict(kind='paper_reference',scale=1.,payload_bits=64,ts=2.,framed=False)]
    for s in (.8,.85,.9,.95,1.):configs.append(dict(kind='framed_1bps',scale=s,payload_bits=256,ts=256/126,framed=True))
    save(OUT/'manifest.json',dict(configs=configs,packets=PACKETS,seed=20370924,source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitations='Legacy sampled colored-noise front end, no explicit 2 Hz receive filter. Fixed PSD versus frequency offset. Initial phase/timing known. CRC and header transmitted, acquisition budget counted but acquisition detection not simulated. No outer error-correcting code in this reference comparison.'))
    results=[]
    for index,c in enumerate(configs):
        path=OUT/f'packet_{index}.json'
        if path.exists():results.append(json.loads(path.read_text()));continue
        row=next(r for r in rows if abs(r['scale']-c['scale'])<1e-8)
        case=dict(amplitude_m_s2=row['amplitude_m_s2'],receiver=asdict(rec),carrier_hz=FC)
        model,det,physical,fit=make_decoder(case,c['ts'],.1,16,3)
        prefix=det.memory_symbols*16
        def one(i,noiseless=False):
            rng=np.random.default_rng(np.random.SeedSequence([20370924,c['payload_bits'],i]))
            payload=rng.integers(0,256,c['payload_bits']//8,dtype=np.uint8)
            data=frame(payload) if c['framed'] else payload
            symbols=encode(data);wave=np.r_[np.ones(prefix,complex),model.waveform(np.r_[symbols,np.zeros(6,int)])]
            if not noiseless:wave+=receiver_noise(rng,len(wave),16/c['ts'],row['amplitude_m_s2'],receiver=physical)
            observed=lfilter(det.taps,[1],wave)[prefix:]
            recovered=decode(det.decode_whitened(observed,6,True)[:len(symbols)])
            recovered_payload=recovered[2:-2] if c['framed'] else recovered
            errors=int(np.unpackbits(recovered_payload^payload).sum())
            accepted=valid(recovered) if c['framed'] else True
            return dict(errors=errors,accepted=accepted,wrong_accepted=bool(errors and accepted),packet_failure=bool(errors or not accepted))
        assert one(0,True)['errors']==0
        progress('packet_screen',index,len(configs))
        with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(one,range(PACKETS)))
        counts=[r['errors'] for r in records];accepted=[r for r in records if r['accepted']]
        drive=duty(c['ts'],.4,.1,inertia=row['inertia_kg_m2'],fc=row['carrier_hz'],drag=.05*row['scale']**3,max_rpm=1800/row['scale'])
        # Stationary random mapped waveform: spectral diagnostic only, no claim of occupied-bandwidth compliance.
        rng=np.random.default_rng(20370925);long=model.waveform(encode(rng.integers(0,256,2048)))
        f,p=welch(long,fs=16/c['ts'],nperseg=8192,return_onesided=False,detrend=False)
        spill=float(p[np.abs(f)>1].sum()/p.sum())
        result=dict(config=c,source=row,packets=PACKETS,records=records,payload_ber=sum(counts)/(PACKETS*c['payload_bits']),packet_aware_ber_ci95=packet_ber_interval(counts,c['payload_bits']),packet_failure_rate=sum(r['packet_failure'] for r in records)/PACKETS,accepted_packets=len(accepted),accepted_payload_ber=sum(r['errors'] for r in accepted)/(len(accepted)*c['payload_bits']) if accepted else None,wrong_accepted_packets=sum(r['wrong_accepted'] for r in records),nominal_payload_bit_s=1. if c['framed'] else 64/60,oracle_correct_packet_payload_bit_s=(1. if c['framed'] else 64/60)*(1-sum(r['packet_failure'] for r in records)/PACKETS),airtime_s=256 if c['framed'] else 60,mechanics=drive,whitener=fit,signal_fraction_outside_2hz=spill,scope='Exploratory paper-compatible acquired-packet screen; not strict-band or acquisition validation. BER counts tentative payloads including rejected frames. CRC rejections also reported.')
        save(path,result);results.append(result)
        print(json.dumps(dict(case=index,diameter=row['diameter_m'],ber=result['payload_ber'],per=result['packet_failure_rate'],rms=drive['worst_valid_message_rms_nm'])),flush=True)
    save(OUT/'summary.json',dict(fields=rows,rate_accounting=rates,packet_results=results))
    lines=['# Paper comparison and refined 0.5 m screen','','The reference rotor is unchanged: 0.50 m diameter and 7.91 kg. The paper field and peak/RMS torque were reproduced numerically. Its 1.333 bit/s mapped rate becomes 1.067 bit/s after the six-symbol trailer; neither includes acquisition or CRC/header.','','## Consistent framing','','With 256 useful bits, 16 header bits, CRC-16, 12 budgeted acquisition symbols and six trailing symbols, the frame has 126 symbols. At the original 2 s/symbol it lasts 252 s and carries 1.01587 useful bit/s. At 256/126 s/symbol it lasts 256 s and carries exactly 1 useful bit/s. This meets airtime, not total startup/settling or demonstrated acquisition latency. The earlier 64/128/256 information-block restriction unnecessarily excluded this packet length.','','## Refined finite-volume field grid','','| Diameter (m) | Mass (kg) | Gaussian benchmark (bit/s) |','|---:|---:|---:|']
    for r in rows:lines.append(f"| {r['diameter_m']:.4f} | {r['mass_kg']:.2f} | {r['gaussian_band_benchmark_bit_s']:.4f} |")
    lines+=['','## Acquired-packet diagnostic','','No outer error-correcting code is used in this baseline reconciliation. Fixed reference receiver PSD versus offset; center tuning does not change noise. New frames transmit and check CRC-16/CCITT-FALSE and a fixed 16-bit header, but reserve rather than simulate acquisition. The front end remains the paper-compatible sampled-noise model, not a new hard 2 Hz receive filter. These results must not be promoted to validated fixed-band links.','','| Case | Diameter (m) | Payload BER | Packet-aware 95% BER interval | Packet failure | RMS torque (Nm) |','|---|---:|---:|---|---:|---:|']
    for r in results:
        lo,hi=r['packet_aware_ber_ci95'];lines.append(f"| {r['config']['kind']} | {r['source']['diameter_m']:.3f} | {r['payload_ber']:.5g} | [{lo:.3g}, {hi:.3g}] | {r['packet_failure_rate']:.3g} | {r['mechanics']['worst_valid_message_rms_nm']:.3f} |")
    lines+=['','BER counts all tentative decoded payloads, including CRC-rejected frames. Correct-packet throughput is an oracle diagnostic, not operational goodput. CRC passing with no observed errors is not an undetected-error guarantee. 200 packets per point only screen candidates; packet-aware confidence bounds govern interpretation. The paper rare-error result is supported by its existing much larger validation, not reproduced statistically by this short diagnostic.','','Next: use variable payload lengths for uncoded/convolutional candidates; retain supported block lengths for each LDPC construction. Apply a consistent finite-band receiver to signal and noise and test detector convergence before selecting hardware. Refine near any observed BER transition; do not infer minimum transmitter size from Gaussian capacity alone.']
    (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    progress('reconciliation_complete',len(configs),len(configs));print('RECONCILIATION_COMPLETE',flush=True)
if __name__=='__main__':main()
