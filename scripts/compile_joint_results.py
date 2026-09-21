"""Compile audited joint-link screens; never label screening winners validated."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/joint_link_study_2026_09_20'
RUN=OUT/'exxact_crossed_v1'
LABEL={'uncoded':'Uncoded','ldpc':'LDPC 1/2','tb_1/2':'Tail-biting 1/2','tb_2/3':'Tail-biting 2/3','tb_3/4':'Tail-biting 3/4'}
rows=[]
for p in RUN.glob('case_*.json'):
    r=json.loads(p.read_text());c=r['config'];case=c['case']
    row=dict(id=c['id'],distance_m=case['distance_m'],receiver=case['receiver_id'],receiver_status=case['receiver_status'],
             resonance_hz=case['receiver']['resonance_hz'],proof_mass_g=1000*case['receiver']['mass_kg'],
             scheme=c['scheme'],information_bits=c['information_bits'],q=c['q'],symbol_s=c['symbol_s'],
             failures=sum(x['failure'] for x in r['records']),packets=r['packets'],ber=r['postdecode_ber'],
             per=r['packet_failure_rate'],per_ci95=r['packet_failure_ci95'],payload_bit_s=r['payload_bit_s'],
             observed_goodput_bit_s=r['observed_correct_payload_bit_s'],packet_minutes=r['packet_duration_s']/60,
             raw_ber=r['raw_ber'],accepted_fraction=r['accepted_packet_fraction'],wrong_accepted=r['wrong_accepted_packets'],
             detection=r['error_detection'],states=r['states'],simulation_s=r['elapsed_s'])
    rows.append(row)
rows.sort(key=lambda x:(x['distance_m'],x['resonance_hz'],x['q'],x['information_bits'],x['scheme']))
assert len(rows)==180 and sum(x['packets'] for x in rows)==9000
(OUT/'compiled_screen.json').write_text(json.dumps(rows,indent=2)+'\n')

def rate(x):return f'{x:.3f}' if x>=.01 else f'{x:.3g}'
def table(data,full=False):
    lines=['| Distance m | Receiver Hz | Code | Payload bits | q | Symbol s | BER after decoding | Failed packets | Payload bit/s | Oracle correct-payload bit/s | Packet min |',
           '|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in data:
        lines.append(f'| {r["distance_m"]:g} | {r["resonance_hz"]:g} | {LABEL[r["scheme"]]} | {r["information_bits"]} | {r["q"]} | {r["symbol_s"]:.4g} | {r["ber"]:.4g} | {r["failures"]}/{r["packets"]} | {rate(r["payload_bit_s"])} | {rate(r["observed_goodput_bit_s"])} | {r["packet_minutes"]:.2f} |')
    return '\n'.join(lines)

summary=[]
for d in (.5,1.,2.):
    group=[r for r in rows if r['distance_m']==d]
    clean=[r for r in group if r['failures']==0]
    bestclean=max(clean,key=lambda r:r['payload_bit_s']) if clean else None
    least=min(group,key=lambda r:(r['per'],-r['payload_bit_s']))
    summary.append(dict(distance_m=d,configurations=len(group),zero_failure_configurations=len(clean),
                        fastest_zero_failure_example=bestclean,lowest_observed_failure_example=least))
(OUT/'compiled_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

report=['# Compiled joint-link screening results','',
'Completed screen: 180 configurations, 50 packets each (9,000 total), run on Exxact and returned locally. An earlier local pilot comprised 150 configurations and 2,400 packets; it is kept separate because payload length and waveform were paired. All returned checkpoint counts and source hashes passed audit. Two noisy 50-packet configurations reproduced identical packet outcomes locally.', '',
'**These are coarse acquired-packet screens, not optimized or reliability-qualified links. No code family or receiver has been selected.**', '',
'## What was compared','',
'The same 7.9135 kg, 0.5 m diameter transmitter and common drive limits at 0.5, 1 and 2 m. Three coupled receiver/carrier choices (50.3, 40 and 25.15 Hz), two waveform settings, two information lengths (64/128 bits), and five schemes: uncoded, short LDPC rate 1/2, and tail-biting convolutional rates 1/2, 2/3 and 3/4. This is 3 × 3 × 2 × 2 × 5 = 180 configurations. All use the seven-tone byte mapper and six-symbol trailer; byte padding is counted.', '',
'The 40/25.15 Hz designs scale receiver dimensions, stiffness and mass together, to 6.164/24.8 g proof masses, and assume preserved Q and readout noise. They are conditional design sensitivities. The original 50.3 Hz reference has 3.1 g proof mass. All three distances admit the same maximum receiver resource envelope; the best entries nevertheless do not use identical physical receivers. Receiver and carrier are paired, and candidate symbol durations depend on receiver benchmark and mechanics. Therefore improvements cannot be attributed to receiver resonance alone. Q/4 variants were tested only in the earlier pilot.', '',
'Code rate and block length are crossed at each receiver/waveform point. However q=0.1 and q=0.25 use different durations: this screen cannot isolate a tone-spacing effect independently of duration. The Gaussian benchmark seeded duration scales; its ratio to tested rate is not a measured coding gap.', '',
'## Distance-level findings','',
'| Distance m | Configurations | Zero-failure configurations | Lowest observed packet failure |',
'|---:|---:|---:|---:|']
for s in summary:
    r=s['lowest_observed_failure_example'];report.append(f'| {s["distance_m"]:g} | {s["configurations"]} | {s["zero_failure_configurations"]} | {r["failures"]}/50 |')
report += ['',
'At 0.5 m all tested points had no errors. Within this grid, the fastest clean entries are uncoded 128-bit packets at 1.135 bit/s and 1.88 minutes, tied across receiver choices. This demonstrates overhead cost at already-clean operating points, not that coding is unnecessary at faster operating points.', '',
'At 1 m, 18 configurations had no packet failures. The fastest such entry is the conditional 25.15 Hz receiver with rate-3/4 tail-biting coding and 128 information bits: 0.274 bit/s, 7.79 minutes per packet. It is an exploratory candidate, not a validated winner.', '',
'At 2 m no configuration was failure-free. The lowest observed failure fraction is 4/50 (8%), for the conditional 25.15 Hz receiver, rate-1/2 tail-biting, 64 information bits: BER 0.01875, payload 0.003282 bit/s, and 325.00 minutes (5.42 hours) per packet. Its pointwise 95% packet-failure interval is 2.22–19.23%. Even this long-packet point misses the illustrative post-decoding BER target of 1e-3. Slow this region down or improve physically supported coupling/noise before calling it a reliable link.', '',
'## Matched code and block-length comparisons','',
'The table below holds receiver/carrier and waveform fixed within each distance: conditional 25.15 Hz receiver, q=0.25. This avoids comparing one code on a different waveform against another. Correct-payload throughput counts a packet as successful only if it is accepted and its complete information payload is correct; it uses simulator knowledge of the transmitted payload and is not operational goodput or an ARQ rate. Without a common CRC, the uncoded and convolutional receivers accept all outputs, including corrupted ones. BER includes tentative outputs of rejected packets as a diagnostic; rejection remains visible in packet failure.', '',
table([r for r in rows if r['resonance_hz']==25.15 and r['q']==.25 and r['distance_m'] in (1.,2.)]), '',
'At 1 m in this matched table all coded entries had zero failures; rate-3/4 tail-biting at 128 bits gives 0.274 bit/s compared with 0.193 bit/s for either rate-1/2 family. At 2 m with 64 bits, LDPC has BER/PER 0.0181/0.14 versus 0.0188/0.08 for rate-1/2 tail-biting; these small-sample counts do not establish a convincing winner. For rate-3/4 tail-biting at 2 m, increasing information length from 64 to 128 bits raises observed packet failure from 0.20 to 0.44. Changing 64 to 128 information bits reduces trailer overhead but lengthens the complete transmission. The observed failure changes must be considered alongside that gain; short screens cannot establish a block-length optimum. Native LDPC syndrome rejection and the convolutional decoder without a CRC have different error-detection behavior, so accepted-but-wrong counts are retained in compiled_screen.json. They are not comparable application-level undetected-error guarantees.', '',
'## Uncertainty and limits','',
'Each configuration has only 50 packets. Zero failures gives a two-sided pointwise 95% exact interval of 0–7.11%, not proof of packet failure <=1e-3. Intervals are not simultaneous across 180 searched configurations. Do not pool different operating points into a single reliability confidence claim. Bit errors are correlated within packets, and 50 × payload_bits must not be treated as that many independent BER trials.', '',
'The runs omit CRC/acquisition overhead, achieved actuator tracking, phase drift, environmental noise and receiver startup/dynamic-range qualification. Analytic drive checks passed; occupied waveform bandwidth and filter-memory robustness still need assessment. Multi-hour packet examples particularly depend on the stationary acquired-phase assumption. Polar coding and alternative mappings are not yet compared.', '',
'## Next calculations indicated by this screen','',
'1. At 0.5 m, refine faster mechanically feasible waveforms before paying for larger reliability runs at uniformly clean, conservative points.',
'2. At 1 m, refine around the promising q=0.25 region with matched code/length comparisons; include CRC-aided polar and common error-detection overhead.',
'3. At 2 m, bracket a slower low-failure operating region using rate-1/2 references, then test whether higher-rate alternatives improve reliable throughput. There is no zero-failure starting point in this crossed screen.',
'4. Cross duration and tone spacing independently in local refinement and revisit Q/readout sensitivity. Do not infer a monotonic preferred code rate or block length from only these points.',
'5. Only after freezing finalists, perform independent high-statistics validation and physical robustness checks.', '',
'## Artifacts','',
'- compiled_screen.json: all 180 flattened results and error-detection diagnostics.',
'- compiled_summary.json: distance summaries and representative candidate identities.',
'- SCREEN_TABLES.md: complete parameter/result table.',
'- screen_tradeoffs.png: rate/failure/latency screening plot.',
'- exxact_crossed_v1/: original manifests, runtime and all per-packet records.',
'- exxact_results_audit.json: returned-data and local/remote replay checks.',
'- Compile command: `C:/ProgramData/anaconda3/python.exe AIP_Advances_Submission/scripts/compile_joint_results.py`.', '',
'Original manuscript and numerical baseline artifacts remain unchanged.']
(OUT/'SCREEN_RESULTS.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
(OUT/'SCREEN_TABLES.md').write_text('# Complete crossed-screen table\n\nAll entries are 50-packet exploratory screens. See SCREEN_RESULTS.md for physical assumptions and uncertainty.\n\n'+table(rows)+'\n',encoding='utf-8')

colors=dict(zip(LABEL,['#555555','#1b9e77','#7570b3','#d95f02','#e7298a']))
fig,axes=plt.subplots(2,3,figsize=(13,7.4),layout='constrained')
for col,d in enumerate((.5,1.,2.)):
    for scheme in LABEL:
        for k,marker in [(64,'o'),(128,'^')]:
            r=[x for x in rows if x['distance_m']==d and x['scheme']==scheme and x['information_bits']==k]
            axes[0,col].scatter([x['payload_bit_s'] for x in r],[x['per'] for x in r],c=colors[scheme],marker=marker,s=36,alpha=.7)
            axes[1,col].scatter([x['packet_minutes'] for x in r],[x['observed_goodput_bit_s'] for x in r],c=colors[scheme],marker=marker,s=36,alpha=.7)
    axes[0,col].set(title=f'{d:g} m',xlabel='Payload rate including trailer (bit/s)',ylim=(-.035,1.04),xscale='log')
    axes[1,col].set(xlabel='Packet duration (minutes)',xscale='log')
    if d>.5:axes[1,col].set_yscale('symlog',linthresh=1e-5)
    axes[1,col].set_ylim(0,1.2*max(x['observed_goodput_bit_s'] for x in rows if x['distance_m']==d))
    for ax in axes[:,col]:ax.grid(alpha=.2)
axes[0,0].set_ylabel('Observed packet failure fraction')
axes[1,0].set_ylabel('Oracle correct-payload rate (bit/s)')
from matplotlib.lines import Line2D
handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=colors[s],label=LABEL[s]) for s in LABEL]
handles += [Line2D([0],[0],marker=m,color='none',markerfacecolor='gray',label=f'{k} information bits') for k,m in [(64,'o'),(128,'^')]]
fig.legend(handles=handles,loc='outside lower center',ncol=4,frameon=False)
fig.suptitle('Three-distance screen: 50 packets per configuration\nAll receiver/waveform candidates shown; zero observed failures is not reliability validation',fontsize=12)
fig.savefig(OUT/'screen_tradeoffs.png',dpi=170)
plt.close(fig)
print(json.dumps(summary,indent=2))