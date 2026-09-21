"""Summarize initial range calculations without promoting screens to results."""
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from range_revision import ROOT,OUT,save
from gravcomm.packet_confidence import packet_ber_interval


def main():
    scenarios=json.loads((OUT/'scenarios.json').read_text())
    cases=scenarios['cases'];screen=json.loads((OUT/'screen.json').read_text())
    all_trials=[];selected={}
    for c in cases:
        trials=[]
        for path in OUT.glob('*_validation.json'):
            r=json.loads(path.read_text())
            if r['case']==c['name']:trials.append((path,r))
        count=len(trials)
        for path,r in trials:
            # Bonferroni across the finite number of candidates validated for
            # THIS scenario. The packet bounds are time-uniform within a trial.
            alpha=.05/count
            lo,hi=packet_ber_interval(r['bit_errors_per_packet'],alpha=alpha)
            status='BER upper bound below target' if hi<=.001 else 'BER lower bound above target' if lo>.001 else 'inconclusive'
            summary={k:r[k] for k in ('case','symbol_s','q','mapped_bit_s','packet_payload_bit_s','packets','bit_errors','errored_packets','observed_ber','mechanics','samples_per_symbol','memory_symbols')}
            summary.update(ber_lower=lo,ber_upper=hi,status=status,confidence_alpha=alpha,validated_candidates_in_scenario=count,
                           source_file=path.name,packet_duration_s=30*r['symbol_s'])
            all_trials.append(summary)
        current=[r for r in all_trials if r['case']==c['name']]
        passed=[r for r in current if r['status']=='BER upper bound below target']
        if passed:selected[c['name']]=max(passed,key=lambda r:r['packet_payload_bit_s'])
        elif current:selected[c['name']]=min(current,key=lambda r:r['ber_upper'])
    save(OUT/'validation_summary.json',dict(trials=all_trials,selected=selected,
        confidence='95% simultaneous over validated candidates within each scenario, not simultaneous over all scenarios. Time-uniform packet-fraction bounds. Failed trials retained.',
        scope='Conditional acquired-packet BER only. Not a global optimum, hardware qualification, or completed sampling/drive validation.'))
    lines=['# Initial range-revision calculations','',
        'These are new calculations for the coauthor revision, separate from the original manuscript and baseline results. The manuscript has not been rewritten.', '',
        '## Methods and status','',
        f'- Screen: {len(screen)} parameter combinations, 32 packets per mechanically feasible point. Screening zeros are not treated as BER evidence.',
        '- Current seven-tone, byte-constrained, full-waveform whitened Viterbi decoder; eight-byte payload and six-symbol trailer; acquired clock/phase and stationary receiver noise.',
        '- Accepted finite-volume source quadrature; 128/256 phase agreement and reproduction of the original 0.5/1/2 m fields. No new FEA.',
        '- Four regression tests cover source similarity, physical carrier/noise translation, baseline message-duty torque, and agreement with the original decoder at its baseline.',
        '- Validation uses separate seeds and packet-aware time-uniform intervals. When two candidates are tested for one scenario, the table uses 97.5% bounds for each (95% simultaneous within that scenario). Bounds are not simultaneous across scenarios.',
        '- Mechanics columns are analytic trajectory/duty screens. New achieved-drive trajectories, environmental noise, acquisition, and hardware phase stability remain to be evaluated.', '',
        '## Desktop source: unchanged receiver and drive','',
        '| Distance (m) | Source A2 (nGal) | Gaussian benchmark (bit/s) | Tested mapped rate (bit/s) | Payload incl. trailer (bit/s) | BER interval | Status |',
        '|---:|---:|---:|---:|---:|---|---|']
    desktop=[c for c in cases if c['name'].startswith('desktop')]
    for c in desktop:
        r=selected.get(c['name'])
        cols=f"{r['mapped_bit_s']:.6g} | {r['packet_payload_bit_s']:.6g} | [{r['ber_lower']:.3g}, {r['ber_upper']:.3g}] | {r['status']}" if r else '— | — | — | not validated'
        lines.append(f"| {c['distance_m']:g} | {c['amplitude_m_s2']/1e-11:.5g} | {c['gaussian_benchmark_bit_s']:.6g} | {cols} |")
    r=selected.get('desktop_2m')
    if r:lines+=['',f"The tested 2 m desktop candidate requires {r['symbol_s']/3600:.2f} hours per symbol and {r['packet_duration_s']/3600:.1f} hours per 64-bit packet including its trailer. Stationary noise and accurate acquired phase over this duration are major untested assumptions; this is not a practical 2 m desktop demonstration."]
    lines+=['','## Enlarged source at 2 m','',
            'Both enlarged scenarios use 63.31 kg moving mass, 1.00 m diameter, and 32 times the desktop polar inertia. Assumed drive ratings are 16 N m peak, 4 N m continuous, and 1.6 kW bidirectional mechanical power; center drag is assumed to be 0.4 N m. These are different drive specifications from the desktop.','',
            '| Scenario | Carrier (Hz) | Receiver resonance (Hz) | Gaussian benchmark (bit/s) | Tested payload rate (bit/s) | BER upper bound | Status |',
            '|---|---:|---:|---:|---:|---:|---|']
    for c in cases:
        if not c['name'].startswith('scale'):continue
        r=selected.get(c['name']);label=c['name'].replace('scale2_2m_','').replace('_',' ')
        tail=f"{r['packet_payload_bit_s']:.5g} | {r['ber_upper']:.5g} | {r['status']}" if r else '— | — | benchmark only'
        lines.append(f"| {label} | {c['carrier_hz']:.3f} | {c['receiver']['resonance_hz']:.3f} | {c['gaussian_benchmark_bit_s']:.6g} | {tail} |")
    lines+=['','Maintaining the original rotation speed increases the centrifugal stress scale by four relative to the desktop at that speed. Halving the speed preserves that scale but moves the carrier to 25.15 Hz. The unchanged resonant receiver performs poorly there; retuning it is a separate design assumption, retaining the assumed proof mass, Q, temperature, and readout noise. Neither case is mechanically qualified.','',
            '## One large rotor versus smaller coherent rotors','',
            'Equal total moving mass (63.31 kg) and equal nominal centrifugal stress scale; actual nonintersecting rotor locations are included. All rotors in an array carry the same waveform and phases are ideally aligned at the receiver. Footprints and drive/support requirements differ.','',
            '| Rotors | Individual mass (kg) | Carrier (Hz) | Total A2 (nGal) | Benchmark: fixed receiver (bit/s) | Benchmark: retuned receiver (bit/s) |',
            '|---:|---:|---:|---:|---:|---:|']
    path=OUT/'coherent_array_benchmarks.json'
    if path.exists():
        for a in json.loads(path.read_text()):
            caps=a['gaussian_benchmarks_bit_s']
            lines.append(f"| {a['rotors']} | {a['individual_mass_kg']:.2f} | {a['carrier_hz']:.2f} | {a['ideal_receiver_aligned_amplitude_m_s2']/1e-11:.4g} | {caps['fixed_50p3Hz_receiver']:.5g} | {caps['retuned_receiver']:.5g} |")
    lines+=['','These are source/noise benchmarks, not array BER results. Smaller rotors reduce total inertia but also reduce coherent gravitational coupling at fixed total mass. The receiver frequency response can reverse the ranking; the retuned column does not hold receiver design fixed. Parallel independent subchannels remain to be simulated.','',
            '## Numerical diagnostics','']
    for path in sorted(OUT.glob('sampling_filter_checks*.json')):
        for r in json.loads(path.read_text()):
            lines.append(f"- {r['case']}, Ts={r['symbol_s']:g} s, q={r['q']:g}: {r['packets']} packets; 3/4-symbol filters give {r['paired_filter_bit_errors_3symbols']}/{r['paired_filter_bit_errors_4symbols']} bit errors with {r['paired_filter_disagreeing_bits']} disagreeing decoded bits. Independent doubled-sampling run: {r['independent_double_sampling_bit_errors']} bit errors. Diagnostic counts are not precise rare-event BER estimates.")
    lines+=['','## Next calculations','',
            '1. Refine duration/spacing around the qualifying candidates; the coarse grid does not identify maximum rates. Include the known 1.75 s / q=0.075 desktop candidate in the refined search.',
            '2. Complete sampling/filter checks for the final chosen candidates, and evaluate achieved rather than prescribed actuator phase.',
            '3. Add receiver environmental-noise and long-term phase/frequency-stability sensitivity; multi-hour symbols make these essential.',
            '4. Reassess enlarged-source mechanical loads and receiver assumptions, then evaluate an explicit parallel-subchannel array.',
            '5. Restructure the main text/appendix and update abstract/conclusions after settling those results.','',
            '## Reproduction','',
            'Use Python with NumPy, SciPy, Matplotlib, and Numba, with one BLAS thread. From the submission folder:', '',
            '```text','python scripts/range_revision.py prepare','python scripts/range_revision.py screen --packets 32 --workers 4',
            'python scripts/range_revision.py validate --max-packets 25000 --workers 4',
            'python scripts/range_array_comparison.py','python scripts/range_revision_checks.py --packets 1000 --workers 4',
            'python tests/test_range_revision.py','python scripts/summarize_range_revision.py','```','',
            'Candidate-specific follow-up commands and seeds are recorded in RUN_COMMANDS.md and individual JSON files. Validation can resume with a larger packet limit. Do not rerun scripts concurrently against the same output file.']
    (OUT/'INITIAL_RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(7.5,4.7))
    ds=[c['distance_m'] for c in desktop];caps=[c['gaussian_benchmark_bit_s'] for c in desktop]
    ax.semilogy(ds,caps,'--',color='#6a7683',label='Gaussian benchmark (not achieved rate)')
    for passed,color,marker,label in [(True,'#176b8a','o','Coarse candidate: BER upper bound ≤ 0.001'),(False,'#be7424','x','Tested candidate: target not established')]:
        pts=[(c['distance_m'],selected[c['name']]['packet_payload_bit_s']) for c in desktop if c['name'] in selected and (selected[c['name']]['status']=='BER upper bound below target')==passed]
        if pts:ax.semilogy(*zip(*pts),linestyle='none',marker=marker,color=color,markersize=7,label=label)
    ax.set(xlabel='Distance from rotor axis (m)',ylabel='Rate (bit/s)',title='Desktop rotor: initial rate–distance calculations',xticks=ds)
    ax.grid(True,which='both',alpha=.17);ax.legend(fontsize=8,loc='upper right')
    fig.text(.12,.015,'Fixed receiver and drive • Conditional acquired packets • No environmental noise\nPayload includes the six-symbol trailer; search and actuator verification are incomplete.',fontsize=8)
    fig.tight_layout(rect=(0,.09,1,1));fig.savefig(OUT/'initial_rate_distance.png',dpi=180);fig.savefig(OUT/'initial_rate_distance.pdf');plt.close(fig)
    print('Wrote INITIAL_RESULTS.md, validation_summary.json, and initial_rate_distance.png/.pdf')

if __name__=='__main__':main()
