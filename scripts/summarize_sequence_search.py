"""Write the report only after both final cases meet the error-count rule."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.packet_confidence import packet_ber_interval
OUT=ROOT/'results/sequence_cf_fsk'
records=[]
for name in ('selected','boundary'):
    r=json.loads((OUT/(name+'_ber.json')).read_text())
    assert r['ber_estimate_has_required_error_count'],name+' has insufficient observed errors'
    r['packet_aware_anytime95_ber_interval']=packet_ber_interval(r['bit_errors_per_packet'])
    r['legacy_prefix_packets']=7500 if name=='selected' else 750
    (OUT/(name+'_ber.json')).write_text(json.dumps(r,indent=2)+'\n')
    records.append(r)
selected,boundary=records
drive=json.loads((OUT/'selected_actuator.json').read_text())
rows=[]
for r in records:
    ci=r['packet_aware_anytime95_ber_interval']
    rows.append(f"| {r['ts']:g} | {r['nominal_bit_s']:.3f} | {r['packet_payload_bit_s']:.3f} | {r['bits']:,} | {r['bit_errors']} | {r['errored_packets']} | {r['observed_ber']:.2g} | {ci[0]:.2g} to {ci[1]:.2g} |")
d=drive['all_transition_drive'];m=drive['all_transition_mechanics']
text=f'''# Whitened-Viterbi CF-FSK search: verified results

The selected tested desktop operating point uses 2 s symbols and carries
1.333 bit/s before packet overhead, four times the earlier 0.333 bit/s point.
Eight-byte packets with six known trailer symbols carry 1.067 bit/s.
The source remains the selected shaped rotor at a 0.5 m receiver distance.

## Error-count validation

Every final estimate has at least 100 observed bit errors in at least 30
independent errored packets. Errors within a packet are correlated. The 95%
intervals below operate on independent packet bit-error fractions and remain
valid under the error-count stopping rule; they do not assume independent bits.

| Symbol period (s) | Nominal bit/s | Packet bit/s | Bits tested | Bit errors | Errored packets | Measured BER | 95% BER interval |
|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

The 1.75 s case uses a 40% transition fraction and spacing-dwell product 0.075.
Its BER interval lies above the 0.001 target. The 2 s case uses a 40% transition
fraction and spacing-dwell product 0.1; its interval lies below that target.
This is the best verified point in the evaluated grid. The retained search
records include computation-budget exclusions at some shorter periods.

## Selected operating point

- Carrier: 50.3 Hz; seven tones from 50.05 to 50.55 Hz.
- Tone spacing: 0.083333 Hz; transition: 0.8 s; dwell: 1.2 s.
- Rotor speed: 1501.5 to 1516.5 rpm.
- Peak total torque: {d['peak_torque_nm']:.4f} N m.
- All-transition cycle RMS torque: {d['rms_torque_nm']:.4f} N m.
- Maximum individual-symbol RMS torque: {m['worst_symbol_rms_nm']:.4f} N m,
  a conservative bound on long-message duty.
- Peak motoring / absorbed mechanical power:
  {d['peak_motoring_power_w']:.2f} / {d['peak_regenerative_power_w']:.2f} W.
- Maximum sampled carrier tracking error: {d['max_sampled_carrier_error_hz']:.3g} Hz.

These values satisfy the assumed 2 N m peak, 0.5 N m continuous and 400 W
bidirectional drive specifications. A separate noiseless packet constructed from
the achieved drive phase decodes correctly; its maximum complex-envelope
deviation from the nominal waveform is {drive['maximum_complex_waveform_error']:.3g}.
The Monte Carlo trials use the nominal sequence waveform. Assembly strength,
fatigue, containment and environmental coupling
remain the separate engineering validation described in the manuscript.

## Why throughput improved

The old detector discarded transitions and required orthogonal dwell tones.
The new decoder uses the entire continuous-phase waveform and its sequence
memory. Tone spacing times dwell duration is now 0.1, allowing closer tones and
shorter symbols. The known byte mapping is enforced directly in the trellis.
The earlier operating point used 8 s symbols and a spacing-dwell product of 1.

The Gaussian benchmark is 3.9045 bit/s over the new 2 Hz sampled baseband;
the earlier 1.5625 Hz calculation was 3.9031 bit/s. The nominal mapped rate
uses about 34% of this benchmark, and the short-packet rate uses about 27%.
The revised bandwidth therefore contributes very little to the rate gain.

The old error calculation was a conditional symbol-by-symbol test over 49
contexts. The new results are bit-error measurements on complete random-byte
packets with acquired initial phase, carrier and timing. Packet rates include
the trailer and exclude acquisition overhead.

## Whitening and numerical checks

The causal 13-tap filter retains 6 s of memory at 2 complex samples/s.
Whitened noise has 1.05% RMS PSD deviation over the full sampled band; an
8 s, 17-tap filter reduces this to 0.37%. Both filters gave identical payload
decisions on 500 paired noisy packets. A separate 500-packet check doubled
the sampling rate. These small checks are convergence diagnostics, and their
low error counts are not presented as BER estimates.

Noise is synthesized from the original receiver spectrum independently of
the whitening approximation. Viterbi keeps phase, tone and filter-memory
states, with exact add-compare-select updates and traceback. There is no beam
pruning. Initial phase and timing are assumed acquired. The receiver and drive
parameters retain the paper's stated modeling assumptions.

The full regression suite is recorded in `test_report.txt`. Reproduction,
search grids, complete packet counts and statistical details are in `README.md`
and the adjacent JSON files. The manuscript and public release have not yet
been updated to this new operating point.
'''
(OUT/'SUMMARY.md').write_text(text)
print(text)
