"""Audit manifest and manuscript-facing summary of completed checks."""
import hashlib,json,platform
from importlib.metadata import version
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/communications_checks'
records=[]
for p in sorted(OUT.glob('fs*.json')):
    r=json.loads(p.read_text()); assert r['complete'], p
    records.append((p.name,r))
assert len(records)==5
drive=json.loads((OUT/'drive_refinement.json').read_text());assert len(drive)==4
energy=json.loads((OUT/'revised_energy.json').read_text())
capacities=json.loads((OUT/'capacity.json').read_text())
lines=['# Communications checks before manuscript revision','',
       'All random-message runs use the selected shaped rotor at 0.5 m, 2 s symbols, '
       '40% transitions, spacing–dwell product 0.1, eight-byte packets and six known trailer '
       'symbols. Initial timing, phase, frequency, tone and prehistory are acquired. '
       'Nominal and trailer-inclusive rates are 1.333 and 1.067 bit/s. No added ECC.', '',
       '## Sampling, whitening and technical-noise checks','',
       '| Sampling (Hz) | FIR memory (s) | Added noise ASD (SI/√Hz) | Bits | Errors / errored packets | BER estimate or upper bound | Packet-aware 95% interval |',
       '|---:|---:|---:|---:|---:|---:|---|']
for name,r in records:
    c=r['config'];lo,hi=r['interval95']
    value=f"{r['observed_ber']:.4g}" if r['adequate_error_count'] else f'< {hi:.4g}'
    lines.append(f"| {c['fs']:g} | {c['memory']:g} | {c['technical_asd']:.1g} | {r['bits']} | {r['bit_errors']} / {r['errored_packets']} | {value} | [{lo:.3g}, {hi:.3g}] |")
lines += ['', 'The error-count stopping criterion is at least 100 bit errors in at least 30 '
          'errored packets. The longer-memory diagnostic has a 6,000-packet cap; '
          'when that criterion is unmet it is reported as an upper bound. Confidence sequences '
          'operate on independent packet error fractions and account for arbitrary within-packet '
          'error clustering. Intervals are individual 95% statements, with no pooled or simultaneous '
          'coverage claim. Configurations use common seed indices and can be correlated with one another.', '',
          'Added technical noise is independent flat Gaussian acceleration noise added in quadrature '
          'to the physical PSD. The whitening filter is redesigned using that total PSD. This assumes '
          'the additional spectrum is known to the receiver. Symbol timing and tone spacing are held fixed; '
          'the waveform has not been reoptimized for either added-noise case. Coherent interference is outside this test. '
          'The corresponding Gaussian benchmarks are 2.9202 and 1.8848 bit/s, compared with 3.9045 bit/s '
          'without added noise. Poor BER at the fixed waveform therefore does not establish a fundamental '
          'rate ceiling under those noise levels.', '',
          '## Achieved-drive sensitivity', '',
          '| Symbol duration (s) | Message | Max envelope error vs nominal | Change under 200→400 Hz drive-grid refinement | Noisy packet decisions changed / paired trials |',
          '|---:|---|---:|---:|---:|']
for r in drive:
    lines.append(f"| {r['symbol_s']} | {r['payload_kind']} | {r['achieved_nominal_max_envelope_difference']:.3g} | {r['coarse_fine_max_envelope_difference']:.3g} | {r['disagreeing_decoded_packets']} / {r['paired_noise_trials']} |")
lines += ['', 'The tests include repeated bytes 251,42 (the worst sustained-duty cycle) and a fixed '
          'random message at both 2 s and 1.75 s. Each pair uses identical physical noise and the '
          'same nominal decoder templates. These are fixed-message sensitivity checks; their '
          'error counts are not estimates of random-message BER. Trailer duty lowers packet RMS '
          'relative to the indefinitely repeated payload. At 1.75 s the drive can follow the waveform '
          'while the worst-message continuous-duty rating remains exceeded.', '',
          '## Energy and centrifugal references','']
d=energy['diagnostics']
lines += [f"The closed all-transition cycle lasts {energy['duration_s']:g} s. Positive motor work is "
          f"{d['motoring_work_j']:.6f} J, absorbed work {d['absorbed_work_j']:.6f} J, and drag heat "
          f"{d['drag_heat_j']:.6f} J. The energy-balance residual including kinetic-energy change is "
          f"{d['energy_balance_error_j']:.3g} J.", '',
          f"At the selected maximum speed of 1516.5 rpm, stored energy is "
          f"{energy['stored_energy_at_max_speed_j']:.3f} J and the scaled 10 mm regional centrifugal "
          f"stress average is {energy['centrifugal_regional_mean_mpa']:.4f} MPa. This remains a centrifugal-only calculation.", '',
          '## Capacity convergence','', '| Distance (m) | Capacity in 2 Hz (bit/s) |', '|---:|---:|']
for r in capacities:
    if r['bandwidth_hz']==2:
        lines.append(f"| {r['distance_m']} | {r['capacity_bit_s']:.9g} |")
lines += ['', f"Largest relative change on doubling the frequency grid: {max(abs(r['relative_grid_change']) for r in capacities):.3g}. "
          'The 0.5 m value is unchanged to displayed precision between 2, 4 and 8 Hz bands.', '',
          '## Reproduction and scope','',
          'Scripts: `check_communications_revision.py`, `check_revision_drive.py`, '
          '`check_revision_capacity.py`, and `summarize_communications_checks.py`. '
          'Run from the submission directory with PYTHONPATH=src in the configured WSL Python '
          'environment, one BLAS thread, and the sequence/test dependencies. BER configurations '
          'are stored in each JSON file; use --fs and --memory accordingly, with --technical-asd '
          'for the noise checks. The 8 s FIR diagnostic uses --max-packets 6000. '
          'The JSON records retain all packet error counts and are resumable.', '',
          'See `test_report.txt` for the regression result and `manifest.json` for SHA-256 hashes '
          'of the checked scripts, core modules and numerical records. The manuscript and public '
          'release remain unchanged during these checks.']
(OUT/'SUMMARY.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
paths=list(OUT.glob('*.json'))
paths=[p for p in paths if p.name!='manifest.json']
paths+=list((ROOT/'src/gravcomm').glob('*.py'))
paths += [ROOT/'scripts'/name for name in ['check_communications_revision.py','check_revision_drive.py',
          'check_revision_capacity.py','summarize_communications_checks.py','desktop_cf_fsk.py','optimize_cf_fsk.py','paper_results.py']]
paths += [ROOT/'fea/results/shaped_waveforms/summary.csv',ROOT/'results/desktop_cf_fsk/summary.json',
          ROOT/'results/sequence_cf_fsk/selected_ber.json']
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
(OUT/'manifest.json').write_text(json.dumps(dict(sha256=manifest,python=platform.python_version(),
    dependencies={name:version(name) for name in ['numpy','scipy','numba','pytest']}),indent=2)+'\n')
print('\n'.join(lines))
