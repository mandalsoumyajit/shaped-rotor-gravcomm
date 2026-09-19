"""Summarize the fixed-duration sweep without treating sparse errors as estimates."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / 'results/spacing_sweep'
rows = sorted((json.loads(p.read_text()) for p in OUT.glob('*_fs8_m3.json')),
              key=lambda r: r['tone_spacing_hz'])
assert len(rows) == 9 and all(r['complete'] for r in rows)
plt.rcParams.update({'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 8,
                     'legend.fontsize': 7, 'mathtext.fontset': 'stix'})
fig, axs = plt.subplots(2, 2, figsize=(7, 4.8), layout='constrained')
x = np.array([r['tone_spacing_hz'] for r in rows])
ax = axs[0, 0]
for r in rows:
    lo, hi = r['packet_aware_anytime95_ber_interval']
    y = r['observed_ber']; xx = r['tone_spacing_hz']
    if r['adequate_error_count']:
        ax.errorbar(xx, y, yerr=[[y-lo], [hi-y]], fmt='o', ms=3, color='#185988', capsize=2)
    else:
        ax.errorbar(xx, hi, yerr=hi*.45, uplims=True, fmt='o', ms=3, color='#b45922')
ax.axhline(1e-3, color='k', ls='--', lw=.8, label='BER target')
ax.set(yscale='log', ylabel='Uncoded bit error rate', title='(a) Viterbi decoding, 1.75 s symbols')
ax.legend(loc='upper left')
ax = axs[0, 1]
ax.plot(x, [max(r['tone_noise_asd'])/r['tone_noise_asd'][3] for r in rows], 'o-', ms=3)
ax.set(ylabel='Outer-tone / resonant noise ASD', title='(b) Receiver noise penalty')
ax = axs[1, 0]
for key, label, style in [('peak_torque_nm', 'Peak', '-'),
                           ('worst_message_rms_nm', 'Worst-message RMS', '--'),
                           ('random_message_rms_nm', 'Random-byte RMS', ':')]:
    ax.plot(x, [r['mechanics'][key] for r in rows], style, label=label)
ax.axhline(2, color='0.4', lw=.7, ls='-.')
ax.axhline(.5, color='0.4', lw=.7, ls='-.')
ax.set(yscale='log', ylabel='Torque (N m)', title='(c) Drive limits: 2 peak, 0.5 continuous')
ax.legend()
ax = axs[1, 1]
ax.plot(x, [r['mechanics']['peak_power_w'] for r in rows], 'o-', ms=3)
ax.axhline(400, color='k', lw=.8, ls='--', label='400 W limit')
ax.set(ylabel='Peak mechanical power magnitude (W)', title='(d) Motoring and absorption requirement')
ax.legend()
for ax in axs.flat:
    ax.set(xscale='log', xlabel='Adjacent-tone spacing (Hz)')
    ax.grid(True, which='major', alpha=.2)
fig.savefig(OUT/'spacing_sweep.png', dpi=220)
fig.savefig(OUT/'spacing_sweep.svg')

fields = ['tone_spacing_hz', 'packets', 'payload_bits', 'bit_errors', 'errored_packets',
          'adequate_error_count', 'observed_ber', 'lower95', 'upper95',
          'peak_torque_nm', 'worst_message_rms_nm', 'random_message_rms_nm', 'peak_power_w']
with (OUT/'spacing_sweep.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fields); writer.writeheader()
    for r in rows:
        d = {k: r[k] for k in fields if k in r}
        d.update(r['mechanics']); d.pop('max_rpm')
        d['lower95'], d['upper95'] = r['packet_aware_anytime95_ber_interval']
        writer.writerow(d)
lines = ['# Fixed-duration CF-FSK tone-spacing sweep', '',
 'Desktop shaped rotor at 0.5 m; 1.75 s symbols; 0.7 s quintic transition; seven tones centered at 50.3 Hz. '
 'Eight payload bytes plus six known trailer symbols. Nominal payload rate 1.524 bit/s; '
 'packet rate including trailer 1.219 bit/s. Acquisition overhead is excluded. No added ECC.', '',
 'Standard Viterbi decoding retains continuous phase and the full finite whitening-filter history, '
 'and enforces the existing byte-to-three-tone mapping. Phase and timing are acquired. '
 'Noise is generated from the physical receiver PSD. Primary runs use 8 complex samples/s and '
 '5.25 s whitening memory. The source trajectory is ideal; actuator requirements are evaluated separately.', '',
 'Each point stops after at least 100 bit errors in at least 30 independent packets, or 10,000 packets. '
 'Intervals are packet-aware, anytime-valid 95% intervals for each point individually. '
 'They are not simultaneous confidence intervals across the sweep. Sparse-error points are reported as upper bounds. '
 'Raw observed fractions remain in the CSV for auditing.', '',
 '| Spacing (Hz) | Bits | Errors / errored packets | BER estimate or 95% upper bound | Worst RMS (N m) | Random RMS (N m) | Peak (N m) | Peak power (W) |',
 '|---:|---:|---:|---:|---:|---:|---:|---:|']
for r in rows:
    m=r['mechanics']; hi=r['packet_aware_anytime95_ber_interval'][1]
    value=f"{r['observed_ber']:.3g}" if r['adequate_error_count'] else f'< {hi:.3g}'
    lines.append(f"| {r['tone_spacing_hz']:.5f} | {r['payload_bits']} | {r['bit_errors']} / {r['errored_packets']} | {value} | {m['worst_message_rms_nm']:.3f} | {m['random_message_rms_nm']:.3f} | {m['peak_torque_nm']:.3f} | {m['peak_power_w']:.1f} |")
lines += ['', 'Torque includes 0.05 N m drag at the carrier. Worst-message RMS is the maximum cycle mean '
          'over all valid byte sequences; random-byte RMS assumes independent uniform bytes. '
          'The continuous rating is 0.5 N m, the peak rating 2 N m, and the bidirectional power rating 400 W. '
          'These duty metrics describe sustained payload transmission; trailer and acquisition duty depend on packet scheduling.', '',
          '## Numerical checks', '', '| Run | Packets | Bit errors | Errored packets |', '|---|---:|---:|---:|']
for p in sorted(OUT.glob('*.json')):
    if '_fs16_' in p.name or '_m4' in p.name:
        r=json.loads(p.read_text())
        assert r['complete']
        lines.append(f"| {p.stem} | {r['packets']} | {r['bit_errors']} | {r['errored_packets']} |")
lines += ['', 'The wider-spacing deterioration persists with doubled sampling and longer whitening memory. '
          'At q=0.1, the first 1000 primary packets contain 15 bit errors, compared with 18 at doubled '
          'sampling and 4 with longer memory. These low counts leave uncertainty in the detailed low-BER curve. '
          'The checks support the broad trend; a precise optimum requires more sampling.', '',
          '## Interpretation', '',
          'The receiver exhibits a useful spacing window: narrow spacings produce poorly separated sequences, '
          'while wider spacings encounter increasing input-referred readout noise. The lowest observed error '
          'fraction occurs at 0.14286 Hz, with only 36 errors in nine packets; its supported statement is '
          'BER < 5.96e-4 at 95% confidence. The sweep does not resolve a precise optimal spacing.', '',
          'At 0.09524 Hz the receiver meets the 1e-3 target with a per-point 95% upper bound of 8.98e-4. '
          'Peak torque (1.570 N m) and power (248 W) fit the drive ratings. Worst-message sustained RMS '
          'torque is 0.583 N m, exceeding the 0.5 N m continuous rating; independent random-byte RMS '
          'is 0.285 N m. At 0.14286 Hz, peak torque is 2.330 N m and worst-message RMS is 0.872 N m. '
          'Consequently, continuous transmission guaranteed for every valid message encounters motor torque '
          'constraints before reaching the apparent receiver optimum at this symbol duration. Traffic with '
          'specified statistical or burst duty can have different thermal requirements.', '',
          'All 75 active regression tests pass. The manuscript and public release have not been changed by this sweep.', '']
(OUT/'SUMMARY.md').write_text('\n'.join(lines))
print('\n'.join(lines))
