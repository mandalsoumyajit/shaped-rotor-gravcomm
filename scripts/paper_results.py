"""Rebuild submission figures and a colored-noise Gaussian rate benchmark.

Reads accepted rounded-rotor results; never runs or substitutes a new FEA.
Run with Ubuntu system Python (-s); geometry plotting also needs Gmsh.
"""
from pathlib import Path
import csv
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gravcomm.receivers import StructuralOscillator


def gaussian_waterfill(psd, bandwidth, signal_power):
    """Midpoint-grid real bandpass capacity (bits/s), one-sided SI PSD.

    Relaxes constant-envelope, actuator and alphabet constraints. Each positive
    frequency bin represents one real bandpass channel. Input power is mean
    squared acceleration, not the peak amplitude squared.
    """
    psd = np.asarray(psd, dtype=float)
    if psd.ndim != 1 or not len(psd) or np.any(~np.isfinite(psd)) or np.any(psd <= 0):
        raise ValueError('PSD must be a positive finite vector')
    if not np.isfinite(bandwidth) or bandwidth <= 0 or not np.isfinite(signal_power) or signal_power < 0:
        raise ValueError('Invalid bandwidth or power')
    if signal_power == 0:
        return 0.0
    df = bandwidth / len(psd)
    lo, hi = float(psd.min()), float(psd.max() + signal_power / bandwidth)
    for _ in range(100):
        level = (lo + hi) / 2
        if np.maximum(level - psd, 0).sum() * df > signal_power:
            hi = level
        else:
            lo = level
    allocated = np.maximum((lo + hi) / 2 - psd, 0)
    np.testing.assert_allclose(allocated.sum() * df, signal_power, rtol=1e-10, atol=0)
    return float(np.log2(1 + allocated / psd).sum() * df)


def main():
    out = ROOT / 'results' / 'paper_revision'
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral'],
                         'mathtext.fontset': 'stix', 'font.size': 7,
                         'axes.labelsize': 7, 'axes.titlesize': 7, 'legend.fontsize': 6.5,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                         'axes.linewidth': .5, 'lines.linewidth': .9,
                         'xtick.major.width': .5, 'ytick.major.width': .5,
                         'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    def save(fig, name):
        fig.tight_layout(pad=.45, w_pad=1.0, h_pad=.6)
        fig.savefig(out / (name + '.pdf'), bbox_inches='tight', pad_inches=.02)
        fig.savefig(out / (name + '.png'), dpi=300, bbox_inches='tight', pad_inches=.02)
        plt.close(fig)

    # Actual accepted CAD curves, not a schematic dumbbell.
    sys.path.insert(0, str(ROOT / 'fea'))
    from design_report import curves
    paths = curves(ROOT / 'fea/results/design_sweep/medium_rounded/h0.013_rpm1800/geometry.brep')
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.1), gridspec_kw={'width_ratios': [1.05, 1.15, 1]})
    for points in paths:
        for ax in (axes[0], axes[2]):
            ax.plot(points[:, 0], points[:, 1], color='#263d4a', lw=.65)
    for x in (-200, 200):
        axes[0].add_patch(Circle((x, 0), 35, color='#bdc7cf', zorder=0))
    axes[0].set(xlim=(-260, 260), ylim=(-260, 260), title='(a) Plan view', xticks=[-200, 0, 200], yticks=[-200, 0, 200])
    axes[0].set_aspect('equal'); axes[0].set_xlabel('x (mm)'); axes[0].set_ylabel('y (mm)')
    axes[1].add_patch(Rectangle((-250, -10), 500, 20, fc='#e4e7e9', ec='#263d4a', lw=.6))
    for x in (-200, 200):
        axes[1].add_patch(Rectangle((x-35, -39.6202), 70, 79.2404, fc='#bdc7cf', ec='#263d4a', lw=.6))
    axes[1].add_patch(Rectangle((-6, -10), 12, 20, fc='white', ec='none'))
    axes[1].set(xlim=(-265, 265), ylim=(-55, 55), title='(b) Axial section, y = 0', xlabel='x (mm)', ylabel='z (mm)', xticks=[-200, 0, 200], yticks=[-40, -20, 0, 20, 40])
    axes[1].text(0, 33, 'Steel: 70 mm diameter\n79.24 mm height', ha='center', va='center', fontsize=6.5)
    axes[1].text(0, -33, 'Aluminium: 20 mm thick', ha='center', va='center', fontsize=6.5)
    axes[2].set(xlim=(20, 34), ylim=(12, 26), title='(c) Hub/spoke fillet', xlabel='x (mm)', ylabel='y (mm)', xticks=[20, 25, 30], yticks=[15, 20, 25])
    axes[2].set_aspect('equal')
    axes[2].annotate('R3', (26, 16), (30, 22), arrowprops={'arrowstyle': '->', 'lw': .5}, fontsize=6.5)
    save(fig, 'geometry')

    wave = np.genfromtxt(ROOT / 'fea/results/shaped_waveforms/medium_rounded_0.5m_waveform.csv', delimiter=',', names=True)
    rows = list(csv.DictReader((ROOT / 'fea/results/shaped_waveforms/summary.csv').open()))
    selected = [r for r in rows if r['case'] == 'medium_rounded']
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.0))
    axes[0].plot(wave['phase_rad'] / np.pi, wave['ac_acceleration_m_s2'] / 1e-11, label='Full AC field')
    axes[0].plot(wave['phase_rad'] / np.pi, wave['second_harmonic_m_s2'] / 1e-11, '--', label='Second harmonic')
    axes[0].set(xlabel=r'Mechanical angle / $\pi$', ylabel='Acceleration (nGal)', title='(a) Waveform at d = 0.5 m', xticks=[0, .5, 1, 1.5, 2])
    axes[0].legend(frameon=False, loc='upper center', fontsize=6)
    distances = np.array([float(r['distance_m']) for r in selected])
    a2 = np.array([float(r['a2_m_s2']) for r in selected])
    q = float(selected[0]['quadrupole_anisotropy_kg_m2'])
    for d, marker in [(.5, 'o'), (1., 's'), (2., '^')]:
        h = np.genfromtxt(ROOT / f'fea/results/shaped_waveforms/medium_rounded_{d:g}m_harmonics.csv', delimiter=',', names=True)
        mask = (h['harmonic'] % 2 == 0) & (h['harmonic'] <= 10)
        axes[1].semilogy(h['harmonic'][mask], h['peak_amplitude_m_s2'][mask] / h['peak_amplitude_m_s2'][1], marker+'-', ms=3, label=f'{d:g} m')
    axes[1].set(xlabel='Mechanical harmonic n', ylabel=r'$A_n/A_2$', title='(b) Harmonic content', xticks=[2,4,6,8,10])
    axes[1].legend(frameon=False)
    axes[2].loglog(distances, a2 / 1e-11, 'o-', ms=3, label='Finite volume')
    axes[2].loglog(distances, 9 * 6.67430e-11 * q / (4 * distances**4) / 1e-11, '--', label='Quadrupole')
    axes[2].set(xlabel='Receiver radius d (m)', ylabel=r'$A_2$ (nGal)', title='(c) Distance dependence')
    axes[2].set_xticks(distances, ['0.5', '1', '2'])
    axes[2].xaxis.set_minor_formatter(NullFormatter())
    axes[2].legend(frameon=False)
    for ax in axes: ax.grid(alpha=.16, lw=.4)
    save(fig, 'field')

    rec = StructuralOscillator()
    bandwidth = 1.5625
    n = 32768
    freq = rec.resonance_hz + bandwidth * ((np.arange(n) + .5) / n - .5)
    psd = rec.acceleration_noise_asd(freq)**2
    capacity = [gaussian_waterfill(psd, bandwidth, a*a/2) for a in a2]
    coarse = gaussian_waterfill(psd.reshape(-1, 2).mean(axis=1), bandwidth, a2[0]**2/2)
    assert abs(coarse / capacity[0] - 1) < 1e-7
    fig, ax = plt.subplots(figsize=(3.4, 1.85))
    for label, values in [('Thermal', rec.thermal_acceleration_asd(freq)), ('Readout, input referred', rec.readout_acceleration_asd(freq)), ('Total', np.sqrt(psd))]:
        ax.semilogy(freq - rec.resonance_hz, values / 1e-11, label=label)
    tones = np.arange(-3, 4) * .15625
    ax.plot(tones, rec.acceleration_noise_asd(tones + rec.resonance_hz) / 1e-11, 'ko', ms=3)
    ax.set(xlabel='Frequency offset from 50.3 Hz (Hz)', ylabel=r'ASD (nGal / $\sqrt{\mathrm{Hz}}$)', ylim=(1, 40))
    ax.legend(frameon=False, ncol=1, loc='lower right'); ax.grid(alpha=.16, lw=.4)
    save(fig, 'receiver')
    # Compact chain: all mathematical explanation belongs in the caption.
    fig, ax = plt.subplots(figsize=(3.4, 1.22))
    ax.set(xlim=(0, 4), ylim=(0, 1.5)); ax.axis('off')
    labels = [['Byte / tone\nmapping', 'Smooth\nfrequency command', 'Constrained\nmotor drive', 'Shaped\nrotor'],
              ['Decoded\nbytes', 'Dwell-bin\nenergy detector', 'Calibration\nand bandpass', 'Gravimeter']]
    w, h = .86, .43
    for row, y in enumerate([1.02, .12]):
        for col, label in enumerate(labels[row]):
            x = col + .05
            ax.add_patch(Rectangle((x, y), w, h, fc='#f2f4f5', ec='#263d4a', lw=.55))
            ax.text(x+w/2, y+h/2, label, ha='center', va='center', fontsize=6.1)
        for col in range(3):
            start, end = ((col+.91, y+h/2), (col+1.05, y+h/2))
            if row: start, end = end, start
            ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=5, lw=.55, color='#263d4a'))
    ax.add_patch(FancyArrowPatch((3.48, 1.02), (3.48, .55), arrowstyle='-|>', mutation_scale=5, lw=.55, color='#263d4a'))
    ax.text(3.35, .78, 'Newtonian\nfield', ha='right', va='center', fontsize=6)
    save(fig, 'modulation_chain')
    record = {'bandwidth_hz': bandwidth, 'grid_points': n,
              'waterfill_coarsening_relative_change': coarse / capacity[0] - 1,
              'rows': [{'distance_m': d, 'a2_m_s2': a, 'gaussian_benchmark_bit_s': c}
                       for d, a, c in zip(distances.tolist(), a2.tolist(), capacity)]}
    (out / 'gaussian_benchmark.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
