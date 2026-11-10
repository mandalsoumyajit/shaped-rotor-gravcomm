"""Figure: torque-limited rate of rigid-rotor modulation compared with CF-FSK.

(a) Rate ceiling versus peak drive torque (relative to the CF-FSK design) for the 24 Hz
    channel, for several lowest information-bearing frequencies f_lo.
(b) Minimum peak torque for 1 bit/s, relative to CF-FSK, versus f_lo, for both carriers.
Curves come from torque_bound.rate (small-noise phase-channel water-filling) and are
cached in results/optimality_curves.json.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mitigation as mg
import torque_bound as tb

HERE = Path(__file__).resolve().parent
CACHE = HERE/'results'/'optimality_curves.json'
FIG = HERE.parent/'modulation_optimality.pdf'
BLUE, ORANGE, GRAY, INK, MUTED = '#2a78d6', '#eb6834', '#7a7974', '#1f1e1c', '#52514e'
RAMP = {0.3: '#86b6ef', 0.1: '#3987e5', 0.02: '#1c5cab', 0.0: '#0d366b'}     # ordinal blue (validated)
DASH = {0.3: (0, (1.2, 1.2)), 0.1: (0, (3, 1.5)), 0.02: (0, (6, 1.5, 1.2, 1.5)), 0.0: '-'}
CARRIERS = {24.: (.8125*mg.A0, '24 Hz (0.5 and 1 m)', BLUE, 'o', '-'),
            18.: (1.0*mg.A0, '18 Hz (2 m)', ORANGE, 's', (0, (4, 1.5)))}


def phase_noise(fc, amp):
    bg = mg.background('baseline', fc); p = amp**2/2
    return (mg.total_input_psd(bg, tb.F) + mg.total_input_psd(bg, -tb.F))/(2*p)


def compute():
    out = dict(torque_ratio=list(np.geomspace(0.005, 2.0, 28)), rate={}, flo=list(np.geomspace(0.003, 0.36, 22)),
               min_torque={})
    nw24 = phase_noise(24., CARRIERS[24.][0])
    for flo in RAMP:
        out['rate'][str(flo)] = [tb.rate(nw24, x*tb.FDOT, flo=flo) for x in out['torque_ratio']]
    for fc, (amp, *_r) in CARRIERS.items():
        nw = phase_noise(fc, amp); rows = []
        for flo in [0.0] + out['flo']:
            try:
                fm = brentq(lambda fd: tb.rate(nw, fd, flo=flo) - 1.0, 1e-4, 1e3, xtol=1e-7)
                rows.append(fm/tb.FDOT)
            except ValueError:
                rows.append(None)                       # 1 bit/s not reachable within the excursion
        out['min_torque'][str(fc)] = rows
        out.setdefault('threshold_flo', {})[str(fc)] = brentq(
            lambda f: tb.rate(nw, 1e3, flo=f) - 1.0, 0.1, 0.5, xtol=1e-5)   # 1 bit/s unreachable above
    CACHE.write_text(json.dumps(out, indent=1) + '\n')
    return out


def main(recompute=False):
    d = compute() if recompute or not CACHE.exists() else json.loads(CACHE.read_text())
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral'], 'mathtext.fontset': 'stix',
                         'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'legend.fontsize': 6,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.linewidth': .5,
                         'lines.linewidth': 1.0, 'xtick.major.width': .5, 'ytick.major.width': .5,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))

    # (a) rate ceiling versus peak torque, 24 Hz channel
    x = np.array(d['torque_ratio'])
    labels = {0.0: 'no limit', 0.02: r'$f_{\rm lo}=0.02$ Hz', 0.1: r'$f_{\rm lo}=0.1$ Hz', 0.3: r'$f_{\rm lo}=0.3$ Hz'}
    for flo in (0.0, 0.02, 0.1, 0.3):
        y = np.array(d['rate'][str(flo)], dtype=float)
        ax[0].semilogx(x, y, color=RAMP[flo], ls=DASH[flo], label=labels[flo])
        ax[0].text(x[-1]*1.06, y[-1] + (0.14 if flo == 0.3 else 0), labels[flo], color=MUTED, fontsize=5.6, va='center')
    ax[0].axhline(1.0, color=GRAY, lw=.5, ls=':', zorder=0)
    ax[0].plot([1.0], [1.0], marker='D', ms=4.5, color=INK, mec='white', mew=.6, zorder=5, ls='none')
    ax[0].annotate('CF-FSK design', (1.0, 1.0), xytext=(0.13, 0.42), fontsize=6, color=INK,
                   arrowprops=dict(arrowstyle='-', lw=.5, color=INK))
    ax[0].set(xlabel='Peak drive torque / CF-FSK design torque', ylabel='Rate ceiling (bit/s)',
              title='(a) Torque-limited rate, 24 Hz channel', xlim=(0.005, 5.5), ylim=(0, 3.6))
    ax[0].set_xticks([0.01, 0.1, 1]); ax[0].set_xticklabels(['0.01', '0.1', '1'])
    ax[0].spines[['top', 'right']].set_visible(False)

    # (b) minimum torque for 1 bit/s versus f_lo
    flo = np.array(d['flo'])
    for fc, (_a, label, c, m, ls) in CARRIERS.items():
        rows = d['min_torque'][str(fc)]
        y = np.array([np.nan if v is None else v for v in rows[1:]], dtype=float)
        ok = np.isfinite(y)
        ax[1].loglog(flo[ok], y[ok], color=c, ls=ls, marker=m, ms=3, mec='white', mew=.4, label=label)
    edge = d['threshold_flo']['24.0']                    # 24 Hz loses 1 bit/s first (18 Hz: 0.34 Hz)
    ax[1].axvspan(edge, 0.4, color='#f0efec', lw=0, zorder=0)
    ax[1].text(edge*1.03, 0.0062, '1 bit/s\nunreach-\nable', color=MUTED, fontsize=5.2, va='bottom')
    ax[1].axhline(1.0, color=GRAY, lw=.5, ls=':')
    ax[1].text(0.0033, 1.12, 'CF-FSK design torque', color=MUTED, fontsize=5.6)
    ax[1].set(xlabel=r'Lowest information-bearing frequency $f_{\rm lo}$ (Hz)',
              ylabel='Minimum peak torque for 1 bit/s\n/ CF-FSK design torque',
              title='(b) Torque floor versus phase-coherence requirement', xlim=(0.003, 0.4), ylim=(0.005, 2))
    ax[1].set_xticks([0.01, 0.1]); ax[1].set_xticklabels(['0.01', '0.1'])
    ax[1].legend(loc='upper left', bbox_to_anchor=(0.0, 0.86), frameon=False)
    ax[1].spines[['top', 'right']].set_visible(False)

    fig.tight_layout(w_pad=2.0)
    fig.savefig(FIG); fig.savefig(HERE/'results'/'modulation_optimality.png', dpi=220)
    print('wrote', FIG)


if __name__ == '__main__':
    import sys
    main(recompute='--recompute' in sys.argv)
