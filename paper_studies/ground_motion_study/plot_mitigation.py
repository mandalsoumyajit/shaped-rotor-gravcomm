"""Ground-motion mitigation figure for the revision (writes ../ground_motion_mitigation.pdf)."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mitigation as mg

HERE = Path(__file__).resolve().parent
A = json.loads((HERE/('results/analysis_hop.json' if __import__('os').environ.get('GC_DESIGN') == 'hop' else 'results/analysis.json')).read_text())
BLUE, ORANGE, AQUA, YELLOW, GRAY = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#7a7974'
STYLE = {'baseline': (BLUE, '-', 'o', 'Two-stage isolation (baseline)'),
         'feedforward': (ORANGE, '--', 's', 'Witness feed-forward (20 dB)'),
         'structural': (AQUA, '-.', '^', 'Low-loss two-stage isolation'),
         'three_stage': (YELLOW, ':', 'D', 'Three-stage isolation'),
         'gradiometer': (AQUA, '-.', '^', 'Gradiometer (40 dB, 0.8 m)')}


RATES = (json.loads((HERE/'results/rate_scaling.json').read_text())
         if mg.DESIGNS is mg.HOP_DESIGNS and (HERE/'results/rate_scaling.json').exists() else None)


def rate_factor(name, carrier, instrument):
    """Information-rate amplitude factor relative to the qualified 18 Hz channel (joint design).
    Thresholds from rate_scaling.py, interpolated linearly in carrier and in log(instrument factor),
    in log(threshold)."""
    pts = {(v['carrier_hz'], v['instrument']): v['threshold_A0'] for v in RATES.values() if v['scenario'] == name}
    by_k = {}
    for f, k in pts:
        by_k.setdefault(k, []).append(f)
    usable = sorted(k for k, fs in by_k.items() if min(fs) - 1e-9 <= carrier <= max(fs) + 1e-9)
    vals = []
    for k in usable:
        fs = sorted(by_k[k])
        vals.append(float(np.interp(carrier, fs, [np.log(pts[(x, k)]) for x in fs])))
    if not (usable[0] - 1e-12 <= instrument <= usable[-1] + 1e-12):
        raise ValueError(f'no rate thresholds bracket instrument factor {instrument} at {carrier} Hz for {name}')
    lt = np.interp(np.log(instrument), np.log(usable), vals)
    return float(np.exp(lt)/RATES['baseline|18|1']['threshold_A0'])


def mass_2m(name, carrier, instrument):
    g = A['grids']['2m|gradiometer' if name == 'gradiometer' else '2m|single']
    s, a = np.array(g['scale']), np.array(g['a2'])
    factor = (rate_factor(name, carrier, instrument) if RATES else
              mg.equivalent_amplitude_factor(mg.background(name, carrier), mg.background('baseline', 18.), instrument))
    target = mg.DESIGNS[2.0][1]*factor
    x = float(np.exp(np.interp(np.log(target), np.log(a), np.log(s))))
    return mg.sized(2.0, x, carrier)


def main():
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral'], 'mathtext.fontset': 'stix',
                         'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'legend.fontsize': 6,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.linewidth': .5,
                         'lines.linewidth': 1.0, 'xtick.major.width': .5, 'ytick.major.width': .5,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.35))

    # (a) residual support acceleration
    f = np.geomspace(5, 40, 300)
    floor = mg.StructuralOscillator().acceleration_noise_asd(50.3)/1e-11
    for name in ('baseline', 'feedforward', 'structural', 'three_stage'):
        c, ls, _, label = STYLE[name]
        ax[0].loglog(f, mg.background(name, 20.).components(f)['support']/1e-11, color=c, ls=ls, label=label)
    ax[0].axvspan(10, 40, color='#f0efec', zorder=0, lw=0)
    ax[0].text(10.6, 400, 'assumed flat\nground ASD', fontsize=5.6, color='#52514e', va='top')
    ax[0].axhline(floor, color=GRAY, lw=.8, ls=(0, (4, 2)))
    ax[0].text(5.4, floor*1.25, 'instrumental floor', color='#52514e', fontsize=6)
    for fc in (18, 24):
        ax[0].axvline(fc, color=GRAY, lw=.4, ls=':')
    ax[0].set(xlabel='Frequency (Hz)', ylabel=r'Support acceleration ASD (nGal/$\sqrt{\mathrm{Hz}}$)',
              title='(a) Residual support motion', xlim=(5, 40), ylim=(1e-3, 1e3))
    ax[0].set_xticks([5, 10, 20, 40]); ax[0].set_xticklabels(['5', '10', '20', '40'])
    ax[0].legend(loc='lower left', frameon=False, handlelength=2.4)

    # (b) 2 m link mass versus carrier
    carriers = np.arange(10, 24.01, .5)
    for name in ('baseline', 'three_stage'):
        c, ls, m, label = STYLE[name]
        for k, lw, tag in ((1., 1.0, r'$\times1$'), (.1, .7, r'$\times0.1$')):
            rows = [mass_2m(name, fc, k) for fc in carriers]
            mass = np.array([r['mass_kg'] for r in rows]); ok = np.array([r['rim_ok'] for r in rows])
            ax[1].plot(carriers[ok], mass[ok], color=c, ls='-' if k == 1 else '--', lw=lw,
                       label=f"{'Two-stage' if name == 'baseline' else 'Three-stage'}, {tag}")
            if (~ok).any():
                ax[1].plot(carriers[~ok], mass[~ok], color=c, ls='-' if k == 1 else '--', lw=lw, alpha=.35)
    m2 = mg.DESIGNS[2.0][3]
    ax[1].plot([18], [m2], 'o', ms=4, mfc='white', mec='#0b0b0b', mew=.8, zorder=5)
    ax[1].annotate('qualified 2 m design' if mg.DESIGNS is mg.HOP_DESIGNS else 'released 2 m design', (18, m2), (12.6, m2 + 45), fontsize=6,
                   arrowprops={'arrowstyle': '-', 'lw': .4})
    ax[1].set(xlabel='Gravitational carrier (Hz)', ylabel='Required rotor mass, 2 m link (kg)',
              title='(b) Carrier choice at 2 m', xlim=(10, 24), ylim=(0, 380))
    ax[1].legend(loc='upper right', frameon=False, fontsize=5.8, handlelength=2.4, ncol=2,
                 title='Isolation, instrumental noise', title_fontsize=5.8, columnspacing=1.0)

    # (c) 2 m link at 18 Hz versus instrumental floor
    k = np.geomspace(.1, 1, 40)
    for name in ('baseline', 'feedforward', 'three_stage', 'gradiometer'):
        c, ls, m, label = STYLE[name]
        ax[2].semilogx(k, [mass_2m(name, 18., x)['mass_kg'] for x in k], color=c, ls=ls, label=label)
    ax[2].set(xlabel='Instrumental noise relative to assumed floor', ylabel='Required rotor mass, 2 m link (kg)',
              title='(c) Sensor improvement at 18 Hz', xlim=(.1, 1), ylim=(0, 260))
    ax[2].set_xticks([.1, .2, .5, 1]); ax[2].set_xticklabels(['0.1', '0.2', '0.5', '1'])
    ax[2].legend(loc='upper left', frameon=False, handlelength=2.4)
    for a in ax:
        a.grid(True, which='major', color='#e4e3df', lw=.4); a.set_axisbelow(True)
        a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.tight_layout(pad=.4, w_pad=1.0)
    fig.savefig(HERE.parent/'ground_motion_mitigation.pdf', bbox_inches='tight', pad_inches=.02)
    fig.savefig(HERE/'results/ground_motion_mitigation.png', dpi=220, bbox_inches='tight', pad_inches=.02)


if __name__ == '__main__':
    main()
