"""Torque-amplitude Pareto front of CF-FSK waveforms from cpm_optimize.py results.

x: threshold amplitude for the required information rate, in dB above the water-filling
   capacity threshold (the strict bound for any signal of that power);
y: peak inertial torque at fixed rotor inertia, relative to the released CF-FSK waveform.
--refine re-evaluates the non-dominated designs with a longer predictor (order 18) and more
symbols, writes results/pareto_front_f<fc>.json, and plots the refined front.
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cpm_optimize as co
import waterfill_tones as wt

HERE = Path(__file__).resolve().parent
RES = HERE/'results'
FIG = HERE.parent/'modulation_pareto.pdf'
BLUE, ORANGE, GRAY, INK, MUTED = '#2a78d6', '#eb6834', '#b9b8b3', '#1f1e1c', '#52514e'
CARRIER_STYLE = {24.: (BLUE, '24 Hz channel (0.5 and 1 m)'), 18.: (ORANGE, '18 Hz channel (2 m)')}


def load(fc):
    rows = [json.loads(l) for l in (RES/f'pareto_f{fc:g}.jsonl').read_text().splitlines()]
    cf = co.make(7, 0.434, co.T_CF, 0.6, 6)
    for r in rows:
        r['torque_ratio'] = r['fdot']/cf.fdot_max
    return rows


def front(rows, bound, key='threshold'):
    pts = sorted((20*np.log10(r[key]/bound), r['torque_ratio'], i) for i, r in enumerate(rows) if r.get(key))
    out, best = [], np.inf
    for x, y, i in pts:
        if y < best*0.995:
            out.append(i); best = y
    return out


def _refine_one(args):
    fc, r = args
    r = dict(r)
    d = co.make(r['M'], r['excursion_req'], r['T'], r['beta'], r['D'])
    r['threshold_refined'] = co.threshold(d, fc, n_sym=15000, p_max=18)
    return r


def refine(fc, rows, bound, workers=18):
    from concurrent.futures import ProcessPoolExecutor
    idx = front(rows, bound)
    cf = dict(M=7, excursion_req=0.434, T=co.T_CF, beta=0.6, D=6)
    jobs = [(fc, rows[i]) for i in idx] + [(fc, cf)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        res = list(pool.map(_refine_one, jobs))
    out, cf_ref = res[:-1], res[-1]['threshold_refined']
    for r in out:
        print(fc, {k: r[k] for k in ('M', 'D', 'beta', 'T', 'excursion', 'torque_ratio', 'threshold', 'threshold_refined')}, flush=True)
    data = dict(carrier_hz=fc, bound_threshold=bound, cffsk_threshold_refined=cf_ref, front=out)
    (RES/f'pareto_front_f{fc:g}.json').write_text(json.dumps(data, indent=1) + '\n')
    return data


def main(do_refine=False):
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral'], 'mathtext.fontset': 'stix',
                         'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'legend.fontsize': 6,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.linewidth': .5,
                         'lines.linewidth': 1.0, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True)
    for ax, fc in zip(axes, (24., 18.)):
        color, label = CARRIER_STYLE[fc]
        rows = load(fc); bound = wt.bound_threshold(fc)
        x = np.array([20*np.log10(r['threshold']/bound) if r['threshold'] else np.nan for r in rows])
        y = np.array([r['torque_ratio'] for r in rows])
        ax.scatter(x, y, s=4, color=GRAY, lw=0, label='all evaluated waveforms', zorder=1)
        path = RES/f'pareto_front_f{fc:g}.json'
        if do_refine:
            ref = refine(fc, rows, bound)
        elif path.exists():
            ref = json.loads(path.read_text())
        else:
            ref = None
        if ref:
            pts = sorted((20*np.log10(r['threshold_refined']/bound), r['torque_ratio']) for r in ref['front'])
            fx, fy, best = [], [], np.inf
            for px, py in pts:                       # re-take the non-dominated set after refinement
                if py < best*0.995:
                    fx.append(px); fy.append(py); best = py
            cfx = 20*np.log10(ref['cffsk_threshold_refined']/bound)
        else:
            idx = front(rows, bound)
            fx, fy = x[idx], y[idx]
            cfx = x[[i for i, r in enumerate(rows) if r['M'] == 7 and r['D'] == 6 and r['beta'] == 0.6
                     and abs(r['excursion_req'] - 0.434) < 1e-6 and abs(r['T'] - co.T_CF) < 1e-6]][0]
        ax.step(fx, fy, where='post', color=color, lw=1.2, label='Pareto front', zorder=3)
        ax.plot(fx, fy, 'o', color=color, ms=3, mec='white', mew=.4, zorder=4)
        ax.plot([cfx], [1.0], 'D', color=INK, ms=4.5, mec='white', mew=.6, zorder=5, label='released CF-FSK')
        ax.axvline(0, color=MUTED, lw=.6, ls='--', zorder=0)
        ax.text(0.08, 0.09, 'capacity bound', rotation=90, color=MUTED, fontsize=5.6, va='bottom')
        ax.set_yscale('log')
        ax.set(xlabel='Threshold amplitude above water-filling bound (dB)', title=label, xlim=(-0.2, 6), ylim=(0.08, 3))
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('Peak torque / CF-FSK torque\n(same rotor inertia)')
    axes[0].legend(loc='upper right', frameon=False, markerscale=2.5)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG); fig.savefig(RES/'modulation_pareto.png', dpi=220)
    print('wrote', FIG)


if __name__ == '__main__':
    main('--refine' in sys.argv)
