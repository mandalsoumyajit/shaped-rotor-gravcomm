"""Figures for the joint (hop-limited) design.

modulation_design.pdf   (a) commanded frequency, released CF-FSK vs joint waveform;
                        (b, c) torque-amplitude Pareto fronts with waveform thresholds and coded links
receiver_profile_joint.pdf  receiver noise with the joint design's seven tones
sequence_detector_joint.pdf iterative receiver schematic
Run with the ground-motion venv:  python figure_scripts/joint_figures.py
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
GMS = PAPER/'ground_motion_study'
sys.path.insert(0, str(GMS))
import mitigation as mg                      # noqa: E402  (puts the release package on the path)
from gravcomm.receivers import StructuralOscillator          # noqa: E402
import cpm_optimize as co                    # noqa: E402
import waterfill_tones as wt                 # noqa: E402

RES = GMS/'results'
BLUE, ORANGE, AQUA, GRAY, INK, MUTED = '#2a78d6', '#eb6834', '#1baf7a', '#b9b8b3', '#1f1e1c', '#52514e'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral'], 'mathtext.fontset': 'stix',
                     'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7, 'legend.fontsize': 6,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.linewidth': .5,
                     'lines.linewidth': 1.0, 'xtick.major.width': .5, 'ytick.major.width': .5,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})

T_CF, DF_CF = 224/162, 0.14464
T_HOP = 224/190; DF_HOP = (3/14)/T_HOP
LINKS = {24.: dict(released=0.8125, joint=0.58), 18.: dict(released=1.0, joint=0.697)}   # qualified A2/A0 (joint = NB-LDPC GF(512), 2026-11-05)


def quintic(u):
    return 10*u**3 - 15*u**4 + 6*u**5


def trajectory(offs, T, beta, df, n=80):
    t, f, prev = [], [], 0
    for j, k in enumerate(offs):
        u = np.linspace(0, 1, n)
        tt = j*T + u*T
        uu = np.minimum(u/beta, 1)
        t.append(tt); f.append(df*(prev + (k - prev)*quintic(uu))); prev = k
    return np.concatenate(t), np.concatenate(f)


def modulation_design():
    fig = plt.figure(figsize=(7.0, 2.75))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1, 1], wspace=.55)
    ax = fig.add_subplot(gs[0])
    rng = np.random.default_rng(4)
    import hop_packets as hp
    tones, prev = [], hp.M//2
    for w in rng.integers(0, 2**hp.BITS, 3):
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    t, f = trajectory(np.array(tones) - 3, T_HOP, 1.0, DF_HOP)
    ax.plot(t, f, color=BLUE, label='hop-limited, no dwell ($\\beta=1$, hops $\\leq3$)')
    cf = rng.integers(-3, 4, 9)
    t2, f2 = trajectory(cf, T_CF, 0.6, DF_CF)
    ax.plot(t2, f2, color=MUTED, lw=.8, ls='--', label='CF-FSK with dwell ($\\beta=0.6$)')
    ax.set(xlabel='Time (s)', ylabel='Commanded offset $f_g-f_c$ (Hz)', xlim=(0, 12), ylim=(-.62, .74),
           title='(a) Frequency trajectories')
    ax.legend(loc='upper left', frameon=False, ncol=1)
    ax.spines[['top', 'right']].set_visible(False)

    hop = {r['carrier_hz']: r for r in json.loads((RES/'hop_threshold.json').read_text())}
    for k, v in json.loads((RES/'finite_length_hop.json').read_text()).items():   # longer simulation (60000 symbols)
        hop[float(k)]['threshold_A0'] = v['threshold']
    for k, fc in enumerate((24., 18.)):
        ax = fig.add_subplot(gs[k + 1])
        color = BLUE if fc == 24. else ORANGE
        bound = wt.bound_threshold(fc)
        rows = [json.loads(l) for l in (RES/f'pareto_f{fc:g}.jsonl').read_text().splitlines()]
        cfd = co.make(7, 0.434, co.T_CF, 0.6, 6)
        x = [20*np.log10(r['threshold']/bound) for r in rows if r['threshold']]
        y = [r['fdot']/cfd.fdot_max for r in rows if r['threshold']]
        ax.scatter(x, y, s=2.5, color=GRAY, lw=0, zorder=1, label='evaluated waveforms')
        ref = json.loads((RES/f'pareto_front_f{fc:g}.json').read_text())
        pts = sorted((20*np.log10(r['threshold_refined']/bound), r['torque_ratio']) for r in ref['front'])
        fx, fy, best = [], [], np.inf
        for px, py in pts:
            if py < best*0.995:
                fx.append(px); fy.append(py); best = py
        ax.step(fx, fy, where='post', color=color, lw=1.1, zorder=3, label='Pareto front')
        cfx = 20*np.log10(ref['cffsk_threshold_refined']/bound)
        hx = 20*np.log10(hop[fc]['threshold_A0']/bound); hy = hop[fc]['fdot']/cfd.fdot_max
        lx = {n: 20*np.log10(a/bound) for n, a in LINKS[fc].items()}
        for (wx, wy, cx, mk, lab) in ((cfx, 1.0, lx['released'], 'D', 'released'), (hx, hy, lx['joint'], '*', 'joint')):
            ax.annotate('', xy=(cx, wy), xytext=(wx, wy),
                        arrowprops=dict(arrowstyle='->', color=MUTED, lw=.6, ls=(0, (2, 1.5))), zorder=2)
            ax.plot([wx], [wy], mk, color=INK, ms=4.5 if mk == 'D' else 7, mec='white', mew=.5, zorder=5)
            ax.plot([cx], [wy], mk, color=color, ms=4.5 if mk == 'D' else 7, mec=INK, mew=.5, zorder=5)
            ax.text((wx + cx)/2, wy*1.09, f'{cx - wx:.1f} dB', ha='center', va='bottom', fontsize=5.8, color=MUTED)
        ax.axvline(0, color=MUTED, lw=.6, ls='--', zorder=0)
        ax.text(.12, .085, 'water-filling bound', rotation=90, color=MUTED, fontsize=5.4, va='bottom')
        ax.set_yscale('log')
        ax.set(xlabel='Amplitude above water-filling bound (dB)', xlim=(-.3, 8.2), ylim=(.08, 3),
               title=f'({"bc"[k]}) {fc:g} Hz: ' + ('0.5 and 1 m links' if fc == 24. else '2 m link'))
        if k == 0:
            ax.set_ylabel('Peak torque relative to CF-FSK\n(same rotor inertia)')
        if k == 1:
            from matplotlib.lines import Line2D
            h = [Line2D([], [], marker='o', ls='', color=GRAY, ms=2.5, label='evaluated waveforms'),
                 Line2D([], [], color=INK, label='Pareto front'),
                 Line2D([], [], marker='D', ls='', color=INK, mec='white', ms=4, label='CF-FSK waveform (threshold)'),
                 Line2D([], [], marker='*', ls='', color=INK, mec='white', ms=6.5, label='hop-limited waveform (threshold)'),
                 Line2D([], [], marker='D', ls='', color=GRAY, mec=INK, ms=4, label='CF-FSK link, single-pass decoding'),
                 Line2D([], [], marker='*', ls='', color=GRAY, mec=INK, ms=6.5, label='hop-limited link, joint decoding')]
            fig.legend(handles=h, loc='lower center', ncol=6, frameon=False, fontsize=6, handletextpad=.3,
                       columnspacing=1.2, bbox_to_anchor=(.5, -.13))
        ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(PAPER/'modulation_design.pdf', bbox_inches='tight', pad_inches=.02)
    fig.savefig(RES/'modulation_design.png', dpi=220, bbox_inches='tight', pad_inches=.02)


def receiver_profile():
    from gravcomm.environmental_noise import SeismicBackground
    rec = StructuralOscillator(); env = SeismicBackground(); nu = np.linspace(-1.25, 1.25, 1001)
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.35))
    th = rec.thermal_acceleration_asd(50.3 + nu); ro = rec.readout_acceleration_asd(50.3 + nu)
    axs[0].plot(nu, th/1e-11, color=AQUA, label='Thermal')
    axs[0].plot(nu, ro/1e-11, color=ORANGE, label='Readout')
    axs[0].plot(nu, np.hypot(th, ro)/1e-11, color=INK, label='Instrument total')
    tones = np.arange(-3, 4)*DF_HOP
    for fc, c in ((18, ORANGE), (24, BLUE)):
        total = np.sqrt(th**2 + ro**2 + env.acceleration_psd(fc + nu))/1e-11
        axs[1].plot(nu, total, color=c, label=f'Total at {fc} Hz')
        t = np.sqrt(rec.acceleration_noise_asd(50.3 + tones)**2 + env.acceleration_psd(fc + tones))/1e-11
        axs[1].plot(tones, t, 'o', ms=3, color=c, mec='white', mew=.4)
    for ax, title in zip(axs, ['(a) Fixed instrumental profile', '(b) Total noise and signaling tones']):
        ax.axvspan(-1, 1, color='#e8eef3', zorder=-1)
        ax.set(xlabel='Offset from tuned carrier (Hz)', ylabel=r'Acceleration ASD (nGal/$\sqrt{\mathrm{Hz}}$)',
               title=title, xlim=(-1.25, 1.25), ylim=(0, 55))
        ax.legend(frameon=False); ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(PAPER/'receiver_profile_joint.pdf', bbox_inches='tight', pad_inches=.02)
    fig.savefig(RES/'receiver_profile_joint.png', dpi=220, bbox_inches='tight', pad_inches=.02)


def detector_schematic():
    fig, ax = plt.subplots(figsize=(7.0, 2.45)); ax.set(xlim=(0, 10), ylim=(0, 3.55)); ax.axis('off')

    def box(x, y, w, h, text, color='#eaf1f7'):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.06', fc=color, ec='#24445a', lw=.8))
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=6.6)

    def arrow(a, b, text=None, off=(0, .1), rad=0.):
        ax.annotate('', xy=b, xytext=a, arrowprops=dict(arrowstyle='->', color='#24445a', lw=.9,
                                                        connectionstyle=f'arc3,rad={rad}'))
        if text:
            ax.text((a[0] + b[0])/2 + off[0], (a[1] + b[1])/2 + off[1], text, ha='center', va='bottom', fontsize=5.8,
                    color=MUTED)
    box(.05, 2.55, 1.75, .9, 'Calibrated complex\nsamples\nCausal SOS filter\n16 Hz → 4 Hz')
    box(2.2, 2.55, 1.75, .9, 'Compact innovations\ntransform\n64-sample predictor')
    box(4.35, 2.55, 2.0, .9, 'Reduced-state graph\nphase + 4-tone history\n(hop-limited tones)')
    box(4.35, .9, 2.0, .9, 'Graph forward–backward\nword posteriors with\nbit priors (labeling)')
    box(7.55, .9, 2.35, .9, 'Deinterleave\nTail-biting log-MAP\nSISO decoder (K5)')
    box(7.55, 2.55, 2.35, .9, '16 iterations\nhard decision\nCRC / payload')
    arrow((1.8, 3.0), (2.2, 3.0)); arrow((3.95, 3.0), (4.35, 3.0))
    arrow((5.35, 2.55), (5.35, 1.8), 'all retained edges', off=(.62, -.05))
    arrow((6.35, 1.5), (7.55, 1.5), 'extrinsic bit LLRs', off=(0, .07))
    arrow((7.55, 1.2), (6.35, 1.2))
    ax.text(6.95, 1.12, 'a-priori bit LLRs\n(interleaved)', ha='center', va='top', fontsize=5.8, color=MUTED)
    arrow((8.72, 1.8), (8.72, 2.55), 'a-posteriori', off=(.45, -.05))
    box(.05, .9, 1.75, .9, 'Declared colored PSD\nindependent noise\nmodel', '#f5f0e6'); arrow((1.8, 1.35), (2.9, 2.55))
    box(2.2, .05, 1.75, .7, 'Candidate means\nsame filter and\ninnovations transform', '#f5f0e6'); arrow((3.95, .4), (4.6, 2.55), rad=-.15)
    ax.text(1.0, .35, 'Known phase, timing\nand initial mean state', ha='center', va='center', fontsize=5.8, color=MUTED)
    fig.savefig(PAPER/'sequence_detector_joint.pdf', bbox_inches='tight', pad_inches=.02)
    fig.savefig(RES/'sequence_detector_joint.png', dpi=220, bbox_inches='tight', pad_inches=.02)


def density_limit(fc, eps=1e-2):
    """Amplitude (A0) at which the frame information density through the detector graph (actual 512-word
    source, reduced-state metric; info_density.py) falls below 256 bits with probability eps (probit in dB)."""
    from statistics import NormalDist
    nd = NormalDist(); o = json.loads((RES/'info_density_hop.json').read_text())
    ks = [k for k in o if k.split(':')[0] == f'{fc:g}']
    x = np.array([20*np.log10(float(k.split(':')[1])) for k in ks])
    z = np.array([nd.inv_cdf(min(max(o[k]['p_below_256_gauss'], 1e-12), 1 - 1e-12)) for k in ks])
    c = np.polyfit(x, z, 1)
    return float(10**((nd.inv_cdf(eps) - c[1])/c[0]/20))


def gap_figure():
    """(a) amplitude ladder from the water-filling bound to the qualified links; (b) measured frame-error
    waterfalls of the NB-LDPC joint design and the superseded binary (K5) joint design, with fits."""
    import nb_waterfall_fit as wf
    g = json.loads((RES/'gap_analysis.json').read_text())
    fl = json.loads((RES/'finite_length_hop.json').read_text())       # 60000-symbol thresholds and NA bounds
    nbf = json.loads((RES/'nb_waterfall_fit2.json').read_text())
    hop_rate = 256/(172*224/190)
    db = lambda a, b: 20*np.log10(a/b)
    rows = []
    for fc in (24., 18.):
        L = g[f'{fc:g}']['ladder']; b = wt.bound_threshold(fc)
        rows.append((f'CF-FSK, single pass, {fc:g} Hz', [('waveform', db(L['A_cm'], b)), ('finite length', db(L['A_cm_na'], L['A_cm'])),
                                              ('bit labeling (single pass)', db(L['A_gmi_na'], L['A_cm_na'])),
                                              ('code (single pass)', db(L['A_op'], L['A_gmi_na']))]))
    for fc in (24., 18.):
        b = wt.bound_threshold(fc, rate=hop_rate); f = fl[f'{fc:g}']; thr, na = f['threshold'], f['normal_approx_eps0.01']
        dens = density_limit(fc); a12 = nbf[f'{fc:g}']['probit2']['amplitude']['fer_1e-2']; a_op = LINKS[fc]['joint']
        rows.append((f'Hop-limited, joint, {fc:g} Hz', [('waveform', db(thr, b)), ('finite length', db(na, thr)),
                                           ('word source + reduced-state metric', db(dens, na)),
                                           ('iterative decoding (FER $10^{-2}$)', db(a12, dens)),
                                           ('qualification margin', db(a_op, a12))]))
    colors = {'waveform': '#9ec5f0', 'finite length': '#c9c8c4', 'bit labeling (single pass)': ORANGE,
              'code (single pass)': '#f2b38f', 'word source + reduced-state metric': AQUA,
              'iterative decoding (FER $10^{-2}$)': BLUE, 'qualification margin': '#e6e1f2'}
    dark = ('bit labeling (single pass)', 'iterative decoding (FER $10^{-2}$)', 'word source + reduced-state metric')
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw=dict(width_ratios=[1.55, 1]))
    ax = axs[0]; seen = set()
    for i, (name, segs) in enumerate(rows[::-1]):
        x = 0
        for lab, w in segs:
            ax.barh(i, w, left=x, height=.6, color=colors[lab], edgecolor='white', lw=.6,
                    hatch='////' if lab == 'qualification margin' else None,
                    label=None if lab in seen else lab)
            if w > .4:
                ax.text(x + w/2, i, f'{w:.1f}', ha='center', va='center', fontsize=5.8,
                        color='white' if lab in dark else INK,
                        bbox=dict(fc=colors[lab], ec='none', pad=.5) if lab == 'qualification margin' else None)
            seen.add(lab); x += w
        ax.text(x + .08, i, f'{x:.1f} dB', va='center', fontsize=6, color=MUTED)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]])
    ax.set(xlabel='Received amplitude above water-filling bound (dB)', xlim=(0, 8.6), title='(a) Distance from the limits')
    ax.legend(frameon=False, fontsize=5.4, loc='upper center', bbox_to_anchor=(.45, -.2), ncol=3, columnspacing=.8)
    ax.spines[['top', 'right']].set_visible(False)
    ax = axs[1]
    k5f = json.loads((RES/'k5_waterfall_fit2.json').read_text())
    qual = json.loads((RES/'qualification_nb.json').read_text())
    for fc, c in ((24., BLUE), (18., ORANGE)):
        pts = nbf[f'{fc:g}']['points']; a = np.array([p['amp'] for p in pts]); fer = np.array([p['fer'] for p in pts])
        n = np.array([p['packets'] for p in pts]); k = np.array([p['failures'] for p in pts])
        xs = np.linspace(20*np.log10(a.min()) - .3, 20*np.log10(LINKS[fc]['joint']) + .6, 200)
        ax.semilogy(10**(xs/20), wf.fer_model('probit2', xs, nbf[f'{fc:g}']['probit2']['params']), color=c, lw=1,
                    label=f'GF(512) LDPC, {fc:g} Hz')
        m = k > 0
        ax.semilogy(a[m], fer[m], 'o', color=c, ms=2.6)
        q = qual[f'qualnb_f{fc:g}_a{LINKS[fc]["joint"]:.4f}']['result']
        ax.semilogy(LINKS[fc]['joint'], q['rejected']/q['packets'], 's', color=c, mfc='white', ms=3.2, mew=.8)
        kp = k5f[f'{fc:g}']; ka = np.array([p['amp'] for p in kp['points']])
        xk = np.linspace(20*np.log10(ka.min()) - .2, 20*np.log10(ka.max()) + .3, 200)
        ax.semilogy(10**(xk/20), wf.fer_model('probit', xk, kp['probit']['params']), color=c, lw=.8, ls='--',
                    label=f'binary conv., {fc:g} Hz')
        ax.semilogy(ka, [max(p['fer'], 1e-5) for p in kp['points']], 'x', color=c, ms=2.6, mew=.6)
    ax.set(xscale='log', xlabel='Received amplitude $A_2/A_0$', ylabel='Frame-error probability', ylim=(1e-4, 1),
           xlim=(0.43, 0.9), title='(b) Packet waterfalls, joint receivers')
    ax.set_xticks([0.45, 0.5, 0.6, 0.7, 0.8, 0.9]); ax.set_xticklabels(['0.45', '0.5', '0.6', '0.7', '0.8', '0.9'])
    ax.minorticks_off()
    ax.legend(frameon=False, fontsize=5.4, loc='upper center', bbox_to_anchor=(.5, -.2), ncol=2, columnspacing=.8)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(PAPER/'modulation_gap.pdf', bbox_inches='tight', pad_inches=.02)
    fig.savefig(RES/'modulation_gap.png', dpi=220, bbox_inches='tight', pad_inches=.02)
    return rows


def link_summary():
    from gravcomm.environmental_noise import SeismicBackground
    d = np.array([.5, 1, 2])
    new = {k: [mg.HOP_DESIGNS[x][i] for x in d] for k, i in (('mass', 3), ('torque', 5), ('power', 6))}
    old = {k: [mg.RELEASED_DESIGNS[x][i] for x in d] for k, i in (('mass', 3), ('torque', 5), ('power', 6))}
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.35))
    ax = axs[0]
    ax.semilogy(d, old['mass'], 'o--', color=GRAY, ms=3, lw=.8, label='CF-FSK, single-pass decoding')
    ax.semilogy(d, new['mass'], 'o-', color=BLUE, ms=3.5, label='hop-limited, joint decoding')
    for x, y in zip(d, new['mass']):
        ax.annotate(f'{y:.2f}', (x, y), xytext=(0, -11), textcoords='offset points', ha='center', fontsize=6)
    ax.set(xlabel='Axis-to-receiver distance (m)', ylabel='Bare rotor mass (kg)', xticks=d, ylim=(3, 450),
           title='(a) Qualified source sizes')
    ax.legend(frameon=False, loc='upper left')
    ax = axs[1]
    ax.semilogy(d, old['torque'], 'o--', color=GRAY, ms=3, lw=.8)
    ax.semilogy(d, np.array(old['power'])/1e3, 's--', color=GRAY, ms=3, lw=.8)
    ax.semilogy(d, new['torque'], 'o-', color=BLUE, ms=3.5, label='Peak torque (N m)')
    ax.semilogy(d, np.array(new['power'])/1e3, 's-', color=ORANGE, ms=3.5, label='Peak power (kW)')
    ax.set(xlabel='Axis-to-receiver distance (m)', ylabel='Peak command (units in legend)', xticks=d,
           title='(b) Inertial drive demand')
    ax.legend(frameon=False, loc='upper left')
    ax = axs[2]
    f = np.geomspace(1, 80, 500); c = SeismicBackground().components(f)
    floor = float(StructuralOscillator().acceleration_noise_asd(50.3)); n = np.sqrt(floor**2 + c['background_upper_psd'])
    ax.axvspan(10, 80, color='#ebebeb', lw=0)
    ax.loglog(f, n/1e-11, color=INK, label='Total')
    ax.loglog(f, np.full_like(f, floor)/1e-11, ls='--', color=MUTED, label='Instrument')
    ax.loglog(f, c['support_acceleration_asd']/1e-11, color=BLUE, label='Support')
    ax.loglog(f, c['newtonian_acceleration_asd']/1e-11, color=AQUA, label='Newtonian')
    for fc in (18, 24):
        nc = SeismicBackground().components(fc)
        ax.plot(fc, float(np.sqrt(floor**2 + nc['background_upper_psd'])/1e-11), 'o', color=INK, ms=3)
    ax.set(xlabel='Tuned carrier frequency (Hz)', ylabel=r'Center ASD (nGal/$\sqrt{\mathrm{Hz}}$)', xlim=(1, 80),
           title='(c) Assumed noise scenario')
    ax.legend(frameon=False, loc='lower left', fontsize=5.6)
    for ax in axs:
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(w_pad=1.2)
    fig.savefig(PAPER/'revised_link_summary_joint.pdf', bbox_inches='tight', pad_inches=.02)
    fig.savefig(RES/'revised_link_summary_joint.png', dpi=220, bbox_inches='tight', pad_inches=.02)


if __name__ == '__main__':
    which = sys.argv[1:] or ['modulation', 'profile', 'schematic', 'summary', 'gap']
    if 'gap' in which:
        gap_figure()
    if 'summary' in which:
        link_summary()
    if 'profile' in which:
        receiver_profile()
    if 'schematic' in which:
        detector_schematic()
    if 'modulation' in which:
        modulation_design()
    print('done')
