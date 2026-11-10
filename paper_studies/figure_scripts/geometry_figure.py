"""Figure 2 for the revision: geometry block of AIP_Advances_submission/scripts/paper_results.py
with the material label spelled 'Aluminum' (Referee 1). Reads the accepted CAD curves only."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

_ROOT = Path(__file__).resolve().parents[2]
AIP = _ROOT if (_ROOT/'paper_studies').is_dir() else _ROOT/'AIP_Advances_submission'   # repo or working tree
sys.path.insert(0, str(AIP/'fea'))
from design_report import curves

OUT = Path(__file__).resolve().parents[1]/'geometry_revised.pdf'


def main():
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral'],
                         'mathtext.fontset': 'stix', 'font.size': 7,
                         'axes.labelsize': 7, 'axes.titlesize': 7, 'legend.fontsize': 6.5,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                         'axes.linewidth': .5, 'lines.linewidth': .9,
                         'xtick.major.width': .5, 'ytick.major.width': .5,
                         'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    paths = curves(AIP/'fea/results/design_sweep/medium_rounded/h0.013_rpm1800/geometry.brep')
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
    axes[1].set(xlim=(-265, 265), ylim=(-55, 55), title='(b) Axial section, y = 0', xlabel='x (mm)', ylabel='z (mm)',
                xticks=[-200, 0, 200], yticks=[-40, -20, 0, 20, 40])
    axes[1].text(0, 33, 'Steel: 70 mm diameter\n79.24 mm height', ha='center', va='center', fontsize=6.5)
    axes[1].text(0, -33, 'Aluminum: 20 mm thick', ha='center', va='center', fontsize=6.5)
    axes[2].set(xlim=(20, 34), ylim=(12, 26), title='(c) Hub/spoke fillet', xlabel='x (mm)', ylabel='y (mm)',
                xticks=[20, 25, 30], yticks=[15, 20, 25])
    axes[2].set_aspect('equal')
    axes[2].annotate('R3', (26, 16), (30, 22), arrowprops={'arrowstyle': '->', 'lw': .5}, fontsize=6.5)
    fig.tight_layout(pad=.45, w_pad=1.0, h_pad=.6)
    fig.savefig(OUT, bbox_inches='tight', pad_inches=.02)


if __name__ == '__main__':
    main()
