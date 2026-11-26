"""Figure 13: dimensioned schematic of the preliminary horizontal-shaft concept, for the 0.5 m-link rotor.

Derived from fea/draw_horizontal_support.py (release 0.5.0), whose paper variant hard-coded the conventional
design (s = 0.9612, 481 mm diameter). Here every rotor dimension follows from the scale s of the qualified
joint design (Sec. 3.2 proportions at R0 = 250 mm):
  outer radius 250 s, rim width 0.08 R, spoke width 30 s, carrier thickness 20 s, hub radius 30 s,
  insert centres at 200 s, insert radius 35 s, boss radius 50 s, steel insert length 79.235 s and
  same-mass tungsten length 0.4361 of it (Appendix A).
Bearings, pedestals, shafts and drive are schematic and unchanged.
Run: python paper_studies/figure_scripts/draw_horizontal_support.py (writes the PDF to paper_studies/,
like the other figure scripts, and a PNG to paper_studies/ground_motion_study/results/).
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Arc

PAPER = Path(__file__).resolve().parents[1]
RES = PAPER/'ground_motion_study'/'results'
S = 0.9010437566657981            # 0.5 m-link joint design (paper_studies/ground_motion_study/results/sizing_nb.json)
R = 250*S; RI = 0.92*R; HALF_T = 10*S; YI = 200*S; RINS = 35*S; RBOSS = 50*S; RHUB = 30*S
L_STEEL = 79.235*S; L_W = 0.4361*L_STEEL
center = 300.

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'none'})
fig, axs = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw={'width_ratios': [1, 1.15]}); fig.patch.set_facecolor('#f7f9fc')
blue = '#377eb8'; dark = '#18354b'; gold = '#b98a36'; gray = '#718096'; base = '#d9e2ec'


def rect(ax, x, y, w, h, c, **kw):
    p = Rectangle((x, y), w, h, facecolor=c, edgecolor=dark, lw=1.2, **kw); ax.add_patch(p); return p


def label(ax, text, xy, at, align='left'):
    ax.annotate(text, xy=xy, xytext=at, ha=align, va='center', fontsize=10, color=dark,
                arrowprops=dict(arrowstyle='-', color=dark, lw=.9), bbox=dict(facecolor='#f7f9fc', edgecolor='none', pad=2))


def dim(ax, a, b, text, offset=8):
    ax.annotate('', a, b, arrowprops=dict(arrowstyle='<->', color=dark, lw=.8))
    ax.text((a[0] + b[0])/2, (a[1] + b[1])/2 + offset, text, ha='center', color=dark, fontsize=9)


ax = axs[0]
ax.add_patch(Circle((0, center), R, fc=blue, ec=dark, lw=1.5)); ax.add_patch(Circle((0, center), RI, fc='#f7f9fc', ec=dark, lw=1))
rect(ax, -RI, center - 15*S, 2*RI, 30*S, blue); rect(ax, -15*S, center - RI, 30*S, 2*RI, blue)
for x in (-YI, YI):
    ax.add_patch(Circle((x, center), RBOSS, fc=blue, ec=dark)); ax.add_patch(Circle((x, center), RINS, fc=gold, ec=dark))
ax.add_patch(Circle((0, center), RHUB, fc=gray, ec=dark)); ax.add_patch(Circle((0, center), 10, fc='white', ec=dark)); ax.plot(0, center, 'o', color=dark, ms=3)
ax.add_patch(Arc((0, center), 160, 160, theta1=35, theta2=125, color=dark, lw=1.5)); ax.annotate('', (-46, 365), (-32, 373), arrowprops=dict(arrowstyle='->', color=dark))
label(ax, 'Two dense inserts', (-YI + 8, center), (-265, 610))
label(ax, 'Aluminum carrier + rim', (125, center + ((R + RI)**2/4 - 125**2)**0.5), (60, 590))
label(ax, 'Hub / shaft axis\n(points toward viewer)', (0, center), (35, 160))
dim(ax, (-R, 30), (R, 30), f'Rotor diameter {2*R:.0f} mm', 12)
ax.set_title('A  Looking along the shaft\nRotor turns in a vertical plane', loc='left', fontsize=13, color=dark, pad=18)
ax.set_xlim(-285, 285); ax.set_ylim(-25, 655)

ax = axs[1]
rect(ax, -205, 0, 390, 18, base)
for x in (-60, 60):                      # stationary pedestals outside the axial swept envelope
    rect(ax, x - 12, 18, 24, 265, base); rect(ax, x - 20, 283, 40, 34, base)
    rect(ax, x - 13, 289, 26, 22, 'white')
    if x < 0:
        ax.plot([x - 11, x - 2], [291, 309], color=blue, lw=2); ax.plot([x + 2, x + 11], [309, 291], color=blue, lw=2)
    else:
        rect(ax, x - 9, 289, 18, 4, blue); rect(ax, x - 9, 307, 18, 4, blue)
rect(ax, -92, 290, 92 - HALF_T, 20, gray); rect(ax, HALF_T, 290, 85 - HALF_T, 20, gray)     # stub shafts
rect(ax, -HALF_T - 5, center - RHUB, 5, 2*RHUB, gray); rect(ax, HALF_T, center - RHUB, 5, 2*RHUB, gray)  # hub flanges
rect(ax, -HALF_T, center - R, 2*HALF_T, 2*R, blue)
for y in (center - YI, center + YI):
    ax.add_patch(Rectangle((-L_STEEL/2, y - RINS), L_STEEL, 2*RINS, fc='none', ec=gold, lw=1.3, ls='--'))
    rect(ax, -L_W/2, y - RINS, L_W, 2*RINS, gold)
rect(ax, -115, 284, 23, 32, gray); rect(ax, -190, 270, 75, 60, base); rect(ax, -181, 18, 55, 252, base)
ax.plot([-202, 105], [300, 300], ls='-.', color=dark, lw=.7)
label(ax, 'Motor + coupling', (-145, 330), (-195, 640))
label(ax, 'Drive-side locating pair\n(back-to-back angular contact)', (-60, 311), (-195, 585))
label(ax, 'Floating roller bearing\n(allows axial expansion)', (60, 310), (85, 420))
label(ax, 'Bolted hub flanges\non two stub shafts', (HALF_T + 5, center + RHUB - 4), (85, 355))
label(ax, 'Tungsten insert', (L_W/2, center + YI), (85, 515))
label(ax, 'Dashed: steel insert\nsame mass, longer', (-L_STEEL/2, center + YI - 10), (-195, 450))
label(ax, 'Stationary pedestals', (-60, 120), (-195, 100))
label(ax, 'Common rigid base', (100, 15), (85, 75))
dim(ax, (-60, 240), (60, 240), '120 mm bearing span', 12)
ax.annotate('Gravity', xy=(155, 175), xytext=(155, 230), ha='center', color=dark, arrowprops=dict(arrowstyle='->', color=dark, lw=1.5))
ax.set_title('B  Edge view of the rotor\nHorizontal shaft supported on both sides', loc='left', fontsize=13, color=dark, pad=18)
ax.set_xlim(-210, 260); ax.set_ylim(-25, 655)
for ax in axs:
    ax.set_aspect('equal'); ax.axis('off')
fig.subplots_adjust(top=.93, bottom=.03, wspace=.15)
for path in (PAPER/'horizontal_shaft_support_paper.pdf', RES/'horizontal_shaft_support_paper.png'):
    fig.savefig(path, dpi=180, facecolor=fig.get_facecolor(), bbox_inches='tight', pad_inches=.04)
print(f'R = {R:.2f} mm, D = {2*R:.1f} mm, steel insert {L_STEEL:.1f} mm, tungsten {L_W:.1f} mm, carrier {2*HALF_T:.2f} mm')
