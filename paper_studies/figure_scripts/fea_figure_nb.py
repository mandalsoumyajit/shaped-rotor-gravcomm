"""FEA geometry and results figure for the NB-LDPC joint design rotors (fea_geometry_results_current.pdf).

Same layout and plotting code as AIP_Advances_submission/scripts/restore_exposition_figures.py, with the
2 m rotor of public_release_v0_5_0/fea/results/three_link_designs_2026_11_05_nb (s = 2.854447700523453).
Panels (a)-(c) read the h0.013/design3_high mesh (model.inp) and integration-point fields (fields.npz),
copied from the WSL run folder (the full solver outputs are not shipped); pass that folder as argv[1].
Panels (d)-(e) use the compiled fine-mesh values (three_link_designs_2026_11_05_nb/RESULTS.md); panel (f)
is the scale-independent s=1, 1800 rpm reference convergence check (unchanged).
Usage: python figure_scripts/fea_figure_nb.py <folder with model.inp and fields.npz>
"""
import ast, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

PAPER = Path(__file__).resolve().parents[1]
AIP = PAPER.parent if (PAPER.parent/'paper_studies').is_dir() else PAPER.parent/'AIP_Advances_submission'   # repo or working tree
module = ast.parse((AIP/'fea/postprocess_carrier.py').read_text())
parser = next(node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == 'read_input')
exec(compile(ast.Module(body=[parser], type_ignores=[]), 'read_input', 'exec'))
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'legend.fontsize': 8, 'pdf.fonttype': 42})
SCALE_2M = 2.854447700523453
case = Path(sys.argv[1])
nodes, elements, materials = read_input(case/'model.inp')
fig, axs = plt.subplots(2, 3, figsize=(10, 6.6), layout='constrained')
segments = []
for eid, conn in elements.items():
    if materials[eid] != 'ALUMINUM':
        continue
    xyz = np.array([nodes[n] for n in conn[:4]])
    if np.min(xyz[:, 2]) <= 0 <= np.max(xyz[:, 2]):
        for i, j in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
            segments.append(xyz[[i, j], :2]*1000)
axs[0, 0].add_collection(LineCollection(segments, colors='#52677a', linewidths=.22))
axs[0, 0].set(xlim=(-125, 125), ylim=(-125, 125), xlabel='x (mm)', ylabel='y (mm)', title='(a) 2 m rotor: hub mesh section')
axs[0, 0].set_aspect('equal')
samples = np.load(case/'fields.npz')['stress']; ss = samples[:, 2:8]
vm = np.sqrt(((ss[:, 0] - ss[:, 1])**2 + (ss[:, 1] - ss[:, 2])**2 + (ss[:, 2] - ss[:, 0])**2)/2 + 3*np.sum(ss[:, 3:]**2, axis=1))/1e6
unique, inv = np.unique(samples[:, 0].astype(int), return_inverse=True); means = np.bincount(inv, weights=vm)/np.bincount(inv)
xy = []; vals = []
for eid, v in zip(unique, means):
    if materials[eid] != 'ALUMINUM':
        continue
    xyz = np.array([nodes[n] for n in elements[eid][:4]])
    if np.min(xyz[:, 2]) <= 0 <= np.max(xyz[:, 2]):
        xy.append(xyz[:, :2].mean(axis=0)); vals.append(v)
xy = np.array(xy)
sc = axs[0, 1].scatter(xy[:, 0], xy[:, 1], c=vals, s=2, cmap='magma', vmin=0, vmax=46, rasterized=True)
axs[0, 1].set(xlabel='x (m)', ylabel='y (m)', title='(b) Centrifugal aluminum stress'); axs[0, 1].set_aspect('equal')
fig.colorbar(sc, ax=axs[0, 1], label='Element-mean von Mises (MPa)', shrink=.8)
uarr = np.load(case/'fields.npz')['u']; coords = np.array([nodes[int(n)] for n in uarr[:, 0]])
zmax = .01*SCALE_2M; sel = np.isclose(coords[:, 2], zmax, atol=1e-7)
sc = axs[0, 2].scatter(coords[sel, 0], coords[sel, 1], c=np.linalg.norm(uarr[sel, 1:], axis=1)*1e6, s=2, cmap='viridis', rasterized=True)
axs[0, 2].set(xlabel='x (m)', ylabel='y (m)', title='(c) Centrifugal displacement'); axs[0, 2].set_aspect('equal')
fig.colorbar(sc, ax=axs[0, 2], label='Displacement magnitude (μm)', shrink=.8)
d = np.array([.5, 1, 2]); ax = axs[1, 0]
ax.plot(d, [8.54, 25.40, 47.85], 'o-', label='Horizontal combined bound')
ax.plot(d, [7.90, 24.16, 45.25], 's--', label='Centrifugal peak')
ax.axhline(120, ls=':', color='k', label='Static screening limit')
ax.set(xlabel='Link distance (m)', ylabel='Aluminum stress (MPa)', title='(d) Scaled steel rotors', xticks=d, ylim=(0, 130)); ax.legend(fontsize=6.7)
ax = axs[1, 1]
for j, v in enumerate([[82.36, 48.22, 27.39], [100.90, 58.57, 32.94], [112.29, 64.23, 35.49]]):
    ax.plot(d, v, 'o-', label=f'Mode {j + 1}' + (' (angular coupled)' if j == 2 else ''))
ax.set(xlabel='Link distance (m)', ylabel='Prestressed frequency (Hz)', title='(e) Center-speed carrier modes', xticks=d); ax.legend(fontsize=7)
ax = axs[1, 2]; ax.plot([18, 13, 10], [50.1718, 58.12655, 58.49185], 'o-')
ax.set(xlabel='Reference nominal mesh size (mm)', ylabel='Peak aluminum stress (MPa)', title='(f) Reference convergence, 1800 rpm', xticks=[10, 13, 18]); ax.invert_xaxis()
for ax in axs[1]:
    ax.grid(alpha=.2)
fig.savefig(PAPER/'fea_geometry_results_current.pdf', bbox_inches='tight', pad_inches=.04)
fig.savefig(PAPER/'ground_motion_study/results/fea_geometry_results_nb.png', dpi=150, bbox_inches='tight', pad_inches=.04)
print('max element-mean Al stress in slice', max(vals))
