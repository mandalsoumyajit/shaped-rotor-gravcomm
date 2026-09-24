"""Final screening report: separate raw screens, convergence, and fresh confirmation."""
from pathlib import Path
import sys,json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from compile_tone_count_study import compile_stage
from tone_count_runtime import OUT
from resume_tone_count_screen import resilient_save as save
from tone_count_study import configurations

def main():
 rows=compile_stage('screen')+compile_stage('extra');confirm=compile_stage('confirm');gates={}
 for f in (OUT/'gates').glob('*_53000.json'):
  x=json.loads(f.read_text());gates[x['config']['id']]=x
 save(OUT/'screen_summary.json',rows);save(OUT/'confirmation_summary.json',confirm)
 selected=json.loads((OUT/'shortlist.json').read_text())
 # Five columns include matched q controls, reusing identical grid points where applicable.
 mapping=[]
 for c in configurations():
  key=(c['tone_count'],c['group_symbols'],c['group_bits'])
  if key not in mapping:mapping.append(key)
 fig,axs=plt.subplots(1,2,figsize=(11,5.3),layout='constrained')
 for ax,fc in zip(axs,[24.,18.]):
  arr=np.full((len(mapping),5),np.nan);labels={}
  for i,(M,L,bits) in enumerate(mapping):
   for j,span in enumerate([.4,.8,1.2,1.6,.2*(M-1)]):
    found=[r for r in rows if r['config']['tone_count']==M and r['config']['group_symbols']==L and r['config']['carrier_hz']==fc and abs(r['config']['span_symbol_product']-span)<1e-10]
    if found and found[0].get('packets')==4:
     x=found[0];arr[i,j]=np.log10(1+x['errors']);labels[i,j]=str(x['errors'])+('*' if not x['numerical_checks'] else '')
  im=ax.imshow(arr,cmap='YlOrRd',vmin=0,vmax=3,aspect='auto')
  for (i,j),v in labels.items():ax.text(j,i,v,ha='center',va='center',fontsize=9,color='black')
  ax.set(xticks=range(5),xticklabels=['0.4','0.8','1.2','1.6','q=0.08'],yticks=range(len(mapping)),yticklabels=[f'{M} tones; {b} bits/{L} symbols' for M,L,b in mapping],xlabel='Span x symbol duration; last column: fixed q',title=f'{fc:g} Hz: payload errors / 896 bits')
 fig.colorbar(im,ax=axs,label='log10(1 + error count)',shrink=.8)
 fig.savefig(OUT/'tone_count_screen.pdf',bbox_inches='tight');fig.savefig(OUT/'tone_count_screen.png',dpi=130,bbox_inches='tight');plt.close(fig)
 lines=['# Tone-count study: results and interpretation','',f'Grid: {sum(r["packets"] for r in rows)} completed packets, {len(rows)} alphabet/mapping/spacing/carrier configurations. Fresh confirmation: {sum(r["packets"] for r in confirm)} packets. Validation: 19 tests passed (15 generalized-alphabet checks and 4 original streaming regressions).','', '## What the comparison establishes','', 'Tone count is a design variable. The three-symbol byte-mapping argument alone cannot justify seven tones as the best mechanical/communications choice. Several alternatives merit independent BER qualification. The existing seven-tone links retain their previously established qualification; this screen does not replace it.','', 'The screen fixes source amplitudes (42.09 nGal at 24 Hz; 51.80 nGal at 18 Hz), 224 payload bits, rate-2/3 code, 224 s frame budget and the instrumental/environmental receiver model. The 24 Hz result serves both 0.5 and 1 m source sizes. Transition fraction is 0.6. Grouping and symbol time change with alphabet size. Acquisition remains an allowance, with known initial phase/timing supplied to the simulation.','', '## Fresh confirmation','', '| Tones | Mapping bits/symbols | Carrier Hz | Tone span Hz | Errors / payload bits | CRC rejected / packets | Peak torque ratio | Detector gate |','|---:|:---:|---:|---:|:---|:---|---:|:---|']
 for x in sorted(confirm,key=lambda x:(x['config']['tone_count'],x['config']['group_symbols'],x['config']['carrier_hz'])):
  c=x['config'];g=gates.get(c['id'],{});status='pass' if g.get('passes_detector_gate') else 'pending / unresolved'
  lines.append(f"| {c['tone_count']} | {c['group_bits']}/{c['group_symbols']} | {c['carrier_hz']:g} | {x['span_hz']:.3f} | {x['errors']} / {x['bits']} | {x['rejections']} / {x['packets']} | {x['torque_ratio_to_published']:.3f} | {status} |")
 lines+=['', 'All tentative payload errors count, including rejected packets. Fresh seeds differ from screening and convergence seeds. Zero errors in 32 packets gives a conservative packet-aware 95% BER upper bound of approximately 0.154 with the saved confidence method, far above 1e-3. These data support candidate selection and expose mechanical tradeoffs; they cannot establish reliability superiority. Individual confidence intervals are in confirmation_summary.json. The table is exploratory and has no simultaneous confidence claim after candidate selection.','', '## Screening and detector limits','', 'The figure labels actual error counts; an asterisk means pruning or another per-run numerical check failed. Absence of an asterisk does not imply history/lookahead convergence; those separate diagnostics are summarized below. The equal-span grid uses span*Ts = 0.4, 0.8, 1.2, 1.6. The final column fixes adjacent spacing times dwell at q=0.08; duplicate points reuse records. This additional control prevents phase-lattice cost from silently eliminating viable eight-tone settings.','', '| Tones | Mapping | Carrier | Selected span*Ts | Pilot errors / 896 | Torque ratio | Convergence status |','|---:|:---:|---:|---:|---:|---:|:---|']
 for c in selected:
  x=next(r for r in rows if r['id']==c['id']);g=gates.get(c['id']);status='pilot only; higher torque than control' if g is None and c['tone_count']==9 else 'not checked' if g is None else 'pass' if g.get('passes_detector_gate') else 'requires receiver refinement'
  lines.append(f"| {c['tone_count']} | {c['group_bits']}/{c['group_symbols']} | {c['carrier_hz']:g} | {c['span_symbol_product']:g} | {x['errors']} | {x['torque_ratio_to_published']:.3f} | {status} |")
 lines+=['','History and lookahead are diagnostic parameters, not channel resources. The bounded screen uses four tone-history symbols and roughly the original lookahead duration in seconds. Wider-history checks use five symbols and up to 200000 states; separate checks double lookahead. Passing requires <1% relative soft-metric change, no changed bit signs or decoded bits, no pruning/missing alternatives, and whitening PSD error <1e-3. One diagnostic packet per point is a screening check; error-containing fresh packets and additional seeds should be replayed before qualification.','', '## Mechanical interpretation','', 'Peak inertial torque is proportional to I*pi*(15/8)*span/(beta*Ts). The ratio column uses the same physical rotor as the published design at each distance. It excludes windage, bearing/electrical losses, drivetrain compliance and acquisition transients. Lower peak torque does not establish lower average electrical consumption. Carrier frequencies and source sizes are unchanged in this study.','', 'The four-tone candidate has Ts=1.067 s, full tone span 0.5625 Hz and peak-torque ratio 0.8402, giving approximately 1.34, 22.06 and 440.87 N m for the existing 0.5, 1 and 2 m rotors. The five-tone candidate has Ts=1.179 s, span 0.6786 Hz and ratio 0.9170. These are commanded-motion comparisons at fixed source size; minimum source amplitudes have not yet been reoptimized.','', '## Scope and next decision','', 'This is a finite alphabet/spacing screen with natural radix mappings, at most four symbols and 12 coded bits per group. It is not a global mapping, alphabet, code-rate, packet-length or transition-fraction optimization. Unused radix words create mapping-dependent redundancy and tone occupancy; this is explicitly part of the compared protocols. Known final padding bits are enforced. Acquisition/flush symbol allowances have different durations as Ts changes and remain within the common 224 s frame budget.','', 'Recommended next step: independently qualify the four- and five-tone candidates, replay error packets with wider receiver settings, and refine their amplitude thresholds at the two normalized channels before changing manuscript design tables. Keep the seven-tone result as the already-qualified benchmark.','', '## Reproduction and progress','', 'Study plan: ../PLAN.md. Immutable job manifests and per-packet JSON are in screen/, extra/, gates/ and confirm/. Code snapshots and SHA256 dependencies are retained. Entry points: tone_count_study.py, tone_count_extra.py, tone_count_gates.py, tone_count_confirm.py, compile_tone_count_study.py and finalize_tone_count_study.py. resume_tone_count_screen.py and robust_tone_host.py only add retry handling for transient Windows/OneDrive checkpoint locks. Existing scientific records and qualified receiver implementation were not modified. Nine-tone points remain at the pilot stage: the unpruned zero-error setting has higher peak torque than the qualified seven-tone control and was not advanced to expensive wider-history diagnostics in this bounded follow-up.']
 (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 print('screen',sum(r['packets'] for r in rows),'confirm',sum(r['packets'] for r in confirm),'gates',len(gates))
if __name__=='__main__':main()
