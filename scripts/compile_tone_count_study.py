"""Compile tone-count screens without converting zero pilot counts into BER claims."""
from pathlib import Path
import sys,json,csv
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from gravcomm.packet_confidence import packet_ber_interval
from tone_count_runtime import OUT,save
from tone_count_study import configurations

def compile_stage(stage):
 groups=defaultdict(list)
 for p in (OUT/stage).glob('*.json'):
  x=json.loads(p.read_text())
  if 'index' in x:groups[x['config']['id']].append(x)
 rows=[]
 for key,data in sorted(groups.items()):
  failures=[x for x in data if 'error' in x];good=[x for x in data if 'error' not in x];c=data[0]['config']
  row=dict(id=key,config=c,packets=len(good),runtime_failures=len(failures))
  if good:
   errors=[x['errors'] for x in good];e=sum(errors);n=len(good);spec=good[0]
   row.update(errors=e,bits=224*n,ber=e/(224*n),ber_interval=packet_ber_interval(errors,224),raw_errors=sum(x['raw_errors'] for x in good),rejections=sum(not x['accepted'] for x in good),wrong_accepted=sum(x['wrong_accepted'] for x in good),symbol_s=spec['symbol_s'],span_hz=spec['tone_span_hz'],torque_per_inertia=spec['peak_torque_per_inertia'],torque_ratio_to_published=spec['peak_torque_per_inertia']/6.161906907745526,pruned_states=sum(x['pruned_states'] for x in good),max_states=max(x['peak_states'] for x in good),max_block_s=max(x['max_block_s'] for x in good),max_latency_s=max(x['latency_budget_s'] for x in good),detector_s_mean=np.mean([x['total_detector_s'] for x in good]),numerical_checks=all(x['passes_detector_gate'] for x in good),convergence_checked=all(x['gate'] for x in good),padding_bits=spec['padding_bits'],acquisition_allowance_s=12*spec['symbol_s'])
  rows.append(row)
 return rows

def main():
 rows=compile_stage('screen')+compile_stage('extra');save(OUT/'screen_summary.json',rows)
 selection=[]
 for M in range(2,10):
  for fc in [24.,18.]:
   # Keep legacy and more efficient seven-tone packing separate.
   for L in ([3,4] if M==7 else [None]):
    options=[r for r in rows if r['config']['tone_count']==M and r['config']['carrier_hz']==fc and (L is None or r['config']['group_symbols']==L) and r.get('packets')==4 and r.get('numerical_checks') and not r['runtime_failures']]
    if options:
     best=min(options,key=lambda r:(r['errors']>0,r['errors'],r['torque_per_inertia']))
     selection.append(best['config'])
 save(OUT/'shortlist.json',selection)
 # Exact result counts; color is logarithmic only for visual dynamic range.
 mapping=[]
 for c in configurations():
  key=(c['tone_count'],c['group_symbols'],c['group_bits'])
  if key not in mapping:mapping.append(key)
 grid=[.4,.8,1.2,1.6];fig,axs=plt.subplots(1,2,figsize=(10,5.2),layout='constrained')
 for ax,fc in zip(axs,[24.,18.]):
  values=np.full((len(mapping),4),np.nan);texts={}
  for i,(M,L,b) in enumerate(mapping):
   for j,g in enumerate(grid):
    found=[r for r in rows if r['config']['tone_count']==M and r['config']['group_symbols']==L and r['config']['carrier_hz']==fc and r['config']['span_symbol_product']==g]
    if found and found[0].get('packets')==4:
     row=found[0];values[i,j]=np.log10(1+row['errors']);texts[i,j]=str(row['errors'])+('*' if not row['numerical_checks'] else '')
  im=ax.imshow(values,cmap='YlOrRd',vmin=0,vmax=3,aspect='auto')
  for (i,j),v in texts.items():ax.text(j,i,v,ha='center',va='center',fontsize=9,color='black')
  ax.set(xticks=range(4),xticklabels=grid,yticks=range(len(mapping)),yticklabels=[f'{M} tones; {b} bits/{L} symbols' for M,L,b in mapping],xlabel='Total tone span × symbol duration',title=f'{fc:g} Hz: payload errors / 896 bits')
 fig.colorbar(im,ax=axs,label='log10(1 + error count)',shrink=.8)
 fig.savefig(OUT/'tone_count_screen.pdf');fig.savefig(OUT/'tone_count_screen.png',dpi=140);plt.close(fig)
 lines=['# Tone-count screening results','',f'Completed {sum(r["packets"] for r in rows)} packet records across {len(rows)} configurations. Four packets/configuration; this is exploratory screening. Asterisks in the figure flag failed numerical/state-budget checks.','', 'The 24 Hz channel represents the qualified 0.5 and 1 m source amplitudes; 18 Hz represents the qualified 2 m source amplitude. Payload 224 bits, rate-2/3 code, frame budget 224 s, fixed receiver and background. Source amplitudes are fixed here; this is not a new minimum-size search.','', '## Lowest-error screened spacing for each mapping and carrier','', '| Tones | Mapping bits/symbols | Carrier Hz | Span × Ts | Span Hz | Errors / bits | Peak torque / published | Rejected | Diagnostic status |','|---:|:---:|---:|---:|---:|:---|---:|---:|:---|']
 for c in selection:
  row=next(r for r in rows if r['id']==c['id']);lines.append(f"| {c['tone_count']} | {c['group_bits']}/{c['group_symbols']} | {c['carrier_hz']:g} | {c['span_symbol_product']:g} | {row['span_hz']:.3f} | {row['errors']} / {row['bits']} | {row['torque_ratio_to_published']:.3f} | {row['rejections']} | Wider-history check pending |")
 lines+=['','These spacings minimize observed error count, breaking zero-error ties by torque. Four packets cannot resolve BER near 1e-3. Packet-aware confidence intervals are in screen_summary.json and are intentionally broad. Selection uses the screen only; later fresh packets use different seeds.','', 'Mapping lengths change symbol time and acquisition/flush allowance duration at fixed frame time. Alphabet size, grouping, spacing, and overhead are therefore reported together. Larger alphabets with pruning or missing alternatives remain numerically unresolved. All candidate waveforms use the same causal filter; tone span beyond 2 Hz incurs its actual attenuation rather than being silently clipped.']
 (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
 print('Completed packets',sum(r['packets'] for r in rows),'shortlisted',len(selection))
if __name__=='__main__':main()
