"""Plot retained capacity calculations and the defined CF-FSK rate/torque family."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
data=json.loads((root/'results/paper_revision/gaussian_benchmark.json').read_text())
cf=json.loads((root/'results/desktop_cf_fsk/summary.json').read_text())
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],
 'mathtext.fontset':'stix','font.size':7,'axes.labelsize':7,'axes.titlesize':7,
 'legend.fontsize':6.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,
 'axes.linewidth':.5,'lines.linewidth':.9,'pdf.fonttype':42})
fig,(a,b)=plt.subplots(1,2,figsize=(7,2.45))
d=np.array([r['distance_m'] for r in data['rows']])
c=np.array([r['gaussian_benchmark_bit_s'] for r in data['rows']])
a.semilogy(d,c,'o-',ms=3,label='Gaussian capacity benchmark')
a.semilogy(.5,1/3,'s',ms=4,color='#c05c20',label='Selected CF-FSK point')
for x,y,label in zip(d,c,['3.903','0.1549','0.000881']):
 a.annotate(label,(x,y),xytext=(5,5),textcoords='offset points',fontsize=6.5)
a.set(xlabel='Receiver radius d (m)',ylabel='Rate (bit/s)',
 title='(a) Capacity and selected data rate',xlim=(.43,2.35),ylim=(3e-4,12),xticks=[.5,1,1.5,2])
a.legend(loc='upper right',frameon=False)
# CAD inertia: verified manuscript value (the plotted expression is analytical).
I=cf["inertia_kg_m2"];rate=np.linspace(.10,.75,250);ts=8/(3*rate)
torque=15*np.pi*I*6/(8*.2*.8*ts**2)
b.plot(rate,torque,label='Peak inertial torque')
b.axhline(2,ls='--',color='.4',lw=.7,label='Assumed peak drive limit')
for t in (16,8,4):
 r=8/(3*t);v=15*np.pi*I*6/(8*.2*.8*t*t)
 b.plot(r,v,'o',ms=3,color='#1f77b4')
 b.annotate(f'{t} s',(r,v),xytext=(-3,7),textcoords='offset points',ha='center',fontsize=6.5)
b.plot(cf["payload_bit_s"],cf["diagnostics"]["peak_torque_nm"],'s',ms=4,color='#c05c20',label='Selected total drive torque')
b.set(xlabel='Mapped payload rate (bit/s)',ylabel='Peak torque (N m)',
 title='(b) Seven-tone CF-FSK mechanical demand',xlim=(.08,.77),ylim=(0,6))
b.legend(loc='upper left',frameon=False)
for ax in (a,b):ax.grid(alpha=.17,lw=.4);ax.tick_params(width=.5,length=2.5)
fig.tight_layout(pad=.6,w_pad=1.5)
for ext in ('pdf','png'):fig.savefig(root/f'results/paper_revision/capacity_rate.{ext}',dpi=220,bbox_inches='tight',pad_inches=.025)
