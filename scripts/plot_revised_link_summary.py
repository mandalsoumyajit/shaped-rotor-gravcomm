from pathlib import Path
import sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.receivers import StructuralOscillator
f=np.geomspace(1,80,500);c=SeismicBackground().components(f);floor=float(StructuralOscillator().acceleration_noise_asd(50.3));n=np.sqrt(floor**2+c['background_upper_psd'])
fig,axs=plt.subplots(1,3,figsize=(10,3.4));plt.rcParams.update({'font.size':9})
d=np.array([.5,1,2]);m=[7.027334750750038,37.69173275922457,227.3153385558173];torque=[1.597850216,26.26002699,524.7263804];power=[.120492,1.980244,29.680143]
axs[0].semilogy(d,m,'o-',color='#20639b');axs[0].set(xlabel='Axis-to-receiver distance (m)',ylabel='Bare rotor mass (kg)',xticks=d,title='(a) Qualified tested source sizes')
for x,y in zip(d,m):axs[0].annotate(f'{y:.2f}',(x,y),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
axs[0].set_ylim(4,450)
axs[1].semilogy(d,torque,'o-',label='Torque (N m)',color='#20639b');axs[1].semilogy(d,power,'s--',label='Power (kW)',color='#c55a11');axs[1].set(xlabel='Axis-to-receiver distance (m)',ylabel='Peak command (units in legend)',xticks=d,title='(b) Inertial drive demand');axs[1].legend(fontsize=8,loc='upper left')
ax=axs[2];ax.axvspan(10,80,color='#ebebeb');ax.loglog(f,n/1e-11,label='Total',color='black');ax.loglog(f,np.full_like(f,floor)/1e-11,label='Instrument',ls='--');ax.loglog(f,c['support_acceleration_asd']/1e-11,label='Support');ax.loglog(f,c['newtonian_acceleration_asd']/1e-11,label='Newtonian');
for fc in [18,24]:
 nc=SeismicBackground().components(fc);y=float(np.sqrt(floor**2+nc['background_upper_psd'])/1e-11);ax.plot(fc,y,'ko',ms=4)
ax.set(xlabel='Tuned carrier frequency (Hz)',ylabel='Center ASD (nGal / sqrt(Hz))',title='(c) Assumed noise scenario',xlim=(1,80));ax.legend(fontsize=7,loc='lower left')
for ax in axs:ax.grid(alpha=.18,which='major');ax.tick_params(labelsize=8)
fig.tight_layout();fig.savefig(ROOT/'figures/revised_link_summary.pdf');fig.savefig(ROOT/'figures/revised_link_summary.png',dpi=180)
