"""Compact vector journal figures from the checked communications records."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter,NullFormatter
from gravcomm.receivers import CARTER_OSCILLATOR as rec
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/paper_revision'
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],'mathtext.fontset':'stix',
 'font.size':7,'axes.labelsize':7,'axes.titlesize':8,'legend.fontsize':6,
 'xtick.labelsize':6.5,'ytick.labelsize':6.5,'lines.linewidth':.9,'axes.linewidth':.6,
 'pdf.fonttype':42,'savefig.pad_inches':.025})
def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.png'),dpi=220,bbox_inches='tight')
    plt.close(fig)

fig,ax=plt.subplots(figsize=(3.35,2.32));ax.set(xlim=(0,3.35),ylim=(0,2.32));ax.axis('off')
def box(x,y,w,h,text,shade=False):
    ax.add_patch(Rectangle((x,y),w,h,facecolor='#eef2f5' if shade else 'white',edgecolor='0.2',lw=.6))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=7)
def arrow(a,b): ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',lw=.7,shrinkA=1,shrinkB=1))
box(.05,1.85,1.53,.4,'Calibrated acceleration\nComplex baseband')
box(.05,1.23,1.53,.36,'FIR whitening')
box(.05,.54,1.53,.46,'Viterbi sequence detector',True)
box(.05,.02,1.53,.3,'Recovered bytes')
arrow((.815,1.85),(.815,1.59));arrow((.815,1.23),(.815,1.0));arrow((.815,.54),(.815,.32))
box(1.92,1.62,1.37,.55,'Noise PSD\nFilter and template design')
box(1.92,.61,1.37,.74,'Acquired initial state\nTiming and carrier\nByte / trailer constraints')
arrow((1.92,1.72),(1.58,1.45));arrow((1.92,1.64),(1.58,.9));arrow((1.92,.82),(1.58,.74))
fig.subplots_adjust(left=0,right=1,bottom=0,top=1);save(fig,'detector')

rows=sorted([json.loads(p.read_text()) for p in (ROOT/'results/spacing_sweep').glob('*_fs8_m3.json')],key=lambda r:r['tone_spacing_hz'])
cap=[r for r in json.loads((ROOT/'results/communications_checks/capacity.json').read_text()) if r['bandwidth_hz']==2]
fig,axes=plt.subplots(1,3,figsize=(7.05,2.48),layout='constrained')
ax=axes[0]
ax.loglog([r['distance_m'] for r in cap],[r['capacity_bit_s'] for r in cap],'o-',ms=3,label='Gaussian benchmark')
ax.plot(.5,4/3,'s',ms=3.4,color='#9b481d',label='Mapped rate')
ax.plot(.5,64/60,'s',ms=3.4,mfc='white',mec='#9b481d',label='With trailer')
ax.set(xlabel='Receiver distance (m)',ylabel='Rate (bit/s)',title='(a) Capacity and verified rate',xlim=(.43,2.2),ylim=(5e-4,8))
ax.set_xticks([.5,1,2],['0.5','1','2']);ax.xaxis.set_minor_formatter(NullFormatter())
ax.legend(loc='lower left',frameon=False,handlelength=1.4,borderpad=.15,labelspacing=.3)
ax=axes[1]
for r in rows:
    lo,hi=r['packet_aware_anytime95_ber_interval'];x=r['tone_spacing_hz'];y=r['observed_ber']
    if r['adequate_error_count']: ax.errorbar(x,y,yerr=[[y-lo],[hi-y]],fmt='o',ms=2.7,color='#19547c',capsize=1.5,lw=.8)
    else: ax.errorbar(x,hi,yerr=hi*.48,uplims=True,fmt='o',ms=2.7,color='#9b481d',capsize=1.5)
ax.axhline(.001,color='0.3',ls='--',lw=.7)
ax.text(.88,.00118,'Target',ha='right',fontsize=6)
ax.set(xscale='log',yscale='log',xlabel='Tone spacing (Hz)',ylabel='Uncoded BER',title='(b) Receiver: 1.75 s symbols',ylim=(8e-5,.3))
ax=axes[2];x=np.array([r['tone_spacing_hz'] for r in rows])
for key,label,style,col in [('peak_torque_nm','Peak','-','#19547c'),('worst_message_rms_nm','Worst-message RMS','--','#9b481d'),('random_message_rms_nm','Uniform-byte RMS',':','#43774d')]:
    ax.loglog(x,[r['mechanics'][key] for r in rows],ls=style,color=col,label=label)
for value,label in [(2,'Peak rating'),(.5,'Continuous rating')]:
    ax.axhline(value,color='0.5',ls='-.',lw=.6)
    ax.text(.88,value*1.08,label,ha='right',fontsize=5.8,color='0.3')
ax.set(xlabel='Tone spacing (Hz)',ylabel='Torque (N m)',title='(c) Drive: same waveforms',ylim=(.11,22))
ax.legend(loc='upper left',frameon=False,handlelength=1.6,labelspacing=.2)
for ax in axes[1:]:
    ax.set_xlim(.043,1.05);ax.set_xticks([.05,.1,.2,.5,1],['0.05','0.1','0.2','0.5','1'])
    ax.xaxis.set_minor_formatter(NullFormatter());ax.axvline(.095238095,color='0.6',ls=':',lw=.6)
for ax in axes:ax.grid(alpha=.13,lw=.4)
save(fig,'communications_performance')

f=np.linspace(49.5,51.1,1500)
fig,ax=plt.subplots(figsize=(3.35,1.85),layout='constrained')
for label,y in [('Thermal',rec.thermal_acceleration_asd(f)),('Readout, input referred',rec.readout_acceleration_asd(f)),('Total',rec.acceleration_noise_asd(f))]:
    ax.semilogy(f-50.3,y/1e-11,label=label)
tones=50.3+np.arange(-3,4)/12
ax.plot(tones-50.3,rec.acceleration_noise_asd(tones)/1e-11,'ko',ms=2.5)
ax.set(xlabel='Frequency offset from resonance (Hz)',ylabel=r'Noise ASD (nGal/$\sqrt{\mathrm{Hz}}$)',ylim=(.5,50))
ax.legend(frameon=False,loc='lower right');ax.grid(alpha=.15,lw=.4);save(fig,'receiver')
