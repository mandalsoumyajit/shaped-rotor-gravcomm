"""Compact conceptual dumbbell diagram; no hardware performance assumptions."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc

out = Path(__file__).resolve().parents[1]/'results/paper_revision'
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'],
    'mathtext.fontset':'stix','font.size':7,'axes.labelsize':7,
    'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5,
    'axes.linewidth':.5,'lines.linewidth':.9,'pdf.fonttype':42})
fig, (ax, wave) = plt.subplots(1,2,figsize=(6.8,2.35),gridspec_kw={'width_ratios':[1.1,1]})
angle=np.pi/5
p=np.array([np.cos(angle),np.sin(angle)])
ax.add_patch(Circle((0,0),1,fill=False,ls='--',lw=.5,color='.55'))
ax.plot([-p[0],p[0]],[-p[1],p[1]],color='.15',lw=1.7)
ax.scatter([p[0],-p[0]],[p[1],-p[1]],s=65,c='#4b6475',zorder=3)
ax.plot(0,0,'+',color='k',ms=6)
ax.plot([-1.15,3.25],[0,0],color='.7',lw=.5)
ax.scatter(3,0,marker='s',s=37,color='#2877a5')
ax.text(p[0]+.12,p[1]+.08,r'$m=M/2$')
ax.text(-p[0]-.04,-p[1]-.26,r'$m=M/2$',ha='center')
ax.text(.37,.40,r'$r$');ax.text(-.38,-.38,r'$r$')
ax.add_patch(Arc((0,0),.8,.8,theta1=0,theta2=36,lw=.6))
ax.text(.50,.11,r'$\phi$')
ax.annotate('',xy=(3,-.48),xytext=(0,-.48),arrowprops={'arrowstyle':'<->','lw':.6})
ax.text(1.5,-.42,r'$d$',ha='center')
ax.text(3,-.22,'Receiver',ha='center');ax.text(3.15,.15,r'$x$')
ax.annotate('',xy=(-.28,.77),xytext=(.25,.78),arrowprops={'arrowstyle':'->','connectionstyle':'arc3,rad=.35','lw':.7})
ax.text(-.05,1.05,r'$\omega$');ax.text(-.12,.13,'Axis',ha='right')
ax.set(xlim=(-1.35,3.5),ylim=(-1.13,1.25));ax.set_aspect('equal');ax.axis('off')
ax.set_title('(a) Idealized equal-mass dumbbell',fontsize=7,pad=2)
phi=np.linspace(0,2*np.pi,1024,endpoint=False)
# Dimensionless r=1, d=3, G=1, total mass=1. Mean removed.
d=3.;a=np.zeros_like(phi)
for sign in (-1,1):
    dx=sign*np.cos(phi)-d; y=sign*np.sin(phi)
    a += .5*dx/(dx*dx+y*y)**1.5
a-=a.mean();quad=9/(4*d**4)
np.testing.assert_allclose(a,np.roll(a,512),rtol=1e-12,atol=1e-14)
wave.plot(phi/(2*np.pi),a/quad,label='Exact point-mass field')
wave.plot(phi/(2*np.pi),-np.cos(2*phi),'--',label='Leading quadrupole')
wave.set(xlim=(0,1),xlabel='Mechanical revolutions',ylabel=r'$[a_x-\overline{a_x}]/A_{2,\mathrm{quad}}$')
wave.set_title('(b) Two signal cycles per revolution ($d/r=3$)',fontsize=7,pad=4)
wave.legend(loc='upper center',bbox_to_anchor=(.5,-.29),ncol=2,frameon=False,handlelength=1.5,fontsize=6)
wave.grid(alpha=.16,lw=.4);wave.tick_params(width=.5,length=2.5)
fig.tight_layout(pad=.45,w_pad=1.3)
fig.savefig(out/'introductory_dumbbell.pdf',bbox_inches='tight')
fig.savefig(out/'introductory_dumbbell.png',dpi=180,bbox_inches='tight')
