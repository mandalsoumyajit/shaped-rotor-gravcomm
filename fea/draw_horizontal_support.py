"""Dimensioned schematic of the preliminary horizontal-shaft concept."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,FancyArrowPatch,Arc
OUT=Path(__file__).resolve().parent/'results/tungsten_comparison_2026_09_20'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,2,figsize=(15,8),gridspec_kw={'width_ratios':[1,1.15]});fig.patch.set_facecolor('#f7f9fc')
blue='#377eb8';dark='#18354b';gold='#b98a36';gray='#718096';base='#d9e2ec'
def rect(ax,x,y,w,h,c,**kw):
 p=Rectangle((x,y),w,h,facecolor=c,edgecolor=dark,lw=1.2,**kw);ax.add_patch(p);return p
def label(ax,text,xy,at,align='left'):
 ax.annotate(text,xy=xy,xytext=at,ha=align,va='center',fontsize=10,color=dark,arrowprops=dict(arrowstyle='-',color=dark,lw=.9),bbox=dict(facecolor='#f7f9fc',edgecolor='none',pad=2))
def dim(ax,a,b,text,offset=8):
 ax.annotate('',a,b,arrowprops=dict(arrowstyle='<->',color=dark,lw=.8));ax.text((a[0]+b[0])/2,(a[1]+b[1])/2+offset,text,ha='center',color=dark,fontsize=9)
R=240.296;center=300.;ri=221.073;s=.9611854132640992
ax=axs[0]
ax.add_patch(Circle((0,center),R,fc=blue,ec=dark,lw=1.5));ax.add_patch(Circle((0,center),ri,fc='#f7f9fc',ec=dark,lw=1))
rect(ax,-ri,center-15*s,2*ri,30*s,blue);rect(ax,-15*s,center-ri,30*s,2*ri,blue)
for x in [-200*s,200*s]:
 ax.add_patch(Circle((x,center),50*s,fc=blue,ec=dark));ax.add_patch(Circle((x,center),35*s,fc=gold,ec=dark))
ax.add_patch(Circle((0,center),30*s,fc=gray,ec=dark));ax.add_patch(Circle((0,center),10,fc='white',ec=dark));ax.plot(0,center,'o',color=dark,ms=3)
ax.add_patch(Arc((0,center),160,160,theta1=35,theta2=125,color=dark,lw=1.5));ax.annotate('',(-46,365),(-32,373),arrowprops=dict(arrowstyle='->',color=dark))
label(ax,'Two dense inserts',(-192,300),(-265,610));label(ax,'Aluminum carrier + rim',(130,480),(60,590));label(ax,'Hub / shaft axis\n(points toward viewer)',(0,300),(35,160))
dim(ax,(-R,30),(R,30),'Rotor diameter 481 mm',12)
ax.set_title('A  Looking along the shaft\nRotor turns in a vertical plane',loc='left',fontsize=13,color=dark,pad=18)
ax.set_xlim(-285,285);ax.set_ylim(-25,655)
ax=axs[1]
rect(ax,-205,0,390,18,base)
# Stationary pedestals are outside the axial swept envelope.
for x in [-60,60]:
 rect(ax,x-12,18,24,265,base);rect(ax,x-20,283,40,34,base)
 rect(ax,x-13,289,26,22,'white')
 if x<0:
  ax.plot([x-11,x-2],[291,309],color=blue,lw=2);ax.plot([x+2,x+11],[309,291],color=blue,lw=2)
 else:
  rect(ax,x-9,289,18,4,blue);rect(ax,x-9,307,18,4,blue)
# Stub shafts, hub face flanges and schematic drive coupling.
rect(ax,-92,290,82,20,gray);rect(ax,10,290,75,20,gray)
rect(ax,-14.61,271.16,5,57.68,gray);rect(ax,9.61,271.16,5,57.68,gray)
rect(ax,-9.61,center-R,19.22,2*R,blue)
for y in [center-192.237,center+192.237]:
 ax.add_patch(Rectangle((-38.08,y-33.64),76.16,67.28,fc='none',ec=gold,lw=1.3,ls='--'))
 rect(ax,-16.61,y-33.64,33.22,67.28,gold)
rect(ax,-115,284,23,32,gray);rect(ax,-190,270,75,60,base);rect(ax,-181,18,55,252,base)
ax.plot([-202,105],[300,300],ls='-.',color=dark,lw=.7)
label(ax,'Motor + coupling',(-145,330),(-195,640))
label(ax,'Drive-side locating pair\n(back-to-back angular contact)',(-60,311),(-195,585))
label(ax,'Floating roller bearing\n(allows axial expansion)',(60,310),(85,420))
label(ax,'Bolted hub flanges\non two stub shafts',(14,325),(85,355))
label(ax,'Tungsten insert', (16.6,495),(85,515))
label(ax,'Dashed: steel insert\nsame mass, longer',(-38,485),(-195,450))
label(ax,'Stationary pedestals',(-60,120),(-195,100))
label(ax,'Common rigid base',(100,15),(85,75))
dim(ax,(-60,240),(60,240),'120 mm bearing span',12)
ax.annotate('Gravity',xy=(155,175),xytext=(155,230),ha='center',color=dark,arrowprops=dict(arrowstyle='->',color=dark,lw=1.5))
ax.set_title('B  Edge view of the rotor\nHorizontal shaft supported on both sides',loc='left',fontsize=13,color=dark,pad=18)
ax.set_xlim(-210,260);ax.set_ylim(-25,655)
for ax in axs:ax.set_aspect('equal');ax.axis('off')
fig.suptitle('Horizontal-shaft support concept',fontsize=21,color=dark,x=.06,ha='left',y=.99)
fig.text(.06,.055,'Illustrated with the 0.5 m-link rotor. Dimensions in mm; supports and drive are schematic.\nViews show different rotor angles (90 degrees apart). Retention details and protective enclosure omitted.',color=dark,fontsize=10)
fig.subplots_adjust(top=.85,bottom=.14,wspace=.15)
for ext in ('png','svg','pdf'):fig.savefig(OUT/f'horizontal_shaft_support.{ext}',dpi=180,facecolor=fig.get_facecolor())
