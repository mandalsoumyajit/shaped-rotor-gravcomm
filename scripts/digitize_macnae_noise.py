"""Reproducible manual pixel digitization of the user-supplied Macnae spectrum."""
from pathlib import Path
import json,csv,shutil,hashlib,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/macnae_noise_digitization_2026_09_21'
# Manually selected centerline pixels. Vertical power-line spikes, annotation
# arrows and the diagonal constant-voltage coil line are excluded.
LOW=[(94,153),(105,175),(116,201),(124,216),(129,186),(133,190),(141,224),(149,253),(157,277),(160,280),(165,253),(169,263),(175,293),(183,319),(191,344),(200,365),(210,388),(220,410),(230,427),(242,440),(254,449),(267,453),(280,454),(291,452),(302,446),(311,427),(318,400),(324,376),(329,367),(335,366),(341,372),(347,389),(353,396),(359,395),(365,387),(371,388),(379,393),(387,397),(393,403),(397,401),(401,383),(404,381),(408,396),(412,414),(419,423),(427,428),(435,439),(443,448),(454,459),(462,471),(474,491),(482,496),(493,510),(503,527),(515,550),(527,575)]
HIGH=[(573,574),(576,549),(582,530),(590,516),(600,507),(611,503),(622,502),(633,504),(645,511),(657,522),(669,536),(681,549),(690,560)]
# Multi-tick least-squares log-axis calibration, original 744 x 810 screenshot.
XT=[(94,-2),(185,-1),(276,0),(367,1),(458,2),(549,3),(640,4)]
YT=[(36,1),(127,0),(218,-1),(309,-2),(400,-3),(491,-4),(582,-5)]
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 parser=argparse.ArgumentParser();parser.add_argument('--source-image',type=Path,required=True);src=parser.parse_args().source_image
 shutil.copy2(src,OUT/'source_screenshot.png')
 xc=np.polyfit(*np.array(XT).T,1);yc=np.polyfit(*np.array(YT).T,1)
 rows=[]
 for name,pts in [('low',LOW),('high',HIGH)]:
  for x,y in pts:
   rows.append(dict(segment=name,x_pixel=x,y_pixel=y,frequency_Hz=10**np.polyval(xc,x),asd_T_sqrtHz=1e-9*10**np.polyval(yc,y),trace_style='dashed' if name=='low' and 200<=x<=318 else 'solid'))
 with (OUT/'anchors.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (OUT/'calibration.json').write_text(json.dumps(dict(x_ticks=XT,y_ticks=YT,x_log10_coefficients=xc.tolist(),y_log10_nT_coefficients=yc.tolist(),pixel_selection='manual centerline',source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),vertical_pixel_uncertainty=3,uncertainty_note='3 pixels is a tracing tolerance only, not environmental or source uncertainty',source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11951293/',original_doi='10.1190/1.1441739'),indent=2))
 fig,ax=plt.subplots(figsize=(9,8));ax.imshow(plt.imread(src));
 for name,pts in [('low',LOW),('high',HIGH)]:
  p=np.array(pts);ax.plot(p[:,0],p[:,1],'-o',ms=2,lw=.9,label=name+' trace')
 ax.set_xlim(70,705);ax.set_ylim(600,120);ax.legend();ax.set_title('Manual trace overlay; interference spikes excluded')
 fig.savefig(OUT/'trace_overlay.png',dpi=150);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
 for seg in ('low','high'):
  a=[r for r in rows if r['segment']==seg];f=np.array([r['frequency_Hz'] for r in a]);b=np.array([r['asd_T_sqrtHz'] for r in a])*1e12
  ax.loglog(f,b,'o-',ms=2,color='navy');ax.fill_between(f,b/10**(3/91),b*10**(3/91),alpha=.2,color='navy')
 ax.set(xlabel='Frequency (Hz)',ylabel='Horizontal magnetic ASD (pT / sqrt(Hz))',title='Macnae historical nominal background: digitized broadband trace')
 ax.axvline(18,color='gray',ls=':',lw=.8);ax.axvline(24,color='gray',ls=':',lw=.8);ax.grid(True,which='both',alpha=.2)
 fig.savefig(OUT/'digitized_spectrum.png',dpi=160);fig.savefig(OUT/'digitized_spectrum.pdf');plt.close(fig)
if __name__=='__main__':main()
