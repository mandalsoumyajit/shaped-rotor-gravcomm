"""Primary comparison: magnetic and gravitational carriers aligned."""
import json
import numpy as np
from scipy.special import erfcinv,logsumexp
from scipy.optimize import brentq
from magnetic_wall_study import Channel,DESIGNS,ROOT,RINFO
from gravcomm.magnetic_link import minimum_capacity_logpower,noise_asd
OUT=ROOT/'results/magnetic_wall_aligned_carriers_2026_09_21'

def bpsk_power(c,sigma,fc,ambient,n=512):
 rate=256/224;alpha=.25;half=rate*(1+alpha)/2
 z=(np.arange(n)+.5)/n*2*half-half
 p=np.ones_like(z)/rate;edge=abs(z)>(1-alpha)*rate/2
 p[edge]=.5/rate*(1+np.cos(np.pi/(alpha*rate)*(abs(z[edge])-(1-alpha)*rate/2)))
 return float(np.log(erfcinv(.002)**2*rate)+logsumexp(np.log(p/p.sum())-c.gain(fc+z,sigma,ambient=ambient)))

def capacity_power(c,sigma,fc,ambient,n=2048):
 edges=np.linspace(fc-1,fc+1,n+1);f=(edges[1:]+edges[:-1])/2
 return minimum_capacity_logpower(c.gain(f,sigma,ambient=ambient),np.diff(edges),RINFO)

def crossing(fun,budget):
 return float(10**brentq(lambda x:fun(10**x)-np.log(budget),-4,12,xtol=1e-8))

def main():
 OUT.mkdir(parents=True,exist_ok=True);rows=[]
 for design in DESIGNS:
  c=Channel(design,design['distance']);fc=18. if design['distance']==2 else 24.
  for pt in (.38,1.,3.14,10.):
   ambient=pt*1e-12
   bf=lambda s:bpsk_power(c,s,fc,ambient)
   cf=lambda s:capacity_power(c,s,fc,ambient)
   bs=crossing(bf,design['power']);cs=crossing(cf,design['power'])
   b0=bf(5.8e7);c0=cf(5.8e7)
   rows.append(dict(distance_m=design['distance'],wall_thickness_m=c.thickness,carrier_Hz=fc,capacity_band_Hz=[fc-1,fc+1],background_pT_sqrtHz=pt,power_ceiling_W=design['power'],bpsk_boundary_S_m=bs,capacity_boundary_S_m=cs,copper_bpsk_W=float(np.exp(b0)),copper_capacity_W=float(np.exp(c0)),quadrature_change_dB=dict(bpsk=float((bpsk_power(c,bs,fc,ambient,n=1024)-bf(bs))*10/np.log(10)),capacity=float((capacity_power(c,cs,fc,ambient,n=4096)-cf(cs))*10/np.log(10)))))
  print(f'Completed {design["distance"]} m at {fc} Hz',flush=True)
 (OUT/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
 lines=['# Aligned-carrier magnetic-wall comparison','',
 'This is the primary comparison requested by the user. Magnetic carrier centers equal the qualified gravitational centers: 24 Hz at 0.5 and 1 m, and 18 Hz at 2 m. Earlier independently optimized sub-Hz magnetic results are supplementary bounds and do not define this comparison.', '',
 'The finite-hardware wall geometry is retained: thicknesses 0.227, 0.547 and 1.202 m. Source outer radii match the gravitational rotor radii; the magnetic source is parallel to the wall with its winding depth included. Noise comprises instrumental, wall thermal and ELF-paper background PSDs. The representative background is 1 pT/sqrt(Hz), bracketed by 0.38 and 3.14 pT/sqrt(Hz), with 10 pT/sqrt(Hz) as stress. These are broadband level assumptions. Neither link receives unmodeled gradiometric cancellation credit.', '',
 'The analytic coherent BPSK benchmark carries 224 payload plus 32 overhead bits in 224 s, roll-off 0.25, and target uncoded BER 1e-3. Its support extends +/-0.7143 Hz about the fixed center. It assumes ideal pre-equalization and whitening. The optimistic Gaussian capacity bound is restricted to a 2 Hz interval about that same center, matching the stated nominal receive bandwidth. This is a hard allowed spectral window for the bound; it does not reproduce the gravitational receiver filter transition band. Water filling may favor one side within this window. Its information-rate target is 1-h2(0.001) bit/s; finite-length coding and acquisition are unqualified.', '',
 'Power ceilings remain numerically equal to gravitational peak mechanical demand (120.49, 1980.24, 29680.14 W), while the magnetic model allows this as average source dissipation. Equal average electrical power remains unresolved. Model scope and omitted driver/thermal constraints follow the finite-hardware study.', '',
 '| Distance (m) | Center (Hz) | Background (pT/sqrt(Hz)) | BPSK crossing (S/m) | Band-limited capacity crossing (S/m) |','|---:|---:|---:|---:|---:|']
 for r in rows:lines.append(f'| {r["distance_m"]:g} | {r["carrier_Hz"]:g} | {r["background_pT_sqrtHz"]:g} | {r["bpsk_boundary_S_m"]:.4g} | {r["capacity_boundary_S_m"]:.4g} |')
 lines+=['','For the representative 1 pT/sqrt(Hz) background, at copper conductivity 5.8e7 S/m:','','| Distance (m) | BPSK required source W | Ideal minimum source W |','|---:|---:|---:|']
 for r in rows:
  if r['background_pT_sqrtHz']==1:lines.append(f'| {r["distance_m"]:g} | {r["copper_bpsk_W"]:.3g} | {r["copper_capacity_W"]:.3g} |')
 lines+=['','Large extrapolated powers diagnose infeasibility within the linear model; they are not physical drive designs. Failure of the capacity bound applies to the specified source, sensor, noise and allowed band. It does not exclude other electromagnetic frequencies, sensor types, geometries or communication pathways.','',f'Maximum absolute change from doubling spectral integration resolution: {max(abs(v) for r in rows for v in r["quadrature_change_dB"].values()):.4g} dB. Earlier finite-wall direct quadrature checks and the fixed-frequency skin-depth audit also apply.','', 'Reproduce: scripts/magnetic_wall_aligned_carriers.py. Raw data: summary.json.']
 (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 print('\n'.join(lines[-12:]))
if __name__=='__main__':main()
