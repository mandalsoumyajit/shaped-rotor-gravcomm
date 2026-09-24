"""Atomic receiver alternatives at aligned magnetic/gravity carriers."""
import json
import numpy as np
from magnetic_wall_study import Channel,DESIGNS,ROOT
from gravcomm.magnetic_wall import wall_channel
from magnetic_wall_aligned_carriers import bpsk_power,capacity_power,crossing
OUT=ROOT/'results/magnetic_wall_atomic_receiver_2026_09_21'
class AtomicChannel(Channel):
 def __init__(self,design,floor):
  super().__init__(design,design['distance'])
  self.floor=floor;self.rx_radius=.005
  H,L,E=wall_channel(self.product,1.,self.magnetic_distance,thickness=self.thickness,coil_radius=self.coil_radius,copper_mass=.1*design['mass'],stand_off=self.z1,receive_radius=self.rx_radius)
  _,_,Er=wall_channel(self.product,1.,self.magnetic_distance,thickness=self.thickness,coil_radius=self.rx_radius,copper_mass=1.,stand_off=self.z2)
  self.logH=np.log(np.maximum(abs(H),1e-300));self.logEr=np.log(np.maximum(Er/self.product,1e-300))
 def gain(self,f,sigma,ambient=1e-12,**kw):
  x=np.log(np.maximum(sigma*f,self.product[0]));h=np.interp(x,np.log(self.product),self.logH)
  e=np.exp(np.interp(x,np.log(self.product),self.logE))*f
  wall=4*1.380649e-23*300*np.exp(np.interp(x,np.log(self.product),self.logEr))/(4*np.pi**2*f)
  return 2*h-np.log(self.copper+e)-np.log(self.floor**2+ambient**2+wall)
def main():
 OUT.mkdir(parents=True,exist_ok=True);rows=[]
 for design in DESIGNS:
  fc=18. if design['distance']==2 else 24.
  for floor,label in [(35e-15,'optimistic compensated SERF'),(5e-12,'Mx benchmark')]:
   c=AtomicChannel(design,floor)
   for pt in (.38,1.,3.14,10.):
    ambient=pt*1e-12
    bs=crossing(lambda s:bpsk_power(c,s,fc,ambient),design['power'])
    cs=crossing(lambda s:capacity_power(c,s,fc,ambient),design['power'])
    rows.append(dict(distance_m=design['distance'],carrier_Hz=fc,sensor=label,instrument_ASD_T_sqrtHz=floor,background_pT_sqrtHz=pt,bpsk_boundary_S_m=bs,capacity_boundary_S_m=cs))
 (OUT/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
 lines=['# Atomic magnetometer alternative','',
 'Aligned carrier centers remain 24, 24 and 18 Hz; the ideal capacity window remains +/-1 Hz. Source geometry and power ceilings are unchanged. The 0.3403 m radius PCB receiving loop is replaced by a compact 5 mm radius disk-average field sensor, 15 mm from the wall face. This is an explicit spatial surrogate for a vapor cell; wall thermal noise is recalculated for this smaller sensing region. It is not a model of a particular packaged instrument.', '',
 'Two narrowband instrument scenarios are used: 35 fT/sqrt(Hz), motivated by a published SERF sensitivity of 33.5 fT/sqrt(Hz) at 10 Hz, and 5 pT/sqrt(Hz), from a NIST Mx instrument benchmark across 1–100 Hz. Constant ASD across the 17–25 Hz comparison bands is an assumption. The SERF case assumes a suitable compensated field and operating environment; it is not demonstrated unshielded performance. Any shielding or feedback that changes the wanted field requires its own transfer model. No gradiometric background suppression is credited. Neither scenario is extrapolated to sub-Hz in this study.', '',
 'The external ELF backgrounds remain 0.38, 1, 3.14 and 10 pT/sqrt(Hz), added in PSD with instrument and wall thermal noise. Atomic sensing removes the induction conversion penalty but does not suppress actual environmental magnetic fluctuations.', '',
 '| Distance (m) | Sensor | Background (pT/sqrt(Hz)) | BPSK boundary (S/m) | Capacity boundary (S/m) |','|---:|---|---:|---:|---:|']
 for r in rows:
  lines.append(f'| {r["distance_m"]:g} | {r["sensor"]} | {r["background_pT_sqrtHz"]:g} | {r["bpsk_boundary_S_m"]:.4g} | {r["capacity_boundary_S_m"]:.4g} |')
 lines+=['','The receiver change includes both noise and sensing geometry. Boundary shifts therefore cannot be attributed to intrinsic sensitivity alone. The same unresolved average-electrical-power accounting and finite-packet limitations as the aligned-carrier study apply.','', 'Sources:', '- https://www.nature.com/articles/s41378-026-01226-z', '- https://tf.nist.gov/ofm/smallclock/Mx_Magnetometer.html', '', 'Reproduce with scripts/magnetic_wall_atomic_receiver.py.']
 (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 print(json.dumps([r for r in rows if r['background_pT_sqrtHz']==1],indent=2))
if __name__=='__main__':main()
