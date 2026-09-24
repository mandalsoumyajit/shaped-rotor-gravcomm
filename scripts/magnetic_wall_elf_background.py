"""ELF-paper environmental levels applied to the finite-hardware wall model."""
from pathlib import Path
import json
import numpy as np
from magnetic_wall_study import Channel, DESIGNS, ROOT
from gravcomm.magnetic_link import noise_asd

OUT=ROOT/'results/magnetic_wall_elf_background_2026_09_21'
LEVELS=[(.01,'quiet sweep endpoint'),(.03,'quiet comparison'),(.1,'sweep'),(.38,'mine-surface lower example'),(1.,'representative'),(3.,'ELF sweep upper endpoint'),(3.14,'mine-surface upper example'),(10.,'tunnel / stress')]

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 rows=[]
 for design in DESIGNS:
  c=Channel(design,design['distance'])
  for pt,label in LEVELS:
   ambient=pt*1e-12
   b=c.crossing('bpsk',ambient=ambient)
   cap=c.crossing('capacity',ambient=ambient)
   lp,fc=c.bpsk(5.8e7,ambient=ambient)
   rows.append(dict(distance_m=design['distance'],wall_thickness_m=c.thickness,power_ceiling_W=design['power'],ambient_pT_sqrtHz=pt,label=label,bpsk_crossing_S_m=b,capacity_crossing_S_m=cap,at_sigma_5_8e7=dict(bpsk_average_W=float(np.exp(lp)),bpsk_carrier_Hz=fc,capacity_minimum_average_W=float(np.exp(c.capacity(5.8e7,ambient=ambient))))))
  print(f"Completed {design['distance']} m",flush=True)
 (OUT/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
 report=['# ELF-background magnetic-wall comparison','',
 'The representative background is 1 pT/sqrt(Hz), bracketed by the ELF paper\'s mine-surface examples of 0.38 and 3.14 pT/sqrt(Hz). The 10 pT/sqrt(Hz) tunnel example is retained as a stress case. The full 0.01–3 pT/sqrt(Hz) ELF sweep is also sampled. These are constant broadband ASD assumptions; the cited 1 kHz observations are not measurements across the entire optimized link band. No gradiometric cancellation is credited to either link. Ambient, instrumental, and wall-thermal PSDs add.', '',
 'Geometry, finite-coil response, noise extrapolation, and power accounting follow ../magnetic_wall_finite_hardware_2026_09_21/RESULTS.md. In particular, the numerical magnetic average-power ceilings equal gravitational peak mechanical powers; equal average electrical consumption is not established. The Gaussian capacity boundary is asymptotic and does not demonstrate a finite-packet implementation.', '',
 '| Distance (m) | Background (pT/sqrt(Hz)) | BPSK boundary (S/m) | Ideal capacity boundary (S/m) |','|---:|---:|---:|---:|']
 for r in rows:
  if r['ambient_pT_sqrtHz'] in (.38,1.,3.14,10.):report.append(f"| {r['distance_m']:g} | {r['ambient_pT_sqrtHz']:g} | {r['bpsk_crossing_S_m']:.4g} | {r['capacity_crossing_S_m']:.4g} |")
 report+=['','The boundaries change little across this environmental bracket because the assumed inductive receiver becomes noisy at low frequencies. This is a property of the specified receiver and broadband assumptions. It does not establish environmental insignificance for a different magnetometer or for an actual low-frequency site spectrum.','', '| Frequency (Hz) | Instrument ASD (pT/sqrt(Hz)) |','|---:|---:|']
 for f in (.1,1.,10.,100.,1000.):report.append(f'| {f:g} | {float(noise_asd(np.array(f),ambient=0))*1e12:.4g} |')
 report+=['','At conductivity 5.8e7 S/m the 2 m geometry still exceeds the BPSK benchmark power budget while passing the optimistic capacity bound. This does not establish that magnetic communication is impossible. Only the three reference gravitational links are qualified.','', 'Source: ELF_Paper/journal/ieee-access/IEEE-Access-ELF.tex, background discussion and PCB receiver example. Reproduce with scripts/magnetic_wall_elf_background.py. The earlier quiet-case map remains labelled as the earlier scenario; this directory contains the revised three-design comparison.']
 (OUT/'RESULTS.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
 print('COMPLETE',flush=True)
if __name__=='__main__':main()
