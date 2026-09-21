"""Map tested received amplitudes to uniform rotor size at all three distances."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import sys,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from scipy.optimize import brentq
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'src')]
from range_revision import MASS,INERTIA,FC,duty,save
from gravcomm import FiniteRotor
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20'

def main():
 fields=json.loads((OUT/'reconciliation_v1/fields.json').read_text());rotor=FiniteRotor.from_npz(ROOT/'fea/results/shaped_waveforms/medium_rounded.npz')
 def one(spec):
  reference,distance=spec;target=reference['amplitude_m_s2'];cache={}
  def amplitude(s):
   key=float(s)
   if key not in cache:cache[key]=s*rotor.harmonic_amplitude(distance/s,samples=64)
   return cache[key]
  s=reference['scale'] if distance==.5 else brentq(lambda x:amplitude(x)-target,.5,min(4.,(distance-.0627)/.25),xtol=1e-8)
  refined=s*rotor.harmonic_amplitude(distance/s,samples=128)
  assert abs(refined/target-1)<1e-7
  protocols=[]
  for name,ts in [('uncoded_256',256/126),('tb_3_4_256',256/162),('tb_2_3_256',256/180),('half_rate_224',224/210)]:
   protocols.append(dict(protocol=name,symbol_s=ts,mechanics=duty(ts,.4,.1,inertia=INERTIA*s**5,fc=FC/s,drag=.05*s**3,max_rpm=1800/s)))
  row=dict(reference_diameter_at_0p5m=reference['diameter_m'],distance_m=distance,scale=s,diameter_m=.5*s,mass_kg=MASS*s**3,inertia_kg_m2=INERTIA*s**5,carrier_hz=FC/s,amplitude_m_s2=refined,amplitude_match_relative=abs(refined/target-1),protocols=protocols,scope='Equal received second-harmonic amplitude under fixed-noise tunable-center model. Not minimum required size or validated BER. Geometry and actuator feasibility remain separate.')
  print(json.dumps({k:row[k] for k in ('reference_diameter_at_0p5m','distance_m','diameter_m','mass_kg')}),flush=True)
  return row
 specs=[(r,d) for r in fields if r['scale'] in (.95,1.) for d in (.5,1.,2.)]
 with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(one,specs))
 save(OUT/'equal_amplitude_size_map.json',rows)
 lines=['# Equal-amplitude transmitter size map','','This maps the two tested 0.5 m-link signal amplitudes to uniform transmitter sizes at longer distances. Under the fixed PSD-versus-offset and tunable-center model, identical prescribed baseband waveforms then have the same channel statistics. It does not establish minimum reliable size, actuator feasibility, or hardware validity.','','| Reference diameter at 0.5 m (m) | Distance (m) | Equal-amplitude diameter (m) | Mass (kg) | Carrier (Hz) | Rate-3/4 RMS torque (Nm) |','|---:|---:|---:|---:|---:|---:|']
 for r in rows:lines.append(f"| {r['reference_diameter_at_0p5m']:.3f} | {r['distance_m']:g} | {r['diameter_m']:.4f} | {r['mass_kg']:.2f} | {r['carrier_hz']:.3f} | {r['protocols'][1]['mechanics']['worst_valid_message_rms_nm']:.2f} |")
 lines+=['','Torque examples use q=0.1, rate 3/4, 256 useful bits, 1 useful bit/s, and assumed drag 0.05 s^3 Nm. Enlarged rotors exceed the reference motor ratings; these are demands, not feasible drive selections. The centrifugal stress similarity at carrier 50.3/s Hz does not validate supports, modulation loads, modal behavior or containment.']
 (OUT/'EQUAL_AMPLITUDE_MAP.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
if __name__=='__main__':main()
