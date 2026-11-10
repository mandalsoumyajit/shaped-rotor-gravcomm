"""Frozen inputs for the hop-limited CPFSK rotor scales (2026-11-05).

The structural pipeline only needs, per link: carrier, commanded gravitational
frequency range, peak angular acceleration (through peak_torque_nm/inertia_kg_m2),
and the geometric scale.  CAD mass and polar inertia scale exactly as s^3 and s^5
from the archived 2026-09-20 frozen records of the same shaped rotor.
"""
import json,math,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OLD=ROOT/'results/three_link_designs_2026_09_20/frozen_designs.json'
TAG=os.environ.get('GRAVCOMM_FEA_TAG','2026_11_05_nb')
OUTS=[ROOT/f'results/three_link_designs_{TAG}',ROOT/f'results/tungsten_comparison_{TAG}']
SYMBOL_S=224/190
HALF_RANGE_HZ=0.545306          # commanded f_g = fc +- HALF_RANGE_HZ
TONE_SPACING_HZ=0.181760         # tone offsets k*spacing, k=-3..3
PEAK_DFG_DT_HZ_S=0.867212        # largest hop (3 tones) over a full-symbol quintic
ALPHA=math.pi*PEAK_DFG_DT_HZ_S   # rad/s^2, alpha = pi * df_g/dt
# Final design points supplied by the user on 2026-11-05 (supersede the earlier _hop scales 0.93416/1.63467/2.96032).
LINKS=[(0.5,24.0,0.9010437566657981),(1.0,24.0,1.5761135310793883),(2.0,18.0,2.854447700523453)]
def main():
 old=json.loads(OLD.read_text());ref=old['channels'][0]['source_designs'][0]
 M0=ref['mass_kg']/ref['scale']**3;I0=ref['inertia_kg_m2']/ref['scale']**5
 for ch in old['channels']:
  x=ch['source_designs'][0];assert abs(x['mass_kg']/x['scale']**3/M0-1)<1e-9 and abs(x['inertia_kg_m2']/x['scale']**5/I0-1)<1e-9
 channels=[]
 for d,fc,s in LINKS:
  I=I0*s**5
  cfg=dict(id=f'hop_cpfsk_f{fc:g}_d{d:g}',scheme='hop-limited continuous-phase FSK, iterative receiver',carrier_hz=fc,
   tone_spacing_hz=TONE_SPACING_HZ,tone_offsets_k=list(range(-3,4)),commanded_half_range_hz=HALF_RANGE_HZ,
   symbol_s=SYMBOL_S,transition_s=SYMBOL_S,transition_profile='quintic over the full symbol, no dwell',max_hop_tones=3,
   peak_dfg_dt_hz_s=PEAK_DFG_DT_HZ_S)
  src=dict(distance_m=d,scale=s,diameter_m=0.5*s,mass_kg=M0*s**3,inertia_kg_m2=I,peak_alpha_rad_s2=ALPHA,peak_torque_nm=I*ALPHA,
   peak_power_w=I*ALPHA*math.pi*(fc+HALF_RANGE_HZ),
   provenance='scale supplied by user 2026-11-05; mass/inertia = archived CAD values scaled by s^3/s^5; peak_torque = inertia*pi*peak_dfg_dt (rigid body, no drag)')
  channels.append(dict(config=cfg,distance_m=[d],source_designs=[src]))
 frozen=dict(selection=f'Hop-limited CPFSK final design points (user, 2026-11-05, tag {TAG}) replacing the 2026-09-20 tb_2/3 designs',
  cad_reference=dict(mass_kg=M0,inertia_kg_m2=I0,source=OLD.relative_to(ROOT.parent).as_posix()),channels=channels)
 for out in OUTS:
  out.mkdir(parents=True,exist_ok=True);(out/'frozen_designs.json').write_text(json.dumps(frozen,indent=2)+'\n')
 print(json.dumps(frozen,indent=1))
if __name__=='__main__':main()
