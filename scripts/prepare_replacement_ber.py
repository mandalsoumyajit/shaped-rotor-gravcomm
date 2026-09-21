import json
from concurrent.futures import ThreadPoolExecutor
from streaming_size_map import one
from streaming_dynamics_runtime import ROOT,save,layout
from compile_streaming_dynamics import mechanical
OUT=ROOT/'results/fixed_rate_sizing_2026_09_20/ber_replacement_v2'
old=json.loads((OUT.parent/'ber_qualification_v1/frozen_designs.json').read_text())['channels'][1]['config']
amps=[1.,1.0625,1.125]
with ThreadPoolExecutor(max_workers=3) as pool:sizes=list(pool.map(one,[(2.,a) for a in amps]))
ts=layout(old['scheme'],old['information_bits'])['symbol_s'];m=mechanical(ts,.6,.08,18.)
channels=[]
for i,(a,size) in enumerate(zip(amps,sizes)):
 c=dict(old,amplitude_scale=a,id=f'tb_2_3_f18_r0.6_q0.08_a{a:g}')
 size.update(stress_ratio=size['scale']*(18+3*.08/(.4*ts))/60,peak_torque_nm=m['peak_torque_nm']*size['inertia_kg_m2'],peak_power_w=m['peak_mechanical_power_w']*size['inertia_kg_m2'])
 assert size['stress_ratio']<=1 and size['radial_gap_m']>=.05
 channels.append(dict(config=c,seed_start=(i+3)*1000000,distance_m=[2.],source_designs=[size]))
save(OUT/'frozen_designs.json',dict(selection='Three larger 2 m candidates frozen together; choose the smallest passing tested candidate; simultaneous confidence across all three.',channels=channels))
print(json.dumps(channels,indent=2))
