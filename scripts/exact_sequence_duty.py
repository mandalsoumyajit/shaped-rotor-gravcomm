import sys,json,itertools
import numpy as np
sys.path.insert(0,'scripts')
from search_sequence_cf_fsk import OUT
from optimize_cf_fsk import I
fc=50.3
rows=[]
for ts,r,q in ((2.,.4,.1),(1.75,.4,.075),(1.75,.4,.1)):
 t=np.linspace(0,ts,20001);u=np.minimum(t/(r*ts),1)
 smooth=10*u**3-15*u**4+6*u**5
 slope=(30*u**2-60*u**3+30*u**4)/(r*ts)
 spacing=q/((1-r)*ts);cost=np.zeros((7,7))
 for a,b in itertools.product(range(7),repeat=2):
  f=fc+spacing*((a-3)+(b-a)*smooth)
  tau=I*np.pi*spacing*(b-a)*slope+.05*f/fc
  cost[a,b]=np.trapz(tau*tau,t)/ts
 words=np.array([[b//49,(b//7)%7,b%7] for b in range(256)])
 edges=np.full((7,7),-np.inf);choice=np.zeros((7,7),int)
 for old in range(7):
  for byte,(a,b,c) in enumerate(words):
   energy=cost[old,a]+cost[a,b]+cost[b,c]
   if energy>edges[old,c]:edges[old,c]=energy;choice[old,c]=byte
 best=(-np.inf,())
 for length in range(1,8):
  for nodes in itertools.combinations(range(7),length):
   for rest in itertools.permutations(nodes[1:]):
    cycle=(nodes[0],)+rest
    avg=sum(edges[a,b] for a,b in zip(cycle,cycle[1:]+cycle[:1]))/(3*length)
    if avg>best[0]:best=(avg,cycle)
 cycle=best[1];payload=[int(choice[a,b]) for a,b in zip(cycle,cycle[1:]+cycle[:1])]
 first=np.bincount(words[:,0],minlength=7)/256
 last=np.bincount(words[:,2],minlength=7)/256
 random=(last@cost@first+np.mean(cost[words[:,0],words[:,1]]+cost[words[:,1],words[:,2]]))/3
 rows.append(dict(ts=ts,r=r,q=q,worst_valid_byte_stream_rms_nm=float(np.sqrt(best[0])),
  uniform_random_byte_stream_rms_nm=float(np.sqrt(random)),worst_repeating_bytes=payload))
(OUT/'exact_message_duty.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
