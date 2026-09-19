import sys,json
sys.path.insert(0,'scripts')
from search_sequence_cf_fsk import evaluate,OUT,trajectory,mechanics
rows=[]
for ts,r,q in ((2.,.4,.1),(2.25,.5,.09375),(2.5,.5,.125)):
 t,f,df,_,_=trajectory(ts,r,q);m=mechanics(ts,r,q,t,f,df)
 print(ts,r,q,m,flush=True)
 if m['feasible']:
  result=evaluate(ts,r,q,200,781731,6.,2.);result['mechanics']=m;rows.append(result)
  print(json.dumps(result),flush=True)
  (OUT/'boundary_refinement.json').write_text(json.dumps(rows,indent=2)+'\n')
