"""Explore below the original search grid; retain computational exclusions."""
import json
from search_sequence_cf_fsk import evaluate,trajectory,mechanics,OUT
rows=[]
for ts in (1.25,1.5,1.75):
    for q in (.015625,.03125,.046875,.0625,.09375):
        r=.5;t,f,df,_,_=trajectory(ts,r,q);mech=mechanics(ts,r,q,t,f,df)
        row=dict(ts=ts,transition_fraction=r,spacing_dwell_product=q,mechanics=mech)
        if not mech['feasible']:row['excluded']='actuator limits'
        else:
            try:row.update(evaluate(ts,r,q,24,981171,6.,2.))
            except ValueError as exc:
                if 'state' not in str(exc):raise
                row['excluded']=str(exc)
        rows.append(row)
        print(json.dumps(row),flush=True)
        (OUT/'fast_refinement.json').write_text(json.dumps(rows,indent=2)+'\n')
