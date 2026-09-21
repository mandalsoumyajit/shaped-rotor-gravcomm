"""Bounded beam graph for long-FIR seven-tone byte-mapped max-log detection.
Approximate search; converge beam width before interpreting code comparisons.
All retained edges are kept for backward max-log metrics, not just final paths.
"""
import numpy as np
from scipy.signal import lfilter

class BeamByteDetector:
    def __init__(self,model,taps,beam_width=512):
        self.model=model;self.taps=np.asarray(taps);self.width=int(beam_width)
        n=model.samples_per_symbol;self.n=n;self.length=len(taps)-1
        if self.length%n or self.length<n:raise ValueError('whole-symbol memory required')
        self.memory=self.length//n;self.modulus=7**(self.memory+1)

        self.phases,inc=model.phase_lattice()
        if self.modulus*self.phases>np.iinfo(np.int64).max:raise ValueError('history too large')
        self.inc=np.array([[inc[a,b] for b in model.alphabet] for a in model.alphabet])
        self.factors=np.exp(2j*np.pi*np.arange(self.phases)/self.phases)
        self.branches=np.array([[model.branch(a,b) for b in model.alphabet] for a in model.alphabet])
        self.current=np.array([[lfilter(self.taps,[1],b) for b in row] for row in self.branches])
        indexes=self.length+np.arange(n)[:,None]-np.arange(self.length)[None,:]
        self.past_matrix=np.where(indexes<=self.length,self.taps[np.minimum(indexes,self.length)],0).T
    def decode(self,observed,tail_symbols=6,forced_bits=None,repair_missing=True):
        blocks=np.asarray(observed).reshape(-1,self.n);payload=len(blocks)-tail_symbols
        if payload<=0 or payload%3:raise ValueError('complete byte triples required')
        history=np.array([sum(3*7**j for j in range(self.memory+1))],dtype=np.int64)
        phase=np.array([0],dtype=np.int64);cost=np.array([0.])
        past=np.ones((1,self.length),complex);graph=[];cost_offset=0.;constraints={}
        for bit,value in (forced_bits or {}).items():
            byte=bit//8;allowed=constraints.get(byte,np.ones(256,bool));allowed &= ((np.arange(256)>>(7-bit%8))&1)==value;constraints[byte]=allowed
        for t,y in enumerate(blocks):
            old=history%7;common=past@self.past_matrix
            template=common[:,None,:]+self.factors[phase,None,None]*self.current[old]
            metric=np.sum(abs(template)**2,axis=2)-2*np.real(np.sum(y.conj()*template,axis=2))
            candidate=np.broadcast_to(np.arange(7),(len(cost),7))
            values=((history//7)%7)[:,None]*49+old[:,None]*7+candidate
            if t>=payload:metric[:,np.arange(7)!=3]=np.inf
            elif t%3==2:metric[values>255]=np.inf
            if t<payload and t//3 in constraints:
                allowed_values=np.flatnonzero(constraints[t//3])
                if t%3==0:allowed=np.isin(candidate,allowed_values//49)
                elif t%3==1:allowed=np.isin(old[:,None]*7+candidate,allowed_values//7)
                else:allowed=np.isin(values,allowed_values)
                metric[~allowed]=np.inf
            next_h=(history[:,None]*7+candidate)%self.modulus
            next_p=(phase[:,None]+self.inc[old])%self.phases
            keys=(next_h*self.phases+next_p).ravel()
            totals=(cost[:,None]+metric).ravel()
            order=np.argsort(totals,kind='stable');order=order[np.isfinite(totals[order])]
            _,first=np.unique(keys[order],return_index=True)
            selected=order[np.sort(first)[:self.width]]
            if not len(selected):raise RuntimeError('empty beam')
            selected_keys=keys[selected]
            sort=np.argsort(selected_keys);sorted_keys=selected_keys[sort]
            where=np.searchsorted(sorted_keys,keys);valid=(where<len(selected_keys))
            valid &= sorted_keys[np.minimum(where,len(selected_keys)-1)]==keys
            dest=np.full(len(keys),-1,dtype=np.int32);dest[valid]=sort[where[valid]]
            source=np.repeat(np.arange(len(cost)),7)
            keep=(dest>=0)&np.isfinite(metric.ravel())
            graph.append((cost.copy(),source[keep],dest[keep],metric.ravel()[keep],values.ravel()[keep],selected//7,selected%7))
            parent=selected//7;choice=selected%7
            newbranch=self.factors[phase[parent],None]*self.branches[old[parent],choice]
            past=np.concatenate((past[parent,self.n:],newbranch),axis=1)
            history=next_h.ravel()[selected];phase=next_p.ravel()[selected]
            cost=totals[selected];cost_offset+=cost.min();cost-=cost.min()
        hard=np.empty(len(blocks),dtype=int);state=int(np.argmin(cost));back=np.zeros(len(cost));llrs=np.empty((payload//3,8));missing=0;missing_indexes=[]
        for t in range(len(graph)-1,-1,-1):
            fwd,src,dst,metric,values,parent,choice=graph[t]
            hard[t]=int(choice[state])-3;state=int(parent[state])
            scores=metric+back[dst]
            if t<payload and t%3==2:
                total=fwd[src]+scores
                for bit in range(8):
                    ones=(values>>(7-bit))&1
                    z=total[ones==0];o=total[ones==1]
                    if not len(z) or not len(o) or not np.isfinite(z.min()) or not np.isfinite(o.min()):
                        missing+=1;missing_indexes.append((t//3)*8+bit)
                    # A missing hypothesis is explicit in diagnostics. Finite
                    # saturation is computational, not a reliability assertion.
                    llrs[t//3,bit]=np.clip((o.min() if len(o) else np.inf)-(z.min() if len(z) else np.inf),-100,100)
            previous=np.full(len(fwd),np.inf);np.minimum.at(previous,src,scores)
            back=previous-np.min(previous)
        repaired=0
        if repair_missing and not forced_bits and missing_indexes:
            vals=(hard[:payload].reshape(-1,3)+3)@np.array([49,7,1])
            bits=np.unpackbits(vals.astype(np.uint8))
            for bit in missing_indexes:
                alternative=self.decode(observed,tail_symbols,{bit:1-int(bits[bit])},False)
                gap=alternative['path_cost']-cost_offset
                llrs.ravel()[bit]=np.clip(gap*(1-2*int(bits[bit])),-100,100);repaired+=1
        return dict(symbols=hard[:payload],llr=llrs.ravel(),missing_bit_hypotheses=missing,repaired_hypotheses=repaired,beam_width=self.width,memory_symbols=self.memory,path_cost=float(cost_offset))



