"""Causal full-history beam search with the matched finite-record noise metric.

No finite-memory state merging or truncated whitening. Beam pruning is the only
search approximation. Timing, carrier/phase and deterministic prehistory are
acquired inputs; this module does not perform acquisition.
"""
import numpy as np
from scipy.linalg import solve_triangular
from scipy.signal import sosfilt, sosfilt_zi

class CausalBeamByteDetector:
    def __init__(self, receiver, symbol_s, amplitude_m_s2, beam_width=256,
                 transition_fraction=.4, q=.1, cached_branches=True):
        if symbol_s<=0 or not 0<transition_fraction<1 or q<=0 or beam_width<1:
            raise ValueError('Invalid waveform or beam parameters')
        self.receiver=receiver;self.symbol_s=float(symbol_s)
        self.amplitude=amplitude_m_s2;self.width=int(beam_width)
        self.r=transition_fraction;self.q=q
        self._whitening_cache=None
        self.cached_branches=cached_branches;self._branch_cache={}

    def decode(self, observed, nbytes, tail_symbols=6, *, initial_previous=0,
               initial_phase_rad=0., initial_filter_state=None,
               sample_offset=0, forced_bits=None, bit_hypotheses=0,
               _save_checkpoints=False, _resume=None):
        """Decode normalized filtered observations; LLR positive means bit zero.

        The initial filter state describes the known *mean*, not random noise.
        If omitted, constant prehistory exp(i*phase) at offset tone zero is
        required. Supply the actual SOS mean state and sample offset when
        continuing an existing stream. Noise covariance is the marginal
        stationary covariance of this record (not conditioned on older noise).
        Missing bit alternatives produce NaN LLRs, never infinite confidence.
        Returned LLRs are max-log over retained complete paths, approximate
        whenever pruning occurred. Forced-bit reruns can assess missing paths.
        """
        if bit_hypotheses<0 or int(bit_hypotheses)!=bit_hypotheses:
            raise ValueError('Nonnegative integer bit quota required')
        if nbytes<1 or int(nbytes)!=nbytes or tail_symbols<0 or int(tail_symbols)!=tail_symbols:
            raise ValueError('Positive byte count and nonnegative integer tail required')
        if initial_previous not in range(-3,4) or sample_offset<0 or int(sample_offset)!=sample_offset:
            raise ValueError('Invalid initial tone or sample offset')
        if initial_filter_state is None and initial_previous!=0:
            raise ValueError('Nonzero initial tone requires explicit filter prehistory')
        rec=self.receiver;fs=rec.input_rate_hz;D=rec.decimation
        total=3*nbytes+tail_symbols
        boundaries=np.ceil(np.arange(total+1)*self.symbol_s*fs-1e-12).astype(int)
        count=len(np.arange((-sample_offset)%D,boundaries[-1],D))
        y=np.asarray(observed,complex)
        if y.shape!=(count,) or np.any(~np.isfinite(y)):raise ValueError('Wrong observation record')
        if self._whitening_cache is None or self._whitening_cache[0]!=count:
            model=rec.finite_record(count,self.amplitude)
            self._whitening_cache=(count,solve_triangular(model.factor,np.eye(count),lower=True))
        W=self._whitening_cache[1]
        wy=W@y
        sos=rec.sos
        zi=sosfilt_zi(sos).astype(complex)*np.exp(1j*initial_phase_rad) if initial_filter_state is None else np.asarray(initial_filter_state,complex)
        if zi.shape!=(len(sos),2):raise ValueError('Wrong SOS mean-state shape')
        state=zi[:,None,:].copy();phase=np.array([initial_phase_rad/(2*np.pi)])
        old=np.array([initial_previous]);paths=np.empty((1,0),np.int8)
        history=np.empty((1,0),complex);cost=np.zeros(1);out_start=0;pruned=0;peak=0;peak_survivors=1
        allowed=np.ones((nbytes,256),bool)
        for bit,value in (forced_bits or {}).items():
            if not 0<=bit<8*nbytes or value not in (0,1):raise ValueError('Invalid forced bit')
            allowed[bit//8]&=((np.arange(256)>>(7-bit%8))&1)==value
        checkpoints={};begin=0
        if _resume is not None:
            begin,state,phase,old,paths,history,cost,out_start,pruned,peak,peak_survivors=_resume
        for t in range(begin,total):
            if _save_checkpoints and t<3*nbytes and t%3==0:
                # Arrays below are replaced, never modified in place. All bit
                # constraints in this byte share precisely this prefix beam.
                checkpoints[t//3]=(t,state,phase,old,paths,history,cost,out_start,pruned,peak,peak_survivors)
            choices=np.array([0]) if t>=3*nbytes else np.arange(-3,4)
            parent=np.repeat(np.arange(len(cost)),len(choices));current=np.tile(choices,len(cost))
            proposed=np.column_stack((paths[parent],current))
            if t<3*nbytes:
                j=t%3;prefix=(proposed[:,-j-1:]+3)@ (7**np.arange(j,-1,-1))
                possible=np.flatnonzero(allowed[t//3])//(7**(2-j))
                valid=np.isin(prefix,possible)
                parent=parent[valid];current=current[valid];proposed=proposed[valid]
            if not len(parent):raise RuntimeError('No valid candidates')
            start,end=boundaries[t:t+2]
            local=np.arange(start,end)/fs-t*self.symbol_s
            tr=self.r*self.symbol_s;u=np.minimum(local/tr,1.)
            integral=tr*(2.5*u**4-3*u**5+u**6)+np.maximum(local-tr,0.)
            spacing=self.q/((1-self.r)*self.symbol_s)
            first=(-(sample_offset+start))%D
            if self.cached_branches:
                key=(t,start,end,first)
                if key not in self._branch_cache:
                    before=np.repeat(np.arange(-3,4),7);after=np.tile(np.arange(-3,4),7)
                    raw=np.exp(2j*np.pi*spacing*(before[:,None]*local+(after-before)[:,None]*integral))
                    branch,branch_state=sosfilt(sos,raw,axis=1,zi=np.zeros((len(sos),49,2),complex))
                    dim=2*len(sos)
                    initial=np.eye(dim,dtype=complex).reshape(dim,len(sos),2).transpose(1,0,2)
                    response,transition=sosfilt(sos,np.zeros((dim,len(local)),complex),axis=1,zi=initial)
                    self._branch_cache[key]=(branch[:,first::D],branch_state.transpose(1,0,2).reshape(49,dim),response[:,first::D],transition.transpose(1,0,2).reshape(dim,dim))
                branch,branch_state,response,transition=self._branch_cache[key]
                flat=state.transpose(1,0,2).reshape(len(cost),-1)
                factor=np.exp(2j*np.pi*phase[parent]);index=(old[parent]+3)*7+current+3
                means=(flat@response)[parent]+factor[:,None]*branch[index]
                flat_next=(flat@transition)[parent]+factor[:,None]*branch_state[index]
                newstate=flat_next.reshape(-1,len(sos),2).transpose(1,0,2)
            else:
                cycles=phase[parent,None]+spacing*(old[parent,None]*local+(current-old[parent])[:,None]*integral)
                raw=np.exp(2j*np.pi*cycles)
                filtered,newstate=sosfilt(sos,raw,axis=1,zi=state[:,parent,:])
                means=filtered[:,first::D]
            endout=out_start+means.shape[1]
            # Apply full dense causal innovations to each candidate's mean history.
            predicted=(history@W[out_start:endout,:out_start].T)[parent] + means@W[out_start:endout,out_start:endout].T
            scores=cost[parent]+np.sum(abs(wy[None,out_start:endout]-predicted)**2,axis=1)
            peak=max(peak,len(scores));order=np.argsort(scores,kind='stable')
            if bit_hypotheses:
                # Keep a common global beam plus the best k paths for every
                # completed bit/value. This preserves alternatives across a
                # packet without merging distinct filter/innovation states.
                keep=np.zeros(len(scores),bool);keep[order[:self.width]]=True
                completed=min((t+1)//3,nbytes)
                if completed:
                    byte_values=((proposed[:,:3*completed].reshape(-1,completed,3)+3)@np.array([49,7,1])).astype(np.uint8)
                    labels=np.unpackbits(byte_values[order],axis=1)
                    for bit in range(labels.shape[1]):
                        for value in (0,1):
                            keep[order[np.flatnonzero(labels[:,bit]==value)[:bit_hypotheses]]]=True
                # Preserve possible prefixes of the next byte, before its bits
                # have labels. This avoids losing a new bit value prematurely.
                if t<3*nbytes and t%3!=2:
                    prefix=(proposed[:,-(t%3+1):]+3)@(7**np.arange(t%3,-1,-1))
                    sorted_prefix=prefix[order]
                    for value in np.unique(prefix):
                        keep[order[np.flatnonzero(sorted_prefix==value)[:bit_hypotheses]]]=True
                selection=order[keep[order]]
            else:selection=order[:self.width]
            peak_survivors=max(peak_survivors,len(selection))
            pruned+=max(0,len(scores)-len(selection))
            history=np.concatenate((history[parent[selection]],means[selection]),axis=1)
            state=newstate[:,selection,:];paths=proposed[selection];cost=scores[selection]
            phase=(phase[parent[selection]]+spacing*self.symbol_s*(self.r*old[parent[selection]]/2+(1-self.r/2)*current[selection]))%1
            old=current[selection];out_start=endout
        values=((paths[:,:3*nbytes].reshape(-1,nbytes,3)+3)@np.array([49,7,1])).astype(np.uint8)
        bits=np.unpackbits(values,axis=1);llr=np.full(8*nbytes,np.nan)
        for b in range(8*nbytes):
            zero=cost[bits[:,b]==0];one=cost[bits[:,b]==1]
            if len(zero) and len(one):llr[b]=one.min()-zero.min()
        best=int(np.argmin(cost))
        return dict(bytes=values[best],symbols=paths[best,:3*nbytes],llr=llr,
                    _checkpoints=checkpoints,
                    path_cost=float(cost[best]),exact_search=pruned==0,
                    pruned_candidates=pruned,peak_candidates=peak,peak_survivors=peak_survivors,
                    bit_hypotheses=bit_hypotheses,global_beam_width=self.width,
                    missing_bit_hypotheses=int(np.isnan(llr).sum()),
                    survivor_bytes=values,survivor_costs=cost,
                    final_mean_filter_state=state[:,best,:],
                    final_phase_rad=float(2*np.pi*phase[best]),final_previous=int(old[best]),
                    next_sample_offset=int(sample_offset+boundaries[-1]))

    def decode_soft(self, observed, nbytes, tail_symbols=6, **initial):
        """Expensive reference soft output using a forced-opposite search per bit.

        Each search keeps the exact causal metric but may prune. Pooling all
        retained candidates makes finite max-log outputs internally consistent.
        No claim of exact LLRs is made unless every search was exhaustive.
        Use to check cheaper soft-output approximations before code comparisons.
        """
        if 'forced_bits' in initial:raise ValueError('Soft reference sets its own bit constraints')
        base=self.decode(observed,nbytes,tail_symbols,**initial)
        bits=np.unpackbits(base['bytes']);all_values=[base['survivor_bytes']]
        all_costs=[base['survivor_costs']];exact=base['exact_search']
        for bit,value in enumerate(bits):
            alternative=self.decode(observed,nbytes,tail_symbols,
                                    forced_bits={bit:1-int(value)},**initial)
            all_values.append(alternative['survivor_bytes']);all_costs.append(alternative['survivor_costs'])
            exact &= alternative['exact_search']
        values=np.concatenate(all_values);costs=np.concatenate(all_costs)
        labels=np.unpackbits(values,axis=1);llr=np.empty(8*nbytes)
        for bit in range(8*nbytes):
            llr[bit]=costs[labels[:,bit]==1].min()-costs[labels[:,bit]==0].min()
        best=int(np.argmin(costs));v=values[best].astype(int)
        return dict(bytes=values[best],symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3,
                    llr=llr,path_cost=float(costs[best]),exact_search=bool(exact),
                    search_runs=1+8*nbytes,missing_bit_hypotheses=0,
                    method='pooled_forced_opposite_full_history_beam')

    def decode_diverse(self, observed, nbytes, tail_symbols=6, *, bit_hypotheses=2, **initial):
        """One forward search with explicit per-bit survivor reservations.

        Retain the best k candidates per completed bit/value plus the global
        beam, and protect unfinished byte prefixes. Quotas prevent vanished
        alternatives, but do not certify that the surviving minimum is global.
        The full causal likelihood is unchanged; no filter state is merged.
        Compare quota/width convergence and forced searches before using LLRs.
        """
        if bit_hypotheses<1:raise ValueError('Positive bit quota required')
        result=self.decode(observed,nbytes,tail_symbols,bit_hypotheses=bit_hypotheses,**initial)
        result['method']='bit_diverse_full_history_beam'
        return result

    def decode_checkpoint_soft(self, observed, nbytes, tail_symbols=6, *, workers=1, parallel_backend="thread", **initial):
        """Same constrained-search reference, sharing prefix beams and covariance.

        For a bit in byte j the constraint first acts at byte j, so unrestricted
        prefix snapshots are reusable without adding a new approximation. Each
        continuation is still beam-pruned. Independent continuations can run on
        a bounded thread pool with numerical-library threads set to one.
        """
        from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
        if 'forced_bits' in initial or 'bit_hypotheses' in initial:
            raise ValueError('Checkpoint reference controls its own constraints')
        if parallel_backend not in ('thread','process'):raise ValueError('Unknown parallel backend')
        if workers<1 or int(workers)!=workers:raise ValueError('Positive integer worker count required')
        base=self.decode(observed,nbytes,tail_symbols,_save_checkpoints=True,**initial)
        bits=np.unpackbits(base['bytes'])
        def continuation(bit):
            return self.decode(observed,nbytes,tail_symbols,forced_bits={bit:1-int(bits[bit])},
                               _resume=base['_checkpoints'][bit//8],**initial)
        # Incremental pooling has bounded output storage, including across threads.
        minima=np.full((len(bits),2),np.inf);best_cost=np.inf;best_bytes=None;exact=True
        def consume(got):
            nonlocal best_cost,best_bytes,exact
            labels=np.unpackbits(got['survivor_bytes'],axis=1);costs=got['survivor_costs']
            for bit in range(len(bits)):
                for value in (0,1):
                    group=costs[labels[:,bit]==value]
                    if len(group):minima[bit,value]=min(minima[bit,value],float(group.min()))
            if got['path_cost']<best_cost:best_cost=got['path_cost'];best_bytes=got['bytes'].copy()
            exact &= got['exact_search']
        consume(base)
        if workers==1:
            for bit in range(len(bits)):consume(continuation(bit))
        elif parallel_backend=='process':
            # A worker receives immutable prefix snapshots once at startup.
            # Spawn callers must use the usual __main__ guard on Windows.
            context=(self,observed,nbytes,tail_symbols,initial,base['_checkpoints'],bits)
            with ProcessPoolExecutor(max_workers=workers,initializer=_checkpoint_worker_init,initargs=(context,)) as pool:
                for start in range(0,len(bits),workers):
                    for got in pool.map(_checkpoint_worker,range(start,min(start+workers,len(bits)))):consume(got)
        else:
            # Submit batches to avoid retaining every completed survivor list.
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for start in range(0,len(bits),workers):
                    for got in pool.map(continuation,range(start,min(start+workers,len(bits)))):consume(got)
        v=best_bytes.astype(int)
        return dict(bytes=best_bytes,symbols=np.column_stack((v//49,v//7%7,v%7)).ravel()-3,
                    llr=minima[:,1]-minima[:,0],path_cost=best_cost,exact_search=bool(exact),
                    missing_bit_hypotheses=0,search_runs=1+len(bits),workers=workers,parallel_backend=parallel_backend,
                    method='checkpoint_shared_forced_opposite_full_history_beam')


def _checkpoint_worker_init(context):
    global _checkpoint_context
    _checkpoint_context=context

def _checkpoint_worker(bit):
    detector,observed,nbytes,tail,initial,snapshots,bits=_checkpoint_context
    return detector.decode(observed,nbytes,tail,forced_bits={bit:1-int(bits[bit])},
                           _resume=snapshots[bit//8],**initial)
