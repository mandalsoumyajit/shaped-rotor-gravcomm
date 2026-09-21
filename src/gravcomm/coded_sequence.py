"""Experimental max-log soft sequence detector and explicit K=7 outer code.

LLR convention: positive favors zero. Channel outputs are max-log bit metrics,
not exact posterior LLRs; outer decoding treats them as additive bit metrics.
No independence of actual channel errors or iterative detection is assumed.
"""
import numpy as np
from numba import njit


def conv_encode(bits):
    """Rate 1/2, octal (171,133); new bit at LSB, six zero tail bits."""
    bits=np.asarray(bits,dtype=np.uint8)
    if bits.ndim!=1 or np.any(bits>1): raise ValueError('binary vector required')
    state=0; output=[]
    for bit in np.r_[bits,np.zeros(6,dtype=np.uint8)]:
        reg=(state<<1)|int(bit)
        output.extend([(reg & 0o171).bit_count()%2,(reg & 0o133).bit_count()%2])
        state=reg & 63
    return np.array(output,dtype=np.uint8)


@njit(cache=True,nogil=True)
def _conv_decode(llr,outputs):
    length=len(llr)//2
    costs=np.full(64,np.inf); costs[0]=0.
    parents=np.empty((length,64),np.int64)
    for t in range(length):
        new=np.full(64,np.inf)
        for s in range(64):
            for b in range(2):
                dest=((s<<1)|b)&63
                val=costs[s]+outputs[s,b,0]*llr[2*t]+outputs[s,b,1]*llr[2*t+1]
                if val<new[dest]:new[dest]=val;parents[t,dest]=s
        costs=new-np.min(new)
    state=0;bits=np.empty(length,np.uint8)
    for t in range(length-1,-1,-1):
        bits[t]=state&1;state=parents[t,state]
    return bits[:-6]


def conv_decode(llr):
    llr=np.asarray(llr,dtype=float)
    if len(llr)%2 or len(llr)<12 or not np.all(np.isfinite(llr)):
        raise ValueError('finite paired metrics including termination required')
    outputs=np.array([[[( ((s<<1)|b)&g).bit_count()%2 for g in (0o171,0o133)] for b in range(2)] for s in range(64)])
    return _conv_decode(llr,outputs)


@njit(cache=True,nogil=True)
def _maxlog(correlation,energy,factors,destinations,edges,phase_indices,byte_values,initial,payload,zero):
    length=correlation.shape[0];states=len(edges)
    forward=np.full((length+1,states),np.inf);forward[0,initial]=0.
    for t in range(length):
        byte_end=t<payload and t%3==2
        for s in range(states):
            if not np.isfinite(forward[t,s]):continue
            for b in range(7):
                if t>=payload and b!=zero:continue
                if byte_end and byte_values[s,b]>255:continue
                e=edges[s,b];d=destinations[s,b]
                metric=energy[e]-2*(factors[phase_indices[s]]*correlation[t,e]).real
                value=forward[t,s]+metric
                if value<forward[t+1,d]:forward[t+1,d]=value
        forward[t+1]-=np.min(forward[t+1])
    backward=np.zeros(states)
    llrs=np.empty((payload//3,8))
    for t in range(length-1,-1,-1):
        new=np.full(states,np.inf); minima=np.full((8,2),np.inf)
        byte_end=t<payload and t%3==2
        for s in range(states):
            for b in range(7):
                if t>=payload and b!=zero:continue
                valuebyte=byte_values[s,b]
                if byte_end and valuebyte>255:continue
                e=edges[s,b];d=destinations[s,b]
                metric=energy[e]-2*(factors[phase_indices[s]]*correlation[t,e]).real
                value=metric+backward[d]
                if value<new[s]:new[s]=value
                if byte_end:
                    total=forward[t,s]+value
                    for k in range(8):
                        bit=(valuebyte>>(7-k))&1
                        if total<minima[k,bit]:minima[k,bit]=total
        if byte_end:
            llrs[t//3]=minima[:,1]-minima[:,0]
        backward=new-np.min(new)
    return llrs.ravel()


class MaxLogByteDetector:
    """Forward/backward min-sum on the existing finite-whitener trellis.

    Exact max-log minima for this trellis, uniform independent byte prior,
    known zero initial state and zero-tone trailer, free final phase.
    Byte padding is not supplied as side information to the detector.
    """
    def __init__(self,decoder):
        self.decoder=decoder
        if decoder.model.alphabet!=(-3,-2,-1,0,1,2,3) or decoder.history_length<2:
            raise ValueError('seven-tone byte history required')
        states=np.arange(decoder.state_count); h=states%decoder.histories;old=h%7
        self.phase_indices=((states//decoder.histories)*decoder.stride-decoder.alpha*(old-3))%decoder.phases
        current=np.arange(7)[None,:]
        self.edges=h[:,None]*7+current
        next_phase=(self.phase_indices[:,None]+decoder.increments[old[:,None],current])%decoder.phases
        reduced=((next_phase+decoder.alpha*(current-3))%decoder.phases)//decoder.stride
        self.destinations=reduced*decoder.histories+(h[:,None]*7+current)%decoder.histories
        self.byte_values=((h//7)%7)[:,None]*49+old[:,None]*7+current

    def decode(self,packet,tail_symbols=6):
        d=self.decoder;blocks=np.asarray(packet).reshape(-1,d.model.samples_per_symbol)
        payload=len(blocks)-tail_symbols
        if payload<=0 or payload%3:raise ValueError('complete byte triples required')
        corr=blocks.conj()@d.templates.T
        return _maxlog(corr,d.energy,d.phase_factors,self.destinations,self.edges,self.phase_indices,
                       self.byte_values,d.initial_history,payload,3)

def ldpc_matrices(k=64):
    """CCSDS 231.1-O-1 (April 2015), Tables 2-1 to 2-6; information lengths 64,128,256.

    Archived primary source: https://ccsds.org/Pubs/231x1o1s.pdf.
    This implements the code only, not a complete CCSDS link protocol.
    """
    specs={
        64: ([[(0,7),2,14,6,None,0,13,0], [6,(0,15),0,1,0,None,0,7],
              [4,1,(0,15),14,11,0,None,3], [0,1,9,(0,13),14,1,0,None]],
             ('0E69166BEF4C0BC2','7766137EBB248418','C480FEB9CD53A713','4EAA22FA465EEA11')),
        128: ([[(0,31),15,25,0,None,20,12,0], [28,(0,30),29,24,0,None,1,20],
               [8,0,(0,28),1,29,0,None,21], [18,30,0,(0,30),25,26,0,None]],
              ('73F5E8390220CE5136ED68E9F39EB162','BAC812C0BCD243794786D9285A09095C',
               '7DF83F76A5FF4C388E6C0D4E025EB712','BAA37B3260CB31C5D0F66A31FAF511BC')),
        256: ([[(0,63),30,50,25,None,43,62,0], [56,(0,61),50,23,0,None,37,26],
               [16,0,(0,55),27,56,0,None,43], [35,56,62,(0,11),58,3,0,None]],
              ('1D21794A22761FAE59945014257E130D74D60540037940142DADEB9CA25EF12E',
               '60E0B6623C5CE5124D2C81ECC7F469AB20678DBFB7523ECE2B54B906A9DBE98C',
               'F6739BCF54273E77167BDA120C6C47744C071EFF5E32A7593138670C095C39B5',
               '28706BD0453002582DAB85F05B9201D08DFDEE2D9D84CA88B371FAE63A4EB07E'))}
    if k not in specs:raise ValueError('Published information lengths: 64,128,256')
    shifts,seeds=specs[k];m=k//4
    h=np.zeros((k,2*k),dtype=np.uint8)
    for a,row in enumerate(shifts):
        for b,values in enumerate(row):
            if values is None:continue
            if isinstance(values,int):values=(values,)
            for shift in values:
                h[m*a+np.arange(m),m*b+(np.arange(m)+shift)%m]^=1
    w=np.zeros((k,k),dtype=np.uint8)
    for a,seed in enumerate(seeds):
        blocks=np.array([int(c) for c in format(int(seed,16),f'0{k}b')],dtype=np.uint8).reshape(4,m)
        for j in range(m):w[a*m+j]=np.roll(blocks,j,axis=1).ravel()
    g=np.c_[np.eye(k,dtype=np.uint8),w]
    if np.any((g@h.T)%2):raise AssertionError('Published generator/parity matrices disagree')
    return g,h


@njit(cache=True,nogil=True)
def _ldpc_decode(llr,neighbors,max_iterations):
    # Flooding sum-product; stable clipping only at tanh/atanh boundaries.
    q=llr[neighbors.ravel()].reshape(neighbors.shape).copy();r=np.zeros_like(q);posterior=llr.copy()
    bits=np.zeros(len(llr),np.uint8)
    for iteration in range(max_iterations):
        for c in range(len(neighbors)):
            for j in range(neighbors.shape[1]):
                product=1.
                for k in range(neighbors.shape[1]):
                    if k!=j:product*=np.tanh(.5*q[c,k])
                product=min(1.-1e-12,max(-1.+1e-12,product))
                r[c,j]=2*np.arctanh(product)
        posterior=llr.copy()
        for c in range(len(neighbors)):
            for j in range(neighbors.shape[1]):posterior[neighbors[c,j]]+=r[c,j]
        bits=(posterior<0).astype(np.uint8)
        valid=True
        for c in range(len(neighbors)):
            parity=0
            for j in range(neighbors.shape[1]):parity^=bits[neighbors[c,j]]
            if parity:valid=False
        if valid:return bits,True,iteration+1
        for c in range(len(neighbors)):
            for j in range(neighbors.shape[1]):q[c,j]=posterior[neighbors[c,j]]-r[c,j]
    return bits,False,max_iterations


class ShortLDPC:
    def __init__(self,k=64):
        self.k=k
        self.g,self.h=ldpc_matrices(k)
        self.neighbors=np.array([np.flatnonzero(row) for row in self.h])

    def encode(self,bits):
        bits=np.asarray(bits,dtype=np.uint8)
        if bits.shape!=(self.k,) or np.any(bits>1):raise ValueError('Wrong information block length')
        return (bits@self.g)%2

    def decode(self,llr,max_iterations=100):
        llr=np.asarray(llr,dtype=float)
        if llr.shape!=(2*self.k,) or not np.all(np.isfinite(llr)):raise ValueError('Wrong channel metric count')
        return _ldpc_decode(llr,self.neighbors,max_iterations)
