"""Anytime confidence for the mean error fraction of independent packets.

For X in [0,1] with mean m, each product of 1+lambda(m)*(X-m) is a
nonnegative mean-one martingale. Fixed mixtures and Ville's inequality give
time-uniform bounds even when stopping depends on accumulated errors.
This is a simple fixed-bet mixture, not the optimized algorithm from
Waudby-Smith and Ramdas, https://arxiv.org/abs/2010.09686 .
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp


def packet_ber_interval(error_counts,bits_per_packet=64,alpha=.05):
    counts=np.asarray(error_counts)
    if counts.ndim!=1 or not len(counts) or np.any(counts<0) or np.any(counts>bits_per_packet):
        raise ValueError('invalid packet error counts')
    values,frequency=np.unique(counts/bits_per_packet,return_counts=True)
    stakes=np.array([.001,.003,.01,.03,.1,.3,.5,.7,.9,.99])
    threshold=np.log(2/alpha)

    def capital(mean,side):
        bets=stakes/mean if side=='lower' else -stakes/(1-mean)
        logs=np.log1p(bets[:,None]*(values[None,:]-mean))@frequency
        return float(logsumexp(logs)-np.log(len(stakes))-threshold)

    eps=1e-14
    lower=0. if np.all(counts==0) else brentq(lambda m:capital(m,'lower'),eps,1-eps,xtol=1e-14)
    upper=1. if np.all(counts==bits_per_packet) else brentq(lambda m:capital(m,'upper'),eps,1-eps,xtol=1e-14)
    return [float(lower),float(upper)]
