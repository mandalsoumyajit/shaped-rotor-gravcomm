"""Iterative (turbo) demodulation and decoding for the CF-FSK links.

Joint decoding of the tone trellis and the convolutional code, as in serially concatenated CPM
(Moqvist & Aulin 2001) and turbo equalization (Koetter, Singer & Tuchler 2004):
  * the released streaming detector (or the hop-limited fork) runs once and its complete trellis
    graph is kept: every retained edge with its whitened-residual metric and word label;
  * a log-domain forward-backward pass over that graph, with a prior on each word's coded bits,
    gives a-posteriori bit LLRs; their extrinsic part is deinterleaved and passed to
  * a soft-in/soft-out (log-MAP, circular tail-biting) decoder of the rate-2/3 code, whose
    extrinsic coded-bit LLRs are interleaved back as the detector's word priors.
The survivor selection of the reduced-state detector is made once, without priors; later passes
reuse its graph. LLR convention: positive favours bit 0. Edge metrics are -log likelihoods.
"""
import math
from collections import deque
import numpy as np
from numba import njit
import mitigation  # noqa: F401  (puts the released gravcomm package on the path)
from gravcomm import convolutional_candidates as cc


# ---------------------------------------------------------------- convolutional code SISO
@njit(cache=True)
def _lse(a, b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    m = max(a, b)
    return m + math.log(math.exp(a - m) + math.exp(b - m))


@njit(cache=True)
def _bcjr_tb(L, outputs, laps=2):
    """Circular log-MAP for the 64-state tail-biting code; L is (K, 2) mother-code LLRs (0 = punctured)."""
    K = L.shape[0]
    g = np.empty((K, 64, 2))
    for t in range(K):
        for s in range(64):
            for b in range(2):
                g[t, s, b] = -(outputs[s, b, 0]*L[t, 0] + outputs[s, b, 1]*L[t, 1])
    alpha = np.zeros((K + 1, 64)); beta = np.zeros((K + 1, 64))
    a = np.zeros(64)
    for lap in range(laps):
        alpha[0] = a
        for t in range(K):
            nxt = np.full(64, -np.inf)
            for s in range(64):
                for b in range(2):
                    d = ((s << 1) | b) & 63
                    nxt[d] = _lse(nxt[d], alpha[t, s] + g[t, s, b])
            m = nxt.max()
            alpha[t + 1] = nxt - m
        a = alpha[K].copy()
    bb = np.zeros(64)
    for lap in range(laps):
        beta[K] = bb
        for t in range(K - 1, -1, -1):
            cur = np.full(64, -np.inf)
            for s in range(64):
                for b in range(2):
                    d = ((s << 1) | b) & 63
                    cur[s] = _lse(cur[s], g[t, s, b] + beta[t + 1, d])
            m = cur.max()
            beta[t] = cur - m
        bb = beta[0].copy()
    info = np.empty(K); app = np.empty((K, 2))
    for t in range(K):
        u0 = -np.inf; u1 = -np.inf; c00 = -np.inf; c01 = -np.inf; c10 = -np.inf; c11 = -np.inf
        for s in range(64):
            for b in range(2):
                d = ((s << 1) | b) & 63
                v = alpha[t, s] + g[t, s, b] + beta[t + 1, d]
                if b == 0:
                    u0 = _lse(u0, v)
                else:
                    u1 = _lse(u1, v)
                if outputs[s, b, 0] == 0:
                    c00 = _lse(c00, v)
                else:
                    c01 = _lse(c01, v)
                if outputs[s, b, 1] == 0:
                    c10 = _lse(c10, v)
                else:
                    c11 = _lse(c11, v)
        info[t] = u0 - u1; app[t, 0] = c00 - c01; app[t, 1] = c10 - c11
    return info, app


def siso_decode(llr, rate='2/3', information_bits=256):
    """Return (info-bit APP LLRs, extrinsic coded-bit LLRs) for the punctured tail-biting code."""
    mask = cc._mask(information_bits, rate)
    L = np.zeros(mask.shape); L[mask] = llr
    info, app = _bcjr_tb(L, cc.OUTPUTS)
    return info, (app - L)[mask]


# ---------------------------------------------------------------- detector graph SISO
class _Recorder(deque):
    """Deque that also keeps every appended graph step (the detector pops consumed steps)."""
    def __init__(self):
        super().__init__(); self.steps = []

    def append(self, x):
        self.steps.append(x); super().append(x)


def record_graph(detector):
    """Attach a recorder to a freshly constructed streaming detector (call before pushing samples)."""
    detector.graph = _Recorder()
    return detector


def pack_graph(detector, group_symbols, ngroups):
    """Flatten the recorded graph into arrays for graph_app."""
    steps = detector.graph.steps
    nodes = [len(s[1]) for s in steps] + [len(detector.cost)]
    eoff = np.cumsum([0] + [len(s[2]) for s in steps])
    src = np.concatenate([s[2] for s in steps]).astype(np.int64)
    dst = np.concatenate([s[3] for s in steps]).astype(np.int64)
    met = -np.concatenate([s[4] for s in steps]).astype(np.float64)
    word = np.concatenate([s[5] for s in steps]).astype(np.int64)
    group = np.array([s[0]//group_symbols if (s[0] % group_symbols == group_symbols - 1 and
                                              s[0] < group_symbols*ngroups) else -1 for s in steps], np.int64)
    noff = np.cumsum([0] + nodes)
    return dict(src=src, dst=dst, met=met, word=word, eoff=eoff, noff=noff, group=group)


@njit(cache=True)
def _graph_app(src, dst, met, word, eoff, noff, group, La, labels):
    """Log-domain forward-backward over the detector graph with word priors; returns APP LLRs (groups, bits)."""
    T = len(eoff) - 1; nb = La.shape[1]
    alpha = np.full(noff[-1], -np.inf); beta = np.full(noff[-1], -np.inf)
    alpha[noff[0]] = 0.0
    prior = np.zeros(len(met))
    for t in range(T):
        g = group[t]
        if g >= 0:
            for e in range(eoff[t], eoff[t + 1]):
                p = 0.0
                for k in range(nb):
                    if labels[word[e], k]:
                        p -= La[g, k]
                prior[e] = p
    for t in range(T):
        for e in range(eoff[t], eoff[t + 1]):
            i = noff[t + 1] + dst[e]
            alpha[i] = _lse(alpha[i], alpha[noff[t] + src[e]] + met[e] + prior[e])
        m = -np.inf
        for i in range(noff[t + 1], noff[t + 2]):
            m = max(m, alpha[i])
        for i in range(noff[t + 1], noff[t + 2]):
            alpha[i] -= m
    for i in range(noff[T], noff[T + 1]):
        beta[i] = 0.0
    for t in range(T - 1, -1, -1):
        for e in range(eoff[t], eoff[t + 1]):
            i = noff[t] + src[e]
            beta[i] = _lse(beta[i], met[e] + prior[e] + beta[noff[t + 1] + dst[e]])
        m = -np.inf
        for i in range(noff[t], noff[t + 1]):
            m = max(m, beta[i])
        for i in range(noff[t], noff[t + 1]):
            beta[i] -= m
    G = La.shape[0]
    l0 = np.full((G, nb), -np.inf); l1 = np.full((G, nb), -np.inf)
    for t in range(T):
        g = group[t]
        if g < 0:
            continue
        for e in range(eoff[t], eoff[t + 1]):
            v = alpha[noff[t] + src[e]] + met[e] + prior[e] + beta[noff[t + 1] + dst[e]]
            for k in range(nb):
                if labels[word[e], k]:
                    l1[g, k] = _lse(l1[g, k], v)
                else:
                    l0[g, k] = _lse(l0[g, k], v)
    return l0 - l1


def graph_app(packed, La, labels):
    return _graph_app(packed['src'], packed['dst'], packed['met'], packed['word'], packed['eoff'],
                      packed['noff'], packed['group'], np.ascontiguousarray(La, dtype=np.float64), labels)


def bit_labels(nwords, nbits):
    return ((np.arange(nwords)[:, None] >> np.arange(nbits - 1, -1, -1)) & 1).astype(np.bool_)


def turbo(packed, permutation, ngroups, nbits, n_coded, iterations=8, information_bits=256, crc_ok=None,
          pad_bits=0, damping=1.0, labels=None, siso=None):
    """Iterate detector-graph APP and code SISO. Returns per-iteration hard info bits (and stop index).

    permutation: coded-bit interleaver as in the transmitter (channel position i carries word[permutation[i]]).
    crc_ok: optional callable(info_bits) for early stopping.
    """
    if labels is None:
        labels = bit_labels(2**nbits, nbits)                  # labels[word] = coded-bit pattern of that word
    if siso is None:
        siso = lambda llr: siso_decode(llr, information_bits=information_bits)
    La_ch = np.zeros(ngroups*nbits)                         # detector prior, channel (interleaved) order
    if pad_bits:
        La_ch[n_coded:] = 30.0                               # known zero pad bits
    history = []
    for it in range(iterations):
        app = graph_app(packed, La_ch.reshape(ngroups, nbits), labels).ravel()
        ext = (app - La_ch)[:n_coded]
        dec_in = np.empty(n_coded); dec_in[permutation] = ext
        info, cext = siso(dec_in)
        hard = (info < 0).astype(np.uint8)
        history.append(hard)
        if crc_ok is not None and crc_ok(hard):
            break
        La_ch[:n_coded] = damping*cext[permutation]
    return history


@njit(cache=True)
def _graph_word_post(src, dst, met, word, eoff, noff, group, G, nwords):
    """Word log-posteriors (no priors) at each word-final step of the detector graph."""
    T = len(eoff) - 1
    alpha = np.full(noff[-1], -np.inf); beta = np.full(noff[-1], -np.inf)
    alpha[noff[0]] = 0.0
    for t in range(T):
        for e in range(eoff[t], eoff[t + 1]):
            i = noff[t + 1] + dst[e]
            alpha[i] = _lse(alpha[i], alpha[noff[t] + src[e]] + met[e])
        m = -np.inf
        for i in range(noff[t + 1], noff[t + 2]):
            m = max(m, alpha[i])
        for i in range(noff[t + 1], noff[t + 2]):
            alpha[i] -= m
    for i in range(noff[T], noff[T + 1]):
        beta[i] = 0.0
    for t in range(T - 1, -1, -1):
        for e in range(eoff[t], eoff[t + 1]):
            i = noff[t] + src[e]
            beta[i] = _lse(beta[i], met[e] + beta[noff[t + 1] + dst[e]])
        m = -np.inf
        for i in range(noff[t], noff[t + 1]):
            m = max(m, beta[i])
        for i in range(noff[t], noff[t + 1]):
            beta[i] -= m
    post = np.full((G, nwords), -np.inf)
    for t in range(T):
        g = group[t]
        if g < 0:
            continue
        for e in range(eoff[t], eoff[t + 1]):
            v = alpha[noff[t] + src[e]] + met[e] + beta[noff[t + 1] + dst[e]]
            post[g, word[e]] = _lse(post[g, word[e]], v)
    for g in range(G):
        m = -np.inf
        for w in range(nwords):
            m = max(m, post[g, w])
        s = 0.0
        for w in range(nwords):
            if post[g, w] > -np.inf:
                s += math.exp(post[g, w] - m)
        for w in range(nwords):
            post[g, w] -= m + math.log(s)
    return post


def word_posteriors(packed, ngroups, nwords):
    return _graph_word_post(packed['src'], packed['dst'], packed['met'], packed['word'], packed['eoff'],
                            packed['noff'], packed['group'], ngroups, nwords)
