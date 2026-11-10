"""Exact word-level APP detection of the released CF-FSK mapping: separates mapping and detector losses.

The released link maps 8 coded bits to 3 base-7 tones (256 of 343 words). Using the same
finite-state auxiliary channel as cpm_rate.info_rate (whitening predictor, exact phase lattice),
this script runs a word-level forward-backward (BCJR) recursion with a uniform prior over the 256
valid words and returns, for one long simulated sequence:
  * the information rate of the actual word input (bits per symbol),
  * exact bit LLRs (log-sum over words) and max-log bit LLRs (max over words),
whose generalized mutual information (GMI) is the bit-interleaved rate with an optimal APP
detector. Comparing with the released detector's GMI (gap_analysis.py) isolates the loss of
its approximations (reduced state, fixed lag, causal mean history) from the mapping loss.
Output: results/word_bcjr.json
"""
import json
from pathlib import Path
import numpy as np
import numba as nb
import cpm_rate as cr
import cpm_optimize as co
from gap_analysis import gmi

OUT = Path(__file__).resolve().parent/'results'/'word_bcjr.json'
WORDS = np.array([[w//49, (w//7) % 7, w % 7] for w in range(256)], np.int64)


@nb.njit(cache=True)
def _word_metrics(C3, E, labels, k, P, M, amp, sig2, words):
    """Metric of every (start state, word) for one word step; C3 holds 3 symbols of correlations.

    State: (psi, s0, s1, s2) with s2 the latest tone; flat index ((psi*M + s0)*M + s1)*M + s2.
    Returns metric[state, word] and end-state index.
    """
    S = P*M**3; W = words.shape[0]
    met = np.empty((S, W)); end = np.empty((S, W), np.int64)
    rot = np.exp(-2j*np.pi*np.arange(P)/P)
    for st in range(S):
        s2 = st % M; s1 = (st//M) % M; s0 = (st//(M*M)) % M; psi = st//(M**3)
        for w in range(W):
            d1, d2, d3 = words[w, 0], words[w, 1], words[w, 2]
            p1 = psi; t1 = ((s0*M + s1)*M + s2)*M + d1
            p2 = (psi + k*labels[s0]) % P; t2 = ((s1*M + s2)*M + d1)*M + d2
            p3 = (p2 + k*labels[s1]) % P; t3 = ((s2*M + d1)*M + d2)*M + d3
            pe = (p3 + k*labels[s2]) % P
            m = 0.0
            m += (2*amp*(rot[p1]*C3[0, t1]).real - amp*amp*E[t1])/sig2
            m += (2*amp*(rot[p2]*C3[1, t2]).real - amp*amp*E[t2])/sig2
            m += (2*amp*(rot[p3]*C3[2, t3]).real - amp*amp*E[t3])/sig2
            met[st, w] = m
            end[st, w] = ((pe*M + d1)*M + d2)*M + d3
    return met, end


@nb.njit(cache=True)
def _lse(a):
    mx = -np.inf
    for x in a:
        if x > mx:
            mx = x
    if mx == -np.inf:
        return mx
    s = 0.0
    for x in a:
        s += np.exp(x - mx)
    return mx + np.log(s)


def run(fc, amp, n_words=3000, seed=11, p_max=12, return_post=False):
    d = cr.Design(**co.CF)
    M = d.M; ns = max(2, int(round(4.0*d.T))); fs = ns/d.T
    labels, k, P = d.lattice()
    r, p = cr.memory(d, ns, p_max, P)
    if r != 2:
        raise ValueError('word recursion assumes a 3-symbol state (r = 2)')
    rng = np.random.default_rng(seed)
    wseq = rng.integers(0, 256, n_words)
    start = M//2
    seq = np.r_[np.full(r + 1, start), WORDS[wseq].ravel()].astype(np.int64)
    hb = float(d.h)*d.beta/2
    x = np.exp(2j*np.pi*(cr.phase_track(d, seq, ns) - hb*d.tones()[start]))
    npsd = co.noise_psd(fc); L = len(x)
    S = fs*npsd(np.fft.fftfreq(L, 1/fs))
    w = np.fft.ifft(np.fft.fft((rng.normal(size=L) + 1j*rng.normal(size=L))/np.sqrt(2))*np.sqrt(S))
    y = amp*x + w
    a, sig2 = cr.predictor(fs*npsd(np.fft.fftfreq(1 << 14, 1/fs)), p)
    yw = cr.whiten(y, a).reshape(len(seq), ns)[r + 1:]
    tup = np.indices((M,)*(r + 2)).reshape(r + 2, -1).T
    tmpl = np.empty((len(tup), ns), complex)
    for t, sy in enumerate(tup):
        ph = cr.phase_track(d, sy, ns)[ns:] - hb*d.tones()[sy[0]]
        tmpl[t] = cr.whiten(np.exp(2j*np.pi*ph), a)[-ns:]
    E = np.sum(np.abs(tmpl)**2, axis=1)
    Cs = yw @ tmpl.conj().T                                  # (symbols, templates)
    Sn = P*M**3; logw = -np.log(256)
    alpha = np.full((n_words + 1, Sn), -np.inf); alpha[0, ((0*M + start)*M + start)*M + start] = 0.0
    logqy = 0.0; logqya = 0.0
    for j in range(n_words):
        met, end = _word_metrics(Cs[3*j:3*j + 3], E, labels, k, P, M, amp, sig2, WORDS)
        logqya += met[_state_index(seq, j, M, labels, k, P), wseq[j]]
        v = (alpha[j][:, None] + met + logw).ravel(); e = end.ravel()
        order = np.argsort(e, kind='stable')
        ue, idx = np.unique(e[order], return_index=True)
        nxt = np.full(Sn, -np.inf); nxt[ue] = np.logaddexp.reduceat(v[order], idx)
        nm = nxt.max(); alpha[j + 1] = nxt - nm; logqy += nm
    logqy += np.log(np.sum(np.exp(alpha[n_words])))
    rate_word_input = (logqya - logqy)/(3*n_words*np.log(2))
    beta = np.zeros(Sn)
    llr_exact = np.empty((n_words, 8)); llr_max = np.empty((n_words, 8)); posts = np.empty((n_words, 256))
    bits = ((wseq[:, None] >> np.arange(7, -1, -1)) & 1).astype(np.int8)
    wbits = ((np.arange(256)[:, None] >> np.arange(7, -1, -1)) & 1).astype(bool)
    for j in range(n_words - 1, -1, -1):
        met, end = _word_metrics(Cs[3*j:3*j + 3], E, labels, k, P, M, amp, sig2, WORDS)   # recomputed
        tot = alpha[j][:, None] + met + beta[end]
        post = np.logaddexp.reduce(tot, axis=0); posts[j] = post - np.logaddexp.reduce(post)
        mx = tot.max(axis=0)
        for bb in range(8):
            llr_exact[j, bb] = np.logaddexp.reduce(post[~wbits[:, bb]]) - np.logaddexp.reduce(post[wbits[:, bb]])
            llr_max[j, bb] = mx[~wbits[:, bb]].max() - mx[wbits[:, bb]].max()
        bv = np.logaddexp.reduce(met + logw + beta[end], axis=1)
        beta = bv - bv.max()
    if return_post:
        return posts, wseq
    ge = gmi(llr_exact.ravel(), bits.ravel()); gm = gmi(llr_max.ravel(), bits.ravel())
    return dict(fc=fc, amp=amp, rate_word_input_bits_per_symbol=rate_word_input, gmi_exact=ge[0], gmi_maxlog=gm[0],
                scale_exact=ge[2], scale_maxlog=gm[2], p=p, n_words=n_words)


def _state_index(seq, j, M, labels, k, P):
    """Trellis state at the start of word j (symbols seq[3j .. 3j+2] are its 3-symbol history)."""
    psi = 0
    for i in range(3*j):
        psi = (psi + k*labels[seq[i]]) % P
    s0, s1, s2 = seq[3*j], seq[3*j + 1], seq[3*j + 2]
    return ((psi*M + s0)*M + s1)*M + s2


def _job(args):
    return run(*args)


def main(workers=5):
    from concurrent.futures import ProcessPoolExecutor
    grids = {24.: (0.40, 0.48, 0.56, 0.64, 0.72), 18.: (0.48, 0.57, 0.66, 0.75, 0.84)}
    jobs = [(fc, a) for fc, amps in grids.items() for a in amps]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        out = list(pool.map(_job, jobs))
    for r in out:
        print(json.dumps(r), flush=True)
    OUT.write_text(json.dumps(out, indent=1, default=float) + chr(10))


if __name__ == '__main__':
    main()
