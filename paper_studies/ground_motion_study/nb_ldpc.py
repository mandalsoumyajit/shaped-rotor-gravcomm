"""Non-binary LDPC codes over GF(2^9) matched to the 512-word hop-limited tone groups.

Each code symbol is one 4-tone group of the hop-limited waveform, so the streaming detector's
word posteriors are the decoder's channel messages directly: no bit labeling and no bit
marginalization (cf. Liva et al. 2016 for non-binary LDPC codes at short block lengths).

Contents: GF(512) arithmetic, progressive-edge-growth (PEG) construction, systematic encoding by
Gaussian elimination, and a probability-domain FFT belief-propagation decoder (Walsh-Hadamard
transform over the additive group of GF(2^9)).
"""
import numpy as np
from numba import njit

P, Q = 9, 512
PRIM = 0x211                                                  # x^9 + x^4 + 1 (primitive)


def _tables():
    exp = np.zeros(2*Q, np.int64); log = np.full(Q, -1, np.int64)
    x = 1
    for i in range(Q - 1):
        exp[i] = x; log[x] = i
        x <<= 1
        if x & Q:
            x ^= PRIM
    exp[Q - 1:2*Q - 2] = exp[:Q - 1]
    return exp, log


EXP, LOG = _tables()
assert len(set(EXP[:Q - 1].tolist())) == Q - 1


def gmul(a, b):
    return 0 if a == 0 or b == 0 else int(EXP[LOG[a] + LOG[b]])


def ginv(a):
    return int(EXP[(Q - 1 - LOG[a]) % (Q - 1)])


# MUL[h] = permutation a -> h*a (row h); used for coefficient permutations of message vectors
MUL = np.zeros((Q, Q), np.int64)
for _h in range(1, Q):
    MUL[_h, 1:] = EXP[LOG[_h] + LOG[np.arange(1, Q)]]


def peg(n, m, dv, seed=0):
    """Progressive edge growth with variable degree dv; returns list of edges (check, var)."""
    rng = np.random.default_rng(seed)
    vnb = [[] for _ in range(n)]; cnb = [[] for _ in range(m)]
    for j in range(n):
        for k in range(dv):
            deg = np.array([len(c) for c in cnb])
            if k == 0:
                cand = np.flatnonzero(deg == deg.min())
            else:
                # BFS from j: find checks not reachable, else the ones reached last
                seen_c = set(vnb[j]); frontier_v = {j}; seen_v = {j}; last = set(vnb[j])
                while True:
                    new_c = set()
                    for v in frontier_v:
                        for c in vnb[v]:
                            for v2 in cnb[c]:
                                if v2 not in seen_v:
                                    seen_v.add(v2)
                                    for c2 in vnb[v2]:
                                        if c2 not in seen_c:
                                            new_c.add(c2)
                    frontier_v = {v2 for c in new_c for v2 in cnb[c]} - {j}
                    if not new_c or len(seen_c | new_c) == m:
                        unreached = [c for c in range(m) if c not in (seen_c | new_c)]
                        pool = unreached if unreached else [c for c in new_c if c not in vnb[j]] or \
                            [c for c in range(m) if c not in vnb[j]]
                        break
                    seen_c |= new_c
                pool = np.array(sorted(pool))
                cand = pool[deg[pool] == deg[pool].min()]
            c = int(rng.choice(cand))
            vnb[j].append(c); cnb[c].append(j)
    return [(c, j) for j in range(n) for c in vnb[j]]


class NBCode:
    """Systematic non-binary LDPC code; H given as edges with GF(512) coefficients."""

    def __init__(self, n, m, edges, coef):
        self.n, self.m = n, m
        self.edges = np.array(edges, np.int64); self.coef = np.array(coef, np.int64)
        H = np.zeros((m, n), np.int64)
        for (c, j), h in zip(edges, coef):
            H[c, j] = h
        self.H = H
        A = H.copy(); piv = []; row = 0
        for col in range(n - 1, -1, -1):                     # pivots from the right: parity at the end
            r = next((i for i in range(row, m) if A[i, col]), None)
            if r is None:
                continue
            A[[row, r]] = A[[r, row]]
            inv = ginv(A[row, col]); A[row] = [gmul(inv, x) for x in A[row]]
            for i in range(m):
                if i != row and A[i, col]:
                    f = A[i, col]; A[i] ^= np.array([gmul(f, x) for x in A[row]])
            piv.append(col); row += 1
            if row == m:
                break
        if row < m:
            raise ValueError('H not full rank')
        self.parity = np.array(piv); self.info = np.array(sorted(set(range(n)) - set(piv)))
        self.k = len(self.info)
        self.A = A                                            # reduced: parity[r] = sum_j A[r, info_j] info_j
        # flattened Tanner graph for the decoder
        order = np.lexsort((self.edges[:, 1], self.edges[:, 0]))
        self.ec = self.edges[order, 0]; self.ev = self.edges[order, 1]; self.eh = self.coef[order]
        self.coff = np.searchsorted(self.ec, np.arange(m + 1))
        vorder = np.argsort(self.ev, kind='stable')
        self.vperm = vorder; self.voff = np.searchsorted(self.ev[vorder], np.arange(n + 1))

    def encode(self, u):
        x = np.zeros(self.n, np.int64); x[self.info] = u
        for r, col in enumerate(self.parity):
            s = 0
            for j in self.info:
                if self.A[r, j]:
                    s ^= gmul(self.A[r, j], int(x[j]))
            x[col] = s
        assert self.syndrome_ok(x)
        return x

    def syndrome_ok(self, x):
        for c in range(self.m):
            s = 0
            for j in np.flatnonzero(self.H[c]):
                s ^= gmul(int(self.H[c, j]), int(x[j]))
            if s:
                return False
        return True

    def decode(self, prior, iterations=30, state=None):
        """prior: (n, Q) probabilities. Returns (posterior, extrinsic, hard, ok); with state (check-to-variable
        messages from a previous call, or True for a fresh start) also returns the updated state (warm start)."""
        c2v = np.full((len(self.ec), Q), 1.0/Q) if (state is None or state is True) else state
        post, ext, ok = _bp(np.ascontiguousarray(prior, np.float64), self.ec, self.ev, self.eh, self.coff,
                            self.vperm, self.voff, MUL, iterations, c2v)
        if state is None:
            return post, ext, post.argmax(1), bool(ok)
        return post, ext, post.argmax(1), bool(ok), c2v


@njit(cache=True)
def _wht(x):
    h = 1; n = x.shape[0]
    while h < n:
        for i in range(0, n, 2*h):
            for j in range(i, i + h):
                a = x[j]; b = x[j + h]
                x[j] = a + b; x[j + h] = a - b
        h *= 2


@njit(cache=True)
def _bp(prior, ec, ev, eh, coff, vperm, voff, MUL, iterations, c2v):
    n, q = prior.shape; E = len(ec); m = len(coff) - 1
    INV = np.zeros_like(MUL)
    for h in range(1, q):
        for a in range(q):
            INV[h, MUL[h, a]] = a
    v2c = np.empty((E, q))
    post = prior.copy(); ok = 0
    F = np.empty((E, q)); fwd = np.empty(q); tmp = np.empty(q)
    for j in range(n):                                       # variable-to-check from prior and incoming c2v
        for k in range(voff[j], voff[j + 1]):
            e = vperm[k]; s = 0.0
            for b in range(q):
                tmp[b] = prior[j, b]
            for k2 in range(voff[j], voff[j + 1]):
                if k2 != k:
                    for b in range(q):
                        tmp[b] *= c2v[vperm[k2], b]
            for b in range(q):
                s += tmp[b]
            for b in range(q):
                v2c[e, b] = tmp[b]/s if s > 0 else prior[j, b]
    for it in range(iterations):
        # check nodes
        for c in range(m):
            for e in range(coff[c], coff[c + 1]):
                h = eh[e]
                for a in range(q):
                    tmp[MUL[h, a]] = v2c[e, a]
                _wht(tmp)
                F[e] = tmp
            for e in range(coff[c], coff[c + 1]):
                for b in range(q):
                    fwd[b] = 1.0
                for e2 in range(coff[c], coff[c + 1]):
                    if e2 != e:
                        for b in range(q):
                            fwd[b] *= F[e2, b]
                _wht(fwd)
                h = eh[e]; s = 0.0
                for a in range(q):
                    v = fwd[MUL[h, a]]
                    v = v if v > 0 else 0.0
                    c2v[e, a] = v; s += v
                if s > 0:
                    for a in range(q):
                        c2v[e, a] /= s
                else:
                    for a in range(q):
                        c2v[e, a] = 1.0/q
        # variable nodes
        for j in range(n):
            for b in range(q):
                fwd[b] = prior[j, b]
            for k in range(voff[j], voff[j + 1]):
                e = vperm[k]
                for b in range(q):
                    fwd[b] *= c2v[e, b]
                s = 0.0
                for b in range(q):
                    s += fwd[b]
                if s > 0:
                    for b in range(q):
                        fwd[b] /= s
                else:
                    for b in range(q):
                        fwd[b] = prior[j, b]
            post[j] = fwd
            for k in range(voff[j], voff[j + 1]):
                e = vperm[k]; s = 0.0
                for b in range(q):
                    tmp[b] = prior[j, b]
                for k2 in range(voff[j], voff[j + 1]):
                    if k2 != k:
                        e2 = vperm[k2]
                        for b in range(q):
                            tmp[b] *= c2v[e2, b]
                for b in range(q):
                    s += tmp[b]
                for b in range(q):
                    v2c[e, b] = tmp[b]/s if s > 0 else prior[j, b]
        # syndrome check on hard decisions
        ok = 1
        for c in range(m):
            s = 0
            for e in range(coff[c], coff[c + 1]):
                j = ev[e]; best = 0
                for b in range(q):
                    if post[j, b] > post[j, best]:
                        best = b
                if best:
                    s ^= MUL[eh[e], best]
            if s:
                ok = 0; break
        if ok:
            break
    ext = np.ones((n, q))
    for j in range(n):
        for k in range(voff[j], voff[j + 1]):
            e = vperm[k]
            for b in range(q):
                ext[j, b] *= c2v[e, b]
        s = 0.0
        for b in range(q):
            s += ext[j, b]
        for b in range(q):
            ext[j, b] = ext[j, b]/s if s > 0 else 1.0/q
    return post, ext, ok


def make_code(n=43, k=29, dv=2, seed=1):
    m = n - k
    edges = peg(n, m, dv, seed)
    rng = np.random.default_rng(1000 + seed)
    coef = rng.integers(1, Q, len(edges))
    return NBCode(n, m, edges, coef)


if __name__ == '__main__':
    import time
    code = make_code()
    print('n, k, edges', code.n, code.k, len(code.edges), 'check degrees', np.bincount(code.edges[:, 0]))
    rng = np.random.default_rng(0)
    # sanity: q-ary symmetric channel with symbol error probability p
    for p in (0.15, 0.25, 0.3):
        fails = 0; t0 = time.time()
        for trial in range(200):
            u = rng.integers(0, Q, code.k); x = code.encode(u)
            y = np.where(rng.random(code.n) < p, rng.integers(0, Q, code.n), x)
            prior = np.full((code.n, Q), p/(Q - 1)); prior[np.arange(code.n), y] = 1 - p
            post, ext, hard, ok = code.decode(prior)
            fails += int(not (ok and np.array_equal(hard, x)))
        print(f'p={p}: FER {fails/200:.3f} ({time.time() - t0:.1f}s)')
