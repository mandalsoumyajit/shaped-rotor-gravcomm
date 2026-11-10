"""General word sets for the hop-limited groups: for each start tone, 512 of its admissible 4-tone
continuations. The released enumerative code (HopCode) takes the lexicographically first 512, which
skews tone use toward the low tones (tone 6 carries 3.4% of symbols, tones 0-2 60%).

Selections:
  'lex'   lexicographically first (reproduces HopCode),
  'dist'  greedy pruning of the most confusable continuation (smallest baseband distance to its
          nearest neighbour) until 512 remain; distance over the 4 group symbols plus one common
          continuation symbol, so the end-phase difference counts.
Word order inside a set (the map from GF(512) symbol to sequence) is irrelevant for the non-binary
code (random coefficients), so sets are stored in lexicographic order.
Usage: python word_tables.py   (writes results/words_<name>.json)
"""
import itertools, json
from fractions import Fraction
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent/'results'
M, D, L, BITS = 7, 3, 4, 9
H, BETA, T = Fraction(3, 14), 1.0, 224/190


def continuations(p):
    seqs = [s for s in itertools.product(range(M), repeat=L)
            if all(abs(a - b) <= D for a, b in zip((p,) + s[:-1], s))]
    return np.array(seqs)                                    # itertools order = lexicographic


class TableHopCode:
    def __init__(self, words):
        self.words = [np.asarray(w, np.int64) for w in words]       # words[p]: (512, L)
        self.lookup = [np.full((M, M**(i + 1)), -1, np.int64) for i in range(L)]
        for p in range(M):
            for idx, s in enumerate(self.words[p]):
                code = 0
                for i in range(L):
                    code = code*M + s[i]
                    self.lookup[i][p, code] = idx if i == L - 1 else 0

    def encode(self, word, start):
        return [int(x) for x in self.words[start][word]]

    def decode(self, seq, start):
        code = 0
        for c in seq:
            code = code*M + c
        return int(self.lookup[L - 1][start, code])


def baseband(p, seqs, fs=16.0, tail=1, tail_tone=M//2):
    """Complex baseband of [p] + seq + tail, sampled over the group and tail symbols, phase 0 at group start."""
    ns = int(round(T*fs)); u = (np.arange(ns) + .5)/ns
    shape = 2.5*u**4 - 3*u**5 + u**6                       # integral of the transition (beta = 1), in units of T
    spacing = float(H)/T; off = lambda a: a - (M - 1)/2
    out = []
    for s in seqs:
        tones = [p] + list(s) + [tail_tone]*tail
        phase0 = 0.0; seg = []
        for a, b in zip(tones[:-1], tones[1:]):
            ph = phase0 + 2*np.pi*spacing*T*(off(a)*u + (off(b) - off(a))*shape)
            seg.append(np.exp(1j*ph))
            phase0 += 2*np.pi*spacing*T*(off(a) + (off(b) - off(a))*0.5)
        out.append(np.concatenate(seg))
    return np.array(out)


def select_dist(p, keep=2**BITS):
    seqs = continuations(p)
    X = baseband(p, seqs)
    G = X @ X.conj().T; nrm = np.real(np.diag(G))
    D2 = nrm[:, None] + nrm[None, :] - 2*np.real(G); np.fill_diagonal(D2, np.inf)
    alive = np.ones(len(seqs), bool)
    while alive.sum() > keep:
        sub = np.where(alive[:, None] & alive[None, :], D2, np.inf)
        i, j = np.unravel_index(np.argmin(sub), sub.shape)
        # drop the one of the pair with the smaller second-nearest distance
        def second(k):
            r = np.sort(sub[k]); return r[1]
        drop = i if second(i) < second(j) else j
        alive[drop] = False
    return seqs[alive], float(np.sqrt(np.min(np.where(alive[:, None] & alive[None, :], D2, np.inf))))


def stats(words):
    allw = np.concatenate(words)
    use = np.bincount(allw.ravel(), minlength=M)/allw.size
    return use


if __name__ == '__main__':
    lex = [continuations(p)[:2**BITS] for p in range(M)]
    dist, dmin = [], []
    for p in range(M):
        w, d = select_dist(p); dist.append(w); dmin.append(d)
        lexd = baseband(p, lex[p]); G = lexd @ lexd.conj().T; n = np.real(np.diag(G))
        D2 = n[:, None] + n[None, :] - 2*np.real(G); np.fill_diagonal(D2, np.inf)
        print(f'start {p}: N={len(continuations(p))}, d_min lex {np.sqrt(D2.min()):.3f}, dist {d:.3f}')
    for name, w in (('lex', lex), ('dist', dist)):
        print(name, 'tone use', np.round(stats(w), 3))
        (OUT/f'words_{name}.json').write_text(json.dumps([x.tolist() for x in w]))
