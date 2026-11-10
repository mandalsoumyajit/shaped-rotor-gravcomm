"""Rank GF(512) LDPC constructions on captured detector word posteriors (surrogate channel).

results/posteriors_hop_f24_a0.550_408.npz holds the non-iterative word posteriors of 408 hop
packets at 24 Hz, 0.55 A0. Each packet (43 consecutive groups) is one surrogate frame; a random
codeword is carried by XOR scrambling (channel symbol = codeword symbol XOR s, with s chosen so the
channel symbol equals the captured transmitted word), so the clustered error statistics of the real
detector are kept. Decoding is BP only (no detector feedback): this ranks codes, it is not a
link performance figure. Usage: python nb_design.py [dv_list] [seeds]
"""
import sys, json, time
from pathlib import Path
import numpy as np
import nb_ldpc as nb

OUT = Path(__file__).resolve().parent/'results'
d = np.load(OUT/'posteriors_hop_f24_a0.550_408.npz')
POSTS = d['posts'].astype(np.float64).reshape(-1, 43, 512); WORDS = d['words'].reshape(-1, 43)


def evaluate(code, frames=None, temper=1.0, seed=0):
    rng = np.random.default_rng(seed)
    idx = np.arange(len(POSTS)) if frames is None else frames
    fails = 0; serr = 0
    a = np.arange(512)
    for f in idx:
        u = rng.integers(0, 512, code.k); x = code.encode(u)
        P = np.exp(temper*(POSTS[f] - POSTS[f].max(1, keepdims=True)))
        prior = np.empty((43, 512))
        for j in range(43):
            s = x[j] ^ WORDS[f, j]
            prior[j] = P[j, a ^ s]
        prior /= prior.sum(1, keepdims=True)
        _, _, hard, ok = code.decode(prior, 50)
        fails += int(not (ok and np.array_equal(hard, x)))
        serr += int(np.sum(prior.argmax(1) != x))
    return fails/len(idx), serr/(len(idx)*43)


def mixed_code(n, k, n3, seed):
    """dv=2 for most symbols, dv=3 for n3 of them (PEG with per-variable degree)."""
    m = n - k
    rng = np.random.default_rng(seed)
    degs = np.r_[np.full(n3, 3), np.full(n - n3, 2)]
    vnb = [[] for _ in range(n)]; cnb = [[] for _ in range(m)]
    edges = []
    # reuse PEG logic by building incrementally with nb.peg-like BFS
    for j in range(n):
        for kk in range(degs[j]):
            deg = np.array([len(c) for c in cnb])
            reach = set(vnb[j]); frontier = set(vnb[j]); dist_ok = None
            while frontier:
                nv = {v for c in frontier for v in cnb[c]}
                nc = {c for v in nv for c in vnb[v]} - reach
                if not nc:
                    break
                if len(reach | nc) == m:
                    dist_ok = [c for c in nc if c not in vnb[j]]
                    break
                reach |= nc; frontier = nc
            unreached = [c for c in range(m) if c not in reach and c not in vnb[j]]
            pool = unreached if unreached else (dist_ok or [c for c in range(m) if c not in vnb[j]])
            pool = np.array(pool); cand = pool[deg[pool] == deg[pool].min()]
            c = int(rng.choice(cand)); vnb[j].append(c); cnb[c].append(j); edges.append((c, j))
    coef = np.random.default_rng(1000 + seed).integers(1, 512, len(edges))
    return nb.NBCode(n, m, edges, coef)


if __name__ == '__main__':
    n3s = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '0').split(',')]
    seeds = range(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
    rows = []
    for n3 in n3s:
        for s in seeds:
            t0 = time.time()
            try:
                code = mixed_code(43, 29, n3, s)
            except ValueError as e:
                print(n3, s, e); continue
            fer, ser = evaluate(code)
            rows.append(dict(n3=n3, seed=s, fer=fer, ser=ser))
            print(f'dv3 symbols {n3:2d} seed {s:2d}: FER {fer:.4f} (raw SER {ser:.3f}) {time.time() - t0:.0f}s', flush=True)
    (OUT/'nb_design.json').write_text(json.dumps(rows, indent=1))
