"""Optimize the bit labeling of the released CF-FSK words for bit-interleaved decoding (BICM).

Exact word posteriors from word_bcjr.run (uniform prior over the 256 words used by the released
mapping) give, for any labeling pi (word -> 8-bit label), exact bit LLRs and hence the
generalized mutual information GMI(pi) per coded bit, without new channel simulations.
pi is improved by a pairwise-swap local search (binary switching, as in labeling design for
BICM; cf. Perotti et al. 2010 for CPM) on a training realization and checked on an
independent realization. The word set (which 256 of the 343 tone triples) is unchanged.
Output: results/label_opt.json (labels and GMI).
"""
import json
from pathlib import Path
import numpy as np
import word_bcjr as wb

OUT = Path(__file__).resolve().parent/'results'/'label_opt.json'
BITS = ((np.arange(256)[:, None] >> np.arange(7, -1, -1)) & 1).astype(bool)


def gmi_of(posts, wtrue, label):
    """GMI per bit for labeling label[word] (exact LLRs from word posteriors, unit scale)."""
    lb = BITS[label]                                     # (words, 8) bit pattern of each word
    p = np.exp(posts)                                    # normalized posteriors
    p1 = p @ lb.astype(float)                            # P(bit = 1) per word step and bit
    p1 = np.clip(p1, 1e-15, 1 - 1e-15)
    t = lb[wtrue]
    q = np.where(t, p1, 1 - p1)                          # posterior of the true bit
    return float(np.mean(1 + np.log2(q)))


def gray(n):
    return n ^ (n >> 1)


def search(posts, wtrue, label, iters=20000, seed=0):
    rng = np.random.default_rng(seed); best = gmi_of(posts, wtrue, label)
    for it in range(iters):
        i, j = rng.integers(0, 256, 2)
        if i == j:
            continue
        label[[i, j]] = label[[j, i]]
        g = gmi_of(posts, wtrue, label)
        if g > best:
            best = g
        else:
            label[[i, j]] = label[[j, i]]
    return label, best


def main():
    fc, amp = 24., 0.48
    tr = wb.run(fc, amp, n_words=3000, seed=21, return_post=True)
    va = wb.run(fc, amp, n_words=3000, seed=22, return_post=True)
    cands = dict(natural=np.arange(256), gray_of_rank=np.array([gray(w) for w in range(256)]))
    out = dict(fc=fc, amp=amp, candidates={})
    for name, lab in cands.items():
        out['candidates'][name] = dict(train=gmi_of(*tr, lab), validate=gmi_of(*va, lab))
        print(name, out['candidates'][name], flush=True)
    best_name = max(cands, key=lambda n: out['candidates'][n]['train'])
    lab, g = search(*tr, cands[best_name].copy())
    out['optimized'] = dict(start=best_name, train=g, validate=gmi_of(*va, lab), label=lab.tolist())
    print('optimized', {k: v for k, v in out['optimized'].items() if k != 'label'}, flush=True)
    OUT.write_text(json.dumps(out) + '\n')


if __name__ == '__main__':
    main()
