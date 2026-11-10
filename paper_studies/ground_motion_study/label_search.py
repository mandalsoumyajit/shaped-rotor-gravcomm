"""Labeling design on the actual detector graphs (joint with the iterative receiver).

For simulated packets of a waveform ('cf' or 'hop') at a given amplitude, the exact word posteriors
of the detector graph (turbo_rx.word_posteriors, no priors) are stored. The labeling (word index ->
coded-bit pattern) is then optimized by pairwise swaps to maximize the bit-level mutual information
of the detector output without priors (the EXIT starting point, which limits convergence of the
iterative receiver, cf. exit_chart.py). Training and validation use different packets.
Output: results/labels_<waveform>.json
Usage: python label_search.py hop 24 0.55
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import json, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver, sampled_cffsk
from gravcomm.streaming_soft import StreamingSoftReceiver
import hop_packets as hp
from hop_receiver import HopLimitedReceiver, sampled_hop
import turbo_rx as tr

RES = Path(__file__).resolve().parent/'results'


def job(args):
    waveform, carrier, amplitude, index = args
    rng = np.random.default_rng(np.random.SeedSequence([20371109, index]))
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    if waveform == 'cf':
        ng, L, nw = 48, 3, 256
        words = rng.integers(0, 256, ng)
        symbols = np.column_stack((words//49, words//7 % 7, words % 7)).ravel() - 3
        ts = 224/(3*ng + 18)
        tx = sampled_cffsk(np.r_[symbols, np.zeros(6, int)], ts, q=.08, transition_fraction=.6)
        nsym = len(symbols) + 6
        det = StreamingSoftReceiver(r, ts, amplitude, ng, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                                    q=.08, transition_fraction=.6)
    else:
        ng, L, nw = 43, hp.L, 512
        words = rng.integers(0, 512, ng)
        tones, prev = [], hp.M//2
        for g in words:
            s = hp.CODE.encode(int(g), prev); tones += s; prev = s[-1]
        ts = 224/(len(tones) + 18); txs = np.r_[tones, np.full(6, hp.M//2)]
        tx = sampled_hop(txs, ts, hp.H, float(hp.BETA), hp.M); nsym = len(txs)
        det = HopLimitedReceiver(r, ts, amplitude, ng, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                                 h=hp.H, beta=hp.BETA, tone_count=hp.M, max_hop=hp.D, group_symbols=L,
                                 group_bits=hp.BITS, pad_bits=0)
    mean = r.stream(1).process(tx)
    noise_rng = np.random.default_rng(np.random.SeedSequence([20371110, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    tr.record_graph(det); position = 0
    for s in range(nsym):
        c = det._geometry(s)[3]; det.push(observed[position:position + c]); position += c
    det.finish()
    return tr.word_posteriors(tr.pack_graph(det, L, ng), ng, nw).astype(np.float32), words


def gmi_of(posts, wtrue, label, nbits):
    lb = tr.bit_labels(len(label), nbits)[label].astype(float)
    p = np.exp(posts.astype(float))
    p1 = np.clip(p @ lb, 1e-15, 1 - 1e-15)
    t = lb[wtrue] > 0
    return float(np.mean(1 + np.log2(np.where(t, p1, 1 - p1))))


def search(posts, wtrue, label, nbits, iters, seed=0):
    rng = np.random.default_rng(seed); best = gmi_of(posts, wtrue, label, nbits); n = len(label)
    for _ in range(iters):
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        label[[i, j]] = label[[j, i]]
        g = gmi_of(posts, wtrue, label, nbits)
        if g > best:
            best = g
        else:
            label[[i, j]] = label[[j, i]]
    return label, best


def bsa(posts, wtrue, label, nbits, iters, seed=0):
    """Binary switching on the expected Hamming cost of word confusions (Zeger-Gersho style).

    C[w, w'] = mean posterior of w' given true word w; cost = sum C * HammingDistance(label w, label w').
    Aggregated statistics make this far less prone to overfitting than direct GMI swaps.
    """
    n = len(label)
    C = np.zeros((n, n)); cnt = np.zeros(n)
    np.add.at(C, wtrue, np.exp(posts.astype(float))); np.add.at(cnt, wtrue, 1)
    C /= np.maximum(cnt, 1)[:, None]; C = C + C.T; np.fill_diagonal(C, 0)
    pop = np.array([bin(i).count('1') for i in range(n)])
    cost = lambda lab: float(np.sum(C*pop[lab[:, None] ^ lab[None, :]]))
    rng = np.random.default_rng(seed); best = cost(label)
    for _ in range(iters):
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        label[[i, j]] = label[[j, i]]
        c = cost(label)
        if c < best:
            best = c
        else:
            label[[i, j]] = label[[j, i]]
    return label, best


def main():
    waveform, carrier, amp = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    npk = int(sys.argv[4]) if len(sys.argv) > 4 else 136
    t0 = time.time()
    cache = RES/f'posteriors_{waveform}_f{carrier:g}_a{amp:.3f}_{npk}.npz'
    if cache.exists():
        z = np.load(cache); posts, words = z['posts'], z['words']
    else:
        with ProcessPoolExecutor(17) as pool:
            res = list(pool.map(job, [(waveform, carrier, amp*mg.A0, i) for i in range(npk)]))
        posts = np.concatenate([r[0] for r in res]); words = np.concatenate([r[1] for r in res])
        np.savez_compressed(cache, posts=posts, words=words)
    ntr = len(posts)*2//3; nbits = 8 if waveform == 'cf' else 9; nw = 2**nbits
    print(f'posteriors for {len(posts)} words ({time.time() - t0:.0f}s)', flush=True)
    starts = {'natural': np.arange(nw), 'gray of index': np.arange(nw) ^ (np.arange(nw) >> 1)}
    if waveform == 'cf':
        p = RES/'label_opt.json'
        if p.exists():
            starts['BICM-optimized'] = np.array(json.loads(p.read_text())['optimized']['label'])
    out = dict(waveform=waveform, carrier=carrier, amp=amp, candidates={})
    for name, lab in starts.items():
        out['candidates'][name] = dict(train=gmi_of(posts[:ntr], words[:ntr], lab, nbits),
                                       validate=gmi_of(posts[ntr:], words[ntr:], lab, nbits))
        print(name, out['candidates'][name], flush=True)
    best_start = max(starts, key=lambda n: out['candidates'][n]['train'])
    results = {}
    lab_b, _c = bsa(posts[:ntr], words[:ntr], starts[best_start].copy(), nbits, iters=30000)
    results['bsa'] = dict(start=best_start, train=gmi_of(posts[:ntr], words[:ntr], lab_b, nbits),
                          validate=gmi_of(posts[ntr:], words[ntr:], lab_b, nbits), label=lab_b.tolist())
    print('bsa', {k: v for k, v in results['bsa'].items() if k != 'label'}, f'{time.time() - t0:.0f}s', flush=True)
    lab_g, g = search(posts[:ntr], words[:ntr], lab_b.copy(), nbits, iters=20000)
    results['bsa+gmi'] = dict(start='bsa', train=g, validate=gmi_of(posts[ntr:], words[ntr:], lab_g, nbits),
                              label=lab_g.tolist())
    print('bsa+gmi', {k: v for k, v in results['bsa+gmi'].items() if k != 'label'}, f'{time.time() - t0:.0f}s', flush=True)
    best = max(results, key=lambda k: results[k]['validate'])
    out['all'] = {k: {kk: vv for kk, vv in v.items() if kk != 'label'} for k, v in results.items()}
    out['optimized'] = dict(results[best], method=best)
    (RES/f'labels_{waveform}.json').write_text(json.dumps(out) + '\n')


if __name__ == '__main__':
    main()
