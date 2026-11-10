"""EXIT curves of the hop-limited detector (joint design) for natural, Gray-of-index and designed labelings.

For simulated packets (fresh seeds, not those used to design the labeling), the detector graph is
recorded once; for each labeling and a-priori mutual information I_A, Gaussian a-priori LLRs are
applied to the 9 labeled bits of each word and the extrinsic output of the graph forward-backward
(turbo_rx.graph_app) is measured (exit_chart.mi_llr). Writes results/exit_hop.json.
Usage: python exit_hop.py 24=0.60,0.65,0.70 18=0.72,0.78,0.84
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import json, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
import hop_packets as hp
from hop_receiver import HopLimitedReceiver, sampled_hop
import turbo_rx as tr
import exit_chart as ec

RES = Path(__file__).resolve().parent/'results'
NB = hp.BITS


def labelings():
    n = 2**NB
    return {'natural': np.arange(n), 'gray of index': np.arange(n) ^ (np.arange(n) >> 1),
            'designed': np.array(json.loads((RES/'labels_hop.json').read_text())['optimized']['label'])}


def job(args):
    carrier, amplitude, index = args
    rng = np.random.default_rng(np.random.SeedSequence([20371115, index]))
    ng = 43; words = rng.integers(0, 2**NB, ng)
    tones, prev = [], hp.M//2
    for g in words:
        s = hp.CODE.encode(int(g), prev); tones += s; prev = s[-1]
    ts = 224/(len(tones) + 18); txs = np.r_[tones, np.full(6, hp.M//2)]
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    mean = r.stream(1).process(sampled_hop(txs, ts, hp.H, float(hp.BETA), hp.M))
    noise_rng = np.random.default_rng(np.random.SeedSequence([20371116, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    det = HopLimitedReceiver(r, ts, amplitude, ng, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                             h=hp.H, beta=hp.BETA, tone_count=hp.M, max_hop=hp.D, group_symbols=hp.L,
                             group_bits=NB, pad_bits=0)
    tr.record_graph(det); position = 0
    for s in range(len(txs)):
        c = det._geometry(s)[3]; det.push(observed[position:position + c]); position += c
    det.finish()
    pk = tr.pack_graph(det, hp.L, ng)
    res = {}
    for name, lab in labelings().items():
        labels = tr.bit_labels(2**NB, NB)[lab]
        bits = labels[words].astype(np.uint8)
        row = []
        for I in ec.IA:
            La = ec.gaussian_prior(bits, I, np.random.default_rng(index*100 + int(I*100)))
            row.append((tr.graph_app(pk, La, labels) - La).ravel().tolist())
        res[name] = (bits.ravel().tolist(), row)
    return res


def curves(carrier, amp, packets=48, workers=17):
    with ProcessPoolExecutor(workers) as pool:
        outs = list(pool.map(job, [(carrier, amp*mg.A0, i) for i in range(packets)]))
    out = {}
    for name in outs[0]:
        bits = np.concatenate([np.array(o[name][0]) for o in outs])
        out[name] = [ec.mi_llr(np.concatenate([np.array(o[name][1][k]) for o in outs]), bits) for k in range(len(ec.IA))]
    return out


def main():
    t0 = time.time(); target = RES/'exit_hop.json'
    out = json.loads(target.read_text()) if target.exists() else dict(IA=ec.IA.tolist(), curves={})
    for arg in sys.argv[1:]:
        fc, amps = arg.split('=')
        for a in amps.split(','):
            key = f'{float(fc):g}_{float(a):.3f}'
            if key in out['curves']:
                continue
            out['curves'][key] = curves(float(fc), float(a))
            print(key, {k: np.round(v, 3).tolist() for k, v in out['curves'][key].items()}, f'{time.time() - t0:.0f}s', flush=True)
            target.write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
