"""Wider-detector replays for the NB-LDPC GF(512) joint receiver.

Re-runs qualification packets of the frozen NB design (qualify_nb.py settings and seeds) with a larger
detector: tone history 5 (state merging over more past symbols), doubled lookahead and a beam wide
enough that no states are pruned. Payload and noise are regenerated from the same seeds, so a replay
sees exactly the noise of the original packet. nb_packets.py is not modified (its hash is frozen by
qualify_nb.py); the detector class used inside nb_packets.run is replaced by a subclass. CPU only.
Usage: python replay_wide_nb.py --case 24=0.58:8702345,8703688 [--history 5] [--lookahead 24]
       [--beam 262144] [--tag nb_h5l24b262k] -> results/replay_<tag>.json
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, json, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import mitigation as mg
import nb_packets as npk
import hop_receiver
import qualify_nb as q

RES = Path(__file__).resolve().parent/'results'
SEEDS = (20371521, 20371522)                                 # qualification attempt 3


def run(args):
    fc, amp, index, history, lookahead, beam, seeds = args
    base = hop_receiver.HopLimitedReceiver

    class Wide(base):
        def __init__(self, *a, **k):
            k.update(history_symbols=history, lookahead_symbols=lookahead, beam_width=beam)
            super().__init__(*a, **k)

    npk.HopLimitedReceiver = Wide
    t0 = time.perf_counter()
    res = npk.run((fc, amp*mg.A0, index, q.ITER, seeds, q.CODE_SEED, q.BP, q.DV, q.WORDS, 'cpu'))
    res.update(seeds=list(seeds), history=history, lookahead=lookahead, beam=beam, carrier=fc, amplitude_A0=amp,
               wall_s=time.perf_counter() - t0)
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='carrier=amplitude:index[,index...]')
    p.add_argument('--history', type=int, default=5)
    p.add_argument('--lookahead', type=int, default=24)
    p.add_argument('--beam', type=int, default=262144)
    p.add_argument('--tag', default=None)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--seeds', default=None, help='payload,noise namespaces (default: qualification attempt 3)')
    a = p.parse_args()
    seeds = tuple(int(x) for x in a.seeds.split(',')) if a.seeds else SEEDS
    jobs = []
    for c in a.case:
        head, idx = c.split(':')
        fc, amp = (float(x) for x in head.split('='))
        jobs += [(fc, amp, int(i), a.history, a.lookahead, a.beam, seeds) for i in idx.split(',')]
    tag = a.tag or f'nb_h{a.history}l{a.lookahead}b{a.beam}'
    out = RES/f'replay_{tag}.json'
    old = json.loads(out.read_text()) if out.exists() else []
    done = {(r['carrier'], r['index']) for r in old}
    jobs = [j for j in jobs if (j[0], j[2]) not in done]
    with ProcessPoolExecutor(a.workers) as pool:
        for r in pool.map(run, jobs):
            old.append(r)
            print(r['carrier'], r['amplitude_A0'], r['index'], 'errors', r['errors'][-1], 'accepted', r['accepted'][-1],
                  'passes', len(r['errors']), 'pruned', r['pruned_states'], f"{r['wall_s']:.0f}s", flush=True)
            out.write_text(json.dumps(old, indent=1) + '\n')


if __name__ == '__main__':
    main()
