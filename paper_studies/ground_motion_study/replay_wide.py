"""Wider-detector replays for the joint (iterative) receiver.

Re-runs packets of the frozen joint design with a larger detector: increased tone history (state
merging over more past symbols), doubled lookahead and/or a wider beam. Payload and noise are
regenerated from the same seeds, so a replay sees exactly the noise of the original packet.
turbo_packets.py is not modified (its hash is frozen by qualify.py); the detector class used inside
turbo_packets.run is replaced by a subclass with the overrides.
Usage: python replay_wide.py --case 24=0.70 --index 8002904 [--index ...] [--history 5] [--lookahead 24]
       [--beam 32768] [--tag h5l24]  -> results/replay_<tag>.json
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, json, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import mitigation as mg
import turbo_packets as tp
import hop_receiver
import qualify as q

RES = Path(__file__).resolve().parent/'results'


def run(args):
    fc, amp, index, history, lookahead, beam = args
    base = hop_receiver.HopLimitedReceiver

    class Wide(base):
        def __init__(self, *a, **k):
            k.update(history_symbols=history, lookahead_symbols=lookahead, beam_width=beam)
            super().__init__(*a, **k)

    tp.HopLimitedReceiver = Wide
    t0 = time.perf_counter()
    res = tp.run(('hop', fc, amp*mg.A0, index, q.ITER, q.SEEDS, q.CODE, q.LABEL))
    res.update(history=history, lookahead=lookahead, beam=beam, carrier=fc, amplitude_A0=amp,
               wall_s=time.perf_counter() - t0)
    return res


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='carrier=amplitude:index[,index...]')
    p.add_argument('--history', type=int, default=5)
    p.add_argument('--lookahead', type=int, default=24)
    p.add_argument('--beam', type=int, default=32768)
    p.add_argument('--tag', default=None)
    p.add_argument('--workers', type=int, default=3)
    a = p.parse_args()
    jobs = []
    for c in a.case:
        head, idx = c.split(':')
        fc, amp = (float(x) for x in head.split('='))
        jobs += [(fc, amp, int(i), a.history, a.lookahead, a.beam) for i in idx.split(',')]
    tag = a.tag or f'h{a.history}l{a.lookahead}b{a.beam}'
    out = RES/f'replay_{tag}.json'
    old = json.loads(out.read_text()) if out.exists() else []
    with ProcessPoolExecutor(a.workers) as pool:
        for r in pool.map(run, jobs):
            old.append(r)
            print(r['carrier'], r['amplitude_A0'], r['index'], 'errors', r['turbo_errors'][-1],
                  'accepted', r['turbo_accepted'][-1], 'pruned', r['pruned_states'], f"{r['wall_s']:.0f}s", flush=True)
            out.write_text(json.dumps(old, indent=1) + '\n')


if __name__ == '__main__':
    main()
