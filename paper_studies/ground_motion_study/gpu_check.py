"""Validate the GPU detector against the CPU detector on identical packets.

For each packet (nb_packets transmitter/channel/seeds) both detectors run on the same samples; the
recorded graphs are compared step by step (edge counts, source/destination/word arrays, edge metrics),
and the full iterative GF(512) receiver is run on both graphs (decoded payload, iterations, symbol
errors). Also reports detection times. Writes results/gpu_check.json.
Usage: python gpu_check.py --carrier 24 --amp 0.50 --indices 7200000-7200015
"""
import argparse, json, time
from pathlib import Path
import numpy as np
import nb_packets as npk

RES = Path(__file__).resolve().parent/'results'


def parse(spec):
    out = []
    for part in spec.split(','):
        a, _, b = part.partition('-')
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return out


def capture(device, carrier, amp, index):
    import mitigation as mg
    return npk.run((carrier, amp*mg.A0, index, 16, (20371103, 20371104), 1, 30, 2, 'lex', device), return_graph=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--carrier', type=float, default=24.0); p.add_argument('--amp', type=float, default=0.50)
    p.add_argument('--indices', default='7200000-7200015')
    a = p.parse_args()
    rows = []
    capture('gpu', a.carrier, a.amp, 7199999)                 # warm-up (CUDA context, kernel compilation)
    for i in parse(a.indices):
        rc, gc = capture('cpu', a.carrier, a.amp, i)
        rg, gg = capture('gpu', a.carrier, a.amp, i)
        same_struct = all(np.array_equal(gc[k], gg[k]) for k in ('src', 'dst', 'word', 'eoff', 'noff', 'group'))
        rel = float(np.max(np.abs(gc['met'] - gg['met'])/np.maximum(np.abs(gc['met']), 1e-12))) if same_struct else None
        same_result = (rc['errors'] == rg['errors'] and rc['accepted'] == rg['accepted']
                       and rc['symbol_errors'] == rg['symbol_errors'])
        row = dict(index=i, graph_identical_structure=same_struct, max_rel_metric_diff=rel,
                   edges=int(len(gc['met'])), edges_gpu=int(len(gg['met'])), result_identical=same_result,
                   cpu_detect_s=rc['detect_s'], gpu_detect_s=rg['detect_s'], cpu_total_s=rc['elapsed_s'],
                   gpu_total_s=rg['elapsed_s'], accepted=rc['accepted'][-1], iterations=len(rc['errors']))
        rows.append(row)
        print(row, flush=True)
    summary = dict(packets=len(rows), structure_identical=sum(r['graph_identical_structure'] for r in rows),
                   results_identical=sum(r['result_identical'] for r in rows),
                   max_rel_metric_diff=max((r['max_rel_metric_diff'] or 0) for r in rows),
                   cpu_detect_s=float(np.median([r['cpu_detect_s'] for r in rows])),
                   gpu_detect_s=float(np.median([r['gpu_detect_s'] for r in rows])))
    print(summary)
    (RES/f'gpu_check_f{a.carrier:g}_a{a.amp:.3f}.json').write_text(json.dumps(dict(summary=summary, rows=rows), indent=1))


if __name__ == '__main__':
    main()
