"""Tally packet screens and estimate the accuracy of equal-whitened-SNR scaling.

For each mitigated case run at two amplitudes, log BER is interpolated linearly in
amplitude to find the amplitude that reproduces the reference (released-isolation)
BER; its ratio to the predicted amplitude measures the scaling error.
"""
import json
from pathlib import Path
import numpy as np

RES = Path(__file__).resolve().parent/'results'


def tally(folder):
    rows = [json.loads(p.read_text()) for p in folder.glob('*.json')]
    n = len(rows); e = sum(r['errors'] for r in rows); f = sum(r['errors'] > 0 or not r['accepted'] for r in rows)
    # Packet-level standard error of BER (errors cluster within packets).
    x = np.array([r['errors']/224 for r in rows])
    return dict(packets=n, errors=e, ber=e/(224*n), ber_se=float(x.std(ddof=1)/np.sqrt(n)),
                packet_failure=f/n, pruned=sum(r['pruned_states'] for r in rows),
                missing=sum(r['missing'] for r in rows), max_whitening=max(r['whitening'] for r in rows))


def main():
    cases = {p.name: tally(p) for p in sorted(RES.glob('*_f*_a*')) if any(p.glob('*.json'))}
    out = dict(cases=cases, equivalence={})
    refs = {'24': 'baseline_f24_a0.6500', '18': 'baseline_f18_a0.8000'}
    pred = {'24': 0.5591, '18': 0.5647}
    for fc, ref in refs.items():
        pts = sorted((float(k.split('_a')[1]), v['ber']) for k, v in cases.items() if k.startswith(f'structural_f{fc}_'))
        if len(pts) == 2 and ref in cases:
            (a1, b1), (a2, b2) = pts
            target = np.log(cases[ref]['ber'])
            a = a1 + (target - np.log(b1))*(a2 - a1)/(np.log(b2) - np.log(b1))
            out['equivalence'][f'structural@{fc}'] = dict(predicted=pred[fc], matched=float(a), ratio=float(a/pred[fc]))
    # Operating-point checks: predicted amplitude, compared with the qualified released channels.
    import mitigation  # noqa: F401  (puts the release on sys.path)
    from gravcomm.packet_confidence import packet_ber_interval
    reference = {'24': dict(ber=3.69e-4, interval=[2.06e-4, 9.93e-4]), '18': dict(ber=3.50e-4, interval=[1.68e-4, 9.93e-4])}
    out['operating_point'] = {}
    for fc in ('24', '18'):
        for p in RES.glob(f'three_stage_f{fc}_a*'):
            rows = [json.loads(q.read_text()) for q in p.glob('*.json')]
            if not rows:
                continue
            lo, hi = packet_ber_interval([r['errors'] for r in rows], 224, .05)
            t = tally(p)
            out['operating_point'][p.name] = dict(t, interval_95=[lo, hi], reference=reference[fc],
                                                 accepted=sum(r['accepted'] for r in rows)/len(rows),
                                                 wrong_accepted=sum(r['accepted'] and r['errors'] > 0 for r in rows))
    (RES/'validation_summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out['operating_point'], indent=1))
    for k, v in cases.items():
        print(f"{k:28s} n={v['packets']:4d} BER={v['ber']:.4f}±{v['ber_se']:.4f} fail={v['packet_failure']:.3f} pruned={v['pruned']} missing={v['missing']} wh={v['max_whitening']:.1e}")
    print(json.dumps(out['equivalence'], indent=1))


if __name__ == '__main__':
    main()
