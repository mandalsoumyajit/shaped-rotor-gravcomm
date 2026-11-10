"""Synchronization loss from paired packet runs (same packets, same noise; acq_packets.py).

For each carrier the receiver with true carrier phase, gain and timing ('known') is compared with the
steady-state tracking receiver ('track') and the cold-start preamble-only receiver ('acq') on identical packet
indices. The frame-error rates are converted into an equivalent amplitude loss with the probit waterfall fit
of the known-parameter receiver (nb_waterfall_fit.json): FER = Phi(-(x - x0)/s), x = 20 log10(A/A0), so a
receiver whose FER is F_r has loss s*(Phi^-1(F_r) - Phi^-1(F_known)) dB. Intervals: paired bootstrap over packets
(the probit slope is held fixed; it is a local conversion, most reliable for small losses).
Output: results/acquisition_loss.json
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import ndtri

RES = Path(__file__).resolve().parent/'results'
CASES = {24: ('0.525', 'known_f24_a0.525', 'track_f24_a0.525', 'acq_f24_a0.525'),
         18: ('0.63', 'known_f18_a0.63', 'track_f18_a0.63', 'acq_f18_a0.63')}


def load(name):
    out = {}
    for p in (RES/name).glob('9*.json'):
        r = json.loads(p.read_text())
        out[r['index']] = (not r['accepted'][-1], int(r['errors'][-1]))
    return out


def main(boot=4000, seed=20371720):
    fit = json.loads((RES/'nb_waterfall_fit.json').read_text())
    rng = np.random.default_rng(seed)
    out = {}
    for fc, (amp, kn, tr, ac) in CASES.items():
        s = fit[str(fc)]['probit']['s_db']
        K = load(kn)
        res = dict(amplitude_A0=float(amp), probit_s_db=s, known={'packets': len(K),
                   'failures': sum(v[0] for v in K.values())})
        for label, name in (('track', tr), ('acq', ac)):
            R = load(name)
            idx = sorted(set(K) & set(R))
            if not idx:
                continue
            fk = np.array([K[i][0] for i in idx]); fr = np.array([R[i][0] for i in idx])
            n = len(idx)
            def loss(a, b):
                pa, pb = (a.sum() + 0.5)/(n + 1), (b.sum() + 0.5)/(n + 1)
                return s*(ndtri(pb) - ndtri(pa))
            bs = []
            for _ in range(boot):
                j = rng.integers(0, n, n)
                bs.append(loss(fk[j], fr[j]))
            res[label] = dict(packets=n, failures_known=int(fk.sum()), failures=int(fr.sum()),
                              fer_known=float(fk.mean()), fer=float(fr.mean()),
                              only_receiver=int((fr & ~fk).sum()), only_known=int((fk & ~fr).sum()),
                              loss_db=float(loss(fk, fr)),
                              loss_db_ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])
            print(fc, label, json.dumps(res[label]), flush=True)
        out[str(fc)] = res
    # synchronization-error statistics of the qualification runs (RMS over packets)
    out['sync_errors'] = {}
    for name in ('acqqual_f24_a0.5800', 'acqqual_f18_a0.6970', 'trackqual_f24_a0.5800', 'trackqual_f18_a0.6970'):
        rows = [json.loads(p.read_text()) for p in (RES/name).glob('9*.json')]
        if not rows:
            continue
        e = {k: np.array([r[k] for r in rows]) for k in ('phase_err_deg', 'timing_err_s', 'amp_err')}
        out['sync_errors'][name] = dict(packets=len(rows),
                                        rejected=sum(not r['accepted'][-1] for r in rows),
                                        **{f'rms_{k}': float(np.sqrt(np.mean(v**2))) for k, v in e.items()},
                                        **{f'max_abs_{k}': float(np.max(np.abs(v))) for k, v in e.items()})
        print(name, json.dumps(out['sync_errors'][name]), flush=True)
    (RES/'acquisition_loss.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
