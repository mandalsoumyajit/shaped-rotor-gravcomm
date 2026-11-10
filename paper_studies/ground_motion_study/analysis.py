"""Equal-whitened-SNR sizing of the mitigation scenarios (writes results/analysis.json).

Received (differential) A2 is tabulated on a grid of uniform scales with the
finite-volume rotor and interpolated in log-log coordinates for inversion.
"""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import mitigation as mg

OUT = Path(__file__).resolve().parent/'results'
NAME = 'analysis_hop.json' if __import__('os').environ.get('GC_DESIGN') == 'hop' else 'analysis.json'
CARRIERS = (10., 12., 15., 18., 24.)
INSTRUMENT = (1., .3, .1)
GRADIOMETER = [('radial', b) for b in (.1, .2, .3, .5)] + [('transverse', b) for b in (.1, .2, .3, .5, .8)]
_rotor = None


def _amp(args):
    global _rotor
    if _rotor is None:
        _rotor = mg.ScaledRotor()
    s, d, b, layout = args
    if layout == 'single':
        return abs(_rotor.coefficient(s, (d, 0, 0), samples=64))
    if layout == 'radial':
        return abs(_rotor.coefficient(s, (d, 0, 0), samples=64) - _rotor.coefficient(s, (d + b, 0, 0), samples=64))
    return abs(_rotor.coefficient(s, (d, b/2, 0), samples=64) - _rotor.coefficient(s, (d, -b/2, 0), samples=64))


def grid(pool, d, layout='single', b=0.):
    smax = (d - .0127)/mg.R0*.999
    s = np.geomspace(.3, smax, 14)
    a = np.array(list(pool.map(_amp, [(x, d, b, layout) for x in s])))
    return s, a


def invert(s, a, target):
    if target > a[-1]:
        return None
    return float(np.exp(np.interp(np.log(target), np.log(a), np.log(s))))


def main():
    OUT.mkdir(exist_ok=True)
    out = dict(noise={}, amplitude={}, sizing={}, gradiometer={}, grids={})
    for name in mg.SCENARIOS:
        for fc in CARRIERS:
            c = mg.background(name, fc).components(np.array([fc]))
            out['noise'][f'{name}@{fc:g}'] = dict(support_ngal=float(c['support'][0]/1e-11),
                                                  newtonian_ngal=float(c['newtonian'][0]/1e-11))
    # Joint (hop) design: amplitude factors are ratios of information-rate thresholds (rate_scaling.py);
    # packet validation showed the equal-whitened-SNR factor to be optimistic for the joint receiver.
    # The SNR factor is kept in 'amplitude_mf' for comparison. Released design: SNR factor (as published).
    rates = json.loads((OUT/'rate_scaling.json').read_text()) if mg.DESIGNS is mg.HOP_DESIGNS else None
    out['amplitude_mf'] = {}
    out['scaling'] = 'information rate' if rates else 'matched-filter SNR'
    for d, (fref, aref, *_) in mg.DESIGNS.items():
        ref = mg.background('baseline', fref)
        for name in mg.SCENARIOS:
            for fc in CARRIERS:
                for k in INSTRUMENT:
                    key = f'{d:g}m|{name}|{fc:g}Hz|inst{k:g}'
                    mf = mg.equivalent_amplitude_factor(mg.background(name, fc), ref, k)
                    out['amplitude_mf'][key] = mf
                    if rates:
                        out['amplitude'][key] = (rates[f'{name}|{fc:g}|{k:g}']['threshold_A0']
                                                 / rates[f'baseline|{fref:g}|1']['threshold_A0'])
                    else:
                        out['amplitude'][key] = mf
    # Field grids depend only on geometry; reuse them from a previous run when present.
    cached = json.loads((OUT/NAME).read_text()) if (OUT/NAME).exists() else {}
    with ProcessPoolExecutor(max_workers=6) as pool:
        for d, (fref, aref, s0, *_) in mg.DESIGNS.items():
            if f'{d:g}m|gradiometer' in cached.get('grids', {}):
                out['grids'][f'{d:g}m|single'] = cached['grids'][f'{d:g}m|single']
                out['grids'][f'{d:g}m|gradiometer'] = cached['grids'][f'{d:g}m|gradiometer']
                out['gradiometer'].update({k: v for k, v in cached['gradiometer'].items() if k.startswith(f'{d:g}m|')})
            else:
                s, a = grid(pool, d)
                out['grids'][f'{d:g}m|single'] = dict(scale=s.tolist(), a2=a.tolist())
                single = _amp((s0, d, 0, 'single'))
                ratios = list(pool.map(_amp, [(s0, d, b, lay) for lay, b in GRADIOMETER]))
                for (lay, b), r in zip(GRADIOMETER, ratios):
                    out['gradiometer'][f'{d:g}m|{lay}|{b:g}'] = r/single
                best = max(zip(ratios, GRADIOMETER))[1]
                sg, ag = grid(pool, d, best[0], best[1])
                out['grids'][f'{d:g}m|gradiometer'] = dict(scale=sg.tolist(), a2=ag.tolist(), layout=best[0], baseline_m=best[1])
            s, a = (np.array(out['grids'][f'{d:g}m|single'][k]) for k in ('scale', 'a2'))
            sg, ag = (np.array(out['grids'][f'{d:g}m|gradiometer'][k]) for k in ('scale', 'a2'))
            for name in mg.SCENARIOS:
                ss, aa = (sg, ag) if name == 'gradiometer' else (s, a)
                for fc in CARRIERS:
                    for k in INSTRUMENT:
                        key = f'{d:g}m|{name}|{fc:g}Hz|inst{k:g}'
                        x = invert(ss, aa, aref*out['amplitude'][key])
                        out['sizing'][key] = None if x is None else mg.sized(d, x, fc)
    (OUT/NAME).write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
