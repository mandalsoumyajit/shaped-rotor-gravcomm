"""Split the coding gap of the released link into block-length, mapping/detector, and code losses.

Inputs: results/llr_capture/*.npz (llr_capture.py) and results/finite_length.json.
For the detector's actual soft outputs (LLR = log P(b=0)/P(b=1), max-log) the generalized mutual
information per coded bit is
    GMI = max_s E[1 - log2(1 + exp(-s (1-2b) LLR))],
the rate an ideal binary code with interleaving achieves on these metrics (bit-interleaved
coded modulation, no iterations). The required GMI is 256/384 = 2/3. Its finite-length version
uses the normal approximation with n = 384 coded bits and the variance of the bit information
density (bits treated as independent after interleaving). Amplitude thresholds are interpolated
in log amplitude. Ladder (amplitude, A0 units):
  A_cm      symbol-level information rate of the waveform, ideal code, infinite length
  A_cm_na   same, 256-bit frame (normal approximation)
  A_gmi     detector soft outputs + ideal binary code, infinite length
  A_gmi_na  same, 384 coded bits
  A_op      released convolutional code (qualified operating point)
"""
import glob, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import norm

RES = Path(__file__).resolve().parent/'results'
N_BITS, K_BITS = 384, 256


def gmi(llr, bits):
    x = np.where(np.isfinite(llr), llr, 0.0)*(1 - 2*bits)
    f = lambda s: -np.mean(1 - np.logaddexp(0, -s*x)/np.log(2))
    s = minimize_scalar(f, bounds=(1e-3, 10), method='bounded').x
    dens = 1 - np.logaddexp(0, -s*x)/np.log(2)
    return float(dens.mean()), float(dens.var()), float(s)


def na_bits(g, v, eps=1e-2, n=N_BITS):
    return n*g - np.sqrt(n*v)*norm.isf(eps) + 0.5*np.log2(n)


def crossing(amps, vals, target):
    la = np.log(amps)
    for i in range(len(vals) - 1):
        if (vals[i] - target)*(vals[i + 1] - target) <= 0 and vals[i] != vals[i + 1]:
            return float(np.exp(la[i] + (target - vals[i])*(la[i + 1] - la[i])/(vals[i + 1] - vals[i])))
    return None


def main():
    fl = json.loads((RES/'finite_length.json').read_text())
    out = {}
    for fc, a_op in (('24', .8125), ('18', 1.0)):
        rows = []
        for path in sorted(glob.glob(str(RES/'llr_capture'/f'f{fc}_a*.npz'))):
            z = np.load(path)
            g, v, s = gmi(z['llr'].ravel(), z['word'].ravel())
            rows.append(dict(amp=float(z['amp_A0']), factor=float(z['amp_factor_of_op']), gmi=g, var=v, scale=s,
                             na_bits=float(na_bits(g, v)), coded_ber=float(np.mean((z['llr'] < 0) != z['word'])),
                             payload_ber=float(z['errors'].sum()/(224*len(z['errors']))),
                             pruned=int(z['pruned'].sum()), missing=int(z['missing'].sum())))
        rows.sort(key=lambda r: r['amp'])
        amps = np.array([r['amp'] for r in rows])
        ladder = dict(A_cm=fl[fc]['threshold'], A_cm_na=fl[fc]['normal_approx_eps0.01'],
                      A_gmi=crossing(amps, [r['gmi'] for r in rows], K_BITS/N_BITS),
                      A_gmi_na=crossing(amps, [r['na_bits'] for r in rows], K_BITS), A_op=a_op)
        db = lambda a, b: None if a is None or b is None else 20*np.log10(a/b)
        split = dict(block_length_symbol_level=db(ladder['A_cm_na'], ladder['A_cm']),
                     mapping_and_detector=db(ladder['A_gmi'], ladder['A_cm']),
                     block_length_bit_level=db(ladder['A_gmi_na'], ladder['A_gmi']),
                     convolutional_code=db(ladder['A_op'], ladder['A_gmi_na']),
                     total=db(ladder['A_op'], ladder['A_cm']))
        out[fc] = dict(points=rows, ladder=ladder, split_db=split)
        print(fc, json.dumps(dict(ladder=ladder, split_db=split), indent=1))
        for r in rows:
            print(f"   x{r['factor']:.2f}: GMI {r['gmi']:.3f}  NA bits {r['na_bits']:.0f}  coded BER {r['coded_ber']:.4f}"
                  f"  payload BER {r['payload_ber']:.2e}  scale {r['scale']:.2f}  pruned {r['pruned']} missing {r['missing']}")
    (RES/'gap_analysis.json').write_text(json.dumps(out, indent=1) + '\n')


if __name__ == '__main__':
    main()
