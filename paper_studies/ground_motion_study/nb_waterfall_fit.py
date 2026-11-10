"""Smooth waterfall of the GF(512) link and threshold / operating-point estimates.

Data: every NB packet folder of the lexicographic-word design (screens nb_f<fc>_a<amp>, seeds
20371103/4, and qualification attempts qualnb_f<fc>_a<amp>, fresh seeds); per amplitude the packet
count N, frame failures F (CRC/syndrome rejected) and payload bit errors E.
Models, in x = 20 log10(A/A0) (dB), fitted by binomial maximum likelihood:
  probit   FER = Q((x - x0)/s)
  cloglog  FER = exp(-exp((x - x0)/s))       (asymmetric alternative; spread = model dependence)
BER(A) = FER(A) * m/224 with m the mean bit errors per failed frame (tested for an amplitude trend).
Outputs: amplitudes at FER 1e-2, 1e-3 and BER 1e-3 with 95% parametric-bootstrap intervals, and the
probability that the qualification rule (time-uniform betting CS, alpha 0.025, checks every 128 packets,
cap 16896) passes at candidate amplitudes, simulated with failures drawn from the fitted FER and
error counts resampled from the observed failed frames. Writes results/nb_waterfall_fit.json.
Usage: python nb_waterfall_fit.py [--sims 300]
"""
import argparse, glob, json, re
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, brentq
import math


def ndtr(z):
    """Standard normal CDF via math.erfc (scipy.special.ndtr crashed after many calls on this
    Python 3.14 / SciPy 1.18 stack)."""
    z = np.asarray(z, dtype=float)
    return np.array([0.5*math.erfc(-v/math.sqrt(2.0)) for v in z.ravel()]).reshape(z.shape)
import mitigation  # noqa: F401  (puts the released gravcomm package on the path)
from gravcomm.packet_confidence import packet_ber_interval

RES = Path(__file__).resolve().parent/'results'
MODELS = {
    'probit': lambda x, x0, s: ndtr(-(x - x0)/s),
    'cloglog': lambda x, x0, s: np.exp(-np.exp((x - x0)/s)),
}
XC = -6.0                                                     # centre (dB) of the quadratic-probit argument


def fer_model(name, x, p):
    """FER at x dB. 'probit2': Q(c0 + c1 z + c2 z^2), z = x - XC (flattening tail allowed)."""
    if name == 'probit2':
        z = x - XC
        return ndtr(-(p[0] + p[1]*z + p[2]*z**2))
    return MODELS[name](x, p[0], abs(p[1]))


def nll_of(name, x, N, F):
    def nll(p):
        q = np.clip(fer_model(name, x, p), 1e-15, 1 - 1e-15)
        if name == 'probit2':                               # monotone over the fitted range plus margin
            zz = np.linspace(x.min() - 1.0, x.max() + 2.0, 31) - XC
            if np.any(p[1] + 2*p[2]*zz <= 0):
                return 1e12
        return -np.sum(F*np.log(q) + (N - F)*np.log(1 - q))
    return nll


def fit_general(name, x, N, F, start=None):
    nll = nll_of(name, x, N, F)
    if start is not None:
        starts = [start]
    elif name == 'probit2':
        starts = [[c0, c1, c2] for c0 in (-2, 0, 2) for c1 in (1, 2, 4) for c2 in (-0.5, 0.0, 0.5)]
    else:
        starts = [[x0, s] for x0 in np.linspace(x.min(), x.max(), 9) for s in (0.1, 0.3, 1.0)]
    best = None
    for s0 in starts:
        r = minimize(nll, s0, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxiter=8000))
        if best is None or r.fun < best.fun:
            best = r
    return np.array(best.x), best.fun


def amp_at_general(name, p, target, lo, hi):
    x = brentq(lambda x: fer_model(name, x, p) - target, lo, hi)
    return 10**(x/20)


def deviance(N, F, exp_f):
    return float(2*np.sum(np.where(F > 0, F*np.log(np.maximum(F, 1e-300)/exp_f), 0) +
                          np.where(N - F > 0, (N - F)*np.log(np.maximum(N - F, 1e-300)/(N - exp_f)), 0)))


DESIGN = 'nb'


def load():
    """Packet rows per (carrier, amplitude). DESIGN 'nb': GF(512) folders; 'k5': K5 joint-receiver folders
    (j_hopk5opt_* screens and qual_hopk5opt_* qualification, fields turbo_errors / turbo_accepted)."""
    rows = {}
    if DESIGN == 'k5':
        pat = re.compile(r'(?:j|qual)_hopk5opt_f(\d+)_a([0-9.]+)$')
        dirs = sorted(glob.glob(str(RES/'j_hopk5opt_f*')) + glob.glob(str(RES/'qual_hopk5opt_f*')))
        ek, ak = 'turbo_errors', 'turbo_accepted'
    else:
        pat = re.compile(r'(?:nb|qualnb)_f(\d+)_a([0-9.]+)$')
        dirs = sorted(glob.glob(str(RES/'nb_f*')) + glob.glob(str(RES/'qualnb_f*')))
        ek, ak = 'errors', 'accepted'
    for d in dirs:
        m = pat.search(Path(d).name)
        if not m:
            continue
        fc, amp = int(m.group(1)), round(float(m.group(2)), 4)
        for f in glob.glob(d + '/[0-9]*.json'):
            r = json.loads(Path(f).read_text())
            if r.get('words', 'lex') != 'lex':
                continue
            rows.setdefault((fc, amp), []).append((not r[ak][-1], r[ek][-1]))
    return rows


def summarize(rows, fc):
    pts = sorted((a, len(v), sum(f for f, _ in v), sum(e for _, e in v), [e for f, e in v if f])
                 for (c, a), v in rows.items() if c == fc)
    return pts


def fit(model, x, N, F):
    f = MODELS[model]
    def nll(p):
        q = np.clip(f(x, p[0], abs(p[1])), 1e-15, 1 - 1e-15)
        return -np.sum(F*np.log(q) + (N - F)*np.log(1 - q))
    best = None
    for x0 in np.linspace(x.min(), x.max(), 9):
        for s in (0.1, 0.3, 1.0):
            r = minimize(nll, [x0, s], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
            if best is None or r.fun < best.fun:
                best = r
    return np.array([best.x[0], abs(best.x[1])]), best.fun


def amp_at(model, p, target):
    f = MODELS[model]
    x = brentq(lambda x: f(x, *p) - target, -20, 5)
    return 10**(x/20)


def pass_probability(model, p, amp, failed_errors, sims, cap=16896, rng=None):
    q = float(MODELS[model](20*np.log10(amp), *p))
    passes, used = 0, []
    for _ in range(sims):
        fails = rng.random(cap) < q
        e = np.zeros(cap, int); e[fails] = rng.choice(failed_errors, fails.sum())
        status = None
        for k in range(128, cap + 1, 128):
            if k < 6144 and e[:k].sum() > 0:
                continue                                       # pass impossible before 6144 packets anyway
            lo, hi = packet_ber_interval(list(e[:k]), 224, 0.025)
            if hi <= 1e-3:
                status = 'pass'; used.append(k); break
            if lo > 1e-3:
                status = 'fail'; break
        passes += status == 'pass'
    return passes/sims, (float(np.median(used)) if used else None), q


def pass_probability_general(name, p, x, failed_errors, sims, rng, cap=16896):
    q = float(fer_model(name, np.array([x]), p)[0])
    passes, used = 0, []
    for _ in range(sims):
        fails = rng.random(cap) < q
        e = np.zeros(cap, int); e[fails] = rng.choice(failed_errors, fails.sum())
        status = None
        for k in range(6144, cap + 1, 128):                   # pass impossible before 6144 packets
            lo, hi = packet_ber_interval(list(e[:k]), 224, 0.025)
            if hi <= 1e-3:
                status = 'pass'; used.append(k); break
            if lo > 1e-3:
                status = 'fail'; break
        passes += status == 'pass'
    return passes/sims, (float(np.median(used)) if used else None), q


def main_v2(a):
    """Per-carrier and pooled fits with the quadratic probit; model-free log-BER interpolation; pass
    probabilities along the pooled quadratic-probit curve. Writes results/nb_waterfall_fit2.json."""
    rows = load(); rng = np.random.default_rng(20371701); out = {}
    data = {}
    for fc in (24, 18):
        pts = summarize(rows, fc)
        data[fc] = dict(amps=np.array([p[0] for p in pts]), N=np.array([p[1] for p in pts]),
                        F=np.array([p[2] for p in pts]), E=np.array([p[3] for p in pts]),
                        failed=np.concatenate([p[4] for p in pts if p[4]]))
    # pooled: 18 Hz amplitudes divided by a fitted carrier scale r (dB offset d = 20 log10 r)
    def pooled_nll(p):
        d, c = p[0], p[1:]
        tot = 0.0
        for fc, off in ((24, 0.0), (18, d)):
            D = data[fc]; x = 20*np.log10(D['amps']) - off
            tot += nll_of('probit2', x, D['N'], D['F'])(c)
        return tot
    best = None
    for d0 in (1.4, 1.6, 1.8):
        for s0 in ([0, 2, 0], [0, 3, 0.3], [1, 4, -0.2]):
            r = minimize(pooled_nll, [d0] + s0, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxiter=10000))
            if best is None or r.fun < best.fun:
                best = r
    dpool, cpool = best.x[0], best.x[1:]
    failed_all = np.concatenate([data[24]['failed'], data[18]['failed']])
    m_all = failed_all.mean()
    for fc in (24, 18):
        D = data[fc]; x = 20*np.log10(D['amps']); m = D['failed'].mean()
        res = {'points': [dict(amp=float(a_), packets=int(n), failures=int(f), fer=f/n, ber=e/(224*n))
                          for a_, n, f, e in zip(D['amps'], D['N'], D['F'], D['E'])],
               'errors_per_failure': float(m)}
        # model-free: log10 BER linear in dB between the bracketing measured points
        ber = D['E']/(224*D['N']); ok = ber > 0
        xs, ls = x[ok], np.log10(ber[ok])
        i = np.flatnonzero((ls[:-1] >= -3) & (ls[1:] < -3))
        if len(i):
            i = i[0]; xi = xs[i] + (-3 - ls[i])*(xs[i + 1] - xs[i])/(ls[i + 1] - ls[i])
            res['ber_1e-3_loglinear'] = float(10**(xi/20))
        for name in ('probit', 'probit2'):
            par, nll = fit_general(name, x, D['N'], D['F'])
            exp_f = D['N']*fer_model(name, x, par)
            lo, hi = x.min() - 1.0, x.max() + 2.0
            targets = {'fer_1e-2': 1e-2, 'fer_1e-3': 1e-3, 'ber_1e-3': 1e-3*224/m, 'ber_1e-4': 1e-4*224/m}
            est = {k: amp_at_general(name, par, t, lo, hi) for k, t in targets.items()}
            boot = {k: [] for k in targets}
            for _ in range(300):
                Fb = rng.binomial(D['N'], fer_model(name, x, par))
                mb = rng.choice(D['failed'], len(D['failed'])).mean()
                pb, _ = fit_general(name, x, D['N'], Fb, start=par)
                for k, t in {'fer_1e-2': 1e-2, 'fer_1e-3': 1e-3, 'ber_1e-3': 1e-3*224/mb, 'ber_1e-4': 1e-4*224/mb}.items():
                    try:
                        boot[k].append(amp_at_general(name, pb, t, lo, hi))
                    except ValueError:
                        pass
            ci = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in boot.items()}
            npar = 2 if name == 'probit' else 3
            res[name] = dict(params=[float(v) for v in par], deviance=deviance(D['N'], D['F'], exp_f),
                             dof=int(len(x) - npar), amplitude=est, ci95=ci, expected_failures=[float(v) for v in exp_f])
            print(f"{fc} Hz {name}: deviance {res[name]['deviance']:.2f} on {len(x) - npar} dof; "
                  + ', '.join(f'{k} {est[k]:.4f} [{ci[k][0]:.4f}, {ci[k][1]:.4f}]' for k in targets), flush=True)
        # pooled curve evaluated at this carrier
        off = 0.0 if fc == 24 else dpool
        xq = x - off
        exp_p = D['N']*fer_model('probit2', xq, cpool)
        tp = {'fer_1e-2': 1e-2, 'fer_1e-3': 1e-3, 'ber_1e-3': 1e-3*224/m_all, 'ber_1e-4': 1e-4*224/m_all}
        estp = {k: 10**((20*np.log10(amp_at_general('probit2', cpool, t, xq.min() - 1, xq.max() + 2)) + off)/20)
                for k, t in tp.items()}
        res['pooled'] = dict(carrier_offset_db=float(dpool), params=[float(v) for v in cpool],
                             deviance=deviance(D['N'], D['F'], exp_p), amplitude=estp,
                             expected_failures=[float(v) for v in exp_p])
        print(f"{fc} Hz pooled probit2 (18/24 scale {10**(dpool/20):.4f}): deviance {res['pooled']['deviance']:.2f}; "
              + ', '.join(f'{k} {estp[k]:.4f}' for k in tp), flush=True)
        if 'ber_1e-3_loglinear' in res:
            print(f"    model-free log-linear BER 1e-3: {res['ber_1e-3_loglinear']:.4f}")
        res['pass_probability'] = {}
        a0 = estp['ber_1e-3']
        for amp in np.round(np.arange(a0, a0*1.10, 0.005*(1.2 if fc == 18 else 1.0)), 4):
            pp, med, q = pass_probability_general('probit2', cpool, 20*np.log10(amp) - off, failed_all, a.sims, rng)
            res['pass_probability'][f'{amp:.4f}'] = dict(p_pass=pp, median_packets=med, fer=q)
            print(f'    pooled A={amp:.4f}: FER {q:.2e}, P(pass) {pp:.2f}, median packets {med}', flush=True)
        out[str(fc)] = res
    (RES/('nb_waterfall_fit2.json' if DESIGN == 'nb' else 'k5_waterfall_fit2.json')).write_text(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--sims', type=int, default=300)
    ap.add_argument('--carriers', type=int, nargs='+', default=[24, 18])
    ap.add_argument('--v2', action='store_true', help='quadratic probit, pooled carriers, model-free check')
    ap.add_argument('--design', default='nb', choices=('nb', 'k5'))
    a = ap.parse_args()
    global DESIGN
    DESIGN = a.design
    if a.v2:
        return main_v2(a)
    rows = load(); rng = np.random.default_rng(20371700); out = {}
    for fc in a.carriers:
        pts = summarize(rows, fc)
        amps = np.array([p[0] for p in pts]); N = np.array([p[1] for p in pts]); F = np.array([p[2] for p in pts])
        E = np.array([p[3] for p in pts]); x = 20*np.log10(amps)
        failed = np.concatenate([p[4] for p in pts if p[4]])
        mpf = [(p[0], float(np.mean(p[4])), len(p[4])) for p in pts if p[4]]
        res = dict(points=[dict(amp=float(p[0]), packets=int(p[1]), failures=int(p[2]), fer=p[2]/p[1],
                                ber=p[3]/(224*p[1])) for p in pts],
                   errors_per_failure=dict(mean=float(failed.mean()), by_amp=mpf))
        m = failed.mean()
        for model in MODELS:
            par, nll = fit(model, x, N, F)
            exp_f = N*MODELS[model](x, *par)
            dev = 2*np.sum(np.where(F > 0, F*np.log(np.maximum(F, 1e-300)/exp_f), 0) +
                           np.where(N - F > 0, (N - F)*np.log(np.maximum(N - F, 1e-300)/(N - exp_f)), 0))
            targets = {'fer_1e-2': 1e-2, 'fer_1e-3': 1e-3, 'ber_1e-3': 1e-3*224/m, 'ber_1e-4': 1e-4*224/m}
            est = {k: amp_at(model, par, t) for k, t in targets.items()}
            boot = {k: [] for k in targets}
            for _ in range(400):
                Fb = rng.binomial(N, MODELS[model](x, *par))
                mb = rng.choice(failed, len(failed)).mean()
                pb, _ = fit(model, x, N, Fb)
                tb = {'fer_1e-2': 1e-2, 'fer_1e-3': 1e-3, 'ber_1e-3': 1e-3*224/mb, 'ber_1e-4': 1e-4*224/mb}
                for k, t in tb.items():
                    try:
                        boot[k].append(amp_at(model, pb, t))
                    except ValueError:
                        pass
            ci = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in boot.items()}
            res[model] = dict(x0_db=float(par[0]), s_db=float(par[1]), nll=float(nll), deviance=float(dev),
                              dof=int(len(x) - 2), amplitude=est, ci95=ci,
                              expected_failures=[float(v) for v in exp_f])
            print(f'{fc} Hz {model}: x0 {par[0]:.3f} dB, s {par[1]:.3f} dB, deviance {dev:.2f} on {len(x) - 2} dof')
            for k in targets:
                print(f'    {k}: {est[k]:.4f} A0  [{ci[k][0]:.4f}, {ci[k][1]:.4f}]')
        # qualification pass probability along the probit fit (cloglog shown for model dependence)
        cands = np.round(np.arange(res['probit']['amplitude']['ber_1e-3'], res['probit']['amplitude']['ber_1e-3']*1.12, 0.01), 3)
        res['pass_probability'] = {}
        for model in MODELS:
            par = np.array([res[model]['x0_db'], res[model]['s_db']])
            for amp in cands:
                pp, med, q = pass_probability(model, par, amp, failed, a.sims, rng=rng)
                res['pass_probability'][f'{model}_{amp:.3f}'] = dict(p_pass=pp, median_packets=med, fer=q)
                print(f'    {model} A={amp:.3f}: FER {q:.2e}, P(pass) {pp:.2f}, median packets {med}')
        out[str(fc)] = res
        for p in pts:
            print(f'    point {p[0]:.3f}: {p[2]}/{p[1]} FER {p[2]/p[1]:.4f} BER {p[3]/(224*p[1]):.2e}')
    path = RES/'nb_waterfall_fit.json'                    # merge per-carrier runs
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update(out); path.write_text(json.dumps(old, indent=1))


if __name__ == '__main__':
    main()
