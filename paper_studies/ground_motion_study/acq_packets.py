"""Packet simulation of the GF(512) joint link with simulated acquisition.

Each packet transmits the previous frame's centre-tone tail (steady), the 12-symbol known preamble
(acquisition.PREAMBLES['maxhop']: 6,3,0,3 x 3, zero net phase), then the frame exactly as nb_packets.run
(172 data symbols + 6 tail symbols), with an unknown carrier phase phi0 ~ U(-pi, pi) and an unknown delay
tau ~ U(-0.25, 0.25) s, through one continuous receiver stream and one stationary noise realisation.
  variant 'acq':   the receiver estimates g = (A/A_nom) e^{i phi} and tau from the preamble (acquisition.Acquirer,
                   exact finite-record likelihood, coarse window +-0.5 s), divides the samples by g_hat,
                   uses amplitude |g_hat| A_nom, timing tau_hat and the preamble-implied filter state;
  variant 'acqcs': 'acq', then up to NSYNC full-frame re-estimates from the decoder's own hard decisions after
                   a CRC failure (code-aided resynchronization; does not recover failed frames);
  variant 'track': steady state on a continuous link (or after one known training frame): the previous frame,
                   with the same g and tau but independent payload and noise, passed its CRC, so its tones are
                   known and give a full-frame ML estimate; it is Fisher-combined with the current preamble estimate;
  variant 'known': identical samples, true g and tau supplied (paired reference).
Detection and decoding are otherwise those of nb_packets.run (CPU detector, 16 passes, 30 BP iterations).
The carrier frequency is assumed known from the clocks (acquisition.freq_requirement).
Usage:
  python acq_packets.py --selftest
  python acq_packets.py --case acq24=24=0.525=acq --case known24=24=0.525=known --packets 2048 --start 9400000
  python acq_packets.py --qualify --case 24=0.58 --case 18=0.697 --start 9500000     (betting CS as qualify_nb)
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, json, math, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from scipy.signal import sosfilt, sosfilt_zi
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
from gravcomm.packet_confidence import packet_ber_interval
import hop_packets as hp
from hop_receiver import HopLimitedReceiver
import turbo_rx as tr
import nb_ldpc as nb
import nb_packets as npk
import acquisition as aq
from acq_receiver import AcqHopReceiver

RES = Path(__file__).resolve().parent/'results'
PRE = aq.PREAMBLES['maxhop']
SEEDS = (20371701, 20371702, 20371703)          # payload, noise, acquisition-offset namespaces
ITER, BP, CODE_SEED = 16, 30, 1
NSYNC = 3                                       # code-aided resynchronisation rounds after a CRC failure
FS = 16.0


def build(carrier, amp_factor, index, seeds, tau_true=None, phi_true=None):
    """Transmit preamble + frame; return receiver outputs and bookkeeping (absolute input index 0 = data start)."""
    K = npk.K
    rng = np.random.default_rng(np.random.SeedSequence([seeds[0], K, index]))
    payload = rng.integers(0, 256, (K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + binascii.crc_hqx(body, 0xffff).to_bytes(2, 'big'), dtype=np.uint8))
    code = npk.get_code(CODE_SEED, 2)
    x = code.encode(npk.info_symbols(info))
    perm = np.random.default_rng(4100 + CODE_SEED).permutation(code.n)
    tones, prev = [], hp.M//2
    for w in x[perm]:
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    ts = 224/(len(tones) + 18)
    assert abs(ts - aq.T) < 1e-12
    orng = np.random.default_rng(np.random.SeedSequence([seeds[2], K, index]))
    tau = orng.uniform(-0.25, 0.25) if tau_true is None else tau_true
    phi = orng.uniform(-np.pi, np.pi) if phi_true is None else phi_true
    g = np.exp(1j*phi)
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    D = r.decimation
    combined = PRE + list(tones) + [hp.M//2]*6
    total = len(tones) + 6
    i_first = -int(math.ceil((len(PRE)*ts + 3.0)*FS/D))*D        # 3 s of centre-tone prehistory, multiple of D
    i_last = int(math.ceil((total*ts + 1.0)*FS))
    t = np.arange(i_first, i_last)/FS
    xin = g*np.exp(2j*np.pi*aq.phase_at(combined, t + len(PRE)*ts - tau))
    mean = r.stream(g).process(xin)                               # output j <-> input i_first + j*D
    nrng = np.random.default_rng(np.random.SeedSequence([seeds[1], K, index]))
    y = mean + r.stationary_noise(len(mean), nrng, amp_factor*mg.A0)
    return dict(info=info, payload=payload, x=x, perm=perm, tones=tones, ts=ts, r=r, D=D, i_first=i_first,
                y=y, tau=tau, g=g, ngroups=len(x), code=code, total=total)


def acquire(b, carrier, amp_factor):
    """ML acquisition on the preamble samples."""
    D, ts = b['D'], b['ts']
    i_s = int(math.ceil(-len(PRE)*ts*FS/D))*D                     # first output-aligned input index in preamble
    acq = aq.Acquirer(carrier, amp_factor, PRE, t_start=i_s/FS + len(PRE)*ts)
    j_s = (i_s - b['i_first'])//D
    return acq.estimate(b['y'][j_s:j_s + acq.n_out])


def detect_decode(b, amp_factor, g_hat, tau_hat, detector=AcqHopReceiver, device='cpu'):
    r, D, ts = b['r'], b['D'], b['ts']
    A_hat = abs(g_hat)*amp_factor*mg.A0
    y = b['y']/g_hat
    start0 = int(np.ceil(tau_hat*FS - 1e-12))
    # unit-gain mean-model state just before input index start0 (prehistory + preamble, delay tau_hat)
    sos = r.sos
    ii = np.arange(b['i_first'], start0)
    xin = np.exp(2j*np.pi*aq.phase_at(PRE, ii/FS + len(PRE)*ts - tau_hat))
    _, zf = sosfilt(sos, xin, zi=sosfilt_zi(sos).astype(complex))
    kw = dict(history_symbols=4, beam_width=32768, lookahead_symbols=12, h=hp.H, beta=hp.BETA, tone_count=hp.M,
              max_hop=hp.D, group_symbols=hp.L, group_bits=hp.BITS, pad_bits=0)
    if device == 'gpu':
        from acq_receiver_gpu import AcqHopReceiverGPU
        from hop_receiver_gpu import GraphGPU
        det = AcqHopReceiverGPU(r, ts, A_hat, b['ngroups'], time_offset=tau_hat, initial_state=zf, **kw)
    elif detector is AcqHopReceiver:
        det = AcqHopReceiver(r, ts, A_hat, b['ngroups'], time_offset=tau_hat, initial_state=zf, **kw)
    else:
        det = detector(r, ts, A_hat, b['ngroups'], **kw)
    if device != 'gpu':
        tr.record_graph(det)
    first_out = (start0 + (-start0) % D - b['i_first'])//D
    pos = first_out
    for symbol in range(b['total']):
        count = det._geometry(symbol)[3]
        det.push(y[pos:pos + count]); pos += count
    stats = det.finish()
    if device == 'gpu':
        app_fn = GraphGPU(det, hp.L, b['ngroups']).word_app
    else:
        packed = tr.pack_graph(det, hp.L, b['ngroups'])
        app_fn = lambda Lp: npk.graph_word_app(packed, Lp)
    code, perm, info, x = b['code'], b['perm'], b['info'], b['x']
    mask = np.zeros((code.n, nb.Q)); mask[code.info[28], (np.arange(nb.Q) & 31) != 0] = npk.FLOOR
    inv = np.argsort(perm); Lp = mask[perm]
    errors, accepted = [], []
    for it in range(ITER):
        app = npk.normalize_log(app_fn(Lp))
        ext = npk.normalize_log(app - Lp)[inv] + mask
        prior = np.exp(ext - ext.max(1, keepdims=True)); prior /= prior.sum(1, keepdims=True)
        post, cext, hard, ok = code.decode(prior, BP)
        bits = npk.symbols_info(hard[code.info])
        errors.append(int(np.sum(bits[16:-16] != info[16:-16])))
        accepted.append(bool(ok and npk.crc_ok(bits)))
        if accepted[-1]:
            break
        Lp = npk.normalize_log(np.log(np.maximum(cext, 1e-300)) + mask)[perm]
    return dict(errors=errors, accepted=accepted, pruned_states=stats['pruned_states'], det=det, hard=hard)


_INNOV = {}


def resync(b, amp_factor, hard, tau0, span=0.25, fisher=False):
    """Code-aided resynchronisation: ML estimate of (g, tau) over preamble + frame, with the data tones taken
    from the decoder's hard codeword decisions (re-encoded exactly as the transmitter does)."""
    r, D, ts = b['r'], b['D'], b['ts']
    tones, prev = [], hp.M//2
    for w in np.asarray(hard)[b['perm']]:
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    seq = PRE + tones + [hp.M//2]*6
    i_s = int(math.ceil(-len(PRE)*ts*FS/D))*D
    n_out = int(((b['total']*ts - aq.GUARD) - i_s/FS)*FS)//D
    key = (b['r'].carrier_hz, amp_factor, n_out)
    if key not in _INNOV:
        _INNOV[key] = r.finite_record(n_out, amp_factor*mg.A0)
    innov = _INNOV[key]
    j_s = (i_s - b['i_first'])//D
    wy = innov.whiten(b['y'][j_s:j_s + n_out])
    pad = aq.Acquirer.PAD
    k = np.arange(-pad, n_out*D)
    t = i_s/FS + k/FS
    cache = {}

    def prof(tau):
        m = r.stream(1.0).process(np.exp(2j*np.pi*aq.phase_at(seq, t + len(PRE)*ts - tau)))[pad//D:]
        wm = innov.whiten(m)
        num = np.vdot(wm, wy); den = np.real(np.vdot(wm, wm)); cache[tau] = num/den
        return -abs(num)**2/den
    grid = np.arange(tau0 - span, tau0 + span + 1e-9, 0.01)
    vals = [prof(x) for x in grid]
    kk = int(np.argmin(vals))
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(prof, bounds=(grid[max(kk - 1, 0)], grid[min(kk + 1, len(grid) - 1)]),
                          method='bounded', options={'xatol': 1e-5})
    tau = float(res.x)
    if tau not in cache:
        prof(tau)
    g = complex(cache[tau])
    if not fisher:
        return g, tau
    model = lambda x: innov.whiten(r.stream(1.0).process(np.exp(2j*np.pi*aq.phase_at(seq, t + len(PRE)*ts - x)))[pad//D:])
    wm = model(tau); dwm = (model(tau + 1e-4) - model(tau - 1e-4))/2e-4
    W = np.column_stack([wm, 1j*wm, g*dwm])
    return g, tau, 2*np.real(W.conj().T @ W)


def preamble_fisher(b, carrier, amp_factor, g, tau):
    D, ts = b['D'], b['ts']
    i_s = int(math.ceil(-len(PRE)*ts*FS/D))*D
    acq = aq.Acquirer(carrier, amp_factor, PRE, t_start=i_s/FS + len(PRE)*ts)
    m = acq.model(tau); dm = (acq.model(tau + 1e-4) - acq.model(tau - 1e-4))/2e-4
    W = np.column_stack([acq.innov.whiten(c) for c in (m, 1j*m, g*dm)])
    return 2*np.real(W.conj().T @ W)


def combine(estimates):
    """Inverse-covariance (Fisher-weighted) combination of (g, tau, J) estimates in (Re g, Im g, tau)."""
    Js = sum(J for _, _, J in estimates)
    v = sum(J @ np.array([g.real, g.imag, tau]) for g, tau, J in estimates)
    th = np.linalg.solve(Js, v)
    return complex(th[0], th[1]), float(th[2])


PREV_SEEDS = (20371711, 20371712, 20371713)     # independent payload/noise of the previous (CRC-verified) frame


def run(job, device='cpu'):
    carrier, amp_factor, index, variant, seeds = job
    t0 = time.perf_counter()
    b = build(carrier, amp_factor, index, seeds)
    if variant in ('acq', 'acqcs'):
        g_hat, tau_hat = acquire(b, carrier, amp_factor)
    elif variant == 'track':
        # Steady-state tracking on a continuous link: the previous frame (same gain and timing, independent
        # payload and noise) passed its CRC, so its tones are known and it gives a full-frame estimate; it is
        # combined with the current preamble estimate by Fisher weighting.
        bp = build(carrier, amp_factor, index, PREV_SEEDS, tau_true=b['tau'], phi_true=float(np.angle(b['g'])))
        gp0, tp0 = acquire(bp, carrier, amp_factor)
        gf, tf, Jf = resync(bp, amp_factor, bp['x'], tp0, fisher=True)
        g1, t1 = acquire(b, carrier, amp_factor)
        g_hat, tau_hat = combine([(gf, tf, Jf), (g1, t1, preamble_fisher(b, carrier, amp_factor, g1, t1))])
    else:
        g_hat, tau_hat = b['g'], b['tau']
    out = detect_decode(b, amp_factor, g_hat, tau_hat, device=device)
    out.pop('det')
    sync = [dict(phase_err_deg=float(np.degrees(np.angle(g_hat/b['g']))), timing_err_s=float(tau_hat - b['tau']),
                 amp_err=float(abs(g_hat) - 1), accepted=out['accepted'][-1], errors=out['errors'][-1])]
    if variant == 'acqcs':
        # On a CRC failure, re-estimate gain and timing over the whole frame from the decoder's decisions and
        # decode again (up to NSYNC rounds). Uses only the receiver's own decisions.
        for _ in range(NSYNC):
            if out['accepted'][-1]:
                break
            g_hat, tau_hat = resync(b, amp_factor, out['hard'], tau_hat)
            out = detect_decode(b, amp_factor, g_hat, tau_hat, device=device); out.pop('det')
            sync.append(dict(phase_err_deg=float(np.degrees(np.angle(g_hat/b['g']))),
                             timing_err_s=float(tau_hat - b['tau']), amp_err=float(abs(g_hat) - 1),
                             accepted=out['accepted'][-1], errors=out['errors'][-1]))
    out['sync_rounds'] = sync
    out['hard'] = None
    out.update(carrier_hz=carrier, amplitude_A0=amp_factor, index=index, variant=variant, seeds=list(seeds),
               device=device,
               tau_true=b['tau'], phi_true=float(np.angle(b['g'])), tau_hat=tau_hat,
               amp_err=float(abs(g_hat) - 1), phase_err_deg=float(np.degrees(np.angle(g_hat/b['g']))),
               timing_err_s=float(tau_hat - b['tau']), elapsed_s=time.perf_counter() - t0)
    return out


def selftest():
    # (1) AcqHopReceiver with zero offset and steady state equals HopLimitedReceiver on the same samples
    b = build(24.0, 0.55, 1, SEEDS, tau_true=0.0, phi_true=0.0)
    kw = dict(history_symbols=4, beam_width=32768, lookahead_symbols=12, h=hp.H, beta=hp.BETA, tone_count=hp.M,
              max_hop=hp.D, group_symbols=hp.L, group_bits=hp.BITS, pad_bits=0)
    dets = [cls(b['r'], b['ts'], 0.55*mg.A0, b['ngroups'], **kw) for cls in (HopLimitedReceiver, AcqHopReceiver)]
    for det in dets:
        pos = (0 - b['i_first'])//b['D']
        for symbol in range(b['total']):
            count = det._geometry(symbol)[3]; det.push(b['y'][pos:pos + count]); pos += count
    print('subclass identical to HopLimitedReceiver at zero offset:',
          bool(np.array_equal(dets[0].cost, dets[1].cost) and np.array_equal(dets[0].hist, dets[1].hist)))
    # (2) noise-free signal with an offset decodes with near-zero path metric when tau_hat = tau
    b = build(24.0, 0.55, 2, SEEDS, tau_true=0.1370, phi_true=1.0)
    b['y'] = b['y'] - (b['r'].stationary_noise(len(b['y']), np.random.default_rng(np.random.SeedSequence([SEEDS[1], npk.K, 2])), 0.55*mg.A0))
    o = detect_decode(b, 0.55, b['g'], b['tau'])
    print('noise-free offset 0.137 s: accepted', o['accepted'][-1], 'errors', o['errors'][-1],
          'best-path whitened metric (0 if offset and initial state are exact): %.3e' % o['det'].path_min)
    # (3) acquisition on one noisy packet
    b = build(24.0, 0.58, 3, SEEDS)
    g_hat, tau_hat = acquire(b, 24.0, 0.58)
    print('acq: true tau %.4f phi %.3f -> est tau %.4f, phase err %.2f deg, amp err %.3f' %
          (b['tau'], np.angle(b['g']), tau_hat, np.degrees(np.angle(g_hat/b['g'])), abs(g_hat) - 1))


def qstatus(folder, start, max_packets):
    rows = {int(p.stem): json.loads(p.read_text()) for p in folder.glob('[0-9]*.json')}
    n = 0
    while start + n in rows:
        n += 1
    out = None
    for k in range(128, n + 1, 128):
        e = [rows[start + i]['errors'][-1] for i in range(k)]
        lo, hi = packet_ber_interval(e, 224, 0.025)
        out = dict(packets=k, errors=int(sum(e)), ber=sum(e)/(224*k), interval=[lo, hi],
                   rejected=sum(not rows[start + i]['accepted'][-1] for i in range(k)))
        if hi <= 1e-3:
            return 'pass', out
        if lo > 1e-3:
            return 'fail', out
    return ('undecided' if n < max_packets else 'exhausted'), out


def run_jobs(jobs, cpu_workers, gpu_workers):
    """(job, tag) pairs on a CPU pool and a GPU pool fed from one queue; yields (result, tag)."""
    from concurrent.futures import wait, FIRST_COMPLETED
    queue = list(jobs)[::-1]; running = {}
    pools = [(d, ProcessPoolExecutor(max_workers=n), n) for d, n in (('cpu', cpu_workers), ('gpu', gpu_workers)) if n]
    try:
        def fill():
            for dev, pool, cap in pools:
                while queue and sum(1 for v in running.values() if v[0] == dev) < cap:
                    job, tag = queue.pop(); running[pool.submit(run, job, dev)] = (dev, tag)
        fill()
        while running:
            done, _ = wait(list(running), return_when=FIRST_COMPLETED)
            for f in done:
                dev, tag = running.pop(f); yield f.result(), tag
            fill()
    finally:
        for _, pool, _ in pools:
            pool.shutdown(wait=True)


def write(res, folder):
    tmp = folder/f"{res['index']}.json.tmp"
    tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--selftest', action='store_true')
    p.add_argument('--qualify', action='store_true')
    p.add_argument('--case', action='append', default=[])
    p.add_argument('--packets', type=int, default=2048)
    p.add_argument('--max-packets', type=int, default=16896)
    p.add_argument('--start', type=int, default=9400000)
    p.add_argument('--workers', type=int, default=18)
    p.add_argument('--batch', type=int, default=512)
    p.add_argument('--gpu-workers', type=int, default=0)
    p.add_argument('--variant', default='acq', help="qualification variant: 'acq' (cold start) or 'track'")
    p.add_argument('--gpucheck', action='store_true', help='CPU vs GPU detector on a few acquisition packets')
    a = p.parse_args()
    if a.selftest:
        return selftest()
    if a.gpucheck:
        rows = []
        for i in range(4):
            job = (24.0, 0.525, 9400000 + i, 'acq', SEEDS)
            c, g = run(job, 'cpu'), run(job, 'gpu')
            same = (c['errors'] == g['errors'] and c['accepted'] == g['accepted'] and c['tau_hat'] == g['tau_hat'])
            rows.append(same)
            print(i, 'cpu', c['errors'], c['accepted'][-1], 'gpu', g['errors'], g['accepted'][-1],
                  'same' if same else 'DIFFERENT', f"cpu {c['elapsed_s']:.1f}s gpu {g['elapsed_s']:.1f}s", flush=True)
        print('all identical:', all(rows))
        return
    if a.qualify:
        cases = []
        for c in a.case:
            fc, amp = (float(v) for v in c.split('='))
            folder = RES/f'{a.variant}qual_f{fc:g}_a{amp:.4f}'; folder.mkdir(parents=True, exist_ok=True)
            cases.append((fc, amp, folder))
        t0 = time.time()
        while True:
            jobs = []
            for fc, amp, folder in cases:
                st, info = qstatus(folder, a.start, a.max_packets)
                print(f'{folder.name}: {st} {info}', flush=True)
                if st != 'undecided':
                    continue
                have = {int(q.stem) for q in folder.glob('[0-9]*.json')}
                todo = [i for i in range(a.start, a.start + a.max_packets) if i not in have][:a.batch]
                jobs += [((fc, amp, i, a.variant, SEEDS), folder) for i in todo]
            if not jobs:
                break
            for res, folder in run_jobs(jobs, a.workers, a.gpu_workers):
                write(res, folder)
            print(f'  batch done ({time.time() - t0:.0f}s)', flush=True)
        summary = {f.name: dict(zip(('status', 'result'), qstatus(f, a.start, a.max_packets))) for _, _, f in cases}
        (RES/f'acquisition_qualification_{a.variant}.json').write_text(json.dumps(summary, indent=1) + '\n')
        print(json.dumps(summary, indent=1))
        return
    jobs = []
    for c in a.case:
        name, fc, amp, variant = c.split('=')
        folder = RES/name; folder.mkdir(parents=True, exist_ok=True)
        for i in range(a.start, a.start + a.packets):
            if not (folder/f'{i}.json').exists():
                jobs.append(((float(fc), float(amp), i, variant, SEEDS), folder))
    print(f'{len(jobs)} packets to run', flush=True)
    t0 = time.time(); done = 0
    for res, folder in run_jobs(jobs, a.workers, a.gpu_workers):
        write(res, folder); done += 1
        if done % 64 == 0:
            print(f'  {done}/{len(jobs)} ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
