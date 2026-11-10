"""Preamble-based acquisition for the hop-limited link: Cramer-Rao bounds and a maximum-likelihood estimator.

Frame (224 s, T = 224/190 s): 172 data symbols, 6 centre-tone tail symbols, 12 acquisition symbols.
The acquisition symbols are a known preamble that follows the previous frame's centre-tone tail, so the
transmitted phase is continuous and the receiver's filter enters the preamble in the steady centre-tone
state. Unknowns: complex gain g = (A/A_nom) e^{i phi} and timing tau (signal delay). The carrier frequency is
taken from the clocks (see freq_requirement); the bound with an unknown frequency offset is also reported.

All quantities use the packet simulations' normalisation: unit nominal signal amplitude, noise PSD divided by
A_nom^2 (CausalReceiver.input_psd). The likelihood is the exact finite-record Gaussian likelihood of the
receiver's output samples (FiniteRecordInnovations: Cholesky whitening of the output covariance).

Usage: python acquisition.py            (bounds for candidate preambles + Monte Carlo of the estimator)
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver
import hop_packets as hp

HERE = Path(__file__).resolve().parent
T = 224/190
FS = 16.0
QUAL = {24.0: 0.58, 18.0: 0.697}            # qualified amplitudes (units of A0)
PREAMBLES = {
    'maxhop': [6, 3, 0, 3]*3,
    'maxhop_mirror': [0, 3, 6, 3]*3,
    'slow': [6, 6, 3, 0, 0, 3]*2,
    'ramp': [4, 5, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3],
    'mixed': [6, 3, 6, 3, 0, 3, 0, 3, 6, 3, 0, 3],
    'centre': [3]*12,
}
SEARCH = 0.5                                 # coarse timing window +-SEARCH s (frame sync assumed to this)
GUARD = 0.6                                  # preamble samples within GUARD s of its end are not used


def phase_at(tones, t, h=hp.H, beta=hp.BETA, M=hp.M):
    """Continuous CPM phase (cycles) of tone sequence `tones` starting at t = 0, preceded by the centre tone
    (phase 0) and followed by its last tone; same pulse as hop_receiver.sampled_hop."""
    offs = np.asarray(tones, float) - (M - 1)/2
    beta = float(beta); spacing = float(h)/T
    before = np.r_[0., offs[:-1]]
    inc = spacing*T*(beta*before/2 + (1 - beta/2)*offs)
    phase0 = np.r_[0, np.cumsum(inc)]
    t = np.asarray(t, float); n = len(offs)
    idx = np.clip(np.floor(t/T).astype(int), 0, n - 1)
    local = t - idx*T; tr = beta*T; u = np.clip(local/tr, 0, 1)
    prim = tr*(2.5*u**4 - 3*u**5 + u**6) + np.maximum(local - tr, 0)
    ph = phase0[idx] + spacing*(before[idx]*local + (offs[idx] - before[idx])*prim)
    ph = np.where(t < 0, 0.0, ph)
    after = t >= n*T
    ph = np.where(after, phase0[n] + spacing*offs[-1]*(t - n*T), ph)
    return ph


class Acquirer:
    """t_start: time of the first observed input sample relative to the nominal preamble start (>= 0).
    The model starts PAD samples earlier in the centre-tone prehistory, so an early signal (tau < 0) enters
    the filter state correctly."""
    PAD = 32                                                      # 2 s of prehistory (multiple of decimation)

    def __init__(self, carrier, amp_factor, preamble, search=SEARCH, guard=GUARD, t_start=0.0):
        self.carrier, self.A = carrier, amp_factor*mg.A0
        self.r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
        self.pre = list(preamble)
        self.D = self.r.decimation
        self.t_start = t_start
        n_in = int(round((len(self.pre)*T - guard - t_start)*FS))
        self.n_out = n_in//self.D                                 # output samples 0..n_out-1 (input index j*D)
        self.n_in = self.n_out*self.D
        self.innov = self.r.finite_record(self.n_out, self.A)
        self.search = search

    def model(self, tau, nu=0.0):
        """Unit-gain receiver output for the preamble delayed by tau, optional frequency offset nu (Hz)."""
        k = np.arange(-self.PAD, self.n_in)
        t = self.t_start + k/FS
        x = np.exp(2j*np.pi*(phase_at(self.pre, t - tau) + nu*t))
        return self.r.stream(1.0).process(x)[self.PAD//self.D:]

    def fisher(self, tau=0.0, with_nu=False, dt=1e-4, dnu=1e-6):
        m = self.model(tau)
        cols = [m, 1j*m, (self.model(tau + dt) - self.model(tau - dt))/(2*dt)]
        if with_nu:
            cols.append((self.model(tau, dnu) - self.model(tau, -dnu))/(2*dnu))
        W = np.column_stack([self.innov.whiten(c) for c in cols])
        return 2*np.real(W.conj().T @ W)

    def bounds(self, with_nu=False):
        C = np.linalg.inv(self.fisher(with_nu=with_nu))
        out = {'amp_rel_std': math.sqrt(C[0, 0]), 'phase_std_deg': math.degrees(math.sqrt(C[1, 1])),
               'timing_std_s': math.sqrt(C[2, 2])}
        if with_nu:
            out['freq_std_Hz'] = math.sqrt(C[3, 3])
            out['phase_drift_std_deg_over_frame'] = math.degrees(2*math.pi*math.sqrt(C[3, 3])*224)
        return out

    def estimate(self, y):
        """ML estimate from the first n_out output samples y: profile out g, search tau on a grid, refine."""
        wy = self.innov.whiten(y[:self.n_out])
        cache = {}

        def prof(tau):
            wm = self.innov.whiten(self.model(tau))
            num = np.vdot(wm, wy); den = np.real(np.vdot(wm, wm))
            cache[tau] = num/den
            return -abs(num)**2/den
        grid = np.arange(-self.search, self.search + 1e-9, 0.025)
        vals = [prof(t) for t in grid]
        k = int(np.argmin(vals))
        lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
        res = minimize_scalar(prof, bounds=(lo, hi), method='bounded', options={'xatol': 1e-5})
        tau = float(res.x); g = cache.get(tau)
        if g is None:
            prof(tau); g = cache[tau]
        return complex(g), tau


def freq_requirement(phase_deg=5.0):
    """Clock offset that keeps the phase drift across one 224 s frame below phase_deg."""
    return phase_deg/360/224


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=2000)
    a = ap.parse_args()
    out = {'T_s': T, 'search_window_s': SEARCH, 'guard_s': GUARD, 'bounds': {}, 'monte_carlo': {},
           'clock_offset_for_5deg_Hz': freq_requirement(), 'clock_ppm_for_5deg': {}}
    for fc, af in QUAL.items():
        out['clock_ppm_for_5deg'][f'{fc:g}'] = freq_requirement()/fc*1e6
        for name, pre in PREAMBLES.items():
            acq = Acquirer(fc, af, pre)
            try:
                b = acq.bounds(); bn = acq.bounds(with_nu=True)
            except np.linalg.LinAlgError:
                b = bn = 'singular'
            out['bounds'][f'{fc:g}|{name}'] = {'known_frequency': b, 'unknown_frequency': bn}
            print(fc, name, b, flush=True)
    for fc, af in QUAL.items():
        acq = Acquirer(fc, af, PREAMBLES['maxhop'])
        b = out['bounds'][f'{fc:g}|maxhop']['known_frequency']
        rng = np.random.default_rng(np.random.SeedSequence([20371601, int(fc)]))
        errs = []
        m0 = acq.model(0.0)
        for i in range(a.trials):
            tau0 = rng.uniform(-0.25, 0.25); phi0 = rng.uniform(-np.pi, np.pi)
            y = np.exp(1j*phi0)*acq.model(tau0) + acq.r.stationary_noise(acq.n_out, rng, acq.A)
            g, tau = acq.estimate(y)
            errs.append((abs(g) - 1, math.degrees(np.angle(g*np.exp(-1j*phi0))), tau - tau0))
        e = np.array(errs)
        out['monte_carlo'][f'{fc:g}'] = {
            'trials': a.trials, 'amp_rel_std': float(e[:, 0].std()), 'amp_bias': float(e[:, 0].mean()),
            'phase_std_deg': float(e[:, 1].std()), 'phase_bias_deg': float(e[:, 1].mean()),
            'timing_std_s': float(e[:, 2].std()), 'timing_bias_s': float(e[:, 2].mean()),
            'max_abs_phase_deg': float(abs(e[:, 1]).max()), 'max_abs_timing_s': float(abs(e[:, 2]).max()),
            'outliers_timing_gt_5sigma': int(np.sum(abs(e[:, 2]) > 5*b['timing_std_s'])),
            'crb': b}
        print(fc, out['monte_carlo'][f'{fc:g}'], flush=True)
    (HERE/'results'/'acquisition_crb.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
