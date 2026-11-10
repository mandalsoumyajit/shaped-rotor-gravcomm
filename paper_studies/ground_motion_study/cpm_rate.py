"""Information rate of torque-limited CF-FSK (continuous-phase modulation) in colored receiver noise.

The rotor's gravitational signal is s(t) = A exp(j phi(t)) at complex baseband, with
instantaneous frequency offset f(t) = dphi/dt / 2pi following the CF-FSK rule: in symbol j
the frequency moves from a_{j-1} df to a_j df over beta T with transition shape s(u), then
dwells. Tones a are equally spaced (integers for odd M, half-integers for even M); the
modulation index is h = df T.  Peak inertial torque is

    tau_pk = pi I_z c D df / (beta T),

with c the peak slope of the shape (15/8 quintic, 3/2 cubic, 1 linear) and D the largest
allowed hop in tone spacings (D = M-1 when every hop is allowed).

The information rate with a stationary Markov input P(a_j | a_{j-1}) is computed by the
simulation-based method of Arnold, Loeliger, Vontobel, Kavcic and Zeng (IEEE Trans. IT 52,
3498, 2006): simulate a long sequence through the true channel (exact colored Gaussian noise)
and evaluate log q(y|a) - log q(y) with a finite-state auxiliary channel q. Here q whitens
the samples with an order-p linear predictor (p <= samples per symbol) and treats the
innovations as white. Its trellis state is (phase lattice point, a_{j-r-1}..a_{j-1}), with r chosen
so the predictor memory fits; phase is
exact because h is rational. The estimate is a lower bound on the true information rate for
that input (auxiliary-channel bound) and becomes exact as the predictor converges.

Units: noise PSDs are normalized to unit signal amplitude, so `amp` enters only as A^2/N.
"""
from dataclasses import dataclass, field
from fractions import Fraction
import numpy as np
import numba as nb
from scipy.linalg import solve_toeplitz

SHAPES = {   # name: (primitive F(u) = int_0^u s, peak slope c)
    'quintic': (lambda u: 2.5*u**4 - 3*u**5 + u**6, 15/8),
    'cubic': (lambda u: u**3 - 0.5*u**4, 1.5),
    'linear': (lambda u: 0.5*u**2, 1.0),
}


@dataclass
class Design:
    M: int                      # number of tones
    h: Fraction                 # modulation index df*T of the (base) tone spacing df (rational)
    beta: float                 # transition fraction of the symbol (0 < beta <= 1)
    T: float                    # symbol duration (s)
    shape: str = 'quintic'
    max_hop: int | None = None  # largest allowed hop in tone indices (None: all hops)
    probs: np.ndarray | None = field(default=None, repr=False)   # Markov matrix P[i,j]=P(a_j|a_{j-1}=i)
    positions: tuple | None = None   # integer tone offsets in units of df (None: equally spaced)
    weights: tuple | None = None     # i.i.d. tone probabilities (None: uniform / maxentropic)

    def __post_init__(self):
        if self.positions is not None:
            self.M = len(self.positions)

    @property
    def D(self):
        return self.M - 1 if self.max_hop is None else min(self.max_hop, self.M - 1)

    @property
    def spacing(self):
        """Base tone spacing df (Hz); tone offsets are tones()*df."""
        return float(self.h)/self.T

    @property
    def max_jump(self):
        """Largest commanded frequency hop, in units of df, over allowed transitions."""
        t = self.tones(); ok = self.transition() > 0
        return float(np.max(np.abs(t[:, None] - t[None, :])[ok]))

    @property
    def fdot_max(self):
        """Peak frequency slew (Hz/s); torque = pi I_z fdot_max."""
        return SHAPES[self.shape][1]*self.max_jump*self.spacing/(self.beta*self.T)

    @property
    def excursion(self):
        return float(np.max(np.abs(self.tones())))*self.spacing

    def tones(self):
        if self.positions is not None:
            return np.asarray(self.positions, float)
        return np.arange(self.M) - (self.M - 1)/2

    def lattice(self):
        """Integer tone labels and the phase lattice (unit = k/P cycles per label)."""
        if self.positions is not None:
            labels = np.asarray(self.positions, np.int64); unit = Fraction(self.h)
        else:
            two = 1 if self.M % 2 else 2
            labels = (self.tones()*two).round().astype(np.int64)
            unit = Fraction(self.h)/two
        return labels, unit.numerator, unit.denominator

    def transition(self):
        if self.probs is not None:
            return np.asarray(self.probs, float)
        i = np.arange(self.M)
        adj = (np.abs(i[:, None] - i[None, :]) <= self.D).astype(float)
        if self.weights is not None:
            w = np.asarray(self.weights, float)*adj
            return w/w.sum(axis=1, keepdims=True)
        if self.D >= self.M - 1:
            return adj/self.M
        w, v = np.linalg.eig(adj)                          # maxentropic chain on the hop graph
        k = np.argmax(w.real); lam = w[k].real; v = np.abs(v[:, k].real)
        return adj*v[None, :]/(lam*v[:, None])

    def entropy_rate(self):
        p = self.transition(); pi = stationary(p)
        with np.errstate(divide='ignore', invalid='ignore'):
            ent = -np.nansum(p*np.where(p > 0, np.log2(p), 0), axis=1)
        return float(pi @ ent)


def stationary(p):
    w, v = np.linalg.eig(p.T)
    s = np.abs(v[:, np.argmin(np.abs(w - 1))].real)
    return s/s.sum()


def phase_track(design, symbols, ns):
    """Baseband phase (cycles) sampled at (n+0.5)/fs, fs = ns/T, for a symbol sequence.

    Symbol 0 starts at phase 0 with previous tone equal to symbols[0]."""
    a = design.tones()[symbols]
    prev = np.r_[a[0], a[:-1]]
    F, _c = SHAPES[design.shape]
    T, b, df = design.T, design.beta, design.spacing
    inc = df*T*(b/2*prev + (1 - b/2)*a)
    start = np.r_[0, np.cumsum(inc)[:-1]]
    tau = (np.arange(ns) + 0.5)/ns*T
    tr = b*T; u = np.minimum(tau/tr, 1)
    prim = tr*F(u) + np.maximum(tau - tr, 0)               # int_0^tau of the step shape
    return (start[:, None] + df*(prev[:, None]*tau[None, :] + (a - prev)[:, None]*prim[None, :])).ravel()


def predictor(psd_samples, p):
    """Order-p forward linear predictor for a discrete complex PSD sampled on an FFT grid."""
    r = np.fft.ifft(psd_samples)
    if p == 0:
        return np.zeros(0, complex), float(r[0].real)
    a = solve_toeplitz((r[:p], r[:p].conj()), r[1:p + 1])
    return a, float((r[0] - a @ r[1:p + 1].conj()).real)


def whiten(x, a):
    """e_n = x_n - sum_k a_k x_{n-k} along the last axis (zero initial conditions)."""
    out = x.copy()
    for k, ak in enumerate(a, start=1):
        out[..., k:] -= ak*x[..., :-k]
    return out


@nb.njit(cache=True)
def _forward(tmpl, yw, labels, k, P, M, r, logp, amp, sig2, start, truth, inc):
    """Forward recursion of the auxiliary trellis; returns (log q(y), log q(y|a)) in nats.

    State: (psi, s_0..s_r) with s_r = a_{j-1}; template index (hist, a_j) = hist*M + a_j.
    """
    n, ns = yw.shape
    H = M**(r + 1); Mr = M**r
    E = np.zeros(tmpl.shape[0])
    for t in range(tmpl.shape[0]):
        for m in range(ns):
            E[t] += tmpl[t, m].real**2 + tmpl[t, m].imag**2
    C = np.empty(tmpl.shape[0], np.complex128)
    alpha = np.full((P, H), -np.inf)
    h0 = 0
    for i in range(r + 1):
        h0 = h0*M + start
    alpha[0, h0] = 0.0
    new = np.empty((P, H))
    rot = np.exp(-2j*np.pi*np.arange(P)/P)
    tmp = np.empty(M)
    logqy = 0.0; logqya = 0.0; psi_true = 0
    for j in range(n):
        for t in range(tmpl.shape[0]):
            acc = 0j
            for m in range(ns):
                acc += np.conj(tmpl[t, m])*yw[j, m]
            C[t] = acc
        for psi2 in range(P):
            for hn in range(H):
                tail = hn//M; i0 = hn % M; sr = tail % M
                lp = logp[sr, i0]
                if lp == -np.inf:
                    new[psi2, hn] = -np.inf; continue
                mx = -np.inf
                for s0 in range(M):
                    psi = (psi2 - k*labels[s0]) % P
                    hist = s0*Mr + tail
                    prev = alpha[psi, hist]
                    if prev == -np.inf:
                        tmp[s0] = -np.inf; continue
                    t = hist*M + i0
                    v = prev + (2*amp*(rot[psi]*C[t]).real - amp*amp*E[t])/sig2 + lp
                    tmp[s0] = v
                    if v > mx:
                        mx = v
                if mx == -np.inf:
                    new[psi2, hn] = -np.inf
                else:
                    sm = 0.0
                    for s0 in range(M):
                        if tmp[s0] > -np.inf:
                            sm += np.exp(tmp[s0] - mx)
                    new[psi2, hn] = mx + np.log(sm)
        norm = -np.inf
        for psi2 in range(P):
            for hn in range(H):
                if new[psi2, hn] > norm:
                    norm = new[psi2, hn]
        for psi2 in range(P):
            for hn in range(H):
                alpha[psi2, hn] = new[psi2, hn] - norm
        logqy += norm
        inc[0, j] = norm
        hist = 0
        for i in range(r + 1):
            hist = hist*M + truth[j + i]
        t = hist*M + truth[j + r + 1]
        g = (2*amp*(rot[psi_true]*C[t]).real - amp*amp*E[t])/sig2
        logqya += g
        inc[1, j] = g
        psi_true = (psi_true + k*labels[truth[j]]) % P
    sm = 0.0
    for psi2 in range(P):
        for hn in range(H):
            if alpha[psi2, hn] > -np.inf:
                sm += np.exp(alpha[psi2, hn])
    return logqy + np.log(sm), logqya


def memory(design, ns, p_target, P=1, max_templates=20000, max_states=40000):
    """Symbols of history r needed for a predictor of order p_target, within template and state budgets.

    Trellis states are P*M^(r+1); templates M^(r+2). The predictor order is reduced if r is capped.
    """
    M = design.M
    r_cap = max(1, int(np.floor(np.log(max_templates)/np.log(M) + 1e-9)) - 2)
    while r_cap > 1 and P*M**(r_cap + 1) > max_states:
        r_cap -= 1
    r = min(max(1, int(np.ceil(p_target/ns))), r_cap)
    return r, min(p_target, r*ns)


def info_rate(design, noise_psd, amp, n_sym=8000, seed=1, fs_target=4.0, p_max=18, return_detail=False,
              return_density=False):
    """Information rate (bit/symbol) of `design` with Markov input in colored Gaussian noise.

    noise_psd(f): complex-baseband noise PSD (per Hz) for unit signal amplitude at offset f (Hz).
    amp: signal amplitude in the same normalization (only amp^2/N matters).
    p_max: predictor order (samples); the trellis carries enough symbol history for it,
    subject to a template budget (see `memory`).
    """
    ns = max(2, int(round(fs_target*design.T)))
    fs = ns/design.T
    labels, k, P = design.lattice()
    r, p = memory(design, ns, p_max, P)
    M = design.M
    hb = float(design.h)*design.beta/2
    trans = design.transition(); logp = np.where(trans > 0, np.log(np.where(trans > 0, trans, 1)), -np.inf)
    rng = np.random.default_rng(seed)
    start = M//2                                          # r+1 known start symbols, then the Markov chain
    seq = np.empty(n_sym + r + 1, np.int64); seq[:r + 1] = start
    cum = np.cumsum(trans, axis=1)
    uu = rng.random(n_sym)
    for j in range(n_sym):
        seq[j + r + 1] = min(np.searchsorted(cum[seq[j + r]], uu[j]), M - 1)
    # remove the constant offset h*beta/2*a_0 of the generator relative to the trellis convention
    x = np.exp(2j*np.pi*(phase_track(design, seq, ns) - hb*design.tones()[start]))
    L = len(x)
    f = np.fft.fftfreq(L, 1/fs)
    S = fs*noise_psd(f)                                   # discrete per-sample PSD
    w = np.fft.ifft(np.fft.fft((rng.normal(size=L) + 1j*rng.normal(size=L))/np.sqrt(2))*np.sqrt(S))
    y = amp*x + w
    a, sig2 = predictor(fs*noise_psd(np.fft.fftfreq(1 << 14, 1/fs)), p)
    yw = np.ascontiguousarray(whiten(y, a).reshape(n_sym + r + 1, ns)[r + 1:])
    # whitened unit-amplitude templates for (s_0..s_r, a_j); the phase at the start of symbol j-r
    # is psi + h (1 - beta/2) s_0, i.e. the generator phase minus h beta/2 s_0
    tup = np.indices((M,)*(r + 2)).reshape(r + 2, -1).T
    tmpl = np.empty((len(tup), ns), complex)
    tones = design.tones()
    for t, sy in enumerate(tup):
        ph = phase_track(design, sy, ns)[ns:] - hb*tones[sy[0]]
        tmpl[t] = whiten(np.exp(2j*np.pi*ph), a)[-ns:]
    inc = np.empty((2, n_sym))
    logqy, logqya = _forward(tmpl, yw, labels, k, P, M, r, logp, amp, sig2, start, seq, inc)
    rate = (logqya - logqy)/(n_sym*np.log(2))
    if return_density:
        # per-symbol information density increments (bits): log q(y_j|a, past) - log q(y_j|past)
        return rate, (inc[1] - inc[0])/np.log(2)
    if return_detail:
        return dict(rate=rate, ns=ns, fs=fs, p=p, r=r, sig2=sig2, P=P, entropy=design.entropy_rate())
    return rate
