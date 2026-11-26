"""Single-mode response to quintic transition sequences (Appendix A, modes and modulation coupling).

Angular acceleration is a train of quintic hops (beta=1, T_s=1.179 s); each mode is a damped
oscillator. Prints the peak response relative to the static (quasistatic) peak for an isolated
hop, alternating hops, a repeating cycle and random hop-limited data (D<=3).
"""
import numpy as np
from scipy.signal import cont2discrete, dlsim
T = 224/190; fs = 4000.0; dt = 1/fs
tau = np.arange(0, 1, dt/T)
pd = 30*tau**2*(1-tau)**2            # d p / d tau for the quintic; acceleration per unit hop
def accel(hops):
    return np.concatenate([h*pd for h in hops])
def amp(a, fn, z):
    wn = 2*np.pi*fn
    A = np.array([[0,1],[-wn**2,-2*z*wn]]); B = np.array([[0],[1]]); C = np.array([[1,0]]); D = np.array([[0]])
    Ad,Bd,Cd,Dd,_ = cont2discrete((A,B,C,D), dt, method='zoh')
    _, y, _ = dlsim((Ad,Bd,Cd,Dd,dt), a[:,None])
    qs = a/wn**2
    return np.max(np.abs(y[:,0]))/np.max(np.abs(qs)), np.max(np.abs(y[:,0]-qs))/np.max(np.abs(qs))
rng = np.random.default_rng(1)
def rand_hops(n):
    s = 3; out = []
    for _ in range(n):
        nxt = rng.integers(max(0,s-3), min(6,s+3)+1); out.append(nxt-s); s = nxt
    return out
cases = {'isolated': [3]+[0]*5, 'alternating +-3': [3,-3]*60, 'alternating +-1': [1,-1]*60,
         'random D<=3': rand_hops(400), 'up-down cycle 3,3,-3,-3': [3,3,-3,-3]*30}
for fn in (27.39, 35.49, 112.29):
    for z in (0.005, 0.02, 0.05):
        r = {k: amp(np.array(accel(v), float), fn, z) for k, v in cases.items()}
        print(f"fn={fn} z={z}: " + "; ".join(f"{k}: peak/static-1={(a-1)*100:+.3f}%  maxdev={b*100:.3f}%" for k,(a,b) in r.items()))
