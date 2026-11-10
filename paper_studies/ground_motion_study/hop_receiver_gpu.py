"""GPU (CuPy) version of the hop-limited streaming detector, for graph-recording receivers.

HopLimitedReceiverGPU reproduces HopLimitedReceiver._step in double precision on the GPU: path
expansion, enumerative-word bookkeeping, filter-state propagation, whitened edge metrics, and
reduced-state merging (lowest score per (history, phase) key, ties to the lower candidate index,
survivors ordered by score, beam truncation) are the same operations. The small per-symbol pieces
(branch waveforms and their filter responses, whitening matrices) are computed on the CPU as before
and copied over. Differences from the CPU detector:
  * it records the detector graph on the GPU (pack_graph_gpu flattens it for the word-level
    forward-backward of nb_packets / turbo_rx) and does not compute the streaming fixed-lag
    hard/LLR outputs, which the graph-based receivers do not use;
  * floating-point summation order differs, so metrics agree to ~1e-12 relative, not bit for bit.
Validation against the CPU detector on identical packets: gpu_check.py.
"""
import numpy as np
import cupy as cp
from scipy.signal import sosfilt
from hop_receiver import HopLimitedReceiver


class HopLimitedReceiverGPU(HopLimitedReceiver):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.state = cp.asarray(self.state)
        self.hist = cp.asarray(self.hist); self.phase = cp.asarray(self.phase)
        self.cost = cp.asarray(self.cost); self.mean_past = cp.asarray(self.mean_past)
        self.inc_g = cp.asarray(self.inc)
        self.G_g = [cp.asarray(g) for g in self.code.G]; self.N_g = [cp.asarray(n) for n in self.code.N]
        self.lookup_g = None if self.table is None else [cp.asarray(l) for l in self.table.lookup]
        self.steps = []                                       # (t, nodes, src, dst, edge, word) on the GPU

    def _step(self, y):
        t = self.symbol; start, end, first, count = self._geometry(t)
        fs = self.receiver.input_rate_hz; Dd = self.receiver.decimation; B, L = self.B, self.L
        local = np.arange(start, end)/fs - t*self.T; tr = self.r*self.T; u = np.minimum(local/tr, 1)
        integral = tr*(2.5*u**4 - 3*u**5 + u**6) + np.maximum(local - tr, 0)
        spacing = float(self.h)/self.T
        before = np.repeat(self.offsets, B); after = np.tile(self.offsets, B)
        raw = np.exp(2j*np.pi*spacing*(before[:, None]*local + (after - before)[:, None]*integral))
        branch, bstate = sosfilt(self.sos, raw, axis=1, zi=np.zeros((len(self.sos), B**2, 2), complex))
        dim = 2*len(self.sos); initial = np.eye(dim, dtype=complex).reshape(dim, len(self.sos), 2).transpose(1, 0, 2)
        response, transition = sosfilt(self.sos, np.zeros((dim, len(local)), complex), axis=1, zi=initial)
        branch = cp.asarray(branch[:, first::Dd]); response = cp.asarray(response[:, first::Dd])
        bstate = cp.asarray(bstate.transpose(1, 0, 2).reshape(B**2, dim))
        transition = cp.asarray(transition.transpose(1, 0, 2).reshape(dim, dim))
        data = t < L*self.nbytes
        S = int(self.cost.shape[0])
        choices = cp.arange(self.M) if data else cp.array([self.zero])
        parent = cp.repeat(cp.arange(S), len(choices)); choice = cp.tile(choices, S)
        old = self.hist % B
        valid = cp.abs(choice - old[parent]) <= self.D
        parent, choice = parent[valid], choice[valid]
        values = cp.zeros(len(choice), cp.int64)
        if data and self.lookup_g is not None:
            pos = t % L
            hpar = self.hist[parent]
            digits = [(hpar//B**j) % B for j in range(pos + 1)]
            start_tone = digits[pos]
            code = cp.zeros(len(choice), cp.int64)
            for i in range(pos):
                code = code*B + digits[pos - 1 - i]
            code = code*B + choice
            values = self.lookup_g[pos][start_tone, code]
            ok = values >= 0
            parent, choice, values = parent[ok], choice[ok], values[ok]
        elif data:
            pos = t % L
            hpar = self.hist[parent]
            digits = [(hpar//B**j) % B for j in range(pos + 1)]
            seq_prev = [digits[pos - i] for i in range(pos + 1)]
            offset = cp.zeros(len(choice), cp.int64)
            for i in range(pos):
                offset += self.G_g[L - 1 - i][seq_prev[i], digits[pos - 1 - i]]
            offset += self.G_g[L - 1 - pos][seq_prev[pos], choice]
            cnt = self.N_g[L - 1 - pos][choice]
            high = cp.minimum(offset + cnt - 1, 2**self.bits - 1)
            stride = 2**self.pad if t//L == self.nbytes - 1 else 1
            ok = ((offset + stride - 1)//stride)*stride <= high
            parent, choice, values = parent[ok], choice[ok], offset[ok]
        index = old[parent]*B + choice
        factor = cp.exp(2j*np.pi*self.phase[parent]/self.phases)
        flat = self.state.transpose(1, 0, 2).reshape(S, dim)
        means = (flat@response)[parent] + factor[:, None]*branch[index]
        Wp, Wc = self.whitener.block_matrices(self.sample, count)
        wy = cp.asarray(Wp@self.obs_past + Wc@y)
        Wp = cp.asarray(Wp); Wc = cp.asarray(Wc)
        prediction = (self.mean_past@Wp.T)[parent] + means@Wc.T
        diff = prediction - wy
        edge = cp.sum(diff.real**2 + diff.imag**2, axis=1)
        scores = self.cost[parent] + edge
        next_hist = (self.hist[parent]*B + choice) % self.modulus
        next_phase = (self.phase[parent] + self.inc_g[old[parent], choice]) % self.phases
        keys = next_hist*self.phases + next_phase
        # numpy: order = argsort(scores, stable); first occurrence of each key in that order; survivors by rank
        o1 = cp.argsort(scores)                               # CuPy sorts are stable (radix / merge)
        o = o1[cp.argsort(keys[o1])]                          # by key, then score, then index
        ks = keys[o]
        isfirst = cp.empty(len(o), cp.bool_); isfirst[0] = True; isfirst[1:] = ks[1:] != ks[:-1]
        reps = o[isfirst]
        rank = cp.empty(len(o1), cp.int64); rank[o1] = cp.arange(len(o1))
        reps = reps[cp.argsort(rank[reps])]                   # survivors in score order
        nuniq = int(reps.shape[0])
        selection = reps[:self.width]
        self.pruned_states += max(0, nuniq - self.width)
        selected_keys = keys[selection]
        sorting = cp.argsort(selected_keys); sorted_keys = selected_keys[sorting]
        nsel = int(selection.shape[0])
        where = cp.searchsorted(sorted_keys, keys); okk = where < nsel
        wc = cp.minimum(where, nsel - 1)
        okk &= sorted_keys[wc] == keys
        destinations = sorting[wc]
        self.steps.append((t, S, parent[okk], destinations[okk], edge[okk], values[okk]))
        self.cost = scores[selection]; self.cost = self.cost - self.cost.min()
        psel = parent[selection]
        updated = (flat@transition)[psel] + factor[selection, None]*bstate[index[selection]]
        self.state = updated.reshape(-1, len(self.sos), 2).transpose(1, 0, 2)
        self.mean_past = cp.ascontiguousarray(
            cp.concatenate((self.mean_past[psel], means[selection]), axis=1)[:, -self.order:])
        self.obs_past = np.r_[self.obs_past, y][-self.order:]
        self.hist = next_hist[selection]; self.phase = next_phase[selection]
        self.symbol += 1; self.sample += count; self.peak_states = max(self.peak_states, nsel)
        self.peak_graph = len(self.steps)
        return []


def pack_graph_gpu(det, group_symbols, ngroups):
    """Same arrays as turbo_rx.pack_graph, from the GPU-recorded graph (returned as numpy)."""
    steps = det.steps
    nodes = [s[1] for s in steps] + [int(det.cost.shape[0])]
    lens = [int(s[2].shape[0]) for s in steps]
    eoff = np.cumsum([0] + lens)
    src = cp.concatenate([s[2] for s in steps]).astype(cp.int64).get()
    dst = cp.concatenate([s[3] for s in steps]).astype(cp.int64).get()
    met = (-cp.concatenate([s[4] for s in steps]).astype(cp.float64)).get()
    word = cp.concatenate([s[5] for s in steps]).astype(cp.int64).get()
    group = np.array([s[0]//group_symbols if (s[0] % group_symbols == group_symbols - 1 and
                                              s[0] < group_symbols*ngroups) else -1 for s in steps], np.int64)
    noff = np.cumsum([0] + nodes)
    return dict(src=src, dst=dst, met=met, word=word, eoff=eoff, noff=noff, group=group)


class GraphGPU:
    """Word-level forward-backward on the GPU-recorded detector graph (same recursion as
    nb_packets._graph_word_app: log domain, per-step max normalization, word log-priors on the
    word-final steps). word_app(Lp) returns the (G, W) word log-APPs as numpy."""

    def __init__(self, det, group_symbols, ngroups, nwords=512):
        import cupyx
        self._smax, self._sadd = cupyx.scatter_max, cupyx.scatter_add
        self.steps = det.steps; self.W = nwords; self.G = ngroups
        self.nodes = [s[1] for s in self.steps] + [int(det.cost.shape[0])]
        self.src = [s[2] for s in self.steps]; self.dst = [s[3] for s in self.steps]
        self.met = [-s[4] for s in self.steps]; self.word = [s[5] for s in self.steps]
        self.group = [s[0]//group_symbols if (s[0] % group_symbols == group_symbols - 1 and
                                             s[0] < group_symbols*ngroups) else -1 for s in self.steps]

    def _lse_into(self, n, idx, v):
        m = cp.full(n, -cp.inf); self._smax(m, idx, v)
        s = cp.zeros(n); self._sadd(s, idx, cp.exp(v - m[idx]))
        with np.errstate(divide='ignore'):
            return m + cp.log(s)

    def word_app(self, Lp):
        Lp = cp.asarray(Lp, dtype=cp.float64)
        T = len(self.steps)
        pri = [Lp[g][self.word[t]] if g >= 0 else None for t, g in enumerate(self.group)]
        edge = [self.met[t] + pri[t] if pri[t] is not None else self.met[t] for t in range(T)]
        alpha = [cp.zeros(1)]
        for t in range(T):
            a = self._lse_into(self.nodes[t + 1], self.dst[t], alpha[t][self.src[t]] + edge[t])
            alpha.append(a - a.max())
        beta = cp.zeros(self.nodes[T]); post = cp.full((self.G, self.W), -cp.inf)
        for t in range(T - 1, -1, -1):
            if self.group[t] >= 0:
                v = alpha[t][self.src[t]] + edge[t] + beta[self.dst[t]]
                post[self.group[t]] = self._lse_into(self.W, self.word[t], v)
            b = self._lse_into(self.nodes[t], self.src[t], edge[t] + beta[self.dst[t]])
            beta = b - b.max()
        return post.get()
