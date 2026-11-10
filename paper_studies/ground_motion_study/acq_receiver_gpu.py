"""GPU twin of acq_receiver.AcqHopReceiver: HopLimitedReceiverGPU with a timing offset and an initial filter state.

_step is HopLimitedReceiverGPU._step (hop_receiver_gpu.py, frozen) with `local` measured from t*T + t0;
_geometry is AcqHopReceiver._geometry. Generated from hop_receiver_gpu.py by a one-line substitution.
"""
import numpy as np
import cupy as cp
from scipy.signal import sosfilt
from hop_receiver_gpu import HopLimitedReceiverGPU
from acq_receiver import AcqHopReceiver


class AcqHopReceiverGPU(HopLimitedReceiverGPU):
    _geometry = AcqHopReceiver._geometry

    def __init__(self, *args, time_offset=0.0, initial_state=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.t0 = float(time_offset)
        if initial_state is not None:
            self.state = cp.asarray(np.asarray(initial_state, complex).reshape(len(self.sos), 1, 2))

    def _step(self, y):
        t = self.symbol; start, end, first, count = self._geometry(t)
        fs = self.receiver.input_rate_hz; Dd = self.receiver.decimation; B, L = self.B, self.L
        local = np.arange(start, end)/fs - t*self.T - self.t0; tr = self.r*self.T; u = np.minimum(local/tr, 1)
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

