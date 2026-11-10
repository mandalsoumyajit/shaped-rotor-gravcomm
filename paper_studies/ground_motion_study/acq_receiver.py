"""HopLimitedReceiver with a timing offset and an initial filter state, for simulated acquisition.

hop_receiver.py is frozen (its hash is recorded in the qualification configs), so the two changes live here:
  * time_offset: symbol t occupies [t*T + time_offset, (t+1)*T + time_offset) on the receiver's sample clock
    (absolute input sample index i <-> time i/fs; the data frame starts nominally at i = 0);
  * initial_state: SOS state of the unit-gain mean model just before the first sample of symbol 0, as computed
    from the known preamble and the estimated timing (default: steady centre tone, as in HopLimitedReceiver).
With time_offset = 0 and the default state, the class is identical to HopLimitedReceiver (_step is the same
code with `local` measured from the shifted symbol start; checked by acq_packets.py --selftest).
"""
import numpy as np
from scipy.signal import sosfilt
from hop_receiver import HopLimitedReceiver


class AcqHopReceiver(HopLimitedReceiver):
    def __init__(self, *args, time_offset=0.0, initial_state=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.t0 = float(time_offset)
        self.path_min = 0.0                       # running sum of the per-step minimum metric (diagnostic only)
        if initial_state is not None:
            self.state = np.asarray(initial_state, complex).reshape(len(self.sos), 1, 2)

    def _geometry(self, t):
        fs = self.receiver.input_rate_hz; D = self.receiver.decimation
        start = int(np.ceil((t*self.T + self.t0)*fs - 1e-12)); end = int(np.ceil(((t + 1)*self.T + self.t0)*fs - 1e-12))
        first = (-start) % D
        return start, end, first, len(range(first, end - start, D))

    def _step(self, y):
        # Identical to HopLimitedReceiver._step except that `local` is measured from t*T + t0.
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
        branch = branch[:, first::Dd]; response = response[:, first::Dd]
        bstate = bstate.transpose(1, 0, 2).reshape(B**2, dim); transition = transition.transpose(1, 0, 2).reshape(dim, dim)
        return self._step_core(t, y, count, branch, bstate, response, transition, dim)

    def _step_core(self, t, y, count, branch, bstate, response, transition, dim):
        B, L = self.B, self.L
        data = t < L*self.nbytes
        choices = np.arange(self.M) if data else np.array([self.zero])
        parent = np.repeat(np.arange(len(self.cost)), len(choices)); choice = np.tile(choices, len(self.cost))
        old = self.hist % B
        valid = np.abs(choice - old[parent]) <= self.D
        parent, choice = parent[valid], choice[valid]
        values = np.zeros(len(choice), np.int64)
        if data and self.table is not None:
            pos = t % L
            digits = [(self.hist[parent]//B**j) % B for j in range(pos + 1)]
            start_tone = digits[pos]
            code = np.zeros(len(choice), np.int64)
            for i in range(pos):
                code = code*B + digits[pos - 1 - i]
            code = code*B + choice
            values = self.table.lookup[pos][start_tone, code]
            ok = values >= 0
            parent, choice, values = parent[ok], choice[ok], values[ok]
        elif data:
            pos = t % L
            digits = [(self.hist[parent]//B**j) % B for j in range(pos + 1)]
            seq_prev = [digits[pos - i] for i in range(pos + 1)]
            offset = np.zeros(len(choice), np.int64)
            for i in range(pos):
                offset += self.code.G[L - 1 - i][seq_prev[i], digits[pos - 1 - i]]
            offset += self.code.G[L - 1 - pos][seq_prev[pos], choice]
            cnt = self.code.N[L - 1 - pos][choice]
            high = np.minimum(offset + cnt - 1, 2**self.bits - 1)
            stride = 2**self.pad if t//L == self.nbytes - 1 else 1
            ok = ((offset + stride - 1)//stride)*stride <= high
            parent, choice, values = parent[ok], choice[ok], offset[ok]
        index = old[parent]*B + choice; factor = np.exp(2j*np.pi*self.phase[parent]/self.phases)
        flat = self.state.transpose(1, 0, 2).reshape(len(self.cost), dim)
        means = (flat@response)[parent] + factor[:, None]*branch[index]
        Wp, Wc = self.whitener.block_matrices(self.sample, count)
        wy = Wp@self.obs_past + Wc@y
        prediction = (self.mean_past@Wp.T)[parent] + means@Wc.T
        edge = np.sum(abs(prediction - wy)**2, axis=1); scores = self.cost[parent] + edge
        next_hist = (self.hist[parent]*B + choice) % self.modulus
        next_phase = (self.phase[parent] + self.inc[old[parent], choice]) % self.phases
        keys = next_hist*self.phases + next_phase
        order = np.argsort(scores, kind='stable'); _, first_key = np.unique(keys[order], return_index=True)
        selection = order[np.sort(first_key)[:self.width]]; selected_keys = keys[selection]
        self.pruned_states += max(0, len(first_key) - self.width)
        sorting = np.argsort(selected_keys); sorted_keys = selected_keys[sorting]
        where = np.searchsorted(sorted_keys, keys); okk = where < len(selection)
        okk &= sorted_keys[np.minimum(where, len(selection) - 1)] == keys
        destinations = sorting[np.minimum(where, len(selection) - 1)]
        self.graph.append((t, self.cost.copy(), parent[okk], destinations[okk], edge[okk], values[okk]))
        self.cost = scores[selection]; self.path_min += min(self.cost); self.cost -= min(self.cost)
        updated = (flat@transition)[parent[selection]] + factor[selection, None]*bstate[index[selection]]
        self.state = updated.reshape(-1, len(self.sos), 2).transpose(1, 0, 2)
        self.mean_past = np.concatenate((self.mean_past[parent[selection]], means[selection]), axis=1)[:, -self.order:]
        self.obs_past = np.r_[self.obs_past, y][-self.order:]
        self.hist = next_hist[selection]; self.phase = next_phase[selection]
        self.symbol += 1; self.sample += count; self.peak_states = max(self.peak_states, len(selection))
        self.peak_graph = max(self.peak_graph, len(self.graph))
        outputs = []
        while self.next_byte < self.nbytes and (t >= L*self.next_byte + L - 1 + self.lag or self.symbol == self.total):
            target = L*self.next_byte + L - 1; back = np.zeros(len(self.cost)); llr = np.full(self.bits, np.nan); hard = None
            last = self.next_byte == self.nbytes - 1
            for step, forward, src, dst, metric, words in reversed(self.graph):
                if step == target:
                    total = forward[src] + metric + back[dst]; hard = int(words[np.argmin(total)])
                    for bit in range(self.bits):
                        if last and bit >= self.bits - self.pad:
                            llr[bit] = np.inf; continue
                        labels = (words >> (self.bits - 1 - bit)) & 1; zero = total[labels == 0]; one = total[labels == 1]
                        if len(zero) and len(one) and np.isfinite(zero.min()) and np.isfinite(one.min()):
                            llr[bit] = one.min() - zero.min()
                    break
                prior = np.full(len(forward), np.inf); np.minimum.at(prior, src, metric + back[dst])
                back = prior - np.min(prior)
            check = llr[:self.bits - self.pad] if (last and self.pad) else llr
            missing = int(np.count_nonzero(~np.isfinite(check))); self.missing += missing
            outputs.append(dict(byte_index=self.next_byte, byte=hard, llr=llr, available_after_symbol=self.symbol,
                                missing_bit_hypotheses=missing))
            self.next_byte += 1
            while self.graph and self.graph[0][0] < L*self.next_byte:
                self.graph.popleft()
        return outputs
