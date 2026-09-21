"""Tail-biting K=7 convolutional reference with published puncturing.

Mother polynomials and masks: ETSI EN 300 744 V1.6.1, section 4.3.3,
Table 2: https://www.etsi.org/deliver/etsi_en/300700_300799/300744/01.06.01_60/en_300744v010601p.pdf
Tail biting is our packet-boundary choice, not the DVB framing standard.
Positive input LLR favors bit zero. Decoding is exact for the sum of supplied
bit metrics; correlated detector metrics need not be true independent LLRs.
"""
import numpy as np
from numba import njit

PATTERNS = {
    '1/2': np.array([[1, 1]], dtype=bool),
    '2/3': np.array([[1, 1], [0, 1]], dtype=bool),
    '3/4': np.array([[1, 1], [0, 1], [1, 0]], dtype=bool),
}
OUTPUTS = np.array([[[((s << 1 | b) & g).bit_count() % 2
                     for g in (0o171, 0o133)] for b in range(2)]
                    for s in range(64)], dtype=np.uint8)


def _mask(information_bits, rate):
    if rate not in PATTERNS:
        raise ValueError('rate must be 1/2, 2/3 or 3/4')
    if not isinstance(information_bits, (int, np.integer)) or information_bits < 1:
        raise ValueError('information_bits must be a positive integer')
    p = PATTERNS[rate]
    return p[np.arange(information_bits) % len(p)]


def metadata(information_bits, rate='1/2'):
    mask = _mask(information_bits, rate)
    n = int(mask.sum())
    return dict(family='tail_biting_convolutional', constraint_length=7,
                generators_octal=['171', '133'], nominal_rate=rate,
                information_bits=int(information_bits), coded_bits=n, n=n,
                effective_rate=information_bits/n, termination_bits=0,
                puncturing_period=len(PATTERNS[rate]),
                puncturing=PATTERNS[rate].astype(int).tolist(),
                decoder='exact_64_start_state_viterbi', crc_bits=0)


def encode(bits, rate='1/2'):
    original = np.asarray(bits)
    if original.ndim != 1 or not np.all((original == 0) | (original == 1)):
        raise ValueError('binary vector required')
    bits = original.astype(np.uint8)
    mask = _mask(len(bits), rate)
    # Most recent history at LSB; cyclic history also handles blocks <6 bits.
    state = sum(int(bits[(-j) % len(bits)]) << (j-1) for j in range(1, 7))
    initial = state
    mother = np.empty((len(bits), 2), dtype=np.uint8)
    for t, bit in enumerate(bits):
        mother[t] = OUTPUTS[state, bit]
        state = ((state << 1) | int(bit)) & 63
    assert state == initial
    return mother[mask]


@njit(cache=True, nogil=True)
def _decode_exact(metrics, outputs):
    length = len(metrics)
    # Precompute once, since 64 candidate starting states use the same metrics.
    branch = np.empty((length, 64, 2))
    for t in range(length):
        for s in range(64):
            for b in range(2):
                branch[t, s, b] = (outputs[s,b,0]*metrics[t,0]
                                   + outputs[s,b,1]*metrics[t,1])
    best = np.inf
    best_start = 0
    for initial in range(64):
        costs = np.full(64, np.inf)
        costs[initial] = 0.
        new = np.empty(64)
        for t in range(length):
            for dest in range(64):
                b = dest & 1
                p = dest >> 1
                v0 = costs[p] + branch[t,p,b]
                v1 = costs[p+32] + branch[t,p+32,b]
                new[dest] = min(v0, v1)
            costs, new = new, costs
        if costs[initial] < best:
            best = costs[initial]
            best_start = initial
    # Only the winning initial state needs a traceback allocation.
    costs = np.full(64, np.inf)
    costs[best_start] = 0.
    new = np.empty(64)
    parents = np.empty((length, 64), np.uint8)
    for t in range(length):
        for dest in range(64):
            b = dest & 1
            p = dest >> 1
            v0 = costs[p] + branch[t,p,b]
            v1 = costs[p+32] + branch[t,p+32,b]
            if v1 < v0:
                new[dest] = v1
                parents[t,dest] = p+32
            else:
                new[dest] = v0
                parents[t,dest] = p
        costs, new = new, costs
    state = best_start
    bits = np.empty(length, np.uint8)
    for t in range(length-1, -1, -1):
        bits[t] = state & 1
        state = parents[t,state]
    return bits


def decode(llr, rate='1/2', information_bits=None):
    """Exact tail-biting minimum additive-metric decoding; no CRC rejection.

    information_bits is required: punctured lengths are not always invertible.
    Puncture phase restarts at every packet, and partial periods are allowed.
    """
    mask = _mask(information_bits, rate)
    llr = np.asarray(llr, dtype=np.float64)
    if llr.ndim != 1 or len(llr) != int(mask.sum()) or not np.all(np.isfinite(llr)):
        raise ValueError('finite one-dimensional LLR vector of coded length required')
    expanded = np.zeros(mask.shape)
    expanded[mask] = llr  # Zero metric for punctured bits.
    return _decode_exact(expanded, OUTPUTS)
