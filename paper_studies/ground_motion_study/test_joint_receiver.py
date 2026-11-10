"""Fast checks of the joint (GF(512) LDPC + hop-limited CF-FSK) link components, acquisition and the
Brack et al. field comparison. Run with pytest from this folder or from the repository root
(tests/test_joint_receiver.py adds this folder to sys.path)."""
import json, math
from fractions import Fraction
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip('numba')
HERE = Path(__file__).resolve().parent


def test_gf512_field_and_code():
    import nb_ldpc as nb
    a = np.arange(1, nb.Q)
    assert all(nb.gmul(int(x), nb.ginv(int(x))) == 1 for x in a[::37])
    assert len(set(nb.EXP[:nb.Q - 1].tolist())) == nb.Q - 1          # x^9 + x^4 + 1 is primitive
    code = nb.make_code(43, 29, 2, 1)
    assert (code.n, code.k) == (43, 29)
    assert np.all(np.bincount(code.edges[:, 1], minlength=43) == 2)  # variable degree 2
    rng = np.random.default_rng(1)
    for _ in range(5):
        x = code.encode(rng.integers(0, nb.Q, code.k))
        assert code.syndrome_ok(x)
        prior = np.full((code.n, nb.Q), 0.2/(nb.Q - 1)); prior[np.arange(code.n), x] = 0.8
        _, _, hard, ok = code.decode(prior, 30)
        assert ok and np.array_equal(hard, x)


def test_hop_code_roundtrip_and_hop_limit():
    import hop_packets as hp
    rng = np.random.default_rng(2)
    for w in rng.integers(0, 2**hp.BITS, 200):
        for start in range(hp.M):
            s = hp.CODE.encode(int(w), start)
            assert len(s) == hp.L and hp.CODE.decode(s, start) == w
            assert np.all(np.abs(np.diff([start] + s)) <= hp.D)


def test_preamble_model_matches_sampled_hop():
    import acquisition as aq, hop_packets as hp
    from hop_receiver import sampled_hop
    tones = aq.PREAMBLES['maxhop'] + [5, 2, 4, 6, 3, 3]
    x = sampled_hop(tones, aq.T, hp.H, float(hp.BETA), hp.M)
    t = np.arange(len(x))/16
    assert np.max(abs(x - np.exp(2j*np.pi*aq.phase_at(tones, t)))) < 1e-12
    assert abs(aq.phase_at(aq.PREAMBLES['maxhop'], [12*aq.T])[0]) < 1e-12         # zero net preamble phase


def _short_frame(tau, groups=3):
    import hop_packets as hp, acq_packets as ap, acquisition as aq
    import mitigation as mg
    from gravcomm.causal_receiver import CausalReceiver
    from scipy.signal import sosfilt, sosfilt_zi
    rng = np.random.default_rng(3)
    tones, prev = [], hp.M//2
    for w in rng.integers(0, 2**hp.BITS, groups):
        s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
    r = CausalReceiver(carrier_hz=24.0, background=mg.background('baseline', 24.0))
    D, T = r.decimation, aq.T
    combined = aq.PREAMBLES['maxhop'] + tones + [hp.M//2]*6
    i_first = -int(math.ceil((12*T + 3)*16/D))*D
    t = np.arange(i_first, int(math.ceil((len(tones) + 7)*T*16)))/16
    y = r.stream(1.0).process(np.exp(2j*np.pi*aq.phase_at(combined, t + 12*T - tau)))
    start0 = int(np.ceil(tau*16 - 1e-12))
    ii = np.arange(i_first, start0)
    _, zf = sosfilt(r.sos, np.exp(2j*np.pi*aq.phase_at(aq.PREAMBLES['maxhop'], ii/16 + 12*T - tau)),
                    zi=sosfilt_zi(r.sos).astype(complex))
    return r, y, i_first, start0, zf, groups, T


def test_acq_receiver_offset_and_state_exact():
    import hop_packets as hp
    from acq_receiver import AcqHopReceiver
    tau = 0.137
    r, y, i_first, start0, zf, groups, T = _short_frame(tau)
    det = AcqHopReceiver(r, T, 0.55*5.180259430854157e-10, groups, time_offset=tau, initial_state=zf,
                         history_symbols=4, beam_width=4096, lookahead_symbols=12, h=hp.H, beta=hp.BETA,
                         tone_count=hp.M, max_hop=hp.D, group_symbols=hp.L, group_bits=hp.BITS, pad_bits=0)
    pos = (start0 + (-start0) % r.decimation - i_first)//r.decimation
    for s in range(det.total):
        c = det._geometry(s)[3]; det.push(y[pos:pos + c]); pos += c
    det.finish()
    assert det.path_min < 1e-12                                      # noise-free: true path has zero metric


def test_brack_mode_shape():
    from brack_comparison import mode_shape, gauss
    x, w = gauss(200, 0, 1.0)
    U = mode_shape(x)
    assert abs(U[0] - 1) < 2e-3 and abs(mode_shape(np.array([0.0]))[0] - 1) < 1e-12
    assert abs(w @ U) < 1e-10                                        # orthogonal to rigid translation
    assert abs(w @ (U*(x - 0.5))) < 1e-10                            # and to rigid rotation


def test_recorded_results_consistent():
    s = json.loads((HERE/'results'/'sizing_nb.json').read_text())
    for name, row in s.items():
        assert abs(row['mass_kg'] - 7.91*row['scale']**3) < 0.01*row['mass_kg']
    q = json.loads((HERE/'results'/'qualification_nb.json').read_text())
    for case, row in q.items():
        assert row['status'] == 'pass' and row['result']['interval'][1] <= 1e-3
    b = json.loads((HERE/'results'/'brack_comparison.json').read_text())
    assert abs(b['at_296mm']['power_W']/3.75e-18 - 1) < 0.02
