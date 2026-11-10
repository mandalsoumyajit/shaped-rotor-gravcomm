"""Packet simulation with the iterative (turbo) receiver, for the released CF-FSK and hop-limited waveforms.

Transmitter, channel and seeds follow run_packets.py / hop_packets.py exactly. For every packet the
released (non-iterative) result and the payload errors after each turbo iteration are recorded,
so both receivers are compared on identical noise. Writes results/<case>/<index>.json.
Usage: python turbo_packets.py --case turbo_cf_f24_a0.70:cf:24:0.70 --packets 64 --iterations 8
       (waveform 'cf' = released CF-FSK, 'hop' = hop-limited candidate; amplitude factor of A0)
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, binascii, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import mitigation as mg
from gravcomm.causal_receiver import CausalReceiver, sampled_cffsk
from gravcomm.streaming_soft import StreamingSoftReceiver
from gravcomm import convolutional_candidates as tb
import hop_packets as hp
from hop_receiver import HopLimitedReceiver, sampled_hop
import turbo_rx as tr
import exit_chart as ec

OUT = Path(__file__).resolve().parent/'results'
K, SCHEME = 256, 'tb_2/3'


def crc_ok(bits):
    data = np.packbits(bits)
    return bool(bytes(data[:2]) == b'\x00\x01' and
                binascii.crc_hqx(bytes(data[:-2]), 0xffff) == int.from_bytes(bytes(data[-2:]), 'big'))


def run(job):
    waveform, carrier, amplitude, index, iterations, seed_ns = job[:6]
    code_name, label_name = (job[6], job[7]) if len(job) > 6 else ('K7 (171,133) released', 'natural (released)')
    code = ec.CODES[code_name]
    if waveform == 'cf':
        label = ec.labelings()[label_name]
    elif label_name.startswith('natural'):
        label = np.arange(2**hp.BITS)
    else:                                                    # labeling designed by label_search.py
        label = np.array(json.loads((OUT/'labels_hop.json').read_text())['optimized']['label'])
    start = time.perf_counter()
    rng = np.random.default_rng(np.random.SeedSequence([seed_ns[0], K, index]))
    payload = rng.integers(0, 256, (K - 32)//8, dtype=np.uint8)
    body = b'\x00\x01' + bytes(payload)
    info = np.unpackbits(np.frombuffer(body + binascii.crc_hqx(body, 0xffff).to_bytes(2, 'big'), dtype=np.uint8))
    word = code.encode(info)
    permutation = np.random.default_rng(4000 + len(word)).permutation(len(word))
    r = CausalReceiver(carrier_hz=carrier, background=mg.background('baseline', carrier))
    if waveform == 'cf':
        packed = np.argsort(label)[np.packbits(word[permutation]).astype(int)]     # byte value v -> word w with label[w] = v
        symbols = np.column_stack((packed//49, packed//7 % 7, packed % 7)).ravel() - 3
        ts = 224/(3*len(packed) + 18)
        tx = sampled_cffsk(np.r_[symbols, np.zeros(6, int)], ts, q=.08, transition_fraction=.6)
        nsym = len(symbols) + 6; L, nbits, ngroups, pad = 3, 8, len(packed), 0
        det = StreamingSoftReceiver(r, ts, amplitude, ngroups, history_symbols=4, beam_width=32768,
                                    lookahead_symbols=12, q=.08, transition_fraction=.6)
    else:
        bits = word[permutation]
        pad = (-len(bits)) % hp.BITS
        values = np.r_[bits, np.zeros(pad, int)].reshape(-1, hp.BITS) @ (1 << np.arange(hp.BITS - 1, -1, -1))
        ranks = np.argsort(label)[values]                    # bit pattern v -> word (rank) w with label[w] = v
        tones, prev = [], hp.M//2
        for w in ranks:
            s = hp.CODE.encode(int(w), prev); tones += s; prev = s[-1]
        ngroups = len(ranks)
        ts = 224/(len(tones) + 18)
        txs = np.r_[tones, np.full(6, hp.M//2)]
        tx = sampled_hop(txs, ts, hp.H, float(hp.BETA), hp.M)
        nsym = len(txs); L, nbits = hp.L, hp.BITS
        # pad bits are known zeros in label space; the turbo receiver uses them as strong priors, so the
        # detector does not prune on them (pad_bits=0)
        det = HopLimitedReceiver(r, ts, amplitude, ngroups, history_symbols=4, beam_width=32768, lookahead_symbols=12,
                                 h=hp.H, beta=hp.BETA, tone_count=hp.M, max_hop=hp.D, group_symbols=L,
                                 group_bits=nbits, pad_bits=0)
    mean = r.stream(1).process(tx)
    noise_rng = np.random.default_rng(np.random.SeedSequence([seed_ns[1], K, index]))
    observed = mean + r.stationary_noise(len(mean), noise_rng, amplitude)
    tr.record_graph(det)
    outputs = []; position = 0
    for symbol in range(nsym):
        count = det._geometry(symbol)[3]
        outputs.extend(det.push(observed[position:position + count])); position += count
    stats = det.finish()
    t_detect = time.perf_counter() - start
    # released (non-iterative) receiver on the same detector run
    channel = np.concatenate([x['llr'] for x in outputs])[:len(word)]
    llr = np.empty(len(word)); llr[permutation] = channel
    released = code_name.endswith('released') and label_name.startswith('natural')
    base = tb.decode(llr, SCHEME[3:], K) if released else info*0
    base_errors = int(np.sum(base[16:-16] != info[16:-16])) if released else -1
    # iterative receiver
    packed_graph = tr.pack_graph(det, L, ngroups)
    labels = tr.bit_labels(2**nbits, nbits)[label]
    hist = tr.turbo(packed_graph, permutation, ngroups, nbits, len(word), iterations=iterations, pad_bits=pad,
                    labels=labels, siso=lambda l: code.siso(l, K))
    errors = [int(np.sum(h[16:-16] != info[16:-16])) for h in hist]
    accepted = [crc_ok(h) for h in hist]
    return dict(code=code_name, labeling=label_name, waveform=waveform, carrier_hz=carrier, amplitude_m_s2=amplitude, index=index,
                released_errors=base_errors, released_accepted=crc_ok(base), turbo_errors=errors,
                turbo_accepted=accepted, pruned_states=stats['pruned_states'],
                missing=stats['missing_bit_hypotheses'], graph_edges=int(len(packed_graph['met'])),
                detect_s=t_detect, elapsed_s=time.perf_counter() - start)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', action='append', required=True, help='name:waveform:carrier:amplitude_factor_of_A0')
    p.add_argument('--packets', type=int, default=64)
    p.add_argument('--seed-start', type=int, default=7000000)
    p.add_argument('--seeds', default='20371103,20371104', help='payload,noise seed namespaces')
    p.add_argument('--iterations', type=int, default=8)
    p.add_argument('--workers', type=int, default=16)
    a = p.parse_args()
    ns = tuple(int(x) for x in a.seeds.split(','))
    jobs = []
    for case in a.case:
        parts = case.split(':'); name, wf, carrier, factor = parts[:4]
        extra = tuple(parts[4:6]) if len(parts) > 4 else ()
        folder = OUT/name; folder.mkdir(parents=True, exist_ok=True)
        for i in range(a.seed_start, a.seed_start + a.packets):
            if not (folder/f'{i}.json').exists():
                jobs.append(((wf, float(carrier), float(factor)*mg.A0, i, a.iterations, ns) + extra, folder))
    print(f'{len(jobs)} packets to run', flush=True)
    t0 = time.time(); done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures = {pool.submit(run, j): folder for j, folder in jobs}
        for f in as_completed(futures):
            res = f.result(); folder = futures[f]
            tmp = folder/f"{res['index']}.json.tmp"
            tmp.write_text(json.dumps(res) + '\n'); os.replace(tmp, folder/f"{res['index']}.json")
            done += 1
            if done % 32 == 0:
                print(f'  {done}/{len(jobs)} ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
