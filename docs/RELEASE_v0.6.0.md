# Release v0.6.0 - joint GF(512) receiver, acquisition and revision studies

This release adds the code and records behind the revised manuscript (Physica Scripta, PHYSSCR-157981).
The v0.5.0 code, results and data archive are unchanged and remain the provenance of the conventional
(CF-FSK, binary convolutional code) design.

## Selected designs

| Distance | Bare mass | Diameter | Carrier | Peak inertial torque | Qualified amplitude |
|---:|---:|---:|---:|---:|---:|
| 0.5 m | 5.79 kg | 0.451 m | 24 Hz | 0.51 N m | 0.580 A0 |
| 1 m | 30.98 kg | 0.788 m | 24 Hz | 8.38 N m | 0.580 A0 |
| 2 m | 184.05 kg | 1.427 m | 18 Hz | 163.2 N m | 0.697 A0 |

Waveform: seven tones, no dwell, hops of at most three tone spacings (hop-limited CF-FSK). Code: (43,29) LDPC over
GF(512) defined on four-tone words, decoded iteratively with the streaming sequence detector (BCJR on the detector
graph, Walsh-Hadamard belief propagation, at most 16 passes). Frame: 224 payload bits, 172 data + 6 tail + 12
preamble symbols, 224 s.

## New material (`paper_studies/`)

`paper_studies/ground_motion_study/` (run from that folder; it puts `src/` and `scripts/` on the path):

- Information rates, bounds and the Pareto search: `cpm_rate.py`, `cpm_optimize.py`, `torque_bound.py`,
  `waterfill_tones.py`, `hop_threshold.py`, `finite_length.py`, `info_density.py`, `gap_analysis.py`,
  `plot_pareto.py`, `plot_optimality.py`.
- Joint receiver: `hop_receiver.py` and `hop_receiver_gpu.py` (waveform, tone mapping, CPU and CuPy detectors;
  frozen, hashed in the qualification configurations), `word_tables.py`, `nb_ldpc.py` (GF(512) code and decoder),
  `nb_packets.py` (packet simulation), `nb_design.py`, `nb_sched.py`, `gpu_check.py`.
- Packet waterfalls and fits: `nb_waterfall_fit.py`; qualification: `qualify_nb.py`; checks, replays and timing:
  `validate_nb.py`, `replay_wide_nb.py`, `timing_nb.py`; drive power: `drive_hop.py`; confinement and capacity:
  `link_limits.py`; sizing: `sizing.py`.
- Superseded joint design with a binary code (kept for the comparison in App. C): `qualify.py`, `hop_packets.py`,
  `turbo_rx.py`, `turbo_packets.py`, `exit_*.py`, `label_*.py`, `llr_capture.py`, `word_bcjr.py`, `replay_wide.py`,
  `validate_hop.py`.
- Ground-motion mitigation: `mitigation.py`, `analysis.py`, `rate_scaling.py`, `run_packets.py`, `summarize.py`,
  `check_weight.py`, `plot_mitigation.py`.
- Acquisition (manuscript App. C.4): `acquisition.py` (Cramer-Rao bounds, maximum-likelihood preamble estimator,
  Monte Carlo), `acq_receiver.py`/`acq_receiver_gpu.py` (detector with timing offset and initial state),
  `acq_packets.py` (packet runs: cold start, code-aided resynchronization, steady-state tracking, paired
  known-parameter reference, sequential qualification), `acq_diagnose.py`, `acq_loss.py`.
- Experimental comparison (App. A.4): `brack_comparison.py` reproduces the rotating-bar experiment of Brack et al.,
  Commun. Phys. 6, 270 (2023); `brack_receiver.py` compares their beam detector with the projected receiver.
- Checks: `test_joint_receiver.py` (also run by `tests/test_paper_studies.py`).

`paper_studies/figure_scripts/`: revised figures. `fea/*_hop.py`, `fea/run_nb_pipeline.sh` and
`scripts/magnetic_wall_{hop,nb}_2026_11_05.py`: finite-element and magnetic-link reruns for the new rotors,
with summaries in `fea/results/*_2026_11_05_nb/` and `results/magnetic_wall_*_2026_11_05_*`.

## Key results reproduced

- Qualification with known synchronization: both carriers pass the betting confidence sequence (alpha = 0.025 per
  channel) at 0.580 and 0.697 A0, 4 rejected frames per channel, accepted payload rate 0.9995 bit/s.
- Acquisition: the 12-symbol preamble alone reaches the Cramer-Rao bound (6 deg, 0.07 s, 11%) but costs
  0.75 dB near threshold on both carriers; with tracking from the previous CRC-verified frame both links pass a fresh qualification
  at the same amplitudes, and paired runs give losses of 0.04 dB (24 Hz) and 0.05 dB (18 Hz).
- Brack et al. comparison: absorbed power 3.78e-18 W against their 3.75e-18 W; distance exponents within 0.03.

## Data

Summary JSON and logs of `paper_studies/ground_motion_study/results/` (files up to 2 MB) are in Git. All packet
records, including the qualification, waterfall, acquisition and replay runs, are in
**gravcomm-v0.6.0-data-supplement.zip** (release asset, also in `datasets/`), with
`DATA_SUPPLEMENT_MANIFEST.sha256` inside. Extract at repository root. The v0.5.0 archive
(`datasets/gravcomm-v0.5.0-data.zip`) is still required for the v0.5.0 workflows and FEA inputs.
Workstation paths in text records were made repository-relative.

v0.5.0 erratum: `results/macnae_noise_digitization_2026_09_21/README.md` in Git is a corrected copy (en dashes, archive
note) of the file in the v0.5.0 data archive, so it does not match `DATA_MANIFEST.sha256`. `scripts/verify_release_v040.py`
now accepts a data file whose Git copy matches `MANIFEST.sha256`; extracting the v0.5.0 archive over the repository
restores the archived copy.

## Requirements

As v0.5.0, plus `numba` for the joint receiver. `cupy` is optional (GPU detector); results record the device used.
GPU and CPU detectors produce identical graphs within 4e-15 relative.

## Limitations

The receiver floor is a projection and the environment a declared scenario. Acquisition is simulated with constant
amplitude, phase and delay, clock-supplied carrier frequency and a known training frame; drift between frames and
clock wander are not modeled. Tungsten variants and mitigation scenarios are not independently BER-qualified.
