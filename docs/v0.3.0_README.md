# Shaped-rotor gravitational communication

Simulation code and generated data for *Near-Field Gravitational Communication
with a Shaped Rotor: Finite-Element Assessment and CF-FSK Modeling*, by
Soumyajit Mandal and Arjuna Madanayake. Paper reproduction version: **v0.3.0**.

The reference source is a 7.913509 kg aluminium carrier with two steel inserts.
The selected seven-tone continuous-frequency frequency-shift keying (CF-FSK)
point uses 2 s symbols and supports 1.333 mapped bit/s at 0.5 m. Eight-byte
packets with six known trailer symbols carry 1.067 bit/s, excluding acquisition.
The finest-sampling confirmation has 144 bit errors in 2,912,000 bits across
31 errored packets: BER 4.95e-5, with a packet-aware 95% interval
[1.59e-5, 1.70e-4]. These are acquired-packet simulations under the stated noise
model. Hardware communication, assembly strength and environmental rejection
remain experimental validation tasks.

## What changed in v0.3.0

- Full-waveform standard Viterbi decoding retains continuous phase and finite
  whitening-filter memory, with exact add-compare-select and traceback.
- The existing byte-to-three-tone mapping is enforced in the trellis.
- Error-count-driven packet trials and anytime confidence intervals account
  for correlated errors within packets.
- The tone-spacing sweep separates receiver discrimination from drive demand.
- Exact worst-valid-message RMS torque and uniform-byte duty replace the
  earlier conservative individual-symbol duty screen in interpretation.
- Sampling, whitening, achieved-drive, technical-noise, energy and capacity
  checks accompany the detector and three-panel performance figures.

The earlier 8 s dwell-bin detector and its 0.333 bit/s result are historical.
Its scripts/data remain for provenance and helper functions; use the current
workflows below for paper claims. Initial sequence-search screens also used
different decoding/termination settings from the final BER confirmations.

## Environment

Python 3.10+; NumPy, SciPy and Matplotlib. Numba accelerates the sequence trellis.
The Monte Carlo worker scripts use Linux `fork`: run them on Linux or WSL.
The numerical library and unit tests also work on Windows.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[sequence,test]'
export PYTHONPATH=src
export OPENBLAS_NUM_THREADS=1
python -m pytest -q tests
python -m unittest discover -s fea -p test_fea.py -v
```

The checked environment is recorded in
`results/communications_checks/manifest.json`. For matching numerical behavior,
use its Python/NumPy/SciPy/Numba versions. On WSL, `python -s` prevents unrelated
user-site packages from overriding the configured environment.

## Fast checks and plots from included results

```sh
python scripts/verify_current_results.py
python scripts/check_revision_capacity.py
python scripts/summarize_communications_checks.py
python scripts/summarize_spacing_sweep.py
python scripts/sequence_paper_figures.py
```

These verify retained packet counts, confidence intervals, source amplitudes,
mechanical references and release-file availability, then regenerate reports
and the current receiver/detector/performance graphics. They do not rerun the
large Monte Carlo experiments. All current manuscript graphics are included.
For geometry/field figure regeneration, first run `scripts/paper_results.py`
and `scripts/introductory_rotor.py`, then `scripts/sequence_paper_figures.py`
to update the communications graphics. The geometry plot requires Gmsh 4.12.1.

## Repeating the communications experiments

```sh
python scripts/check_communications_revision.py --fs 4 --memory 6
python scripts/check_communications_revision.py --fs 8 --memory 6
python scripts/check_communications_revision.py --fs 4 --memory 8 --max-packets 6000
python scripts/check_communications_revision.py --technical-asd 1e-10
python scripts/check_communications_revision.py --technical-asd 2e-10
python scripts/sweep_sequence_spacing.py
python scripts/sweep_sequence_spacing.py --q .1 .5 1 --memory-symbols 4 --max-packets 1000
python scripts/sweep_sequence_spacing.py --q .1 .5 1 --sample-rate 16 --max-packets 1000
python scripts/check_revision_drive.py
python scripts/exact_sequence_duty.py
```

Run from the repository root. The confirmation driver resumes checkpoints;
both it and the sweep retain completed results. The sweep restarts incomplete
points. To regenerate from zero, move the relevant result JSON files to a
separate backup folder first. Keep included records for comparison.
Confirmation runs stop at at least 100 bit errors in at least 30 errored packets.
The longer-memory diagnostic has a 6000-packet cap; its sparse count supports
an upper bound. All seeds and per-packet counts are retained. Sampling/filter
variants share seed indices and must not be pooled as independent replications.
Runtime ranges from seconds for high-error noise cases to many minutes for
the primary confirmations and longer FIR trellises.

The default sweep uses 1.75 s symbols, 40% transitions and nine spacings.
Its receiver simulation uses ideal source trajectories; mechanical demands
are evaluated separately. The achieved-drive script performs fixed-message
paired sensitivity checks and energy accounting. It does not estimate
random-message BER from those fixed messages.

## Data map

| Result | Location |
|---|---|
| Current BER confirmations and checks | `results/communications_checks/` |
| Fixed-duration spacing sweep | `results/spacing_sweep/*_fs8_m3.json` |
| Earlier search stages and exact byte duty | `results/sequence_cf_fsk/` |
| Current paper figures | `results/paper_revision/` |
| Historical dwell detector | `results/desktop_cf_fsk/` |
| Rounded designs, mass, inertia, stress windows | `fea/results/design_sweep/` |
| Finite-volume waveforms and selected quadrature | `fea/results/shaped_waveforms/` |
| Original-carrier load/mode checks | `fea/results/stress_modes/` |
| Rotating-bar benchmark | `fea/results/benchmark/` |

`results/communications_checks/SUMMARY.md` is the current numerical reference.
Older `sequence_cf_fsk/SUMMARY.md` records the first 2 Hz confirmation; the
current paper uses the subsequent 8 Hz result. The historical water-filling
figure generated by `paper_results.py` uses a 1.5625 Hz band; current capacities
are in `communications_checks/capacity.json` and use a 2 Hz band.

The selected quadrature loads with
`gravcomm.FiniteRotor.from_npz('fea/results/shaped_waveforms/medium_rounded.npz')`.
Coordinates and masses are SI. Rotation is about z and the radial receiver at
(d,0,0) measures x acceleration. Source power is the second harmonic's peak
amplitude squared divided by two. Higher harmonics are tabulated separately.

## FEA regeneration

The accepted workflow uses Gmsh 4.12.1 and CalculiX 2.21 on Ubuntu/WSL.
Install `calculix-ccx` and the specified Gmsh Python module; check actual
versions. From `fea/`, run:

```sh
python run_calculix.py
python shaped_rotor.py
python summarize.py
python postprocess_carrier.py
python design_sweep.py
python refine_designs.py
python shape_pair.py
python filleted_designs.py
python static_regions.py
python shaped_waveforms.py
python design_report.py
```

The selected BREP and metadata are included. Large solver outputs and meshes
are regenerated; full runs can require several GB. Fixed-bore, perfectly
bonded, linear-elastic models quantify centrifugal loads and carrier modes.
Regional stress convergence supports geometry comparison; local failure,
retention, fatigue, shaft/bearings and containment require assembly analysis.

## Receiver and statistical scope

The structural-damping receiver has a 3.1 g proof mass, 50.3 Hz resonance,
Q=637000, temperature 300 K and assumed 100 fm/sqrt(Hz) displacement readout.
Frequency-dependent noise is input-referred once. Initial phase, timing,
frequency and tone history are acquired inputs. Physical oscillator startup,
coherent technical interference and hardware synchronization are outside the
stationary packet model. The calibrated-baseband model samples a finite band;
sample-rate and whitening checks quantify numerical sensitivity.

Standard Viterbi finds the best path for the finite-FIR squared-error metric.
Residual whitening error remains in the simulated physical noise. Individual
95% confidence sequences describe each run; the report makes no simultaneous
coverage claim over the search grid. The best verified feasible point leaves
the global maximum rate unresolved.

## License, provenance and citation

Original code, documentation and generated data are released under the MIT
License. Third-party software retains its own licenses and is not bundled.
The manuscript, reference papers and private editorial files are excluded.
Both authors developed the concept; SM performed the simulations and wrote
the paper. No funding or conflicts of interest were reported.

`MANIFEST.sha256` records the exported files; regeneration can change hashes.
The nested communications manifest records the checked computational inputs.
`docs/model_notes.md` contains historical modeling notes. Cite version v0.3.0
using `CITATION.cff` and the versioned GitHub release:
https://github.com/mandalsoumyajit/shaped-rotor-gravcomm/releases/tag/v0.3.0
