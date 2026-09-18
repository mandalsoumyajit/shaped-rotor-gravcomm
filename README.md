# Shaped-rotor gravitational communication

Reproduction code and model-generated data for *Near-Field Gravitational
Communication with a Shaped Rotor: Finite-Element Assessment and CF-FSK Modeling*,
by Soumyajit Mandal and Arjuna Madanayake (2026 draft).

The default source is a 7.913509 kg rotor with an aluminium carrier and two
steel inserts. It is the highest-signal member of three completed rounded
designs at fixed mass, not a global optimum. The numerical CF-FSK point gives
0.333 bit/s at 0.5 m under an idealized frequency-dependent receiver model.
This is not a demonstrated link or a certified mechanical design.

## License and provenance

Original code, documentation and model-generated data in this repository
are released under the MIT License; see `LICENSE`. Third-party software
(CalculiX, Gmsh, NumPy, SciPy and Matplotlib) retains its own licenses and is
not bundled. The manuscript, third-party papers and earlier private drafts
are not included in this code repository. No funding or conflicts of interest
were reported by the authors. Both authors developed the concept; SM performed
the simulations and wrote the paper.

This package contains accepted result summaries, waveform samples, the selected
source quadrature, CAD, and simulation code. Large solver `.dat`, `.frd` and
mesh files are regenerated rather than committed. `MANIFEST.sha256` records
the original exported files; rerunning simulations may change output hashes.
`docs/model_notes.md` documents conventions and historical diagnostics;
the workflows below identify the current paper results.

## Environment

Python 3.10+ with NumPy, SciPy and Matplotlib runs the analytical models,
actuator, detector and unit tests. For a Python virtual environment:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -e .
PYTHONPATH=src python -m unittest discover -s tests -v
```

On Windows, activate `.venv/Scripts/Activate.ps1`, use `python` and set
`$env:PYTHONPATH = "$PWD/src"` before testing. FEA is tested on Ubuntu under
WSL with CalculiX 2.21 and Gmsh 4.12.1. Install the OS packages for NumPy,
SciPy, Matplotlib and CalculiX (`calculix-ccx`); make the Gmsh 4.12.1 Python
module available to that same Python. `python3 -s` avoids mixing user-site
packages into the Ubuntu system environment. Check the actual Gmsh and
CalculiX versions rather than assuming the distribution provides them.

The paper plots use Gmsh to read the supplied BREP; they do not need to run
CalculiX. The FEA scripts use a temporary Linux directory for solver I/O and
copy results back. Full regeneration can require several GB and substantial
CPU time. Solver/Gmsh version changes may alter mesh numbering and local peaks.

## Fast reproduction from included data

From the repository root in the configured Ubuntu environment:

```sh
PYTHONPATH=src python3 -s -m unittest discover -s tests -v
python3 -s -m unittest discover -s fea -p test_fea.py -v
python3 -s scripts/paper_results.py
python3 -s scripts/desktop_cf_fsk.py --plot-only
```

This runs 57 model tests and 5 FEA tests, then regenerates the paper's compact
geometry, field, receiver, modulation-chain and actuator figures. The separate
Gaussian water-filling benchmark is an optimistic comparison, not CF-FSK rate.
Run `python3 -s scripts/desktop_cf_fsk.py` to regenerate the symbol-period
screen, constrained all-transition trajectory, conditional Monte Carlo results,
and technical-noise sweep. Seeds are fixed in the script. Its millions of
independent conditional decisions are not a continuous noisy packet experiment.

## Data map

| Result | Location |
|---|---|
| Three rounded designs and mass/inertia | `fea/results/design_sweep/filleted_screen.csv` |
| Stress windows, volumes and convergence | `fea/results/design_sweep/averaged_stress_*.csv` |
| Selected CAD and exact properties | `fea/results/design_sweep/medium_rounded/h0.013_rpm1800/` |
| Finite-source waveform and spectra | `fea/results/shaped_waveforms/` |
| Load/mode checks on original carrier | `fea/results/stress_modes/` |
| Bar verification | `fea/results/benchmark/` |
| CF-FSK assumptions, trials and trajectory | `results/desktop_cf_fsk/` |
| Current paper plots and Gaussian benchmark | `results/paper_revision/` |

The selected NPZ may be loaded using
`gravcomm.FiniteRotor.from_npz('fea/results/shaped_waveforms/medium_rounded.npz')`.
It stores positive quadrature masses and SI coordinates for both materials.
Rotation is about z; the radial receiver is at (d,0,0) measuring along x.
Mass is not renormalized to conceal the faceted-boundary error.

## Full FEA and field regeneration

Run these commands from `fea/`, in order:

```sh
python3 -s run_calculix.py
python3 -s shaped_rotor.py
python3 -s summarize.py
python3 -s postprocess_carrier.py
python3 -s design_sweep.py
python3 -s refine_designs.py
python3 -s shape_pair.py
python3 -s filleted_designs.py
python3 -s static_regions.py
python3 -s shaped_waveforms.py
python3 -s design_report.py
```

`run_calculix.py` is a bar benchmark; `shaped_rotor.py` is the original
unfilleted carrier. `filleted_designs.py` generates the three accepted rounded
variants. It only reuses a result when metadata and the necessary solver files
exist and the geometry matches. `--geometry-only` checks fixed-mass designs
without solving. `static_regions.py` averages integration-point von Mises
stress over fixed 5 and 10 mm radius aluminium windows. See
`fea/STRESS_MODES_AND_DESIGN.md` for convergence limitations and distinctions
between original and selected geometry.

Re-run the two paper scripts from the repository root after regeneration.
Older helpers `reproduce_figures.py`, the `gravcomm-figures` entry point and
`operating_envelope.py` are historical sensitivity/illustration workflows;
their deployment classes or binary illustration are not the current paper's
default rotor or seven-tone operating point.

## Interpretation limits

- Fixed bore, bonded inserts, linear elastic materials; no actual retention,
  rotor-bearing, fatigue, gyroscopic or containment qualification.
- Region-average stress convergence does not prove local peak convergence.
- Signal power is A2 squared/2, not waveform peak squared/2; higher harmonics
  are calculated separately and excluded from the decoded channel.
- Noise is a one-sided frequency-dependent PSD. Resonance gain is not a second
  reduction of input thermal noise, and the linewidth is not a hard band limit.
- Receiver timing and carrier acquisition, calibration and startup rejection
  are assumed. Environmental coherent coupling and hardware BER are untested.

## Versioned source

Repository: https://github.com/mandalsoumyajit/shaped-rotor-gravcomm

Paper reproduction version: `v0.2.0`. The release includes the source, documented
assumptions and accepted numerical data. Cite this version rather than an
unspecified moving branch. A DOI archive may be added for long-term preservation;
no DOI has yet been assigned.
