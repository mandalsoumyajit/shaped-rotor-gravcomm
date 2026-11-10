# Shaped-rotor gravitational communication

**v0.6.0 - joint GF(512) receiver, simulated acquisition and Physica Scripta revision studies**

Numerical code for a proposed rotating-source Newtonian gravitational communication link. The current study sizes steel-insert transmitters at 0.5, 1 and 2 m for a nominal payload rate of 1 bit/s, with fixed instrumental noise/bandwidth and a tunable center frequency. Results are conditional simulations, not a hardware demonstration.

| Distance | Bare mass | Diameter | Carrier | Peak inertial torque |
|---:|---:|---:|---:|---:|
| 0.5 m | 5.79 kg | 0.451 m | 24 Hz | 0.51 N m |
| 1 m | 30.98 kg | 0.788 m | 24 Hz | 8.38 N m |
| 2 m | 184.05 kg | 1.427 m | 18 Hz | 163.2 N m |

These designs use a hop-limited seven-tone waveform without dwell and a (43,29) LDPC code over GF(512) defined on four-tone words, decoded iteratively with the streaming sequence detector. Independent packet-aware qualification supports post-decoding BER <=1e-3 at 0.580 A0 (24 Hz) and 0.697 A0 (18 Hz), with 224 s frames and 0.9995 bit/s accepted throughput. With simulated preamble acquisition and tracking from the previous CRC-verified frame, the same amplitudes pass a fresh qualification (paired loss 0.04-0.05 dB).

## Added in v0.6.0

See `docs/RELEASE_v0.6.0.md`. The revision studies are in `paper_studies/` (information rates and Pareto search, the joint receiver with CPU and GPU detectors, waterfalls, qualification, replays and timing, acquisition, ground-motion mitigation, the comparison with Brack et al. 2023, and figure scripts), with the finite-element and magnetic-link reruns for the new rotors in `fea/*_hop.py` and `scripts/magnetic_wall_*_2026_11_05.py`. Packet records are in the supplementary data archive `gravcomm-v0.6.0-data-supplement.zip`. The conventional-design results of v0.5.0 (7.03, 37.69 and 227.32 kg) are retained unchanged for comparison.

## Added in v0.5.0

See `docs/RELEASE_v0.5.0.md` for new data, limitations and entry points. This release adds the tone-count study, planar-wall magnetic comparison, atomic receiver alternatives, historical environmental spectrum and array illustration. The conventional seven-tone design results below are unchanged.

## Retained from v0.4.0

- Causal sixth-order receive filtering, decimation and aliased-noise covariance.
- Explicit environmental spectrum with declared isolation and high-frequency extrapolation.
- Bounded-state streaming soft sequence detection with history/lookahead diagnostics.
- Published convolutional mother polynomials and puncturing, tail biting, interleaving and CRC framing; short LDPC candidates retained for comparison.
- Independent qualification records, packet-aware confidence sequences and source/amplitude sizing.
- Three-distance centrifugal/gravity/modulation FEA, modal coupling and refinement checks.
- Same-mass tungsten comparison and horizontal-shaft support drawing.

The original v0.3.0 workflows remain for provenance. Their uncoded 50.3 Hz operating point and ideal-band capacity ratios are not the current link results. See `docs/v0.3.0_README.md` for the historical workflow.

## Install and test

Python 3.10+ with NumPy, SciPy, Matplotlib and Numba. Linux/WSL is recommended for parallel studies and CalculiX. Numerical-library threads should be limited when using multiple workers.

```sh
python -m pip install -e '.[sequence,test]'
export PYTHONPATH=src
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python -m pytest -q tests   # includes tests/test_paper_studies.py (needs numba)
python -m unittest discover -s fea -p test_fea.py -v
```

NumPy 2.4 removed `np.trapz`; `tests/conftest.py` restores it for the test session, and the historical scripts that use it (left byte-identical because their hashes are recorded) should be run with NumPy < 2.4. The Gmsh Python module and CalculiX executable are additional FEA dependencies. The retained mechanical runs used Gmsh 4.12.1 and CalculiX 2.21. Source checksums preserve exact bytes (`.gitattributes` disables line-ending conversion).

## Results and data asset

Summary JSON and reports are included in Git. Download **gravcomm-v0.5.0-data.zip** from the v0.5.0 release and extract at repository root before auditing full packet records or rerunning FEA. The archive includes 36,992 original qualification packets plus 26,240 four-tone packets across the initial and replacement studies, exploratory/diagnostic records, frozen configurations, mechanical summaries, and the two original reference input decks needed by the scaled steel workflow. `DATA_MANIFEST.sha256` checks the extracted data.

Large solver result fields, meshes and logs are regenerated, not shipped. Serialized workstation folder paths have been made repository-relative; numerical inputs and original scientific source hashes are preserved. Run commands from repository root unless stated otherwise. Historical scripts whose data are not included require their original inputs and are not current-release entry points.

```sh
python scripts/verify_release_v040.py
python scripts/compile_ber_qualification.py
python scripts/compile_replacement_ber.py
python scripts/plot_revised_link_summary.py
```

The two qualification compilers check source hashes, all packets and the first decisive stopping prefix. They reproduce reports without new Monte Carlo sampling. They can take several minutes. Do not pool screening, convergence replays or shared-channel samples into qualification counts.

To run new independent qualifications, preserve the released data and use a separate results directory/configuration with fresh seeds. The existing drivers resume included checkpoints rather than manufacture extra independent evidence:

```sh
python scripts/qualify_streaming_ber.py --workers 8 --batch 128
python scripts/qualify_replacement_ber.py --workers 8 --batch 128
```

See `scripts/streaming_dynamics_runtime.py` for exact framing and receiver settings. The final code uses 256 information bits, 224 payload bits, rate 2/3 and 384 coded bits; 12 acquisition symbols are budgeted; acquisition is simulated for the v0.6.0 joint design (`paper_studies/ground_motion_study/acquisition.py`, `acq_packets.py`).

## Mechanical reproduction

After extracting the data asset, run the following in a Linux/WSL environment with Gmsh and CalculiX:

```sh
python fea/three_link_fea.py --prepare --workers 4
python fea/compile_three_link_fea.py --mesh .018
python fea/compile_three_link_fea.py --mesh .013
python fea/three_link_modal_response.py --mesh .018
python fea/three_link_modal_response.py --mesh .013
python fea/refine_three_link_reference.py
python fea/assembly_screen.py
python fea/report_three_link_fea.py

python fea/tungsten_link_fea.py --prepare --workers 4
python fea/compile_tungsten_link_fea.py --mesh .018
python fea/compile_tungsten_link_fea.py --mesh .013
python fea/tungsten_modal_response.py --mesh .018
python fea/tungsten_modal_response.py --mesh .013
python fea/refine_tungsten_reference.py
python fea/tungsten_field_comparison.py
python fea/draw_horizontal_support.py
python fea/report_tungsten_comparison.py
```

Included summary-only solver folders are evidence, not completed reusable solver caches. Use `scripts/prepare_fea_rerun.py` first to remove only summary markers whose required raw solver files are absent, or run in a separate clean results tree. This prevents accidentally skipping calculations because a public summary exists. Each four-worker run can require many GB of memory and disk. The reports contain explicit attachment and modal limitations.

## Scope

The instrument profile is translated unchanged with carrier; physical tuning, acquisition and site noise are unvalidated. The environmental model is a declared scenario, not measured noise. The receiver is reduced-state, with numerical convergence checks, not an exact unlimited-memory decoder. Mechanical models fix the bore and bond inserts; actual retention, contact, fatigue and drivetrain dynamics remain unresolved. Commanded power excludes drag and electrical losses. Tungsten reduces protrusion; no CFD drag savings or new tungsten BER qualification is asserted.

Original code and generated data retain the MIT License. Third-party dependencies retain their licenses. Manuscript drafts, coauthor correspondence and private editorial files are not included. Cite this software version with `CITATION.cff`; the paper itself is maintained separately.

The complete numerical archives are also included in `datasets/gravcomm-v0.5.0-data.zip` and `datasets/gravcomm-v0.6.0-data-supplement.zip` so that repository-based archival services preserve the data with the code.
