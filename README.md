# Shaped-rotor gravitational communication

**v0.5.0 - publication archive with tone-count and electromagnetic comparisons**

Numerical code for a proposed rotating-source Newtonian gravitational communication link. The current study sizes steel-insert transmitters at 0.5, 1 and 2 m for a nominal payload rate of 1 bit/s, with fixed instrumental noise/bandwidth and a tunable center frequency. Results are conditional simulations, not a hardware demonstration.

| Distance | Bare mass | Diameter | Carrier | Peak inertial torque |
|---:|---:|---:|---:|---:|
| 0.5 m | 7.03 kg | 0.481 m | 24 Hz | 1.60 N m |
| 1 m | 37.69 kg | 0.841 m | 24 Hz | 26.26 N m |
| 2 m | 227.32 kg | 1.531 m | 18 Hz | 524.73 N m |

Independent packet-aware qualification supports post-decoding BER <=1e-3 for the selected channels, with 224 s frame budgets and approximately 0.99 bit/s empirical accepted throughput. The 0.5 and 1 m designs share the same normalized channel and qualification sample. The 2 m replacement was qualified in a separate statistical family; no joint 95% claim across both studies is made.

## Added in v0.5.0

See `docs/RELEASE_v0.5.0.md` for new data, limitations and entry points. This release adds the tone-count study, planar-wall magnetic comparison, atomic receiver alternatives, historical environmental spectrum and array illustration. The seven-tone gravitational design results below are unchanged.

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
python -m pytest -q tests
python -m unittest discover -s fea -p test_fea.py -v
```

The Gmsh Python module and CalculiX executable are additional FEA dependencies. The retained mechanical runs used Gmsh 4.12.1 and CalculiX 2.21. Source checksums preserve exact bytes (`.gitattributes` disables line-ending conversion).

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

See `scripts/streaming_dynamics_runtime.py` for exact framing and receiver settings. The final code uses 256 information bits, 224 payload bits, rate 2/3 and 384 coded bits; 12 acquisition symbols are budgeted but acquisition is not simulated.

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

The complete numerical archive is also included in `datasets/gravcomm-v0.5.0-data.zip` so that repository-based archival services preserve the data with the code.
