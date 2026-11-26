# Release v0.6.1 - final-revision studies

This release adds the small analyses and figure updates behind the final revision of the manuscript
(Physica Scripta, PHYSSCR-157981). The code, results and data archives of v0.6.0 are unchanged; no packet
simulation was rerun, and no new data archive is needed.

## New material

`paper_studies/ground_motion_study/` (run from that folder):

- `calibration_tolerance.py` (manuscript Sec. 2.4 and App. B.1): accuracy needed when tracking the receiver
  resonance. A calibration with the wrong resonance frequency f0 or quality factor Q multiplies the signal and
  every noise term by the same response ratio; the script computes the resulting loss of effective
  signal-to-noise ratio for a statistic matched to the transmitted waveform under the assumed noise spectrum.
  Result (`results/calibration_tolerance.json`): a 0.05 dB loss needs an f0 error of 0.17 mHz (24 Hz) or
  0.16 mHz (18 Hz), about 7e-6 of the carrier; Q must be wrong by a factor of about 80 to cost the same.
- `torsional_mode.py` (App. A.8): the 24.17 Hz drivetrain torsional mode of the 2 m rotor near its 18 Hz carrier.
  Steady excitation is amplified 2.1-2.4 times across the tone band; commanded transitions excite the mode
  negligibly (overshoot < 0.22%); one newton-metre of steady torque twists the rotor by 1.6e-6 rad
  (`results/torsional_mode.json`).
- `floor_transfer.py` (Sec. 5.6): transfer of a force at the transmitter mount through a shared concrete floor slab
  (infinite Kirchhoff plate), which sets the tolerable twice-per-revolution force (prints its results).
- `modal_sequences.py` (App. A.6): single-mode response to isolated, alternating, cyclic and random hop sequences
  (prints its results).

`paper_studies/figure_scripts/`:

- `draw_horizontal_support.py` (Figure 13): the horizontal-shaft support concept redrawn for the qualified
  0.5 m-link rotor (451 mm diameter, 71 mm steel and 31 mm same-mass tungsten inserts). The v0.5.0 drawing
  (`fea/draw_horizontal_support.py`, unchanged) shows the conventional 481 mm rotor.
- `receiver_schematic.tex` (Figure 2): panel (a) marks the 50 Hz resonance as that of the reference sensor, tuned
  to 18 or 24 Hz in the links; panel (b) now shows the nonbinary LDPC decoder exchanging word probabilities with
  the sequence detector, as in the v0.6.0 receiver, and states that acquisition is simulated separately.

## Reproduction

`calibration_tolerance.py` imports the v0.6.0 waveform code (`hop_packets.py`, `hop_receiver.py`) and noise model
(`mitigation.py`); the other scripts need only NumPy, SciPy and Matplotlib. All run in seconds to minutes.

## Requirements and limitations

As v0.6.0. The calibration tolerance is a linear template-mismatch calculation, not a packet simulation; the
torsional-mode check is a single-mode model of the screen in App. A.8, without motor or coupling inertia.
