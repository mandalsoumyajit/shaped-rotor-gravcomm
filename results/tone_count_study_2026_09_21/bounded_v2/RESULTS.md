# Tone-count study: results and interpretation

Grid: 320 completed packets, 80 alphabet/mapping/spacing/carrier configurations. Fresh confirmation: 352 packets. Validation: 19 tests passed (15 generalized-alphabet checks and 4 original streaming regressions).

## What the comparison establishes

Tone count is a design variable. The three-symbol byte-mapping argument alone cannot justify seven tones as the best mechanical/communications choice. Several alternatives merit independent BER qualification. The existing seven-tone links retain their previously established qualification; this screen does not replace it.

The screen fixes source amplitudes (42.09 nGal at 24 Hz; 51.80 nGal at 18 Hz), 224 payload bits, rate-2/3 code, 224 s frame budget and the instrumental/environmental receiver model. The 24 Hz result serves both 0.5 and 1 m source sizes. Transition fraction is 0.6. Grouping and symbol time change with alphabet size. Acquisition remains an allowance, with known initial phase/timing supplied to the simulation.

## Fresh confirmation

| Tones | Mapping bits/symbols | Carrier Hz | Tone span Hz | Errors / payload bits | CRC rejected / packets | Peak torque ratio | Detector gate |
|---:|:---:|---:|---:|:---|:---|---:|:---|
| 3 | 3/2 | 18 | 0.979 | 0 / 7168 | 0 / 32 | 1.907 | pass |
| 4 | 2/1 | 18 | 0.562 | 0 / 7168 | 0 / 32 | 0.840 | pass |
| 4 | 2/1 | 24 | 0.562 | 0 / 7168 | 0 / 32 | 0.840 | pass |
| 5 | 9/4 | 18 | 0.679 | 0 / 7168 | 0 / 32 | 0.917 | pass |
| 5 | 9/4 | 24 | 0.679 | 0 / 7168 | 0 / 32 | 0.917 | pass |
| 7 | 8/3 | 18 | 0.868 | 0 / 7168 | 0 / 32 | 1.000 | pass |
| 7 | 8/3 | 24 | 0.868 | 2 / 7168 | 1 / 32 | 1.000 | pass |
| 7 | 11/4 | 18 | 0.846 | 0 / 7168 | 0 / 32 | 0.951 | pass |
| 7 | 11/4 | 24 | 0.846 | 0 / 7168 | 0 / 32 | 0.951 | pass |
| 8 | 3/1 | 18 | 0.913 | 0 / 7168 | 0 / 32 | 0.948 | pass |
| 8 | 3/1 | 24 | 0.913 | 0 / 7168 | 0 / 32 | 0.948 | pass |

All tentative payload errors count, including rejected packets. Fresh seeds differ from screening and convergence seeds. Zero errors in 32 packets gives a conservative packet-aware 95% BER upper bound of approximately 0.154 with the saved confidence method, far above 1e-3. These data support candidate selection and expose mechanical tradeoffs; they cannot establish reliability superiority. Individual confidence intervals are in confirmation_summary.json. The table is exploratory and has no simultaneous confidence claim after candidate selection.

## Screening and detector limits

The figure labels actual error counts; an asterisk means pruning or another per-run numerical check failed. Absence of an asterisk does not imply history/lookahead convergence; those separate diagnostics are summarized below. The equal-span grid uses span*Ts = 0.4, 0.8, 1.2, 1.6. The final column fixes adjacent spacing times dwell at q=0.08; duplicate points reuse records. This additional control prevents phase-lattice cost from silently eliminating viable eight-tone settings.

| Tones | Mapping | Carrier | Selected span*Ts | Pilot errors / 896 | Torque ratio | Convergence status |
|---:|:---:|---:|---:|---:|---:|:---|
| 2 | 1/1 | 24 | 0.4 | 0 | 2.053 | requires receiver refinement |
| 2 | 1/1 | 18 | 0.4 | 0 | 2.053 | requires receiver refinement |
| 3 | 3/2 | 24 | 0.8 | 0 | 1.907 | requires receiver refinement |
| 3 | 3/2 | 18 | 0.8 | 0 | 1.907 | pass |
| 4 | 2/1 | 24 | 0.6 | 0 | 0.840 | pass |
| 4 | 2/1 | 18 | 0.6 | 0 | 0.840 | pass |
| 5 | 9/4 | 24 | 0.8 | 0 | 0.917 | pass |
| 5 | 9/4 | 18 | 0.8 | 0 | 0.917 | pass |
| 6 | 5/2 | 24 | 0.8 | 0 | 0.752 | requires receiver refinement |
| 6 | 5/2 | 18 | 0.8 | 0 | 0.752 | requires receiver refinement |
| 7 | 8/3 | 24 | 1.2 | 0 | 1.000 | pass |
| 7 | 11/4 | 24 | 1.2 | 0 | 0.951 | pass |
| 7 | 8/3 | 18 | 1.2 | 0 | 1.000 | pass |
| 7 | 11/4 | 18 | 1.2 | 0 | 0.951 | pass |
| 8 | 3/1 | 24 | 1.4 | 0 | 0.948 | pass |
| 8 | 3/1 | 18 | 1.4 | 0 | 0.948 | pass |
| 9 | 6/2 | 24 | 1.6 | 0 | 1.083 | pilot only; higher torque than control |
| 9 | 6/2 | 18 | 1.6 | 0 | 1.083 | pilot only; higher torque than control |

History and lookahead are diagnostic parameters, not channel resources. The bounded screen uses four tone-history symbols and roughly the original lookahead duration in seconds. Wider-history checks use five symbols and up to 200000 states; separate checks double lookahead. Passing requires <1% relative soft-metric change, no changed bit signs or decoded bits, no pruning/missing alternatives, and whitening PSD error <1e-3. One diagnostic packet per point is a screening check; error-containing fresh packets and additional seeds should be replayed before qualification.

## Mechanical interpretation

Peak inertial torque is proportional to I*pi*(15/8)*span/(beta*Ts). The ratio column uses the same physical rotor as the published design at each distance. It excludes windage, bearing/electrical losses, drivetrain compliance and acquisition transients. Lower peak torque does not establish lower average electrical consumption. Carrier frequencies and source sizes are unchanged in this study.

The four-tone candidate has Ts=1.067 s, full tone span 0.5625 Hz and peak-torque ratio 0.8402, giving approximately 1.34, 22.06 and 440.87 N m for the existing 0.5, 1 and 2 m rotors. The five-tone candidate has Ts=1.179 s, span 0.6786 Hz and ratio 0.9170. These are commanded-motion comparisons at fixed source size; minimum source amplitudes have not yet been reoptimized.

## Scope and next decision

This is a finite alphabet/spacing screen with natural radix mappings, at most four symbols and 12 coded bits per group. It is not a global mapping, alphabet, code-rate, packet-length or transition-fraction optimization. Unused radix words create mapping-dependent redundancy and tone occupancy; this is explicitly part of the compared protocols. Known final padding bits are enforced. Acquisition/flush symbol allowances have different durations as Ts changes and remain within the common 224 s frame budget.

Recommended next step: independently qualify the four- and five-tone candidates, replay error packets with wider receiver settings, and refine their amplitude thresholds at the two normalized channels before changing manuscript design tables. Keep the seven-tone result as the already-qualified benchmark.

## Reproduction and progress

Study plan: ../PLAN.md. Immutable job manifests and per-packet JSON are in screen/, extra/, gates/ and confirm/. Code snapshots and SHA256 dependencies are retained. Entry points: tone_count_study.py, tone_count_extra.py, tone_count_gates.py, tone_count_confirm.py, compile_tone_count_study.py and finalize_tone_count_study.py. resume_tone_count_screen.py and robust_tone_host.py only add retry handling for transient Windows/OneDrive checkpoint locks. Existing scientific records and qualified receiver implementation were not modified. Nine-tone points remain at the pilot stage: the unpruned zero-error setting has higher peak torque than the qualified seven-tone control and was not advanced to expensive wider-history diagnostics in this bounded follow-up.
