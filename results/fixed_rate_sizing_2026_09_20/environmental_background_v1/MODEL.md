# Assumed environmental background for carrier selection

Adopt NHNM plus an explicit isolation model as a provisional seismic-background baseline; NLNM is the quiet sensitivity case. This is a reproducible design scenario, not a measured laboratory spectrum or a bound on all disturbances. Instrumental sensitivity and 2 Hz receive-filter bandwidth remain fixed; background noise is evaluated at absolute physical frequency.

## Model and provenance

Peterson ground-acceleration PSD: 10^((A+B log10(1/f))/10) in (m/s^2)^2/Hz, using the existing tabulated coefficients in src/gravcomm/noise.py. Model data end at 10 Hz. From 10 to 100 Hz we explicitly hold ground-acceleration ASD at its 10 Hz value; this is an assumed engineering continuation, not published Peterson data or a universal upper bound. The scenario domain is 0.1-100 Hz. Sources: https://doi.org/10.3133/ofr93322 and https://search.r-project.org/CRAN/refmans/IRISSeismic/html/noiseModels.html.

Convert vertical ground acceleration to displacement by dividing ASD by (2*pi*f)^2. For a single horizontal test mass in an isotropic Rayleigh field, gravitational ASD is (2*pi*G*rho*gamma/sqrt(2))*exp(-2*pi*f*h/c_R)*displacement_ASD, following Harms, Eq. (109) of the 2015 review: https://link.springer.com/article/10.1007/lrr-2015-3. We assume rho=2000 kg/m^3, gamma=0.8 and h=0; these are declared ground/site assumptions, not fitted measurements. Rayleigh speed 200 m/s is immaterial at h=0. This does not treat all seismic wave types as measured Rayleigh waves.

For support motion, use an explicit hypothetical two-stage horizontal isolation chain: each stage has resonance 0.3 Hz and damping ratio 0.05. Its base-to-mass transfer is (1+2*i*zeta*r)/(1-r^2+2*i*zeta*r), r=f/f_iso; multiply the two transfers. This cascade is an idealized feed-forward stage model, not a validated coupled suspension. Assume horizontal ground ASD equals the vertical seismic envelope for this scenario. Tilt, actuator/control noise and the coupled six-degree-of-freedom isolation dynamics require a later hardware budget. These assumptions must accompany any numerical operating-point result.

Because support and Newtonian disturbances can share a seismic origin, use the conservative coherent ASD sum before squaring: S_background=(N_Newtonian+N_support)^2. Add independent instrumental PSD afterward. This is a worst-phase envelope, not a claim of statistically independent seismic paths.

## Calculated levels at the tuned carrier center

| Seismic case | Frequency (Hz) | Newtonian ASD (nGal/sqrt Hz) | Support ASD | Total with instrumental floor |
|---|---:|---:|---:|---:|
| NLNM | 0.1 | 0.0007804 | 821.6 | 821.6 |
| NLNM | 0.2 | 0.002648 | 2.828e+04 | 2.828e+04 |
| NLNM | 1 | 5.752e-06 | 5.196 | 7.316 |
| NLNM | 2 | 1.267e-06 | 0.3226 | 5.159 |
| NLNM | 5 | 2.223e-07 | 0.0228 | 5.149 |
| NLNM | 10 | 4.784e-08 | 0.003912 | 5.149 |
| NLNM | 20 | 1.196e-08 | 0.0009163 | 5.149 |
| NLNM | 50.3 | 1.891e-09 | 0.0001421 | 5.149 |
| NHNM | 0.1 | 0.1951 | 2.054e+05 | 2.054e+05 |
| NHNM | 0.2 | 0.3919 | 4.185e+06 | 4.185e+06 |
| NHNM | 1 | 0.001727 | 1560 | 1560 |
| NHNM | 2 | 0.0005267 | 134.1 | 134.2 |
| NHNM | 5 | 0.000704 | 72.21 | 72.4 |
| NHNM | 10 | 0.0003198 | 26.15 | 26.65 |
| NHNM | 20 | 7.994e-05 | 6.124 | 8.001 |
| NHNM | 50.3 | 1.264e-05 | 0.9499 | 5.236 |

All entries above 10 Hz use the stated extrapolation. The fixed instrumental center floor is 5.1494 nGal/sqrt(Hz).

## Implications and next calculation

The Rayleigh gravitational component alone is small relative to the current instrumental floor at frequencies of a few hertz and above in these scenarios. Residual support motion is the dominant frequency-dependent penalty of the assumed terrestrial receiver. Do not assert that Newtonian noise necessarily sets an optimum in the tens-of-hertz band. A frequency optimum must follow the actual noise budget, rather than an imposed low-frequency penalty.

Varying isolation resonance to 0.15, 0.3 and 0.6 Hz is a hardware sensitivity study. The crossings below compare background with the instrumental center floor only; they are not operating-frequency optima or BER thresholds:

- NLNM, isolation resonance 0.15 Hz: background remains below center floor above 0.6345999009812868 Hz on the evaluated grid.
- NLNM, isolation resonance 0.3 Hz: background remains below center floor above 1.0046157902783952 Hz on the evaluated grid.
- NLNM, isolation resonance 0.6 Hz: background remains below center floor above 1.9431216918284249 Hz on the evaluated grid.
- NHNM, isolation resonance 0.15 Hz: background remains below center floor above 10.889300933334333 Hz on the evaluated grid.
- NHNM, isolation resonance 0.3 Hz: background remains below center floor above 21.802183971859467 Hz on the evaluated grid.
- NHNM, isolation resonance 0.6 Hz: background remains below center floor above 43.65158322401661 Hz on the evaluated grid.

The optional CausalReceiver(background=SeismicBackground(), carrier_hz=fc) interface now uses S_total(nu;fc)=S_instrument_reference(50.3+nu)+S_background(fc+nu) through the same causal filter and alias sum, preserving the fixed instrumental profile. It validates the entire acquisition-frequency interval and rejects out-of-domain carriers. Four environmental/integration tests pass. Evaluate full waveform likelihoods and power at candidate carriers; do not use a center-only SNR to choose a carrier. Equal-amplitude simulation reuse is now conditional on carrier and environmental scenario.

The present 16 complex samples/s envelope interface spans offsets +-8 Hz. At fc<=8.1 Hz it reaches below this scenario domain, and for fc<=8 Hz crosses zero physical frequency. Do not substitute abs(f) or clip silently: low-carrier packet simulations require a consistent real-passband/acquisition model or a separately validated narrower sampled-envelope interface. A 2 Hz modulation band also requires a positive lower band edge. The scalar spectrum calculation here is valid independently of that unresolved interface.

Atmospheric gravity, nearby moving objects, coherent lines and technical/control disturbances are not quantified by this seismic scenario. It is a defensible assumed baseline, not complete environmental qualification. Add explicit spectra/coherent disturbance cases when data are available; do not credit chopping with rejecting genuine gravity noise. No code ranking, chosen carrier or post-decoding BER claim follows from these levels.
