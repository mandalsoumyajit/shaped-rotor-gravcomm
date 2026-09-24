# Digitized Macnae nominal magnetic-background spectrum

The supplied screenshot was manually traced with 69 centerline anchors, calibrated against multiple logarithmic ticks on each axis. The image axis is horizontal magnetic-field amplitude spectral density in nT/sqrt(Hz). CSV values and the callable model use T/sqrt(Hz). For simulation, PSD = ASD squared with a one-sided convention; the schematic does not explicitly document its sidedness.

## Deliverables

- anchors.csv: original pixel coordinates, calibrated frequency/ASD, segment and solid/dashed metadata.
- calibration.json: tick calibration, screenshot hash, provenance and tracing tolerance.
- source_screenshot.png: supplied image, preserved for audit.
- trace_overlay.png: selected points superposed on the source drawing.
- digitized_spectrum.png / .pdf: standalone scientific plot.
- ../../src/gravcomm/macnae_noise.py: callable log-log interpolation model.
- ../../scripts/digitize_macnae_noise.py: reproducible digitization/plot script.

The interpolation is piecewise linear in log(f), log(ASD), equivalent to a power law between adjacent anchors. It preserves the visible geomagnetic features, cavity-resonance structure, asynchronous-grid-noise bump, and high-frequency broad hump. Values along dashed source portions are included and marked in the CSV. The separate diagonal constant-voltage coil-sensor line is excluded because it describes instrumental conversion, not environmental noise.

The trace ends below 1 kHz and resumes above 1 kHz; the unshown minimum is left undefined. The model raises for that gap and outside the trace domain; unsupported='nan' returns NaNs instead. It does not silently bridge, clamp, or extrapolate. Narrow power-line harmonics and VLF lines are omitted: their linewidths and integrated powers cannot be inferred reliably from this schematic. Specify those separately from measured RMS amplitudes and bandwidths. No gradiometric cancellation is included.

## Usage

```python
from gravcomm.macnae_noise import macnae_asd
ambient_asd = macnae_asd(frequency_Hz)  # T/sqrt(Hz)
noise_psd = instrument_asd**2 + ambient_asd**2 + wall_noise_psd
```

An optional positive `scale` multiplies ASD for explicitly declared scenario sensitivity. It is not an inferred confidence interval. No link results or manuscript content have been changed by this digitization. The existing magnetic Channel.gain function accepts a frequency-matched array for its ambient argument, enabling direct use over supported bands.

| Frequency (Hz) | Approximate ASD (pT/sqrt(Hz)) |
|---:|---:|
| 0.01 | 520 |
| 0.1 | 6.6 |
| 1 | 0.26 |
| 10 | 1.4 |
| 18 | 0.99 |
| 24 | 1.6 |
| 100 | 0.19 |

The shaded plotting band represents only +/-3 vertical pixels (about +/-8% in amplitude). Horizontal selection tolerance of 3 pixels corresponds to about 8% in frequency. These tolerances exclude uncertainty in the original curve, environment, epoch, calibration and interpretation. Values should be quoted approximately, especially near narrow bumps. All three tests pass: anchor/unit fidelity and power-law interpolation; gap/extrapolation handling; amplitude-to-PSD scaling and invalid-input rejection. Visual overlay was inspected.

## Provenance and interpretation

User-supplied reproduction of Figure 1 in Zhou et al., Measurement of Ambient Magnetic Field Noise for Through-the-Earth (TTE) Communications and Historical Comparisons, IEEE Transactions on Electromagnetic Compatibility 66(3), 720â€“727 (2024), DOI 10.1109/TEMC.2024.3354735:
https://pmc.ncbi.nlm.nih.gov/articles/PMC11951293/

That article attributes the drawing to J. C. Macnae, Y. Lamontagne and G. F. West, Noise processing techniques for time-domain EM systems, Geophysics 49(7), 934â€“948 (1984), DOI 10.1190/1.1441739. Original full text was not retrieved; attribution is verified through the reproducing article.

Zhou et al. explicitly caution that the original plot's acquisition and construction methods were unspecified and its numerical values are not definitive. They also report much higher recent mine noise than historical comparisons. Consequently this model is a historical nominal spectral-shape scenario, not a calibrated contemporary site spectrum. Measured-site or higher-noise sensitivity cases remain appropriate companions. The sub-Hz increase is represented, but the figure does not establish temporal statistics or spatial covariance for a gradiometer model.

Public archive note: the source screenshot and overlay listed above document the original workflow and are excluded from this distribution. Supply the cited source image separately using --source-image to regenerate the overlay.
