"""Input-referred oscillator noise; all ASDs are one-sided SI quantities.

Carter et al., Scientific Reports 14, 17775 (2024), Eq. (3): structural
damping. https://doi.org/10.1038/s41598-024-68623-0
The nominal 100 fm/sqrt(Hz) readout below is an idealized flat approximation,
NOT digitized measured noise. Environmental/control noise is not included.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
import numpy as np

from .constants import K_B


@dataclass(frozen=True)
class StructuralOscillator:
    mass_kg: float = 0.0031
    resonance_hz: float = 50.3
    quality_factor: float = 637000.0
    temperature_k: float = 300.0
    readout_displacement_asd: float = 100e-15
    # Conservative modeling domain below the reported 209 Hz higher mode.
    minimum_frequency_hz: float = 0.2
    maximum_frequency_hz: float = 200.0

    def __post_init__(self):
        for name in ("mass_kg", "resonance_hz", "quality_factor", "temperature_k",
                     "minimum_frequency_hz", "maximum_frequency_hz"):
            x = getattr(self, name)
            if not math.isfinite(x) or x <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not math.isfinite(self.readout_displacement_asd) or self.readout_displacement_asd < 0:
            raise ValueError("readout_displacement_asd must be finite and nonnegative")
        if not self.minimum_frequency_hz < self.resonance_hz < self.maximum_frequency_hz:
            raise ValueError("resonance must be inside the modeling domain")

    def _frequencies(self, frequency_hz):
        f = np.asarray(frequency_hz, dtype=float)
        if np.any(~np.isfinite(f)) or np.any(f < self.minimum_frequency_hz) or np.any(f > self.maximum_frequency_hz):
            raise ValueError("frequency is outside the oscillator modeling domain")
        return f

    def displacement_per_acceleration(self, frequency_hz):
        """|x/a| in s^2 for structural (constant loss-angle) damping.

        The structural term is omega0^2/Q, not omega*omega0/Q (viscous).
        """
        f = self._frequencies(frequency_hz)
        return 1.0 / ((2*math.pi)**2 * np.hypot(self.resonance_hz**2-f*f,
                                               self.resonance_hz**2/self.quality_factor))

    def thermal_acceleration_asd(self, frequency_hz):
        f = self._frequencies(frequency_hz)
        return np.sqrt(8*math.pi*K_B*self.temperature_k*self.resonance_hz**2 /
                       (self.mass_kg*f*self.quality_factor))

    def readout_acceleration_asd(self, frequency_hz):
        return self.readout_displacement_asd / self.displacement_per_acceleration(frequency_hz)

    def acceleration_noise_asd(self, frequency_hz):
        return np.hypot(self.thermal_acceleration_asd(frequency_hz),
                        self.readout_acceleration_asd(frequency_hz))

    def band_noise_upper_asd(self, carrier_hz: float, bandwidth_hz: float) -> float:
        """Conservative constant ASD for [fc-B/2, fc+B/2], not a linewidth cut.

        Thermal PSD decreases with f. Flat-readout input PSD is maximal at
        one endpoint. Sum their separate maxima to bound total PSD everywhere
        in the band, even when those maxima occur at opposite endpoints.
        Used only to form a conservative WHITE-NOISE channel benchmark; the
        result is not the exact colored-noise capacity or measured sensitivity.
        """
        if not math.isfinite(bandwidth_hz) or bandwidth_hz <= 0:
            raise ValueError("bandwidth_hz must be finite and positive")
        edges = self._frequencies([carrier_hz-bandwidth_hz/2, carrier_hz+bandwidth_hz/2])
        thermal = float(self.thermal_acceleration_asd(edges[0]))
        readout = float(np.max(self.readout_acceleration_asd(edges)))
        return math.hypot(thermal, readout)


CARTER_OSCILLATOR = StructuralOscillator()
