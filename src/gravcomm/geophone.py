"""Geophone receiver model and two-element gradiometer.

Models a moving-coil velocity geophone (default: I/O Sensor SM-24) as an in-band
receiver for the 40-265 Hz carriers, and a differential pair as a gravity
gradiometer (Experiment D of the experimental-validation roadmap).

A geophone is a velocity sensor: its response is flat in *velocity* above the
suspension resonance, so its *acceleration* noise floor rises with frequency.
The dominant self-noise terms are coil Johnson noise (electrical) and suspension
Brownian noise (mechanical); a preamp voltage-noise term is included optionally.

Caveat on double counting: the Johnson term here uses the coil resistance and is
treated as independent of the mechanical Brownian term parameterized by a
mechanical Q.  Real geophone damping is part electrical (already in Johnson) and
part mechanical; the split is approximate and the floor should be *measured*
(Experiment A) before being quoted in the manuscript.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import K_B
from .link_budget import amplitude_for_snr, solve_monotonic_range
from .rotor import TwoMassRotor


@dataclass(frozen=True)
class Geophone:
    """Moving-coil velocity geophone.

    ``sensitivity`` is the open-circuit transduction constant (V per m/s),
    ``damping`` the total damping ratio in the deployed (shunted) condition,
    ``mechanical_q`` the suspension quality factor used for Brownian noise only.
    """

    natural_hz: float = 10.0
    damping: float = 0.6
    sensitivity: float = 28.8
    moving_mass: float = 11e-3
    coil_resistance: float = 375.0
    mechanical_q: float = 20.0
    amp_voltage_asd: float = 0.0
    temperature: float = 293.0
    label: str = "SM-24"

    @property
    def angular_natural(self) -> float:
        return 2.0 * math.pi * self.natural_hz

    def velocity_response(self, frequency_hz: float) -> float:
        """Magnitude of output voltage per ground velocity (V per m/s).

        High-pass velocity response: -> sensitivity for f >> f0, ~f^2 below f0.
        """

        if frequency_hz <= 0:
            raise ValueError("frequency_hz must be positive")
        ratio = frequency_hz / self.natural_hz
        numerator = ratio * ratio
        denominator = math.sqrt((1.0 - ratio * ratio) ** 2 + (2.0 * self.damping * ratio) ** 2)
        return self.sensitivity * numerator / denominator

    def johnson_acceleration_asd(self, frequency_hz: float) -> float:
        """Coil Johnson noise referred to input acceleration (m/s^2/sqrt(Hz))."""

        voltage = math.sqrt(4.0 * K_B * self.temperature * self.coil_resistance)
        velocity = voltage / self.velocity_response(frequency_hz)
        return velocity * 2.0 * math.pi * frequency_hz

    def amplifier_acceleration_asd(self, frequency_hz: float) -> float:
        """Preamp voltage noise referred to input acceleration."""

        if self.amp_voltage_asd <= 0:
            return 0.0
        velocity = self.amp_voltage_asd / self.velocity_response(frequency_hz)
        return velocity * 2.0 * math.pi * frequency_hz

    def brownian_acceleration_asd(self) -> float:
        """Suspension Brownian noise referred to acceleration (white).

        S_a = 4 k_B T omega_0 / (m Q_mech).
        """

        return math.sqrt(
            4.0 * K_B * self.temperature * self.angular_natural
            / (self.moving_mass * self.mechanical_q)
        )

    def acceleration_noise_asd(self, frequency_hz: float) -> float:
        """Total self-noise referred to input acceleration (m/s^2/sqrt(Hz))."""

        johnson = self.johnson_acceleration_asd(frequency_hz)
        brownian = self.brownian_acceleration_asd()
        amp = self.amplifier_acceleration_asd(frequency_hz)
        return math.sqrt(johnson * johnson + brownian * brownian + amp * amp)


SM24 = Geophone()


@dataclass(frozen=True)
class Gradiometer:
    """Two co-aligned geophones at fixed baseline, outputs subtracted.

    The baseline points along the line to the source, so the differential
    output approximates the gravity gradient.  ``cmrr_db`` is the common-mode
    rejection of correlated (e.g. floor vibration) input after calibration.
    """

    element: Geophone
    baseline: float
    cmrr_db: float = 40.0

    def __post_init__(self) -> None:
        if self.baseline <= 0:
            raise ValueError("baseline must be positive")

    def differential_signal_amplitude(self, rotor: TwoMassRotor, distance: float) -> float:
        """Second-harmonic differential acceleration, subtracting complex phasors."""

        near = distance - 0.5 * self.baseline
        far = distance + 0.5 * self.baseline
        if near <= rotor.max_arm:
            raise ValueError("near element too close to the rotor")
        return abs(rotor.harmonic_coefficient(near) - rotor.harmonic_coefficient(far))

    def device_noise_asd(self, frequency_hz: float) -> float:
        """Differential device noise: two independent elements add in quadrature."""

        return self.element.acceleration_noise_asd(frequency_hz) * math.sqrt(2.0)

    def effective_noise_asd(self, frequency_hz: float, ambient_asd: float = 0.0) -> float:
        """Differential device noise plus CMRR-suppressed common-mode ambient."""

        residual_common = ambient_asd * 10.0 ** (-self.cmrr_db / 20.0)
        device = self.device_noise_asd(frequency_hz)
        return math.hypot(device, residual_common)


def solve_geophone_range(
    rotor: TwoMassRotor,
    geophone: Geophone,
    carrier_hz: float,
    bandwidth: float,
    snr: float = 10.0,
) -> float:
    """Range where the single-geophone link reaches a target SNR."""

    noise = geophone.acceleration_noise_asd(carrier_hz)
    target = amplitude_for_snr(snr, noise, bandwidth)
    return solve_monotonic_range(
        lambda d: rotor.harmonic_amplitude(d),
        target,
        lower=rotor.max_arm * 1.001,
    )


def solve_gradiometer_range(
    rotor: TwoMassRotor,
    gradiometer: Gradiometer,
    carrier_hz: float,
    bandwidth: float,
    snr: float = 10.0,
    ambient_asd: float = 0.0,
) -> float:
    """Range where the gradiometer link reaches a target SNR.

    Uses the differential signal and the CMRR-limited effective noise.  The
    differential signal falls off one power of distance faster than the single
    sensor, so the range is shorter but common-mode disturbances are rejected.
    """

    noise = gradiometer.effective_noise_asd(carrier_hz, ambient_asd=ambient_asd)
    target = amplitude_for_snr(snr, noise, bandwidth)
    lower = (rotor.max_arm + 0.5 * gradiometer.baseline) * 1.001
    return solve_monotonic_range(
        lambda d: gradiometer.differential_signal_amplitude(rotor, d),
        target,
        lower=lower,
    )


def noise_floor_table(geophone: Geophone, carriers_hz: tuple[float, ...]) -> dict[float, float]:
    """Acceleration noise floor at each carrier, for the measured Table I tier."""

    return {f: geophone.acceleration_noise_asd(f) for f in carriers_hz}
