"""Torsion-balance receiver model for the low-frequency validation experiment.

This module predicts the signal, noise, and SNR for detecting a slowly rotating
quadrupole source with a Cavendish-style torsion balance (Experiment E1 of the
experimental-validation roadmap).  The source rotates about a vertical axis
(parallel to the fibre) at horizontal distance ``d`` from the balance; the
time-varying horizontal gravitational pull torques the boom about the fibre.
A balanced quadrupole produces its leading torque at twice the rotation rate.

Two regimes matter, and the model exposes both:

* **Thermal-limited.**  The Brownian torque noise passes through the same
  transfer function as the signal, so the SNR is independent of where the
  carrier sits relative to resonance; tuning to resonance does not by itself
  improve a thermal-limited measurement (the signal torque is set by geometry,
  not rotation rate).
* **Readout-limited.**  When flat angle-readout noise dominates, the signal is
  boosted by Q on resonance while the readout noise is not, so resonance lifts
  the signal above the readout floor and helps.

This is the bench analogue of the manuscript's "on-resonance optomechanical"
receiver tier, and the distinction is the scientific content of the experiment.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import G, K_B
from .link_budget import snr_power
from .rotor import TwoMassRotor

# Shear moduli (Pa) for common torsion-fibre materials.
SHEAR_MODULUS = {
    "tungsten": 1.61e11,
    "fused_silica": 3.1e10,
    "steel": 7.9e10,
}


@dataclass(frozen=True)
class TorsionBalance:
    """Symmetric two-test-mass torsion balance on a vertical fibre.

    ``test_mass`` is each of the two test masses, placed at ``+/-test_arm`` along
    the horizontal boom.  ``resonance_hz`` and ``quality_factor`` describe the
    torsional mode; ``readout_asd`` is the flat angle-readout noise floor of the
    optical lever (rad/sqrt(Hz)).
    """

    test_mass: float
    test_arm: float
    resonance_hz: float
    quality_factor: float
    readout_asd: float
    boom_inertia: float = 0.0
    temperature: float = 293.0

    def __post_init__(self) -> None:
        for name in ("test_mass", "test_arm", "resonance_hz", "quality_factor", "readout_asd"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")

    @property
    def polar_inertia(self) -> float:
        """Rotational inertia about the fibre axis (kg m^2)."""

        return 2.0 * self.test_mass * self.test_arm * self.test_arm + self.boom_inertia

    @property
    def angular_resonance(self) -> float:
        return 2.0 * math.pi * self.resonance_hz

    @property
    def stiffness(self) -> float:
        """Torsional spring constant kappa = I_p omega_0^2 (N m/rad)."""

        return self.polar_inertia * self.angular_resonance ** 2

    @property
    def thermal_torque_asd(self) -> float:
        """Suspension Brownian torque noise ASD (N m/sqrt(Hz)).

        S_tau = 4 k_B T kappa / (omega_0 Q), the fluctuation-dissipation result
        for a damped torsional oscillator.
        """

        return math.sqrt(
            4.0 * K_B * self.temperature * self.stiffness
            / (self.angular_resonance * self.quality_factor)
        )

    def transfer_magnitude(self, angular_frequency: float | np.ndarray) -> np.ndarray:
        """Torque-to-angle transfer |H(omega)| (rad per N m).

        |H| = 1 / (kappa * sqrt((1 - r^2)^2 + (r/Q)^2)), r = omega/omega_0.
        |H(0)| = 1/kappa; |H(omega_0)| = Q/kappa.
        """

        r = np.asarray(angular_frequency, dtype=float) / self.angular_resonance
        denom = self.stiffness * np.sqrt((1.0 - r * r) ** 2 + (r / self.quality_factor) ** 2)
        return 1.0 / denom

    def fibre_radius(self, length: float, material: str = "tungsten") -> float:
        """Radius (m) of a cylindrical fibre giving the required stiffness.

        kappa = pi G_shear r^4 / (2 L)  ->  r = (2 L kappa / (pi G_shear))^{1/4}.
        """

        if length <= 0:
            raise ValueError("length must be positive")
        shear = SHEAR_MODULUS[material]
        return (2.0 * length * self.stiffness / (math.pi * shear)) ** 0.25


@dataclass(frozen=True)
class TorsionSource:
    """Rotating point-mass source seen by the torsion balance.

    The source rotor turns about a vertical axis located at horizontal distance
    ``distance`` from the fibre, in a plane offset vertically by ``height_offset``
    from the boom.  ``rotation_hz`` is the mechanical rotation rate; the leading
    quadrupole torque appears at ``carrier_hz = 2 * rotation_hz``.
    """

    rotor: TwoMassRotor
    distance: float
    rotation_hz: float
    height_offset: float = 0.0

    def __post_init__(self) -> None:
        if self.distance <= self.rotor.max_arm:
            raise ValueError("distance must exceed the source arm to avoid contact")
        if self.rotation_hz <= 0:
            raise ValueError("rotation_hz must be positive")

    @property
    def carrier_hz(self) -> float:
        return 2.0 * self.rotation_hz


def torque_waveform(
    balance: TorsionBalance,
    source: TorsionSource,
    samples: int = 4096,
) -> tuple[np.ndarray, np.ndarray]:
    """One cycle of source phase and the mean-subtracted fibre torque (N m).

    Exact point-mass sum: for test masses at (+/-test_arm, 0, 0) and source
    masses sweeping a horizontal circle of radius |arm| about (distance, 0,
    height_offset), the torque about the vertical fibre axis is sum of
    x_t * F_y over all test-mass/source-mass pairs.
    """

    phase = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    torque = np.zeros_like(phase)
    for x_t in (balance.test_arm, -balance.test_arm):
        for mass, arm in zip(source.rotor.masses, source.rotor.arms):
            if mass == 0.0:
                continue
            x_s = source.distance + arm * np.cos(phase)
            y_s = arm * np.sin(phase)
            dx = x_s - x_t
            dy = y_s
            dz = source.height_offset
            radius = np.sqrt(dx * dx + dy * dy + dz * dz)
            force_y = G * balance.test_mass * mass * dy / radius ** 3
            torque += x_t * force_y
    return phase, torque - float(np.mean(torque))


def harmonic_torque_amplitude(
    balance: TorsionBalance,
    source: TorsionSource,
    harmonic: int = 2,
    samples: int = 4096,
) -> float:
    """Amplitude (N m) of the torque component at the given phase harmonic."""

    if harmonic < 1:
        raise ValueError("harmonic must be >= 1")
    phase, torque = torque_waveform(balance, source, samples=samples)
    cos_proj = np.sum(torque * np.cos(harmonic * phase))
    sin_proj = np.sum(torque * np.sin(harmonic * phase))
    return float(2.0 / samples * math.hypot(cos_proj, sin_proj))


def peak_torque(balance: TorsionBalance, source: TorsionSource, samples: int = 4096) -> float:
    """Peak absolute AC torque on the fibre over one cycle (N m)."""

    _, torque = torque_waveform(balance, source, samples=samples)
    return float(np.max(np.abs(torque)))


def signal_angle_amplitude(
    balance: TorsionBalance,
    source: TorsionSource,
    harmonic: int = 2,
    samples: int = 4096,
) -> float:
    """Steady-state boom-angle amplitude (rad) at the carrier harmonic."""

    torque_amp = harmonic_torque_amplitude(balance, source, harmonic=harmonic, samples=samples)
    angular_carrier = 2.0 * math.pi * (harmonic * source.rotation_hz)
    return torque_amp * float(balance.transfer_magnitude(angular_carrier))


def angle_noise_asd(
    balance: TorsionBalance,
    source: TorsionSource,
    harmonic: int = 2,
) -> dict[str, float]:
    """Angle-noise ASD (rad/sqrt(Hz)) at the carrier: thermal, readout, total."""

    angular_carrier = 2.0 * math.pi * (harmonic * source.rotation_hz)
    transfer = float(balance.transfer_magnitude(angular_carrier))
    thermal = balance.thermal_torque_asd * transfer
    readout = balance.readout_asd
    total = math.hypot(thermal, readout)
    return {"thermal": thermal, "readout": readout, "total": total}


def signal_to_noise(
    balance: TorsionBalance,
    source: TorsionSource,
    bandwidth: float,
    harmonic: int = 2,
    samples: int = 4096,
) -> float:
    """Power SNR of the lock-in detected carrier in a given bandwidth."""

    signal = signal_angle_amplitude(balance, source, harmonic=harmonic, samples=samples)
    noise = angle_noise_asd(balance, source, harmonic=harmonic)["total"]
    return float(snr_power(signal, noise, bandwidth))


def integration_time_for_snr(
    balance: TorsionBalance,
    source: TorsionSource,
    target_snr: float = 10.0,
    harmonic: int = 2,
    samples: int = 4096,
) -> float:
    """Approximate integration time (s) to reach a target SNR (t ~ 1/B)."""

    if target_snr <= 0:
        raise ValueError("target_snr must be positive")
    signal = signal_angle_amplitude(balance, source, harmonic=harmonic, samples=samples)
    noise = angle_noise_asd(balance, source, harmonic=harmonic)["total"]
    bandwidth = (signal * signal / 2.0) / (noise * noise * target_snr)
    if bandwidth <= 0:
        return math.inf
    return 1.0 / bandwidth


@dataclass(frozen=True)
class TorsionDesignPoint:
    carrier_hz: float
    on_resonance: bool
    peak_torque: float
    carrier_torque: float
    signal_angle: float
    thermal_angle_asd: float
    readout_angle_asd: float
    snr_unit_bandwidth: float
    time_to_snr10: float


def design_point(
    balance: TorsionBalance,
    source: TorsionSource,
    bandwidth: float = 1.0,
    harmonic: int = 2,
    samples: int = 4096,
) -> TorsionDesignPoint:
    """Summarize a balance/source pair for a design-point table."""

    noise = angle_noise_asd(balance, source, harmonic=harmonic)
    near_resonance = abs(source.carrier_hz - balance.resonance_hz) < 0.1 * balance.resonance_hz
    return TorsionDesignPoint(
        carrier_hz=source.carrier_hz,
        on_resonance=near_resonance,
        peak_torque=peak_torque(balance, source, samples=samples),
        carrier_torque=harmonic_torque_amplitude(balance, source, harmonic=harmonic, samples=samples),
        signal_angle=signal_angle_amplitude(balance, source, harmonic=harmonic, samples=samples),
        thermal_angle_asd=noise["thermal"],
        readout_angle_asd=noise["readout"],
        snr_unit_bandwidth=signal_to_noise(balance, source, bandwidth, harmonic=harmonic, samples=samples),
        time_to_snr10=integration_time_for_snr(balance, source, 10.0, harmonic=harmonic, samples=samples),
    )
