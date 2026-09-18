"""Newtonian field models for rotating point-mass transmitters."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import DEFAULT_QUADRUPOLE_COEFFICIENT, G


@dataclass(frozen=True)
class TwoMassRotor:
    """Two point masses on one rigid rotating line.

    Arms are signed distances from the rotation axis.  Opposite signs place the
    two masses on opposite sides of the axis at phase zero.
    """

    masses: tuple[float, float]
    arms: tuple[float, float]
    label: str = "two-mass rotor"

    @classmethod
    def balanced_quadrupole(
        cls,
        total_mass: float,
        long_arm: float,
        mass_ratio: float = 1.0,
        label: str = "balanced quadrupole",
    ) -> "TwoMassRotor":
        """Build a balanced quadrupole.

        `mass_ratio` is heavy/light.  The light mass sits at `+long_arm`; the
        heavy mass sits at `-long_arm/mass_ratio`, so `sum(m_i r_i) = 0`.
        `mass_ratio=1` gives the equal-mass dumbbell.
        Balance removes the dipole multipole, not all odd temporal harmonics:
        only the equal-mass, equal-arm case has half-turn symmetry.
        """

        if total_mass <= 0:
            raise ValueError("total_mass must be positive")
        if long_arm <= 0:
            raise ValueError("long_arm must be positive")
        if mass_ratio < 1 or not np.isfinite(mass_ratio):
            raise ValueError("mass_ratio (heavy/light) must be finite and >= 1")

        light = total_mass / (mass_ratio + 1.0)
        heavy = total_mass - light
        short_arm = long_arm / mass_ratio
        return cls((heavy, light), (-short_arm, long_arm), label=label)

    @classmethod
    def asymmetric_dipole(
        cls,
        mass: float,
        arm: float,
        label: str = "asymmetric dipole",
    ) -> "TwoMassRotor":
        """Build a one-sided point-mass dipole model."""

        if mass <= 0:
            raise ValueError("mass must be positive")
        if arm <= 0:
            raise ValueError("arm must be positive")
        return cls((mass, 0.0), (arm, 0.0), label=label)

    @property
    def total_mass(self) -> float:
        return float(sum(self.masses))

    @property
    def max_arm(self) -> float:
        return float(max(abs(a) for a in self.arms))

    @property
    def dipole_moment(self) -> float:
        return float(sum(m * a for m, a in zip(self.masses, self.arms)))

    @property
    def quadrupole_moment(self) -> float:
        return float(sum(m * a * a for m, a in zip(self.masses, self.arms)))

    def acceleration_component(
        self,
        distance: float,
        phase: np.ndarray | float,
        component: str = "x",
    ) -> np.ndarray:
        """Acceleration at sensor `(distance, 0)` from the rotating masses.

        The rotation plane is the x-z plane.  `component` may be `"x"` or `"z"`.
        """

        if distance <= 0:
            raise ValueError("distance must be positive")
        phase_arr = np.asarray(phase, dtype=float)
        acc = np.zeros_like(phase_arr, dtype=float)
        sensor_x = distance

        for mass, arm in zip(self.masses, self.arms):
            if mass == 0.0:
                continue
            x = arm * np.cos(phase_arr)
            z = arm * np.sin(phase_arr)
            dx = x - sensor_x
            dz = z
            radius = np.sqrt(dx * dx + dz * dz)
            if component == "x":
                acc += G * mass * dx / radius**3
            elif component == "z":
                acc += G * mass * dz / radius**3
            else:
                raise ValueError("component must be 'x' or 'z'")

        return acc

    def ac_signal(
        self,
        distance: float,
        samples: int = 4096,
        component: str = "x",
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return one cycle of phase and mean-subtracted acceleration."""

        phase = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
        signal = self.acceleration_component(distance, phase, component)
        return phase, signal - float(np.mean(signal))

    def peak_ac_amplitude(
        self,
        distance: float,
        samples: int = 4096,
        component: str = "x",
    ) -> float:
        """Peak absolute value of the AC acceleration over one rotation."""

        _, ac = self.ac_signal(distance, samples=samples, component=component)
        return float(np.max(np.abs(ac)))

    def harmonic_coefficient(
        self, distance: float, harmonic: int = 2, samples: int = 4096,
        component: str = "x",
    ) -> complex:
        """Complex peak phasor: a_n(phi) = Re[c_n exp(i n phi)].

        Uniform periodic quadrature of the exact Newtonian waveform. The DC
        term is removed before projection. This retains phase for differential
        receivers; subtract phasors before taking their magnitude.
        """
        if not isinstance(harmonic, (int, np.integer)) or harmonic < 1:
            raise ValueError("harmonic must be a positive integer")
        if not isinstance(samples, (int, np.integer)) or samples <= 2 * harmonic:
            raise ValueError("samples must be an integer greater than 2*harmonic")
        if not np.isfinite(distance) or distance <= self.max_arm:
            raise ValueError("harmonic evaluation requires distance > max_arm")
        phase, ac = self.ac_signal(distance, samples=samples, component=component)
        return complex(2.0 * np.mean(ac * np.exp(-1j * harmonic * phase)))

    def harmonic_amplitude(
        self, distance: float, harmonic: int = 2, samples: int = 4096,
        component: str = "x",
    ) -> float:
        """Peak amplitude of ONE harmonic; its mean-square power is A_n^2/2."""
        return abs(self.harmonic_coefficient(distance, harmonic, samples, component))

    def rms_ac_acceleration(self, distance: float, samples: int = 4096, component: str = "x") -> float:
        """RMS of ALL AC harmonics, not the selected carrier amplitude."""
        _, ac = self.ac_signal(distance, samples=samples, component=component)
        return float(np.sqrt(np.mean(ac * ac)))


def far_field_dipole_amplitude(dipole_moment: float, distance: np.ndarray | float) -> np.ndarray:
    """Peak acceleration from the leading far-field dipole term."""

    d = np.asarray(distance, dtype=float)
    return 2.0 * G * abs(dipole_moment) / d**3


def far_field_quadrupole_amplitude(
    quadrupole_moment: float,
    distance: np.ndarray | float,
    coefficient: float = DEFAULT_QUADRUPOLE_COEFFICIENT,
) -> np.ndarray:
    """Peak acceleration from the draft's leading far-field quadrupole term."""

    d = np.asarray(distance, dtype=float)
    return coefficient * G * abs(quadrupole_moment) / d**4
