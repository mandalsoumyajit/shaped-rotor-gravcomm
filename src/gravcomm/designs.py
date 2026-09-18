"""Receiver tiers and representative transmitter deployment classes."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import MICROGAL, NANOG
from .rotor import TwoMassRotor


@dataclass(frozen=True)
class ReceiverTier:
    key: str
    label: str
    noise_asd: float
    assumption: str


@dataclass(frozen=True)
class Deployment:
    key: str
    label: str
    total_mass: float
    long_arm: float
    tip_speed: float
    steady_power: float
    mass_ratio: float = 1.0

    @property
    def signal_frequency(self) -> float:
        """Quadrupole signal frequency, `2 omega / 2 pi = v / (pi L)`."""

        return self.tip_speed / (math.pi * self.long_arm)

    @property
    def default_bandwidth(self) -> float:
        """Legacy B=fc/2 benchmark ASSUMPTION, not an achievable bandwidth.

        Retained solely for like-for-like comparisons of the flat-noise
        reference curves. No actuator or hardware noise support is implied.
        The separate resonant model uses explicitly selected bandwidths.
        """

        return 0.5 * self.signal_frequency

    def rotor(self) -> TwoMassRotor:
        return TwoMassRotor.balanced_quadrupole(
            self.total_mass,
            self.long_arm,
            mass_ratio=self.mass_ratio,
            label=self.label,
        )


RECEIVERS = (
    ReceiverTier("mems", "MEMS floor extrapolation", 0.1 * MICROGAL,
                 "Hypothetical flat noise: Gao measured this floor at 0.14 Hz only; "
                 "reported bandwidth 108 Hz. Not a validated carrier-band specification."),
    ReceiverTier("opto_bb", "Carter broadband reference", 0.1 * NANOG,
                 "0.1 ng/sqrt(Hz), NOT nGal. Flat reference only; actual spectrum "
                 "is frequency dependent and this approximation is not validated over each signal band."),
)
# No frequency-independent resonant tier: see receivers.CARTER_OSCILLATOR.

DEPLOYMENTS = (
    Deployment("bench", "Bench", 50.0, 0.3, 250.0, 100.0),
    Deployment("field", "Field", 500.0, 0.5, 250.0, 100.0),
    Deployment("trailer", "Trailer", 5_000.0, 1.0, 250.0, 100.0),
    Deployment("fixed", "Fixed installation", 50_000.0, 2.0, 250.0, 200.0),
    Deployment("heavy", "Heavy fixed (CFRP)", 200_000.0, 5.0, 1000.0, 500.0),
)

# Values printed in the MERCon 2026 submission table.  They are retained as
# reference data so the AIP exact-model results can be compared explicitly.
PREVIOUS_MERCON_RANGE_TABLE_M = {
    "bench": {"mems": 0.4, "opto_bb": 1.3, "opto_res": 2.4},
    "field": {"mems": 0.7, "opto_bb": 3.0, "opto_res": 5.4},
    "trailer": {"mems": 1.9, "opto_bb": 7.6, "opto_res": 14.0},
    "fixed": {"mems": 4.7, "opto_bb": 19.0, "opto_res": 34.0},
    "heavy": {"mems": 16.0, "opto_bb": 48.0, "opto_res": 85.0},
}
