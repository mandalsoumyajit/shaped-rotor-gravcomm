"""Published acceleration-noise reference levels."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import MICROGAL, NANOG


@dataclass(frozen=True)
class NoiseReference:
    key: str
    label: str
    noise_asd: float
    source: str
    frequency_hz: float | None
    note: str = ""


# Peterson's empirical models are piecewise A + B log10(T), where T is period
# in seconds and the result is acceleration PSD in dB relative to
# (m s^-2)^2 / Hz.  The standard models cover periods down to 0.1 s (10 Hz).
_PETERSON_NLNM = np.array(
    [
        [0.10, -162.36, 5.64],
        [0.17, -166.70, 0.00],
        [0.40, -170.00, -8.30],
        [0.80, -166.40, 28.90],
        [1.24, -168.60, 52.48],
        [2.40, -159.98, 29.81],
        [4.30, -141.10, 0.00],
        [5.00, -71.36, -99.77],
        [6.00, -97.26, -66.49],
        [10.00, -132.18, -31.57],
        [12.00, -205.27, 36.16],
        [15.60, -37.65, -104.33],
        [21.90, -114.37, -47.10],
        [31.60, -160.58, -16.28],
        [45.00, -187.50, 0.00],
        [70.00, -216.47, 15.70],
        [101.00, -185.00, 0.00],
        [154.00, -168.34, -7.61],
        [328.00, -217.43, 11.90],
        [600.00, -258.28, 26.60],
        [10000.00, -346.88, 48.75],
        [100000.00, -346.88, 48.75],
    ],
    dtype=float,
)

_PETERSON_NHNM = np.array(
    [
        [0.10, -108.73, -17.23],
        [0.22, -150.34, -80.50],
        [0.32, -122.31, -23.87],
        [0.80, -116.85, 32.51],
        [3.80, -108.48, 18.08],
        [4.60, -74.66, -32.95],
        [6.30, 0.66, -127.18],
        [7.90, -93.37, -22.42],
        [15.40, 73.54, -162.98],
        [20.00, -151.52, 10.01],
        [354.80, -206.66, 31.63],
        [10000.00, -206.66, 31.63],
        [100000.00, -206.66, 31.63],
    ],
    dtype=float,
)


def peterson_acceleration_psd_db(frequency_hz: float, model: str) -> float:
    """Peterson NLNM/NHNM acceleration PSD in dB at a frequency."""

    if frequency_hz <= 0:
        raise ValueError("frequency_hz must be positive")
    period = 1.0 / frequency_hz
    data = {"NLNM": _PETERSON_NLNM, "NHNM": _PETERSON_NHNM}[model]
    if period < data[0, 0] or period > data[-1, 0]:
        raise ValueError(f"Peterson {model} is defined for {1.0 / data[-1, 0]:g}-{1.0 / data[0, 0]:g} Hz")
    idx = np.searchsorted(data[:, 0], period, side="right") - 1
    idx = max(0, min(idx, len(data) - 1))
    _, a, b = data[idx]
    return float(a + b * math.log10(period))


def peterson_acceleration_asd(frequency_hz: float, model: str) -> float:
    """Peterson NLNM/NHNM acceleration ASD in m s^-2 / sqrt(Hz)."""

    psd_db = peterson_acceleration_psd_db(frequency_hz, model)
    return math.sqrt(10.0 ** (psd_db / 10.0))


PETERSON_REFERENCE_FREQUENCY_HZ = 10.0

PUBLISHED_NOISE_REFERENCES = (
    NoiseReference(
        "carter_broadband",
        "Carter optomech.",
        0.1 * NANOG,
        "Carter et al. 2024",
        None,
        "0.1 ng/sqrt(Hz) = 98.0665 nGal/sqrt(Hz); flat reference, not a measured full-band spectrum.",
    ),
    NoiseReference(
        "gao_mems",
        "Gao MEMS",
        0.1 * MICROGAL,
        "Gao et al. 2026",
        0.14,
        "Device floor at 0.14 Hz only; not established at modeled carriers (reported bandwidth 108 Hz).",
    ),
    NoiseReference(
        "peterson_nlnm_10hz",
        "Peterson NLNM (10 Hz)",
        peterson_acceleration_asd(PETERSON_REFERENCE_FREQUENCY_HZ, "NLNM"),
        "Peterson 1993",
        PETERSON_REFERENCE_FREQUENCY_HZ,
        "Ground acceleration at 10 Hz, not Newtonian noise or a carrier-band site spectrum. Assumed unit coupling only.",
    ),
    NoiseReference(
        "peterson_nhnm_10hz",
        "Peterson NHNM (10 Hz)",
        peterson_acceleration_asd(PETERSON_REFERENCE_FREQUENCY_HZ, "NHNM"),
        "Peterson 1993",
        PETERSON_REFERENCE_FREQUENCY_HZ,
        "Ground acceleration at 10 Hz, not Newtonian noise or a carrier-band site spectrum. Assumed unit coupling only.",
    ),
)
