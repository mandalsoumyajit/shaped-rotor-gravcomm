"""SNR, capacity, and range solvers."""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np

from .rotor import TwoMassRotor


def _positive(value: float, name: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def snr_power(amplitude: np.ndarray | float, noise_asd: float, bandwidth: float) -> np.ndarray:
    """Power SNR for a sinusoidal peak amplitude in white acceleration noise."""

    _positive(noise_asd, "noise_asd")
    _positive(bandwidth, "bandwidth")
    amp = np.asarray(amplitude, dtype=float)
    return (amp * amp / 2.0) / (noise_asd * noise_asd * bandwidth)


def capacity_bps(amplitude: np.ndarray | float, noise_asd: float, bandwidth: float) -> np.ndarray:
    """Unconstrained white Gaussian-channel benchmark, NOT a rotor data rate.

    bandwidth is an explicit positive-frequency channel width; noise_asd is
    one-sided input ASD. A single harmonic's mean-square power is A^2/2.
    This does not impose constant-envelope, torque, coding or acquisition limits.
    """

    return bandwidth * np.log1p(snr_power(amplitude, noise_asd, bandwidth)) / math.log(2.0)


def amplitude_for_snr(snr: float, noise_asd: float, bandwidth: float) -> float:
    """Peak sinusoid amplitude required to reach a target power SNR."""

    _positive(snr, "snr")
    _positive(noise_asd, "noise_asd")
    _positive(bandwidth, "bandwidth")
    return math.sqrt(2.0 * snr) * noise_asd * math.sqrt(bandwidth)


def solve_monotonic_range(
    amplitude_at_distance: Callable[[float], float],
    target_amplitude: float,
    lower: float,
    upper: float = 1e4,
    iterations: int = 80,
) -> float:
    """Solve `amplitude(distance) = target` for a decreasing positive curve."""

    _positive(target_amplitude, "target_amplitude")
    _positive(lower, "lower")
    if not math.isfinite(upper) or upper <= lower:
        raise ValueError("upper must be finite and greater than lower")
    lo = lower
    hi = upper
    while amplitude_at_distance(hi) > target_amplitude:
        hi *= 2.0
        if hi > 1e9:
            raise RuntimeError("range search failed to bracket target")

    if amplitude_at_distance(lo) < target_amplitude:
        raise ValueError("target is unattainable at the minimum allowed distance")

    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if amplitude_at_distance(mid) > target_amplitude:
            lo = mid
        else:
            hi = mid
    return hi


def solve_rotor_range_for_snr(
    rotor: TwoMassRotor,
    noise_asd: float,
    bandwidth: float,
    snr: float = 10.0,
    samples: int = 4096,
    harmonic: int = 2,
) -> float:
    """Range where the selected Fourier harmonic reaches a target power SNR."""

    target = amplitude_for_snr(snr, noise_asd, bandwidth)
    lower = rotor.max_arm * 1.001
    return solve_monotonic_range(
        lambda d: rotor.harmonic_amplitude(d, harmonic=harmonic, samples=samples),
        target,
        lower=lower,
    )


def solve_rotor_range_for_capacity(
    rotor: TwoMassRotor,
    noise_asd: float,
    bandwidth: float,
    target_capacity: float,
    samples: int = 4096,
    harmonic: int = 2,
) -> float:
    """Range where exact point-mass amplitude reaches target Shannon capacity."""

    _positive(target_capacity, "target_capacity")
    _positive(bandwidth, "bandwidth")
    required_snr = math.expm1(math.log(2.0) * target_capacity / bandwidth)
    return solve_rotor_range_for_snr(
        rotor,
        noise_asd=noise_asd,
        bandwidth=bandwidth,
        snr=required_snr,
        samples=samples,
        harmonic=harmonic,
    )


def multicarrier_capacity_bps(amplitudes, noise_asds, bandwidths) -> float:
    """Sum white Gaussian benchmarks for disjoint bands; ASDs are SQUARED.

    Each amplitude is the peak of one independently modulated carrier. This
    does not treat the harmonics of one rotor as independent data streams.
    """
    a, n, b = (np.asarray(x, dtype=float) for x in (amplitudes, noise_asds, bandwidths))
    if a.ndim != 1 or a.shape != n.shape or a.shape != b.shape or a.size == 0:
        raise ValueError("provide nonempty equal-length one-dimensional arrays")
    if not all(np.all(np.isfinite(x)) for x in (a, n, b)) or np.any(n <= 0) or np.any(b <= 0):
        raise ValueError("inputs must be finite; noise and bandwidth must be positive")
    return float(np.sum(b * np.log1p(a*a/(2*n*n*b)) / math.log(2.0)))


def coherent_tone_power_snr(amplitude, noise_asd: float, integration_time):
    """Known-frequency sinusoid matched-filter SNR, one-sided white ASD.

    For an integer number of cycles, var(A_hat)=ASD^2/T. Thus the coherent
    power SNR is A^2*T/ASD^2, equal to snr_power at ENBW=1/(2T).
    Channel bandwidth B and this estimator ENBW are NOT interchangeable.
    """
    _positive(noise_asd, "noise_asd")
    t = np.asarray(integration_time, dtype=float)
    if np.any(~np.isfinite(t)) or np.any(t <= 0):
        raise ValueError("integration_time must be finite and positive")
    return np.asarray(amplitude, dtype=float)**2 * t / noise_asd**2
