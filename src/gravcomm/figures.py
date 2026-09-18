"""Command-line reproduction of manuscript figures and numerical tables."""

from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch

from .constants import G, MICROGAL, MU0, NANOGAL
from .designs import DEPLOYMENTS, PREVIOUS_MERCON_RANGE_TABLE_M, RECEIVERS, Deployment
from .link_budget import (
    capacity_bps,
    solve_rotor_range_for_capacity,
    solve_rotor_range_for_snr,
    snr_power,
    coherent_tone_power_snr,
)
from .mechanics import rotor_mechanical_state, energy_loss_power, modulation_step
from .noise import PUBLISHED_NOISE_REFERENCES
from .receivers import CARTER_OSCILLATOR, StructuralOscillator
from .plotting import apply_paper_style
from .rotor import TwoMassRotor, far_field_dipole_amplitude, far_field_quadrupole_amplitude


def _save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.pdf")
    fig.savefig(outdir / f"{stem}.png", dpi=200)
    plt.close(fig)


def _exact_capacity_curve(deployment: Deployment, distances: np.ndarray, noise_asd: float) -> np.ndarray:
    rotor = deployment.rotor()
    values = np.full_like(distances, np.nan, dtype=float)
    mask = distances > 1.02 * rotor.max_arm
    amplitudes = np.array([rotor.harmonic_amplitude(d, samples=2048) for d in distances[mask]])
    values[mask] = capacity_bps(amplitudes, noise_asd, deployment.default_bandwidth)
    return values


def _exact_amplitude_curve(deployment: Deployment, distances: np.ndarray) -> np.ndarray:
    rotor = deployment.rotor()
    values = np.full_like(distances, np.nan, dtype=float)
    mask = distances > 1.02 * rotor.max_arm
    values[mask] = np.array([rotor.harmonic_amplitude(d, samples=2048) for d in distances[mask]])
    return values


def figure_1(outdir: Path) -> None:
    """Balanced quadrupole geometry and exact near-field time trace."""

    rotor = TwoMassRotor.balanced_quadrupole(
        total_mass=600.0,
        long_arm=1.0,
        mass_ratio=5.0,
        label="Fig. 1 balanced rotor",
    )
    m1, m2 = rotor.masses
    l1, l2 = abs(rotor.arms[0]), abs(rotor.arms[1])
    d_sensor = 1.5

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(7.0, 2.5),
        gridspec_kw={"width_ratios": [1.0, 1.45]},
    )

    ax1.set_xlim(-1.12, 1.85)
    ax1.set_ylim(-1.05, 1.05)
    ax1.set_aspect("equal")
    ax1.axis("off")

    theta = np.linspace(0, 2.0 * np.pi, 240)
    ax1.plot(l1 * np.cos(theta), l1 * np.sin(theta), "k--", lw=0.4, alpha=0.45)
    ax1.plot(l2 * np.cos(theta), l2 * np.sin(theta), "k--", lw=0.4, alpha=0.35)

    phi = np.deg2rad(25.0)
    heavy = (-l1 * np.cos(phi), -l1 * np.sin(phi))
    light = (l2 * np.cos(phi), l2 * np.sin(phi))
    ax1.plot([heavy[0], light[0]], [heavy[1], light[1]], "k-", lw=1.2)
    ax1.add_patch(Circle(heavy, 0.086, color="#333333", zorder=5))
    ax1.add_patch(Circle(light, 0.050, color="#888888", zorder=5))
    ax1.text(heavy[0] - 0.12, heavy[1] - 0.23, r"$M_1$", fontsize=9, ha="center")
    ax1.text(light[0] + 0.04, light[1] + 0.13, r"$M_2$", fontsize=9, ha="center")
    ax1.plot(0, 0, marker="+", color="k", markersize=10, mew=1.2, zorder=6)
    ax1.text(0.05, 0.13, "axis", fontsize=7)

    mid1 = (0.5 * heavy[0], 0.5 * heavy[1])
    mid2 = (0.5 * light[0], 0.5 * light[1])
    ax1.annotate("", xy=mid1, xytext=(mid1[0] - 0.05, mid1[1] - 0.23),
                 arrowprops=dict(arrowstyle="-", color="dimgray", lw=0.4))
    ax1.text(mid1[0] - 0.09, mid1[1] - 0.35, r"$L_1$", color="dimgray", ha="center")
    ax1.annotate("", xy=mid2, xytext=(mid2[0] + 0.02, mid2[1] + 0.18),
                 arrowprops=dict(arrowstyle="-", color="dimgray", lw=0.4))
    ax1.text(mid2[0] + 0.02, mid2[1] + 0.27, r"$L_2$", color="dimgray", ha="center")

    omega_r = 0.19
    ax1.add_patch(Arc((0, 0), 2 * omega_r, 2 * omega_r, theta1=110, theta2=170, color="k", lw=0.8))
    th = np.deg2rad(168.0)
    ax1.annotate(
        "",
        xy=(omega_r * np.cos(th), omega_r * np.sin(th)),
        xytext=(omega_r * np.cos(th + 0.05), omega_r * np.sin(th + 0.05)),
        arrowprops=dict(arrowstyle="-|>", color="k", lw=0.8),
    )
    ax1.text(-0.08, omega_r + 0.08, r"$\omega$", fontsize=9, ha="center")

    ax1.plot(d_sensor, 0, marker="s", color="#1f77b4", markersize=10, mec="k", mew=0.6)
    ax1.text(d_sensor, -0.20, "gravimeter", fontsize=7, ha="center", color="#1f77b4")
    ax1.annotate(
        "",
        xy=(d_sensor, 0.0),
        xytext=(0.0, 0.0),
        arrowprops=dict(arrowstyle="<->", color="#1f77b4", lw=0.5, shrinkA=6, shrinkB=8),
    )
    ax1.text(d_sensor / 2.0, 0.08, r"$d$", fontsize=9, color="#1f77b4", ha="center")
    ax1.text(-1.05, 0.92, "(a)", fontsize=10, fontweight="bold")

    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel(r"$\Delta g_x$ ($\mu$Gal)")
    times = np.linspace(0.0, 4.0, 4000, endpoint=False)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for freq, color in zip([0.5, 1.0, 3.0], colors):
        phase = 2.0 * np.pi * freq * times
        signal = rotor.acceleration_component(d_sensor, phase)
        ax2.plot(times, (signal - np.mean(signal)) / MICROGAL, color=color, lw=0.9, label=f"{freq:.1f} Hz")

    noise = RECEIVERS[0].noise_asd
    ax2.axhline(0, color="k", lw=0.3, alpha=0.5)
    ax2.axhspan(-noise / MICROGAL, noise / MICROGAL, color="gray", alpha=0.20,
                label=r"Assumed noise (1 Hz BW)")
    ax2.set_xlim(0, 4)
    ax2.set_title(
        rf"$M_1={m1:.0f}$ kg @ $L_1={l1:.1f}$ m, "
        rf"$M_2={m2:.0f}$ kg @ $L_2={l2:.1f}$ m; $d={d_sensor:.1f}$ m",
        fontsize=7.0,
    )
    ax2.legend(loc="lower right", frameon=False, ncol=2, columnspacing=0.7, handlelength=1.3, fontsize=6.2)
    ax2.text(0.02, 0.86, "(b)", transform=ax2.transAxes, fontsize=10, fontweight="bold")
    fig.tight_layout(pad=0.4)
    _save(fig, outdir, "fig1_concept")


def figure_2_quadrupole_validity(outdir: Path) -> None:
    """Exact Newtonian field vs. leading quadrupole approximation."""

    fig, ax = plt.subplots(figsize=(3.45, 2.25))
    distance_over_l = np.logspace(np.log10(1.08), 2.05, 180)
    cases = [
        (1.0, "equal arms\n$M_1=M_2$", "#1f77b4"),
        (5.0, "5:1 mass ratio\n(Fig. 1)", "#ff7f0e"),
        (20.0, "20:1 mass ratio", "#2ca02c"),
    ]

    for mass_ratio, label, color in cases:
        rotor = TwoMassRotor.balanced_quadrupole(
            total_mass=1.0,
            long_arm=1.0,
            mass_ratio=mass_ratio,
            label=label,
        )
        exact = np.array([rotor.harmonic_amplitude(float(d), samples=4096) for d in distance_over_l])
        approx = far_field_quadrupole_amplitude(rotor.quadrupole_moment, distance_over_l)
        ax.semilogx(distance_over_l, exact / approx, color=color, lw=1.1, label=label)

    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.75)
    ax.axvspan(1.08, 3.0, color="#eeeeee", alpha=0.65, lw=0)
    ax.text(1.18, 3.2, "near field", fontsize=7, color="#555555")
    ax.text(24.0, 1.08, "quadrupole limit", fontsize=7, color="#333333")
    ax.set_xlabel(r"Distance ratio $d/L_\mathrm{max}$")
    ax.set_ylabel(r"Exact $A_2$ / leading quadrupole")
    ax.set_xlim(1.08, 110.0)
    ax.set_ylim(0.0, 4.0)
    ax.grid(True, which="both", alpha=0.25, lw=0.3)
    ax.legend(loc="upper right", frameon=False, fontsize=6.1, handlelength=1.5)
    fig.tight_layout(pad=0.3)
    _save(fig, outdir, "fig2_quadrupole_validity")


def figure_2(outdir: Path) -> None:
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.0, 2.25))

    distances = np.logspace(-0.55, 2.25, 260)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b"]
    for deployment, color in zip(DEPLOYMENTS, colors):
        amplitudes = _exact_amplitude_curve(deployment, distances)
        ax_a.loglog(distances, amplitudes / MICROGAL, color=color, label=deployment.label)

    ax_a.axhline(RECEIVERS[0].noise_asd / MICROGAL, color="red", ls="--", lw=0.7, alpha=0.8,
                 label=r"MEMS extrapolation ($0.1\,\mu$Gal/$\sqrt{\rm Hz}$)")
    ax_a.axhline(RECEIVERS[1].noise_asd / MICROGAL, color="purple", ls=":", lw=0.8, alpha=0.8,
                 label=r"Carter flat reference ($98.1\,$nGal/$\sqrt{\rm Hz}$)")
    ax_a.set_xlabel("Distance $d$ (m)")
    ax_a.set_ylabel(r"Carrier amplitude $A_2$ ($\mu$Gal)")
    ax_a.set_xlim(0.28, 180)
    ax_a.set_ylim(1e-6, 1e3)
    ax_a.grid(True, which="both", alpha=0.25, lw=0.3)
    ax_a.legend(loc="lower left", frameon=False, fontsize=5.8, handlelength=1.6)
    ax_a.text(0.03, 0.95, "(a)", transform=ax_a.transAxes, fontsize=9, fontweight="bold")

    cap_distances = np.logspace(-0.5, 1.9, 180)
    for deployment, color in zip(DEPLOYMENTS[1:4], ["#ff7f0e", "#2ca02c", "#9467bd"]):
        c_mems = _exact_capacity_curve(deployment, cap_distances, RECEIVERS[0].noise_asd)
        c_opto = _exact_capacity_curve(deployment, cap_distances, RECEIVERS[1].noise_asd)
        ax_b.loglog(cap_distances, c_mems, color=color, lw=1.0, label=f"{deployment.label} (MEMS)")
        ax_b.loglog(cap_distances, c_opto, color=color, ls="--", lw=1.0, alpha=0.65)

    ax_b.plot(0.7, 1.0 / 60.0, "*", color="black", markersize=10, mew=0.5, mec="white", zorder=10)
    ax_b.annotate(
        "Groszek 2025\n(1 bit/min)",
        xy=(0.7, 1.0 / 60.0),
        xytext=(1.25, 4e-3),
        fontsize=6.5,
        arrowprops=dict(arrowstyle="-", color="k", lw=0.4),
    )
    ax_b.set_xlabel("Distance $d$ (m)")
    ax_b.set_ylabel("Shannon benchmark (bits/s)")
    ax_b.set_xlim(0.3, 80)
    ax_b.set_ylim(1e-3, 1e5)
    ax_b.set_title(r"Flat-noise assumptions; $B=f_c/2$. Not achieved rates.", fontsize=7)
    ax_b.grid(True, which="both", alpha=0.25, lw=0.3)
    ax_b.legend(loc="upper right", frameon=False, fontsize=5.8, handlelength=1.7)
    ax_b.text(0.03, 0.95, "(b)", transform=ax_b.transAxes, fontsize=9, fontweight="bold")
    fig.tight_layout(pad=0.35)
    _save(fig, outdir, "fig2_signal_capacity")


def _block(ax: plt.Axes, xy: tuple[float, float], width: float, height: float, text: str, color: str) -> None:
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.035,rounding_size=0.03",
        linewidth=0.8,
        edgecolor="#222222",
        facecolor=color,
    )
    ax.add_patch(box)
    ax.text(xy[0] + width / 2.0, xy[1] + height / 2.0, text, ha="center", va="center", fontsize=7.0)


def _arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], label: str | None = None) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.7,
            color="#222222",
            shrinkA=2,
            shrinkB=2,
        )
    )
    if label:
        ax.text(
            0.5 * (start[0] + end[0]),
            0.5 * (start[1] + end[1]) + 0.075,
            label,
            ha="center",
            va="bottom",
            fontsize=6.1,
            color="#333333",
        )


def figure_3(outdir: Path) -> None:
    """Modulation and demodulation signal chain."""

    fig, ax = plt.subplots(figsize=(3.45, 1.5))
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(1.15, 2.68)
    ax.axis("off")

    w = 0.78
    h = 0.48
    top_y = 2.12
    bottom_y = 1.24
    xs = [0.12, 1.10, 2.08, 3.06]

    top_labels = [
        "bits\n$b_n$",
        "$M$-ary\nCP/CF-FSK\nmapper",
        "frequency\ncommand\n$\\omega_c+\\Delta\\omega m(t)$",
        "motor +\nbearing servo",
    ]
    bottom_labels = [
        "matched filters\nor MLSE",
        "gravimeter\n$a_e(t)$",
        "selected\nfield harmonic\n$2\\phi(t)$",
        "balanced\nquadrupole\n$\\phi=\\int\\omega dt$",
    ]
    top_colors = ["#f3f3f3", "#d9ecff", "#d9ecff", "#e6f2d8"]
    bottom_colors = ["#f3f3f3", "#eadffc", "#fff7da", "#fff1c7"]

    for x, label, color in zip(xs, top_labels, top_colors):
        _block(ax, (x, top_y), w, h, label, color)
    for x, label, color in zip(xs, bottom_labels, bottom_colors):
        _block(ax, (x, bottom_y), w, h, label, color)

    for idx in range(len(xs) - 1):
        _arrow(ax, (xs[idx] + w, top_y + h / 2.0), (xs[idx + 1], top_y + h / 2.0))
    _arrow(ax, (xs[-1] + w / 2.0, top_y), (xs[-1] + w / 2.0, bottom_y + h))
    for idx in range(len(xs) - 1, 0, -1):
        _arrow(ax, (xs[idx], bottom_y + h / 2.0), (xs[idx - 1] + w, bottom_y + h / 2.0))

    # Mathematical explanation belongs in the manuscript caption, not below
    # the diagram. Crop the canvas around the two rows of blocks.
    fig.tight_layout(pad=0.2)
    _save(fig, outdir, "fig3_modulation")


def figure_4(outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(3.5, 2.05))

    br = 1.5
    volume = 1e-3
    magnetic_mass = 7.5
    magnetic_moment = br * volume / MU0
    b_amp = MU0 * magnetic_moment / (4.0 * np.pi)
    b_noise = 0.1e-15

    grav_mass = 7.5
    grav_arm = 0.5
    grav_amp = TwoMassRotor.asymmetric_dipole(grav_mass, grav_arm).harmonic_amplitude(1.0, harmonic=1)

    freq = np.logspace(-1, 4, 240)
    # Known-frequency matched filter over one full cycle: T=1/f,
    # coherent POWER SNR=A^2*T/ASD^2. These are hypothetical flat noise curves.
    b_snr = coherent_tone_power_snr(b_amp, b_noise, 1/freq)
    g_snr_mems = coherent_tone_power_snr(grav_amp, RECEIVERS[0].noise_asd, 1/freq)
    g_snr_opto = coherent_tone_power_snr(grav_amp, RECEIVERS[1].noise_asd, 1/freq)

    ax.loglog(freq, b_snr, color="#ff7f0e", label=f"B-field: 1-L NdFeB ({magnetic_mass:.1f} kg)")
    ax.loglog(freq, g_snr_mems, color="#1f77b4", label=f"g-field: {grav_mass:.1f} kg, $L$={grav_arm} m, MEMS")
    ax.loglog(freq, g_snr_opto, color="#1f77b4", ls="--",
              label=f"g-field: {grav_mass:.1f} kg, $L$={grav_arm} m, optomech.")
    ax.axhline(1.0, color="k", ls=":", lw=0.5)
    ax.set_xlabel("Frequency $f$ (Hz)")
    ax.set_ylabel("Single-period coherent power SNR")
    ax.set_xlim(0.1, 1e4)
    ax.set_ylim(1e-5, 1e27)
    ax.set_title("Equal source mass, $d=1$ m; assumed flat noise", fontsize=7)
    ax.grid(True, which="both", alpha=0.25, lw=0.3)
    ax.legend(loc="center right", frameon=False, fontsize=6.3, handlelength=1.7)
    fig.tight_layout(pad=0.3)
    _save(fig, outdir, "fig4_compare")


def figure_5(outdir: Path) -> None:
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.0, 2.45), gridspec_kw={"width_ratios": [1.25, 1.0]})

    distances = np.logspace(-0.4, 2.35, 220)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b"]
    for deployment, color in zip(DEPLOYMENTS, colors):
        cap = _exact_capacity_curve(deployment, distances, RECEIVERS[0].noise_asd)
        ax_a.loglog(distances, cap, color=color, label=deployment.label)

    ax_a.axhline(1.0, color="k", ls="--", lw=0.6, alpha=0.7)
    ax_a.axhline(0.01, color="k", ls=":", lw=0.6, alpha=0.7)
    ax_a.plot(0.7, 1.0 / 60.0, "*", color="black", markersize=9, mew=0.5, mec="white", zorder=10)
    ax_a.set_xlabel("Distance $d$ (m)")
    ax_a.set_ylabel("Shannon benchmark (bits/s)")
    ax_a.set_title(r"MEMS floor extrapolated; assumed $B=f_c/2$", fontsize=7)
    ax_a.set_xlim(0.4, 220)
    ax_a.set_ylim(1e-4, 1e6)
    ax_a.grid(True, which="both", alpha=0.25, lw=0.3)
    ax_a.legend(loc="upper right", frameon=False, fontsize=6.0, handlelength=1.6)
    ax_a.text(0.03, 0.95, "(a)", transform=ax_a.transAxes, fontsize=9, fontweight="bold")

    targets = (1.0, 0.01)
    x = np.arange(len(DEPLOYMENTS))
    width = 0.30
    offsets = np.linspace(-width/2, width/2, len(RECEIVERS))
    receiver_colors = ["#d62728", "#9467bd", "#2ca02c"]
    for ridx, receiver in enumerate(RECEIVERS):
        ranges_1 = []
        ranges_low = []
        for deployment in DEPLOYMENTS:
            rotor = deployment.rotor()
            ranges_1.append(
                solve_rotor_range_for_capacity(
                    rotor,
                    receiver.noise_asd,
                    deployment.default_bandwidth,
                    targets[0],
                    samples=2048,
                )
            )
            ranges_low.append(
                solve_rotor_range_for_capacity(
                    rotor,
                    receiver.noise_asd,
                    deployment.default_bandwidth,
                    targets[1],
                    samples=2048,
                )
            )
        xpos = x + offsets[ridx]
        ax_b.bar(xpos, ranges_low, width=width, color=receiver_colors[ridx], alpha=0.25, edgecolor="none")
        ax_b.bar(xpos, ranges_1, width=width, color=receiver_colors[ridx], alpha=0.90, label=receiver.label)

    ax_b.axhline(10.0, color="k", ls=":", lw=0.6, alpha=0.6)
    ax_b.axhline(100.0, color="k", ls=":", lw=0.6, alpha=0.6)
    ax_b.set_ylabel("Flat-noise benchmark range (m)")
    ax_b.set_yscale("log")
    ax_b.set_ylim(0.2, 300)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(["Bench", "Field", "Trailer", "Fixed", "Heavy"], rotation=25, ha="right")
    ax_b.grid(True, axis="y", which="both", alpha=0.25, lw=0.3)
    ax_b.legend(loc="upper right", frameon=False, fontsize=6.0, handlelength=1.4)
    ax_b.text(0.03, 0.95, "(b)", transform=ax_b.transAxes, fontsize=9, fontweight="bold")
    fig.tight_layout(pad=0.35)
    _save(fig, outdir, "fig5_tradeoff")


def _scaled_deployment(
    base: Deployment,
    *,
    total_mass: float | None = None,
    long_arm: float | None = None,
    tip_speed: float | None = None,
) -> Deployment:
    return Deployment(
        key=base.key,
        label=base.label,
        total_mass=base.total_mass if total_mass is None else total_mass,
        long_arm=base.long_arm if long_arm is None else long_arm,
        tip_speed=base.tip_speed if tip_speed is None else tip_speed,
        steady_power=base.steady_power,
        mass_ratio=base.mass_ratio,
    )


def _capacity_range(deployment: Deployment, noise_asd: float, target_capacity: float) -> float:
    return solve_rotor_range_for_capacity(
        deployment.rotor(),
        noise_asd,
        deployment.default_bandwidth,
        target_capacity,
        samples=4096,
    )


def sensitivity_rows() -> tuple[float, list[dict[str, float | str]]]:
    """Exact-model range sensitivity about a representative design point."""

    base = DEPLOYMENTS[3]
    receiver = RECEIVERS[1]
    target_capacity = 1.0
    baseline_range = _capacity_range(base, receiver.noise_asd, target_capacity)

    cases = [
        (
            "Receiver noise ASD",
            "10x noisier",
            _capacity_range(base, receiver.noise_asd * 10.0, target_capacity),
            "10x quieter",
            _capacity_range(base, receiver.noise_asd * 0.1, target_capacity),
        ),
        (
            "Total rotor mass",
            "0.1x mass",
            _capacity_range(_scaled_deployment(base, total_mass=0.1 * base.total_mass), receiver.noise_asd, target_capacity),
            "10x mass",
            _capacity_range(_scaled_deployment(base, total_mass=10.0 * base.total_mass), receiver.noise_asd, target_capacity),
        ),
        (
            "Arm length",
            "0.5x L",
            _capacity_range(_scaled_deployment(base, long_arm=0.5 * base.long_arm), receiver.noise_asd, target_capacity),
            "2x L",
            _capacity_range(_scaled_deployment(base, long_arm=2.0 * base.long_arm), receiver.noise_asd, target_capacity),
        ),
        (
            "Tip speed",
            "0.25x v",
            _capacity_range(_scaled_deployment(base, tip_speed=0.25 * base.tip_speed), receiver.noise_asd, target_capacity),
            "4x v",
            _capacity_range(_scaled_deployment(base, tip_speed=4.0 * base.tip_speed), receiver.noise_asd, target_capacity),
        ),
        (
            "Target bit rate",
            "100 bit/s",
            _capacity_range(base, receiver.noise_asd, 100.0),
            "0.01 bit/s",
            _capacity_range(base, receiver.noise_asd, 0.01),
        ),
    ]

    rows: list[dict[str, float | str]] = []
    for parameter, low_label, low_range, high_label, high_range in cases:
        rows.append(
            {
                "parameter": parameter,
                "low_label": low_label,
                "low_range_m": low_range,
                "low_multiplier": low_range / baseline_range,
                "high_label": high_label,
                "high_range_m": high_range,
                "high_multiplier": high_range / baseline_range,
            }
        )
    return baseline_range, rows


def figure_6(outdir: Path) -> None:
    baseline_range, rows = sensitivity_rows()

    fig, ax = plt.subplots(figsize=(3.45, 2.55))
    y = np.arange(len(rows))
    low = np.array([float(row["low_multiplier"]) for row in rows])
    high = np.array([float(row["high_multiplier"]) for row in rows])
    labels = [str(row["parameter"]) for row in rows]

    ax.hlines(y, low, high, color="#666666", lw=1.2, alpha=0.85)
    ax.plot(low, y, "o", color="#b24a4a", ms=4.2)
    ax.plot(high, y, "o", color="#2f6fba", ms=4.2)
    ax.axvline(1.0, color="k", ls="--", lw=0.7, alpha=0.65)
    ax.set_xscale("log")
    ax.set_xlim(0.35, 2.35)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_ylim(len(rows)-0.4, -0.6)
    ax.set_xlabel("Range multiplier vs. baseline")
    ax.set_title(rf"Flat Carter reference, 1 bit/s benchmark: {baseline_range:.1f} m", fontsize=7.2)
    ax.grid(True, axis="x", which="both", alpha=0.25, lw=0.3)

    for idx, row in enumerate(rows):
        ax.text(low[idx], idx + 0.22, str(row["low_label"]), ha="center", va="center", fontsize=5.7, color="#7a3030")
        ax.text(high[idx], idx - 0.22, str(row["high_label"]), ha="center", va="center", fontsize=5.7, color="#244f86")

    fig.tight_layout(pad=0.35)
    _save(fig, outdir, "fig6_sensitivity")


def noise_floor_rows(target_capacity: float = 1.0) -> list[dict[str, float | str]]:
    """Range impact of published acceleration-noise reference levels."""

    deployments = (DEPLOYMENTS[2], DEPLOYMENTS[3], DEPLOYMENTS[4])
    rows: list[dict[str, float | str]] = []
    for reference in PUBLISHED_NOISE_REFERENCES:
        for deployment in deployments:
            rows.append(
                {
                    "reference": reference.label,
                    "reference_key": reference.key,
                    "source": reference.source,
                    "reference_frequency_hz": "" if reference.frequency_hz is None else reference.frequency_hz,
                    "deployment": deployment.label,
                    "target_capacity_bps": target_capacity,
                    "noise_microgal_per_sqrt_hz": reference.noise_asd / MICROGAL,
                    "range_m": _capacity_range(deployment, reference.noise_asd, target_capacity),
                    "note": reference.note,
                }
            )
    return rows


def figure_7(outdir: Path) -> None:
    """Range penalty from deployed acceleration-noise floors."""

    target_capacity = 1.0
    noise_microgal = np.logspace(-3, 4, 65)
    deployments = (DEPLOYMENTS[2], DEPLOYMENTS[3], DEPLOYMENTS[4])
    colors = {
        "trailer": "#4c78a8",
        "fixed": "#54a24b",
        "heavy": "#b279a2",
    }

    fig, ax = plt.subplots(figsize=(3.45, 2.55))
    for deployment in deployments:
        ranges = [
            _capacity_range(deployment, value * MICROGAL, target_capacity)
            for value in noise_microgal
        ]
        ax.loglog(
            noise_microgal,
            ranges,
            lw=1.35,
            color=colors[deployment.key],
            label=deployment.label.replace(" installation", ""),
        )

    for reference in PUBLISHED_NOISE_REFERENCES:
        x = reference.noise_asd / MICROGAL
        ax.axvline(x, color="#777777", lw=0.45, alpha=0.45)

    # Carter and MEMS reference values nearly coincide after correcting units.
    label_y = [0.16, 0.55, 0.40, 0.60]
    for reference, ypos in zip(PUBLISHED_NOISE_REFERENCES, label_y):
        x = reference.noise_asd / MICROGAL
        ha = "left"
        text_x = x * 1.05
        if reference.key in ("peterson_nhnm_10hz", "gao_mems"):
            ha = "right"
            text_x = x * 0.95
        ax.text(
            text_x,
            ypos,
            reference.label,
            transform=ax.get_xaxis_transform(),
            rotation=90,
            ha=ha,
            va="bottom",
            fontsize=5.1,
            color="#444444",
        )

    ax.axhline(10.0, color="#666666", ls=":", lw=0.65, alpha=0.8)
    ax.axhline(100.0, color="#666666", ls=":", lw=0.65, alpha=0.8)
    ax.text(0.00115, 10.4, "10 m", fontsize=5.4, color="#555555")
    ax.text(0.00115, 104.0, "100 m", fontsize=5.4, color="#555555")

    ax.set_xlim(7e-4, 1.1e4)
    ax.set_ylim(0.8, 300.0)
    ax.set_xlabel(r"Acceleration-noise reference ($\mu$Gal/$\sqrt{\mathrm{Hz}}$)")
    ax.set_ylabel("1 bit/s benchmark range (m)")
    ax.grid(True, which="both", alpha=0.23, lw=0.3)
    ax.legend(loc="lower left", frameon=False, fontsize=6.0, handlelength=1.4)
    fig.tight_layout(pad=0.35)
    _save(fig, outdir, "fig7_noise_floor")


def resonant_operating_point(
    deployment: Deployment, bandwidth_hz: float = 1.0,
    model: StructuralOscillator = CARTER_OSCILLATOR,
    tune_receiver: bool = False,
) -> dict:
    """Tune source frequency toward model resonance, optionally redesign f0.

    Leave B/2 received-frequency headroom under the assumed tip-speed ceiling.
    This is a conservative frequency-excursion allowance, not a servo design.
    tune_receiver=True explicitly sets receiver f0 to the selected carrier.
    Its m, Q and readout then remain hypothetical design assumptions, not
    reported performance of a redesigned or cryogenic Carter sensor.
    """
    if not np.isfinite(bandwidth_hz) or bandwidth_hz <= 0:
        raise ValueError("bandwidth must be finite and positive")
    carrier = min(model.resonance_hz, deployment.signal_frequency-bandwidth_hz/2)
    original_resonance = model.resonance_hz
    if tune_receiver:
        model = replace(model, resonance_hz=carrier)
    band_asd = model.band_noise_upper_asd(carrier, bandwidth_hz)
    tip_speed = np.pi * deployment.long_arm * carrier
    return {
        "class": deployment.label,
        "receiver_model": "Structural oscillator; independent thermal and flat displacement-readout noise",
        "carrier_hz": carrier,
        "receiver_resonance_hz": model.resonance_hz,
        "starting_receiver_resonance_hz": original_resonance,
        "receiver_tuning_requested": tune_receiver,
        "temperature_k": model.temperature_k,
        "mass_kg": model.mass_kg,
        "quality_factor": model.quality_factor,
        "readout_displacement_asd_si": model.readout_displacement_asd,
        "on_resonance": bool(np.isclose(carrier, model.resonance_hz, rtol=0, atol=1e-10)),
        "selected_tip_speed_m_per_s": tip_speed,
        "assumed_tip_speed_ceiling_m_per_s": deployment.tip_speed,
        "channel_bandwidth_hz": bandwidth_hz,
        "carrier_thermal_asd_si": float(model.thermal_acceleration_asd(carrier)),
        "carrier_readout_asd_si": float(model.readout_acceleration_asd(carrier)),
        "carrier_total_asd_si": float(model.acceleration_noise_asd(carrier)),
        "band_noise_upper_asd_si": band_asd,
        "noise_status": "Analytical model, not measured spectrum; m/Q/readout held fixed under tuning/cooling; no environmental/control noise",
    }


def resonant_rows() -> list[dict]:
    rows = []
    for dep in DEPLOYMENTS:
        for bandwidth in (1.0, 0.01, 0.0001):
            op = resonant_operating_point(dep, bandwidth)
            asd = op["band_noise_upper_asd_si"]
            # SNR 10 benchmark gives the SAME spectral efficiency in each band.
            # Do not demand 1 bit/s from a 0.0001 Hz band and overflow 2^(R/B).
            rows.append({**op,
                         "snr": 10.0,
                         "white_noise_benchmark_bps": bandwidth*np.log2(11),
                         "range_snr10_m": solve_rotor_range_for_snr(dep.rotor(), asd, bandwidth),
                         "rate_status": "Conservative white-noise replacement benchmark; NOT rotor throughput"})
    return rows


def receiver_design_rows() -> list[dict]:
    """Explicit design hypotheses: custom resonance and cooling are allowed.

    This sweep holds m, Q, readout and mode-domain assumptions fixed. It
    calculates potential improvements, not fabrication/cryogenic performance.
    """
    scenarios = [("reference_300K", 300., False), ("matched_300K", 300., True),
                 ("matched_77K", 77., True), ("matched_4K", 4., True)]
    rows = []
    for dep in DEPLOYMENTS:
        for name, temperature, tune in scenarios:
            model = replace(CARTER_OSCILLATOR, temperature_k=temperature)
            for bandwidth in (1., .01, .0001):
                op = resonant_operating_point(dep, bandwidth, model=model, tune_receiver=tune)
                rows.append({"scenario": name, **op, "snr": 10.,
                             "white_noise_benchmark_bps": bandwidth*np.log2(11),
                             "range_snr10_m": solve_rotor_range_for_snr(dep.rotor(), op['band_noise_upper_asd_si'], bandwidth),
                             "rate_status": "Conservative white-noise benchmark; no actuator or hardware validation"})
    return rows


def figure_10_design_sweep(outdir: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 2.8))
    f0s = np.geomspace(1, 100, 150)
    for temperature in (300., 77., 4.):
        floors = [float(replace(CARTER_OSCILLATOR, temperature_k=temperature,
                                resonance_hz=f0).thermal_acceleration_asd(f0))/NANOGAL for f0 in f0s]
        ax1.loglog(f0s, floors, label=f"{temperature:g} K")
    ax1.set_xlabel("Designed resonance $f_0$ (Hz)")
    ax1.set_ylabel(r"Resonant thermal ASD (nGal/$\sqrt{\mathrm{Hz}}$)")
    ax1.set_title("Thermal limit with fixed mass and Q", fontsize=8)
    ax1.legend(fontsize=7)
    ax1.grid(True, which="both", alpha=.25)
    bands = np.geomspace(.0001, 1, 150)
    for temperature in (300., 77., 4.):
        model = replace(CARTER_OSCILLATOR, temperature_k=temperature)
        bounds = [model.band_noise_upper_asd(model.resonance_hz, b)/NANOGAL for b in bands]
        ax2.loglog(bands, bounds, label=f"{temperature:g} K")
    ax2.set_xlabel("Signal bandwidth (Hz)")
    ax2.set_ylabel(r"Band noise bound (nGal/$\sqrt{\mathrm{Hz}}$)")
    ax2.set_title("50.3 Hz example; fixed 100 fm readout", fontsize=8)
    ax2.legend(fontsize=7)
    ax2.grid(True, which="both", alpha=.25)
    fig.tight_layout(pad=.4)
    _save(fig, outdir, "fig10_receiver_design_sweep")


def figure_8_receiver(outdir: Path) -> None:
    model = CARTER_OSCILLATOR
    frequencies = np.unique(np.r_[np.geomspace(.2, 200, 700), np.linspace(48, 52, 401), model.resonance_hz])
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6))
    for ax in axes:
        ax.plot(frequencies, model.thermal_acceleration_asd(frequencies)/NANOGAL, label="Thermal (300 K)")
        ax.plot(frequencies, model.readout_acceleration_asd(frequencies)/NANOGAL, label="Readout (assumed 100 fm)")
        ax.plot(frequencies, model.acceleration_noise_asd(frequencies)/NANOGAL, color="k", label="Total model")
        ax.set_yscale("log")
        ax.set_xlabel("Frequency (Hz)")
        ax.grid(True, which="both", alpha=.25)
    axes[0].set_xscale("log")
    axes[0].set_xlim(.2, 200)
    axes[0].set_ylim(1, 2e5)
    axes[0].set_ylabel(r"Input acceleration ASD (nGal/$\sqrt{\mathrm{Hz}}$)")
    axes[0].set_title("Carter oscillator parameters; analytical model", fontsize=8)
    axes[0].legend(fontsize=6)
    axes[1].set_xlim(49.3, 51.3)
    axes[1].set_ylim(1, 100)
    op = model.band_noise_upper_asd(model.resonance_hz, 1)
    axes[1].axhline(op/NANOGAL, ls="--", color="purple", label="1 Hz band noise bound")
    axes[1].axvspan(49.8, 50.8, alpha=.1, color="purple")
    axes[1].set_title("Resonance and a 1 Hz signal band", fontsize=8)
    axes[1].legend(fontsize=6, loc="upper center")
    fig.tight_layout(pad=.5)
    _save(fig, outdir, "fig8_receiver_model")


def figure_9_resonant(outdir: Path) -> None:
    """Keep the explicit resonant calculation separate from flat references."""
    fig, ax = plt.subplots(figsize=(5.5, 3.1))
    rows = resonant_rows()
    for index, band in enumerate((1., .01, .0001)):
        values = [row["range_snr10_m"] for row in rows if row["channel_bandwidth_hz"] == band]
        ax.bar(np.arange(5)+(index-1)*.25, values, width=.25,
               label=f"B={band:g} Hz; benchmark {band*np.log2(11):.3g} bit/s")
    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(["Bench", "Field", "Trailer", "Fixed*", "Heavy"])
    ax.set_yscale("log")
    ax.set_ylabel("Range at SNR 10 (m)")
    ax.set_title("Oscillator model; conservative band-noise bound", fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=.25)
    fig.text(.5, .01, "*Reference receiver only: Fixed is off resonance. Custom tuning is in receiver_design_sweep.csv.", ha="center", fontsize=6)
    fig.tight_layout(rect=(0,.04,1,1))
    _save(fig, outdir, "fig9_resonant_ranges")


def write_tables(results_dir: Path) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)

    rows = resonant_rows()
    with (results_dir / "resonant_receiver_table.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    rows = receiver_design_rows()
    with (results_dir / "receiver_design_sweep.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    with (results_dir / "design_table.csv").open("w", newline="") as f:
        fieldnames = [
            "class",
            "total_mass_kg",
            "long_arm_m",
            "tip_speed_m_per_s",
            "signal_frequency_hz",
            "steady_power_w",
            "power_status",
            "receiver",
            "noise_asd_si",
            "noise_assumption",
            "harmonic",
            "noise_bandwidth_hz",
            "range_exact_snr10_m",
            "range_previous_mercon_m",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for deployment in DEPLOYMENTS:
            rotor = deployment.rotor()
            for receiver in RECEIVERS:
                exact = solve_rotor_range_for_snr(
                    rotor,
                    receiver.noise_asd,
                    bandwidth=1.0,
                    snr=10.0,
                    samples=4096,
                )
                writer.writerow(
                    {
                        "class": deployment.label,
                        "total_mass_kg": deployment.total_mass,
                        "long_arm_m": deployment.long_arm,
                        "tip_speed_m_per_s": deployment.tip_speed,
                        "signal_frequency_hz": deployment.signal_frequency,
                        "steady_power_w": deployment.steady_power,
                        "power_status": "Legacy assumed budget; not predicted loss power",
                        "receiver": receiver.label,
                        "noise_asd_si": receiver.noise_asd,
                        "noise_assumption": receiver.assumption,
                        "harmonic": 2,
                        "noise_bandwidth_hz": 1.0,
                        "range_exact_snr10_m": exact,
                        "range_previous_mercon_m": PREVIOUS_MERCON_RANGE_TABLE_M[deployment.key][receiver.key],
                    }
                )

    with (results_dir / "range_capacity_table.csv").open("w", newline="") as f:
        fieldnames = ["class", "receiver", "target_capacity_bps", "range_m", "carrier_hz",
                      "assumed_channel_bandwidth_hz", "noise_asd_si", "noise_assumption", "rate_status"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for deployment in DEPLOYMENTS:
            rotor = deployment.rotor()
            for receiver in RECEIVERS:
                for target in (1.0, 0.01):
                    writer.writerow(
                        {
                            "class": deployment.label,
                            "receiver": receiver.label,
                            "target_capacity_bps": target,
                            "carrier_hz": deployment.signal_frequency,
                            "assumed_channel_bandwidth_hz": deployment.default_bandwidth,
                            "noise_asd_si": receiver.noise_asd,
                            "noise_assumption": receiver.assumption,
                            "rate_status": "Unconstrained white Gaussian benchmark; no actuator or device validation",
                            "range_m": solve_rotor_range_for_capacity(
                                rotor,
                                receiver.noise_asd,
                                deployment.default_bandwidth,
                                target,
                                samples=4096,
                            ),
                        }
                    )

    baseline_range, rows = sensitivity_rows()
    with (results_dir / "sensitivity_table.csv").open("w", newline="") as f:
        fieldnames = [
            "baseline_class",
            "baseline_receiver",
            "baseline_target_capacity_bps",
            "baseline_range_m",
            "parameter",
            "low_label",
            "low_range_m",
            "low_multiplier",
            "high_label",
            "high_range_m",
            "high_multiplier",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "baseline_class": DEPLOYMENTS[3].label,
                    "baseline_receiver": RECEIVERS[1].label,
                    "baseline_target_capacity_bps": 1.0,
                    "baseline_range_m": baseline_range,
                    **row,
                }
            )

    with (results_dir / "noise_floor_table.csv").open("w", newline="") as f:
        fieldnames = [
            "reference",
            "reference_key",
            "source",
            "reference_frequency_hz",
            "deployment",
            "target_capacity_bps",
            "noise_microgal_per_sqrt_hz",
            "range_m",
            "note",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in noise_floor_rows():
            writer.writerow(row)

    with (results_dir / "modulation_torque_table.csv").open("w", newline="") as f:
        fieldnames = [
            "class",
            "moment_of_inertia_kg_m2",
            "signal_frequency_hz",
            "torque_per_signal_hz_per_second_nm",
            "power_at_start_for_1hz_1s_step_w",
            "peak_power_for_1hz_1s_step_w",
            "step_energy_j",
            "upward_step_within_assumed_tip_limit",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for deployment in DEPLOYMENTS:
            rotor = deployment.rotor()
            inertia = rotor.quadrupole_moment
            step = modulation_step(inertia, deployment.signal_frequency, 1.0, 1.0)
            writer.writerow(
                {
                    "class": deployment.label,
                    "moment_of_inertia_kg_m2": inertia,
                    "signal_frequency_hz": deployment.signal_frequency,
                    "torque_per_signal_hz_per_second_nm": np.pi * inertia,
                    "power_at_start_for_1hz_1s_step_w": step["initial_mechanical_power_w"],
                    "peak_power_for_1hz_1s_step_w": step["peak_absolute_mechanical_power_w"],
                    "step_energy_j": step["energy_change_j"],
                    "upward_step_within_assumed_tip_limit": False,
                }
            )

    with (results_dir / "rotor_mechanics_table.csv").open("w", newline="") as f:
        fieldnames = [
            "class",
            "material",
            "angular_speed_rad_per_s",
            "moment_of_inertia_kg_m2",
            "kinetic_energy_j",
            "kinetic_energy_kwh",
            "specific_energy_kj_per_kg",
            "peak_retention_force_n",
            "required_tie_area_m2",
            "equivalent_tie_diameter_m",
            "illustrative_loss_fraction_per_hour",
            "illustrative_loss_power_w",
            "legacy_assumed_total_power_w",
            "max_loss_fraction_per_hour_if_entire_budget_used",
            "model_limitations",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for deployment in DEPLOYMENTS:
            state = rotor_mechanical_state(deployment)
            writer.writerow(
                {
                    "class": deployment.label,
                    "material": state.material.name,
                    "angular_speed_rad_per_s": state.angular_speed,
                    "moment_of_inertia_kg_m2": state.moment_of_inertia,
                    "kinetic_energy_j": state.kinetic_energy,
                    "kinetic_energy_kwh": state.kinetic_energy / 3.6e6,
                    "specific_energy_kj_per_kg": state.specific_energy / 1.0e3,
                    "peak_retention_force_n": state.peak_retention_force,
                    "required_tie_area_m2": state.required_tie_area,
                    "equivalent_tie_diameter_m": state.equivalent_tie_diameter,
                    "illustrative_loss_fraction_per_hour": 0.01,
                    "illustrative_loss_power_w": energy_loss_power(state.kinetic_energy, 0.01),
                    "legacy_assumed_total_power_w": deployment.steady_power,
                    "max_loss_fraction_per_hour_if_entire_budget_used": deployment.steady_power*3600/state.kinetic_energy,
                    "model_limitations": "Tip-mass-only inertia and retention; no arm self-load, attachments, rotor dynamics or containment. 1%/h is an example, not a loss prediction.",
                }
            )

    with (results_dir / "array_subcarrier_scaling.csv").open("w", newline="") as f:
        fieldnames = [
            "num_subcarriers",
            "same_range_per_actuator_torque_multiplier",
            "same_range_total_torque_multiplier",
            "same_range_total_mass_multiplier",
            "fixed_mass_per_actuator_torque_multiplier",
            "fixed_mass_total_torque_multiplier",
            "fixed_mass_range_multiplier",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for n in (1, 2, 4, 8, 16, 32, 64, 100):
            writer.writerow(
                {
                    "num_subcarriers": n,
                    "same_range_per_actuator_torque_multiplier": n ** (-2.5),
                    "same_range_total_torque_multiplier": n ** (-1.5),
                    "same_range_total_mass_multiplier": n ** 0.5,
                    "fixed_mass_per_actuator_torque_multiplier": n ** (-3.0),
                    "fixed_mass_total_torque_multiplier": n ** (-2.0),
                    "fixed_mass_range_multiplier": n ** (-0.125),
                }
            )

    with (results_dir / "bench_demonstrator_table.csv").open("w", newline="") as f:
        fieldnames = [
            "distance_m",
            "peak_acceleration_microgal",
            "carrier_amplitude_microgal",
            "total_ac_rms_microgal",
            "noise_assumption",
            "snr_1hz_band",
            "bandwidth_for_snr10_hz",
            "one_over_bandwidth_seconds",
            "coherent_time_for_snr10_seconds",
            "coherent_time_convention",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        bench = DEPLOYMENTS[0]
        rotor = bench.rotor()
        noise_asd = RECEIVERS[0].noise_asd
        target_snr = 10.0
        for distance in (0.5, 0.7, 1.0):
            amplitude = rotor.harmonic_amplitude(distance, samples=4096)
            bandwidth_for_snr10 = amplitude * amplitude / (2.0 * target_snr * noise_asd * noise_asd)
            writer.writerow(
                {
                    "distance_m": distance,
                    "peak_acceleration_microgal": rotor.peak_ac_amplitude(distance)/MICROGAL,
                    "carrier_amplitude_microgal": amplitude / MICROGAL,
                    "total_ac_rms_microgal": rotor.rms_ac_acceleration(distance)/MICROGAL,
                    "noise_assumption": RECEIVERS[0].assumption,
                    "snr_1hz_band": float(snr_power(amplitude, noise_asd, 1.0)),
                    "bandwidth_for_snr10_hz": bandwidth_for_snr10,
                    "one_over_bandwidth_seconds": 1.0 / bandwidth_for_snr10,
                    "coherent_time_for_snr10_seconds": np.ceil(target_snr*noise_asd**2/amplitude**2 * bench.signal_frequency)/bench.signal_frequency,
                    "coherent_time_convention": "Known-frequency tone; rounded up to an integer number of carrier cycles",
                }
            )


def reproduce(output_dir: Path, results_dir: Path) -> None:
    apply_paper_style()
    figure_1(output_dir)
    figure_2_quadrupole_validity(output_dir)
    figure_2(output_dir)
    figure_3(output_dir)
    figure_4(output_dir)
    figure_5(output_dir)
    figure_6(output_dir)
    figure_7(output_dir)
    figure_8_receiver(output_dir)
    figure_9_resonant(output_dir)
    figure_10_design_sweep(output_dir)
    write_tables(results_dir)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    # Never silently replace paper figures while its captions still describe
    # the old peak-power/receiver assumptions. Output location is cwd-independent.
    root = Path(__file__).resolve().parents[2]
    corrected = root / "results" / "corrected_model"
    parser.add_argument("--output-dir", type=Path, default=corrected / "Figures", help="figure output directory")
    parser.add_argument("--results-dir", type=Path, default=corrected, help="CSV output directory")
    args = parser.parse_args(argv)
    reproduce(args.output_dir, args.results_dir)
    print(f"Wrote figures to {args.output_dir}")
    print(f"Wrote tables to {args.results_dir}")


if __name__ == "__main__":
    main()
