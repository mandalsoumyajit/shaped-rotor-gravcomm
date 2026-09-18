"""Reproducible illustrative drive comparison; ratings are NOT hardware claims."""
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gravcomm.actuator import Actuator, simulate_actuator
from gravcomm.rotor import TwoMassRotor


def main():
    out = ROOT / 'results' / 'corrected_model' / 'actuator'
    out.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0, 4, 4001)
    # Smooth raised-cosine excursion: +0.4 Hz in 1 s, hold, then return.
    fc = np.full_like(t, 35.)
    dfdt = np.zeros_like(t)
    for start, sign in ((.25, 1), (2.25, -1)):
        u = np.clip(t-start, 0, 1)
        fc += sign*.2*(1-np.cos(np.pi*u))
        dfdt += sign*.2*np.pi*np.sin(np.pi*u)
    drive = Actuator(inertia_kg_m2=2e5, peak_torque_nm=250000.,
                     continuous_torque_nm=150000., response_time_s=0.,
                     motoring_power_w=25e6, regenerative_power_w=20e6,
                     max_speed_rad_s=125., viscous_drag_nm_s=10.)
    command = drive.inertia_kg_m2*np.pi*dfdt + drive.viscous_drag_nm_s*np.pi*fc
    rotor = TwoMassRotor.balanced_quadrupole(50000., 2.)
    cases = [('instantaneous', 0., 0.), ('filtered_50ms', .05, 0.),
             ('filtered_200ms', .2, 0.), ('filtered_50ms_feedback', .05, 8e5)]
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    axes[0].plot(t, fc, 'k--', label='Target')
    axes[1].plot(t, command/1e3, 'k--', label='Feedforward demand')
    summaries = []
    for name, response, gain in cases:
        a = replace(drive, response_time_s=response)
        trace = simulate_actuator(a, t, command, np.pi*fc[0],
                                  initial_demand_nm=command[0], target_carrier_hz=fc,
                                  speed_gain_nm_s=gain, max_step_s=.0005)
        signal = trace.acceleration(rotor, 10.)
        with (out / f'{name}.csv').open('w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['time_s', 'target_carrier_hz', 'actual_carrier_hz', 'phase_rad',
                             'feedforward_torque_nm', 'filtered_demand_nm', 'actual_torque_nm',
                             'mechanical_power_w', 'harmonic2_acceleration_m_s2'])
            writer.writerows(zip(t, fc, trace.carrier_hz(), trace.phase_rad, command,
                                 trace.filtered_demand_nm, trace.actual_torque_nm,
                                 trace.mechanical_power_w, signal))
        summaries.append(dict(case=name, **trace.diagnostics))
        axes[0].plot(t, trace.carrier_hz(), label=name)
        axes[1].plot(t, trace.actual_torque_nm/1e3)
        axes[2].plot(t, trace.mechanical_power_w/1e6)
    with (out / 'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    (out / 'assumptions.json').write_text(json.dumps(dict(
        description='Illustrative assumed drive, not a selected or validated motor. No BER claim.',
        base_drive=asdict(drive), cases=cases, distance_m=10., harmonic=2,
        source_total_mass_kg=50000., source_arm_m=2., max_step_s=.0005,
        target='35 Hz base; raised-cosine +0.4 Hz at 0.25-1.25 s, return at 2.25-3.25 s',
        speed_gain_units='N m per rad/s', initial_demand='steady drag compensation'), indent=2))
    axes[0].set_ylabel('Carrier (Hz)')
    axes[1].set_ylabel('Torque (kN m)')
    axes[2].set_ylabel('Mechanical power (MW)')
    axes[2].set_xlabel('Time (s)')
    axes[0].legend(fontsize=8, ncol=2)
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.25)
    fig.suptitle('Illustrative constrained actuator: assumed ratings, no demonstrated bit rate')
    fig.tight_layout()
    fig.savefig(out / 'actuator_comparison.png', dpi=160)
    fig.savefig(out / 'actuator_comparison.pdf')
    plt.close(fig)
    print(out)


if __name__ == '__main__':
    main()
