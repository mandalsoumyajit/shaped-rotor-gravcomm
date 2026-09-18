"""Rigid-rotor actuator dynamics; all drive ratings are explicit design inputs.

The low-pass state represents torque demand before an ideal instantaneous
torque/power limiter. This is a reduced drive model, not an electrical motor
model: envelope clipping can change torque faster than the demand filter.
RMS torque is a cycle-duty proxy, not a winding-temperature prediction.
"""
from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class Actuator:
    inertia_kg_m2: float
    peak_torque_nm: float
    continuous_torque_nm: float
    response_time_s: float
    motoring_power_w: float
    regenerative_power_w: float
    max_speed_rad_s: float
    viscous_drag_nm_s: float = 0.0

    def __post_init__(self):
        for name in ("inertia_kg_m2", "peak_torque_nm", "continuous_torque_nm",
                     "motoring_power_w", "max_speed_rad_s"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("response_time_s", "regenerative_power_w", "viscous_drag_nm_s"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.continuous_torque_nm > self.peak_torque_nm:
            raise ValueError("continuous torque cannot exceed peak torque")

    def limited_torque(self, demand, omega):
        """Four-quadrant, symmetric torque limit with asymmetric power ratings.

        Constant torque below base speed, constant mechanical power above it.
        No explicit voltage/back-EMF calculation; supply a measured envelope
        via a subclass if this approximation is insufficient. Regenerative
        rating is mechanical absorption capacity, not recovered grid power.
        """
        limit = self.peak_torque_nm
        if omega != 0:
            power = self.motoring_power_w if demand * omega >= 0 else self.regenerative_power_w
            limit = min(limit, power / abs(omega))
        return float(np.clip(demand, -limit, limit))


@dataclass
class ActuatorTrace:
    time_s: np.ndarray
    phase_rad: np.ndarray
    speed_rad_s: np.ndarray
    command_torque_nm: np.ndarray
    filtered_demand_nm: np.ndarray
    actual_torque_nm: np.ndarray
    mechanical_power_w: np.ndarray
    diagnostics: dict

    def carrier_hz(self, harmonic=2):
        _harmonic(harmonic)
        return harmonic * self.speed_rad_s / (2 * np.pi)

    def acceleration(self, rotor, distance, harmonic=2, component="x"):
        """Selected AC harmonic using achieved phase and exact spatial phasor.

        harmonic=None returns the full Newtonian field INCLUDING static field.
        No receiver transfer function or noise is implicitly applied.
        """
        if not math.isfinite(distance) or distance <= rotor.max_arm:
            raise ValueError("receiver must lie outside the rotor")
        if harmonic is None:
            return rotor.acceleration_component(distance, self.phase_rad, component)
        _harmonic(harmonic)
        c = rotor.harmonic_coefficient(distance, harmonic, component=component)
        return np.real(c * np.exp(1j * harmonic * self.phase_rad))


def _harmonic(n):
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("harmonic must be a positive integer")


def simulate_actuator(actuator, time_s, torque_command_nm, initial_speed_rad_s,
                      *, initial_phase_rad=0., initial_demand_nm=0.,
                      target_carrier_hz=None, speed_gain_nm_s=0., harmonic=2,
                      max_step_s=None):
    """Integrate a linearly interpolated torque command using subdivided RK4.

    Optional proportional speed feedback is added BEFORE demand filtering.
    The caller supplies feedforward torque, including desired drag compensation.
    target_carrier_hz is linearly interpolated; no hidden controller or derivative
    estimate is supplied. Zero response time gives ideal instantaneous demand.

    State equations: phi'=omega, I*omega'=torque-b*omega,
    demand'=(command+Kp*(omega_target-omega)-demand)/response_time.
    Positive/negative work, drag heat and torque squared are integrated with
    the same RK4 stages. Power/torque bounds hold at every evaluation. Speed
    and cycle-RMS limits are assessed, not imposed by unphysical speed clipping.
    Refine max_step_s for convergence near clipping boundaries. Output samples
    must resolve the carrier if the acceleration waveform will be processed.
    """
    _harmonic(harmonic)
    t = np.asarray(time_s, dtype=float)
    command = np.asarray(torque_command_nm, dtype=float)
    if (t.ndim != 1 or len(t) < 2 or command.shape != t.shape
            or not np.all(np.isfinite(t)) or not np.all(np.isfinite(command))
            or np.any(np.diff(t) <= 0)):
        raise ValueError("time and command must be finite equal-length vectors with increasing time")
    if not all(math.isfinite(x) for x in (initial_speed_rad_s, initial_phase_rad,
                                         initial_demand_nm, speed_gain_nm_s)) or speed_gain_nm_s < 0:
        raise ValueError("initial states must be finite; speed gain must be finite and nonnegative")
    target = None if target_carrier_hz is None else np.asarray(target_carrier_hz, dtype=float)
    if target is not None and (target.shape != t.shape or not np.all(np.isfinite(target))):
        raise ValueError("target carrier must be a finite vector matching time")
    if speed_gain_nm_s and target is None:
        raise ValueError("speed feedback requires a target carrier")
    step = float(np.min(np.diff(t))) if max_step_s is None else max_step_s
    if not math.isfinite(step) or step <= 0:
        raise ValueError("max_step_s must be finite and positive")
    # Resolve both actuator and closed-loop dynamics, including large P gains.
    if actuator.response_time_s:
        step = min(step, actuator.response_time_s / 10)
    rate = (speed_gain_nm_s + actuator.viscous_drag_nm_s) / actuator.inertia_kg_m2
    if rate:
        step = min(step, .1 / rate)
    # phi, omega, filtered demand, motoring work, absorbed work, drag heat, int(torque^2)
    y = np.array([initial_phase_rad, initial_speed_rad_s, initial_demand_nm, 0., 0., 0., 0.])
    states = np.empty((len(t), 7))
    demands = np.empty(len(t))
    torques = np.empty(len(t))
    peak_torque = peak_motoring = peak_regen = peak_speed = 0.

    def evaluate(now, state):
        nonlocal peak_torque, peak_motoring, peak_regen, peak_speed
        omega = state[1]
        cmd = float(np.interp(now, t, command))
        if target is not None:
            cmd += speed_gain_nm_s * (2*np.pi*np.interp(now, t, target)/harmonic - omega)
        demand = state[2] if actuator.response_time_s else cmd
        torque = actuator.limited_torque(demand, omega)
        power = torque * omega
        peak_torque = max(peak_torque, abs(torque))
        peak_motoring = max(peak_motoring, power)
        peak_regen = max(peak_regen, -power)
        peak_speed = max(peak_speed, abs(omega))
        derivative = np.array([omega, (torque-actuator.viscous_drag_nm_s*omega)/actuator.inertia_kg_m2,
                               (cmd-demand)/actuator.response_time_s if actuator.response_time_s else 0.,
                               max(power, 0.), max(-power, 0.), actuator.viscous_drag_nm_s*omega**2,
                               torque**2])
        return derivative, demand, torque

    states[0] = y
    _, demands[0], torques[0] = evaluate(t[0], y)
    for i in range(1, len(t)):
        count = max(1, math.ceil((t[i]-t[i-1])/step))
        h = (t[i]-t[i-1])/count
        for j in range(count):
            now = t[i-1]+j*h
            k1 = evaluate(now, y)[0]
            k2 = evaluate(now+h/2, y+h*k1/2)[0]
            k3 = evaluate(now+h/2, y+h*k2/2)[0]
            k4 = evaluate(now+h, y+h*k3)[0]
            y = y+h*(k1+2*k2+2*k3+k4)/6
            if not np.all(np.isfinite(y)):
                raise ValueError("nonfinite integration state; check scales and integration step")
        states[i] = y
        _, demands[i], torques[i] = evaluate(t[i], y)
    kinetic_change = .5*actuator.inertia_kg_m2*(y[1]**2-initial_speed_rad_s**2)
    rms = math.sqrt(max(0., y[6]/(t[-1]-t[0])))
    diagnostics = dict(peak_torque_nm=peak_torque, rms_torque_nm=rms,
                       peak_motoring_power_w=peak_motoring, peak_regenerative_power_w=peak_regen,
                       motoring_work_j=y[3], absorbed_work_j=y[4], drag_heat_j=y[5],
                       kinetic_energy_change_j=kinetic_change,
                       energy_balance_error_j=y[3]-y[4]-y[5]-kinetic_change,
                       peak_speed_rad_s=peak_speed,
                       speed_limit_exceeded=bool(peak_speed > actuator.max_speed_rad_s),
                       continuous_torque_exceeded=bool(rms > actuator.continuous_torque_nm))
    if target is not None:
        error = harmonic*states[:, 1]/(2*np.pi)-target
        diagnostics["max_sampled_carrier_error_hz"] = float(np.max(abs(error)))
        diagnostics["final_carrier_error_hz"] = float(error[-1])
    return ActuatorTrace(t.copy(), states[:, 0], states[:, 1], command.copy(), demands,
                         torques, torques*states[:, 1], diagnostics)
