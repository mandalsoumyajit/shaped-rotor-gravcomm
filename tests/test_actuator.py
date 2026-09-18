"""Analytic dynamics, four-quadrant limits, conservation and convergence."""
from dataclasses import replace
import unittest
import numpy as np
from gravcomm.actuator import Actuator, simulate_actuator
from gravcomm.rotor import TwoMassRotor


class ActuatorTest(unittest.TestCase):
    def setUp(self):
        self.a = Actuator(2., 100., 80., 0., 1e6, 1e6, 100.)
        self.t = np.linspace(0, 2, 201)

    def test_constant_torque_phase_and_work(self):
        r = simulate_actuator(self.a, self.t, np.full(201, 4.), 3.)
        np.testing.assert_allclose(r.speed_rad_s, 3+2*self.t, atol=1e-12)
        np.testing.assert_allclose(r.phase_rad, 3*self.t+self.t**2, atol=1e-12)
        self.assertAlmostEqual(r.diagnostics['rms_torque_nm'], 4.)
        self.assertAlmostEqual(r.diagnostics['motoring_work_j'], 40.)
        self.assertAlmostEqual(r.diagnostics['energy_balance_error_j'], 0., places=10)

    def test_filtered_torque_analytic_phase(self):
        tau = .2
        a = replace(self.a, response_time_s=tau)
        r = simulate_actuator(a, self.t, np.full(201, 4.), 3.)
        decay = np.exp(-self.t/tau)
        np.testing.assert_allclose(r.actual_torque_nm, 4*(1-decay), atol=1e-7)
        np.testing.assert_allclose(r.speed_rad_s, 3+2*(self.t-tau*(1-decay)), atol=1e-8)
        np.testing.assert_allclose(r.phase_rad, 3*self.t+self.t**2-2*tau*self.t+2*tau**2*(1-decay), atol=1e-8)

    def test_power_limited_acceleration_and_braking(self):
        a = replace(self.a, motoring_power_w=10., regenerative_power_w=6.)
        for command, power in ((100., 10.), (-100., -6.)):
            r = simulate_actuator(a, self.t, np.full(201, command), 10.)
            np.testing.assert_allclose(r.speed_rad_s, np.sqrt(100+power*self.t), atol=1e-10)
            np.testing.assert_allclose(r.mechanical_power_w, power, atol=1e-12)
            self.assertAlmostEqual(r.diagnostics['energy_balance_error_j'], 0., places=9)

    def test_four_quadrants_and_zero_regen(self):
        a = replace(self.a, motoring_power_w=20., regenerative_power_w=10.)
        for omega in (-10., 10.):
            for command in (-500., 500.):
                power = a.limited_torque(command, omega)*omega
                self.assertAlmostEqual(power, 20. if omega*command > 0 else -10.)
        self.assertEqual(a.limited_torque(500, 0), 100.)
        self.assertEqual(replace(a, regenerative_power_w=0).limited_torque(-10, 10), 0.)

    def test_drag_coastdown_and_energy(self):
        r = simulate_actuator(replace(self.a, viscous_drag_nm_s=1.), self.t, np.zeros(201), 10.)
        np.testing.assert_allclose(r.speed_rad_s, 10*np.exp(-self.t/2), atol=1e-10)
        self.assertAlmostEqual(r.diagnostics['drag_heat_j'], 100*(1-np.exp(-2)), places=8)
        self.assertAlmostEqual(r.diagnostics['energy_balance_error_j'], 0., places=8)

    def test_speed_feedback_closed_form(self):
        target = np.full(201, 10/np.pi)
        r = simulate_actuator(self.a, self.t, np.zeros(201), 5.,
                              target_carrier_hz=target, speed_gain_nm_s=4.)
        np.testing.assert_allclose(r.speed_rad_s, 10-5*np.exp(-2*self.t), atol=1e-8)

    def test_violations_are_reported_without_speed_clamping(self):
        a = replace(self.a, continuous_torque_nm=1., max_speed_rad_s=4.)
        r = simulate_actuator(a, self.t, np.full(201, 200.), 3.)
        self.assertTrue(r.diagnostics['speed_limit_exceeded'])
        self.assertTrue(r.diagnostics['continuous_torque_exceeded'])
        self.assertGreater(r.speed_rad_s[-1], 4.)
        self.assertLessEqual(r.diagnostics['peak_torque_nm'], 100.)

    def test_signal_uses_achieved_phase_and_phasor_sign(self):
        rotor = TwoMassRotor.balanced_quadrupole(50, .3)
        r = simulate_actuator(self.a, self.t, np.full(201, 4.), 3., initial_phase_rad=.4)
        c = rotor.harmonic_coefficient(.5, 2)
        phi = .4+3*self.t+self.t**2
        np.testing.assert_allclose(r.acceleration(rotor, .5), c.real*np.cos(2*phi)-c.imag*np.sin(2*phi), atol=1e-20)
        np.testing.assert_allclose(r.acceleration(rotor, .5, harmonic=None), rotor.acceleration_component(.5, phi), atol=1e-20)

    def test_clipped_filtered_cycle_convergence(self):
        a = replace(self.a, response_time_s=.05, motoring_power_w=20., regenerative_power_w=12.)
        command = 200*np.sin(2*np.pi*self.t)
        coarse = simulate_actuator(a, self.t, command, 10., max_step_s=.001)
        fine = simulate_actuator(a, self.t, command, 10., max_step_s=.0002)
        np.testing.assert_allclose(coarse.speed_rad_s, fine.speed_rad_s, atol=3e-5, rtol=0)
        self.assertLess(abs(fine.diagnostics['energy_balance_error_j']), 1e-6)

    def test_invalid_inputs(self):
        for kwargs in ({'inertia_kg_m2': 0}, {'response_time_s': -1}, {'regenerative_power_w': float('nan')},
                       {'continuous_torque_nm': 101}):
            with self.assertRaises(ValueError):
                replace(self.a, **kwargs)
        with self.assertRaises(ValueError):
            simulate_actuator(self.a, [0, 0], [1, 1], 3.)
        with self.assertRaises(ValueError):
            simulate_actuator(self.a, [0, 1], [1, 1], 3., speed_gain_nm_s=1.)


if __name__ == '__main__':
    unittest.main()
