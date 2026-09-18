"""Independent physics checks, not just round trips through model helpers."""
import math
import unittest
from dataclasses import replace

import numpy as np

from gravcomm.constants import G, K_B, NANOG, NANOGAL
from gravcomm.designs import DEPLOYMENTS, RECEIVERS
from gravcomm.figures import resonant_operating_point
from gravcomm.geophone import Gradiometer, SM24
from gravcomm.link_budget import (
    capacity_bps, coherent_tone_power_snr, multicarrier_capacity_bps,
    snr_power, solve_monotonic_range, solve_rotor_range_for_snr,
    solve_rotor_range_for_capacity,
)
from gravcomm.mechanics import energy_loss_power, modulation_step
from gravcomm.receivers import CARTER_OSCILLATOR, StructuralOscillator
from gravcomm.rotor import TwoMassRotor


class HarmonicsTest(unittest.TestCase):
    def test_equal_mass_against_independent_gauss_quadrature(self):
        # Closed-form equal-dumbbell radial field, independently integrated
        # using nonuniform Gauss-Legendre nodes (production uses uniform DFT).
        m, L, d = 50., .3, .5
        nodes, weights = np.polynomial.legendre.leggauss(256)
        phi = math.pi * (nodes+1)
        c = np.cos(phi)
        exact = -G*m/2 * ((d-L*c)/(d*d+L*L-2*d*L*c)**1.5
                         +(d+L*c)/(d*d+L*L+2*d*L*c)**1.5)
        a2 = abs(np.sum(weights*exact*np.cos(2*phi)))
        rotor = TwoMassRotor.balanced_quadrupole(m, L)
        np.testing.assert_allclose(rotor.harmonic_amplitude(d), a2, rtol=1e-11, atol=0)
        self.assertAlmostEqual(a2/1e-8, 1.4890001001410673, places=9)
        self.assertAlmostEqual(float(snr_power(a2, 1e-9, 1)), 110.856064910005, places=6)

    def test_parseval_and_peak_is_not_carrier(self):
        rotor = DEPLOYMENTS[0].rotor()
        _, ac = rotor.ac_signal(.5, samples=8192)
        coeff = 2*np.fft.rfft(ac)/len(ac)
        power = np.sum(abs(coeff[1:-1])**2)/2 + abs(coeff[-1])**2/4
        np.testing.assert_allclose(power, np.mean(ac**2), rtol=1e-13, atol=0)
        self.assertGreater(rotor.peak_ac_amplitude(.5), 1.7*rotor.harmonic_amplitude(.5))
        np.testing.assert_allclose(rotor.rms_ac_acceleration(.5)**2, power, rtol=1e-12, atol=0)

    def test_balanced_is_not_necessarily_half_turn_symmetric(self):
        symmetric = TwoMassRotor.balanced_quadrupole(600, 1)
        unequal = TwoMassRotor.balanced_quadrupole(600, 1, 5)
        self.assertEqual(unequal.dipole_moment, 0)
        self.assertLess(symmetric.harmonic_amplitude(1.5, 1)/symmetric.harmonic_amplitude(1.5, 2), 1e-12)
        self.assertAlmostEqual(unequal.harmonic_amplitude(1.5, 1)/1e-8, .22371514882047808, places=9)

    def test_quadrupole_radial_and_transverse_asymptotes(self):
        rotor = TwoMassRotor.balanced_quadrupole(50, .3)
        d = 30.
        scale = G*rotor.quadrupole_moment/d**4
        self.assertAlmostEqual(rotor.harmonic_amplitude(d)/scale, 9/4, delta=.0004)
        self.assertAlmostEqual(rotor.harmonic_amplitude(d, component="z")/scale, 3/2, delta=.0003)

    def test_convergence_near_rotor(self):
        for ratio in (1., 5.):
            rotor = TwoMassRotor.balanced_quadrupole(600, 1, ratio)
            for d in (1.02, 1.5, 10.):
                np.testing.assert_allclose(rotor.harmonic_amplitude(d, samples=2048),
                                           rotor.harmonic_amplitude(d, samples=32768), rtol=2e-10, atol=0)

    def test_gradiometer_matches_difference_waveform(self):
        rotor = TwoMassRotor.balanced_quadrupole(600, 1, 5)
        grad = Gradiometer(SM24, .3)
        phi, near = rotor.ac_signal(1.5-.15)
        _, far = rotor.ac_signal(1.5+.15)
        differential = near-far
        expected = 2*math.hypot(np.mean(differential*np.cos(2*phi)),
                               np.mean(differential*np.sin(2*phi)))
        np.testing.assert_allclose(grad.differential_signal_amplitude(rotor, 1.5), expected, rtol=1e-12, atol=0)

    def test_range_solves_carrier_not_peak(self):
        rotor = DEPLOYMENTS[0].rotor()
        d = solve_rotor_range_for_snr(rotor, 1e-9, 1)
        self.assertAlmostEqual(float(snr_power(rotor.harmonic_amplitude(d), 1e-9, 1)), 10, places=10)
        self.assertGreater(float(snr_power(rotor.peak_ac_amplitude(d), 1e-9, 1)), 10)

    def test_unattainable_range_rejected(self):
        with self.assertRaises(ValueError):
            solve_monotonic_range(lambda d: 1/d**4, 2, lower=1)


class ReceiverTest(unittest.TestCase):
    def test_nanog_is_not_nanogal(self):
        self.assertAlmostEqual(.1*NANOG/NANOGAL, 98.0665)
        self.assertAlmostEqual(RECEIVERS[1].noise_asd/NANOGAL, 98.0665)
        self.assertEqual(len(RECEIVERS), 2)  # No universal constant resonant tier.

    def test_resonance_against_fluctuation_dissipation(self):
        r = CARTER_OSCILLATOR
        # At resonance both structural and viscous models give sqrt(4kT*b)/m,
        # where b=m*omega0/Q. Independent dimensional formulation.
        b = r.mass_kg*2*math.pi*r.resonance_hz/r.quality_factor
        expected = math.sqrt(4*K_B*r.temperature_k*b)/r.mass_kg
        np.testing.assert_allclose(r.thermal_acceleration_asd(r.resonance_hz), expected, rtol=1e-14, atol=0)
        self.assertAlmostEqual(expected/NANOGAL, 5.14938813486474, places=10)
        np.testing.assert_allclose(r.displacement_per_acceleration(r.resonance_hz),
                                   r.quality_factor/(2*math.pi*r.resonance_hz)**2, rtol=1e-14, atol=0)

    def test_resonance_reduces_readout_not_thermal_by_Q(self):
        r = CARTER_OSCILLATOR
        self.assertLess(r.readout_acceleration_asd(r.resonance_hz), r.readout_acceleration_asd(40)/10000)
        np.testing.assert_allclose(r.thermal_acceleration_asd(40)/r.thermal_acceleration_asd(50.3),
                                   math.sqrt(50.3/40), rtol=1e-14)
        self.assertGreaterEqual(r.acceleration_noise_asd(50.3), r.thermal_acceleration_asd(50.3))

    def test_band_bound_covers_spectrum_and_allows_more_than_linewidth(self):
        r = CARTER_OSCILLATOR
        for fc, B in ((50.3, 1.), (50.3, .01), (39., 2.)):
            f = np.linspace(fc-B/2, fc+B/2, 10001)
            self.assertGreaterEqual(r.band_noise_upper_asd(fc, B), np.max(r.acceleration_noise_asd(f)))
        self.assertGreater(r.band_noise_upper_asd(50.3, 1), r.band_noise_upper_asd(50.3, .01))
        with self.assertRaises(ValueError):
            r.acceleration_noise_asd(265)
        with self.assertRaises(ValueError):
            r.band_noise_upper_asd(50.3, 110)

    def test_noise_zero_readout_limit(self):
        r = StructuralOscillator(readout_displacement_asd=0)
        np.testing.assert_allclose(r.acceleration_noise_asd([40., 50.3, 60.]),
                                   r.thermal_acceleration_asd([40., 50.3, 60.]), rtol=1e-14)

    def test_tuning_does_not_exceed_source_speed(self):
        for dep in DEPLOYMENTS:
            op = resonant_operating_point(dep, 1)
            self.assertLessEqual(op['carrier_hz']+.5, dep.signal_frequency+1e-12)
            self.assertLessEqual(op['selected_tip_speed_m_per_s'], dep.tip_speed)
            self.assertEqual(op['on_resonance'], dep.key != 'fixed')

    def test_custom_receiver_tunes_to_fixed_class(self):
        dep = DEPLOYMENTS[3]
        op = resonant_operating_point(dep, 1., tune_receiver=True)
        self.assertTrue(op['on_resonance'])
        self.assertAlmostEqual(op['receiver_resonance_hz'], dep.signal_frequency-.5)
        self.assertEqual(op['quality_factor'], CARTER_OSCILLATOR.quality_factor)
        reference = resonant_operating_point(dep, 1.)
        self.assertLess(op['carrier_total_asd_si'], reference['carrier_total_asd_si'])

    def test_cooling_and_lower_resonance_scaling(self):
        r = CARTER_OSCILLATOR
        base = r.thermal_acceleration_asd(r.resonance_hz)
        for temperature in (77., 4.):
            cold = replace(r, temperature_k=temperature)
            np.testing.assert_allclose(cold.thermal_acceleration_asd(cold.resonance_hz)/base,
                                       math.sqrt(temperature/300), rtol=1e-14)
            # Cooling has NOT silently reduced optical readout noise.
            self.assertEqual(cold.readout_acceleration_asd(40), r.readout_acceleration_asd(40))
        lower = replace(r, resonance_hz=20.)
        np.testing.assert_allclose(lower.thermal_acceleration_asd(20)/base,
                                   math.sqrt(20/50.3), rtol=1e-14)


class LinkAndEnergyTest(unittest.TestCase):
    def test_multicarrier_matches_two_independent_white_channels(self):
        # Noise unit changes must cancel from power SNR; a missing ASD square fails.
        expected = 3*math.log2(1+4/(2*.25*3)) + 2*math.log2(1+9/(2*4*2))
        self.assertAlmostEqual(multicarrier_capacity_bps([2e-9,3e-9],[.5e-9,2e-9],[3,2]), expected, places=12)

    def test_low_snr_capacity_numerically_stable(self):
        self.assertAlmostEqual(float(capacity_bps(1e-12, 1, 1))/ (0.5e-24/math.log(2)), 1, places=12)
        d = solve_rotor_range_for_capacity(DEPLOYMENTS[0].rotor(), 1e-9, 10, .01)
        self.assertAlmostEqual(float(capacity_bps(DEPLOYMENTS[0].rotor().harmonic_amplitude(d), 1e-9, 10)), .01, places=12)

    def test_single_cycle_against_filter_weight_energy(self):
        # White one-sided PSD N^2: var(output)=N^2/2 * integral(w^2 dt).
        # w=2cos(2pi*f*t)/T integrates a known-phase tone to its peak amplitude.
        T, A, noise = .01, 3e-9, 2e-9
        t = np.linspace(0,T,10000,endpoint=False)
        weights = 2/T*np.cos(2*math.pi*t/T)
        variance = noise**2/2 * np.mean(weights**2)*T
        self.assertAlmostEqual(float(coherent_tone_power_snr(A,noise,T)), A*A/variance, places=14)
        self.assertAlmostEqual(float(coherent_tone_power_snr(A,noise,2*T))/float(coherent_tone_power_snr(A,noise,T)), 2)

    def test_energy_loss_units(self):
        self.assertAlmostEqual(energy_loss_power(1.5625e6,.01), 4.340277777777778)
        self.assertAlmostEqual(energy_loss_power(1.5625e9,.01), 4340.277777777778)
        self.assertAlmostEqual(energy_loss_power(1e11,.01), 277777.7777777778)
        self.assertEqual(energy_loss_power(1e11,0), 0)

    def test_modulation_work_energy_identity_acceleration_and_braking(self):
        I, fc, T = 2e5, 125/math.pi, 1.
        for delta in (1., -1.):
            result = modulation_step(I,fc,delta,T)
            omega1, omega2 = math.pi*fc, math.pi*(fc+delta)
            work = result['torque_nm']*.5*(omega1+omega2)*T
            self.assertAlmostEqual(work/result['energy_change_j'], 1, places=12)
            self.assertAlmostEqual(result['torque_nm']/(math.pi*I*delta), 1, places=12)
        self.assertAlmostEqual(modulation_step(I,fc,1,T)['initial_mechanical_power_w'], 78539816.33974482, places=5)


if __name__ == '__main__':
    unittest.main()
