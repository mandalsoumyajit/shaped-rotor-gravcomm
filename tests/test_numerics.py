import math
import unittest

from gravcomm.constants import MICROGAL
from gravcomm.designs import DEPLOYMENTS
from gravcomm.geophone import SM24, Geophone, Gradiometer
from gravcomm.link_budget import amplitude_for_snr, snr_power
from gravcomm.mechanics import (
    ALUMINUM,
    CarrierPlate,
    cantilever_bending_stiffness,
    disk_rotor_mechanical_state,
    rotor_mechanical_state,
    whirl_critical_hz,
)
from gravcomm.noise import PUBLISHED_NOISE_REFERENCES, peterson_acceleration_asd
from gravcomm.rotor import TwoMassRotor, far_field_dipole_amplitude
from gravcomm.torsion import (
    TorsionBalance,
    TorsionSource,
    signal_to_noise,
)


class NumericsTest(unittest.TestCase):
    def test_balanced_quadrupole_has_zero_dipole(self):
        rotor = TwoMassRotor.balanced_quadrupole(600.0, 1.0, mass_ratio=5.0)
        self.assertAlmostEqual(rotor.dipole_moment, 0.0, places=12)
        self.assertAlmostEqual(rotor.quadrupole_moment, 120.0, places=12)

    def test_snr_threshold_round_trip(self):
        noise = 0.1 * MICROGAL
        amplitude = amplitude_for_snr(10.0, noise, 1.0)
        self.assertAlmostEqual(float(snr_power(amplitude, noise, 1.0)), 10.0)

    def test_dipole_far_field_formula(self):
        self.assertTrue(math.isclose(far_field_dipole_amplitude(3.75, 1.0), 2.0 * 6.6743e-11 * 3.75))

    def test_bench_rotor_mechanics(self):
        state = rotor_mechanical_state(DEPLOYMENTS[0])
        self.assertAlmostEqual(state.kinetic_energy, 1.5625e6)
        self.assertAlmostEqual(state.peak_retention_force, 5.208333333333333e6)
        self.assertAlmostEqual(state.required_tie_area, 0.010416666666666666)

    def test_published_noise_references_are_in_ascending_order(self):
        floors = [reference.noise_asd for reference in PUBLISHED_NOISE_REFERENCES]
        self.assertEqual(floors, sorted(floors))

    def test_peterson_nlnm_10hz(self):
        self.assertAlmostEqual(peterson_acceleration_asd(10.0, "NLNM") / MICROGAL, 0.3981071705534972)

    def test_geophone_flat_velocity_band(self):
        # Velocity response approaches the transduction constant well above f0.
        self.assertAlmostEqual(SM24.velocity_response(2000.0), SM24.sensitivity, places=2)

    def test_geophone_acceleration_floor_rises_with_frequency(self):
        # A velocity sensor's acceleration noise floor increases with frequency.
        low = SM24.acceleration_noise_asd(40.0)
        high = SM24.acceleration_noise_asd(240.0)
        self.assertGreater(high, low)

    def test_gradiometer_signal_below_single_and_rejects_common_mode(self):
        rotor = DEPLOYMENTS[3].rotor()  # fixed installation
        grad = Gradiometer(SM24, baseline=0.3, cmrr_db=40.0)
        single = rotor.peak_ac_amplitude(8.0)
        differential = grad.differential_signal_amplitude(rotor, 8.0)
        self.assertLess(differential, single)
        # 40 dB CMRR suppresses a common-mode ambient by 100x in amplitude.
        ambient = 10.0 * MICROGAL
        quiet = grad.effective_noise_asd(40.0, ambient_asd=0.0)
        noisy = grad.effective_noise_asd(40.0, ambient_asd=ambient)
        self.assertLess(noisy, math.hypot(quiet, ambient))

    def test_torsion_thermal_limited_snr_independent_of_resonance(self):
        # With negligible readout noise the link is thermal-limited, and the
        # transfer function cancels: SNR is the same on and off resonance.
        balance = TorsionBalance(
            test_mass=0.02, test_arm=0.05, resonance_hz=1.0,
            quality_factor=50.0, readout_asd=1e-30,
        )
        rotor = TwoMassRotor.balanced_quadrupole(4.0, 0.1)
        on_res = TorsionSource(rotor, distance=0.3, rotation_hz=0.5)   # carrier 1 Hz = f0
        off_res = TorsionSource(rotor, distance=0.3, rotation_hz=0.2)  # carrier 0.4 Hz
        snr_on = signal_to_noise(balance, on_res, 1.0)
        snr_off = signal_to_noise(balance, off_res, 1.0)
        self.assertAlmostEqual(snr_on, snr_off, delta=0.01 * snr_on)

    def test_torsion_contact_geometry_rejected(self):
        rotor = TwoMassRotor.balanced_quadrupole(4.0, 0.2)
        with self.assertRaises(ValueError):
            TorsionSource(rotor, distance=0.15, rotation_hz=0.5)

    def test_spoked_plate_lighter_than_solid(self):
        plate = CarrierPlate(0.25, 0.020, 0.030, 4, 0.025, 0.030, ALUMINUM)
        self.assertLess(plate.fill_factor, 1.0)
        self.assertGreater(plate.polar_inertia, 0.0)

    def test_disk_carrier_adds_inertia_and_preserves_insert_inertia(self):
        bench = DEPLOYMENTS[0]
        plate = CarrierPlate(0.25, 0.020, 0.030, 4, 0.025, 0.030, ALUMINUM)
        state = disk_rotor_mechanical_state(bench, plate)
        # The signal-bearing inertia is unchanged from the masses-only rotor.
        self.assertAlmostEqual(state.insert_inertia, bench.rotor().quadrupole_moment)
        # This mechanics helper adds carrier inertia; it does not compute its field.
        self.assertGreater(state.modulation_torque_factor, 1.0)
        self.assertAlmostEqual(
            state.total_inertia, state.insert_inertia + state.plate_inertia
        )

    def test_disk_raises_whirl_critical_speed(self):
        # A short stiff stub shaft has a higher critical speed than a slender arm.
        # rod arm: 10 mm dia (5 mm radius), 0.20 m, 2 kg tip mass.
        rod = whirl_critical_hz(cantilever_bending_stiffness(200e9, 0.005, 0.20), 2.0)
        # disk on a 12 mm dia (6 mm radius) stub, 0.05 m overhang, 7.6 kg rotor.
        disk = whirl_critical_hz(cantilever_bending_stiffness(200e9, 0.006, 0.05), 7.6)
        self.assertGreater(disk, rod)


if __name__ == "__main__":
    unittest.main()
