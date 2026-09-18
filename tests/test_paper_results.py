import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location('paper_results', Path(__file__).resolve().parents[1] / 'scripts/paper_results.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class GaussianBenchmarkTests(unittest.TestCase):
    def test_white_noise_matches_shannon(self):
        for power in [0., 1e-22, 1e-19]:
            actual = mod.gaussian_waterfill(np.full(1000, 3e-21), 1.5625, power)
            expected = 1.5625 * np.log2(1 + power / (3e-21 * 1.5625))
            self.assertAlmostEqual(actual, expected, places=12)

    def test_high_noise_half_band_is_unused(self):
        # Unit bandwidth: allocate all power 0.5 to the half with PSD=1.
        self.assertAlmostEqual(mod.gaussian_waterfill([1., 100.], 1., .5), .5)

    def test_invalid_inputs(self):
        for psd, bw, power in [([], 1, 1), ([0], 1, 1), ([1], 0, 1), ([1], 1, -1)]:
            with self.assertRaises(ValueError):
                mod.gaussian_waterfill(psd, bw, power)
