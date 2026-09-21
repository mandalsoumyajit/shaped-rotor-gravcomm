import itertools
import unittest
import numpy as np
from gravcomm.convolutional_candidates import encode, decode, metadata, PATTERNS


class TailBitingTests(unittest.TestCase):
    def test_noiseless_nonzero_and_arbitrary_lengths(self):
        rng = np.random.default_rng(279)
        for rate in PATTERNS:
            for k in (1, 2, 5, 6, 7, 63, 64, 65, 128, 256):
                for bits in (np.ones(k, dtype=np.uint8), rng.integers(0, 2, k, dtype=np.uint8)):
                    coded = encode(bits, rate)
                    got = decode(20*(1-2*coded.astype(float)), rate, k)
                    np.testing.assert_array_equal(got, bits)
                    self.assertEqual(len(coded), metadata(k, rate)['n'])

    def test_independent_cyclic_polynomial_and_puncturing(self):
        bits = np.array([1,0,0,1,1,0,1,0,1,1,1], dtype=np.uint8)
        # Bit j of polynomial multiplies input delayed j samples.
        mother = np.array([[sum(int(bits[(t-j)%len(bits)]) for j in range(7)
                                if (g >> j)&1)%2 for g in (0o171, 0o133)]
                           for t in range(len(bits))])
        self.assertEqual(mother.reshape(-1).tolist(), encode(bits).tolist())
        self.assertEqual(mother.reshape(-1)[np.resize([1,1,0,1],22).astype(bool)].tolist(),
                         encode(bits,'2/3').tolist())
        self.assertEqual(mother.reshape(-1)[np.resize([1,1,0,1,1,0],22).astype(bool)].tolist(),
                         encode(bits,'3/4').tolist())

    def test_exact_against_exhaustive_codebook(self):
        rng = np.random.default_rng(808)
        for rate in PATTERNS:
            for k in (3, 5, 8):
                words = np.array(list(itertools.product((0,1), repeat=k)), dtype=np.uint8)
                codebook = np.array([encode(w, rate) for w in words])
                for _ in range(8):
                    llr = rng.normal(size=codebook.shape[1])
                    scores = codebook @ llr
                    result = decode(llr, rate, k)
                    self.assertAlmostEqual(float(encode(result, rate) @ llr), float(scores.min()), places=11)

    def test_validation_and_partial_period_rate(self):
        self.assertEqual(metadata(64, '3/4')['n'], 86)
        self.assertEqual(metadata(128, '3/4')['n'], 171)
        self.assertEqual(metadata(256, '3/4')['n'], 342)
        for bad in ([2], [0.2], [[0,1]], []):
            with self.assertRaises(ValueError): encode(bad)
        with self.assertRaises(ValueError): decode([np.nan]*128, information_bits=64)
        with self.assertRaises(ValueError): decode([1]*127, information_bits=64)
        with self.assertRaises(ValueError): encode([0,1], '5/6')


if __name__ == '__main__': unittest.main()
