"""Independent checks for the manuscript desktop operating point."""
import math
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from desktop_cf_fsk import encode,decode,all_transitions,trajectory,monte_carlo,FC


class DesktopCFFSKTest(unittest.TestCase):
    def test_actual_byte_mapping_and_invalid_word(self):
        np.testing.assert_array_equal(decode(encode(np.arange(256))),np.arange(256))
        with self.assertRaises(ValueError):decode([3,3,3])

    def test_all_hops_and_polynomial_integrals(self):
        symbols,previous=all_transitions()
        self.assertEqual(len(set(zip(np.r_[previous,symbols[:-1]],symbols))),49)
        ts=8.;t,f,df,ddf,phase=trajectory(np.array([3]),-3,ts,fs=1000)
        integrate=getattr(np,'trapezoid',np.trapz)
        self.assertAlmostEqual(integrate(df,t),6/6.4,places=10)
        self.assertAlmostEqual(integrate(f-FC,t),phase[-1]/(2*np.pi),places=9)
        self.assertAlmostEqual(integrate(df**2,t),(6/6.4)**2/1.6*10/7,places=10)
        self.assertAlmostEqual(max(df),1.875*(6/6.4)/1.6,places=12)
        self.assertEqual(df[0],0);self.assertEqual(df[-1],0)
        self.assertEqual(ddf[0],0);self.assertEqual(ddf[-1],0)

    def test_noncoherent_noise_normalization_against_exact_white_result(self):
        # Orthogonal equal-energy tones: exact noncoherent M-FSK symbol error.
        gamma=10.;M=7
        exact=sum((-1)**(k+1)*math.comb(M-1,k)/(k+1)*math.exp(-k*gamma/(k+1)) for k in range(1,M))
        result=monte_carlo(np.sqrt(2*gamma)*np.eye(M,dtype=complex),2*np.eye(M),np.arange(M)-3,30000,619)
        # E_s/N0=gamma: each noise quadrature has variance one.
        standard_error=np.sqrt(exact*(1-exact)/result['total_trials'])
        self.assertLess(abs(result['mean_ser']-exact),5*standard_error)


if __name__=='__main__':unittest.main()
