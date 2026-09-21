import unittest
import numpy as np
from numpy.testing import assert_allclose
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.noise import peterson_acceleration_asd
from gravcomm.constants import G
class EnvironmentalNoiseTests(unittest.TestCase):
 def test_units_and_rayleigh_conversion(self):
  c=SeismicBackground().components(np.array([1.,10.]))
  assert_allclose(c['ground_acceleration_asd'],[peterson_acceleration_asd(f,'NHNM') for f in (1.,10.)])
  assert_allclose(c['newtonian_acceleration_asd'],2*np.pi*G*2000*.8/np.sqrt(2)*c['ground_displacement_asd'])
  assert_allclose(c['background_upper_psd'],(c['newtonian_acceleration_asd']+c['support_acceleration_asd'])**2)
 def test_extension_and_domain_are_explicit(self):
  with self.assertRaises(ValueError):SeismicBackground(high_frequency_extension='reject').components(11.)
  with self.assertRaises(ValueError):SeismicBackground().components(.01)
  c=SeismicBackground().components(np.array([10.,20.]))
  assert_allclose(c['ground_acceleration_asd'][0],c['ground_acceleration_asd'][1])
  assert_allclose(c['newtonian_acceleration_asd'][0],4*c['newtonian_acceleration_asd'][1])
  self.assertEqual(c['extrapolated'].tolist(),[False,True])
 def test_isolator_resonance_and_height(self):
  c=SeismicBackground().components(.3)
  expected=(1+.1**2)/(.1**2)
  assert_allclose(abs(c['support_transfer']),expected)
  ground=SeismicBackground().components(1.)['newtonian_acceleration_asd']
  raised=SeismicBackground(sensor_height_m=1).components(1.)['newtonian_acceleration_asd']
  assert_allclose(raised/ground,np.exp(-2*np.pi/200))
 def test_causal_receiver_uses_absolute_background_frequency(self):
  from gravcomm.causal_receiver import CausalReceiver
  background=SeismicBackground();plain=CausalReceiver()
  low=CausalReceiver(carrier_hz=20,background=background)
  high=CausalReceiver(carrier_hz=50.3,background=background)
  f=np.array([-1.,0.,1.]);amplitude=5e-10
  assert_allclose(low.input_psd(f,amplitude)-plain.input_psd(f,amplitude),32*background.acceleration_psd(20+f)/amplitude**2)
  self.assertGreater(low.input_psd([0],amplitude)[0],high.input_psd([0],amplitude)[0])
  model=low.finite_record(128,amplitude)
  assert_allclose(model.factor@model.factor.conj().T,model.covariance,atol=1e-10)
  with self.assertRaises(ValueError):CausalReceiver(carrier_hz=5,background=background)
 def test_vectorized_peterson_matches_scalar_reference(self):
  frequencies=np.geomspace(.1,100,501)
  for model in ('NLNM','NHNM'):
   expected=np.array([peterson_acceleration_asd(float(min(f,10)),model) for f in frequencies])
   got=SeismicBackground(seismic_model=model).components(frequencies)['ground_acceleration_asd']
   assert_allclose(got,expected,rtol=2e-14,atol=0)
if __name__=='__main__':unittest.main()
