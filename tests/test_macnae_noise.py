import numpy as np
import pytest
from gravcomm.macnae_noise import macnae_asd,_anchors,DEFAULT_DATA

def test_units_anchors_and_power_law():
 for f,a in _anchors(str(DEFAULT_DATA)):
  np.testing.assert_allclose(macnae_asd(f),a,rtol=1e-12)
  np.testing.assert_allclose(macnae_asd(np.sqrt(f[:-1]*f[1:])),np.sqrt(a[:-1]*a[1:]),rtol=1e-12)
 assert 1e-10<float(macnae_asd(.01))<1e-9

def test_gap_and_no_extrapolation():
 for f in (.001,1000.,1e6):
  with pytest.raises(ValueError):macnae_asd(f)
 assert np.isnan(macnae_asd(1000.,unsupported='nan'))

def test_asd_scale_and_invalid_input():
 np.testing.assert_allclose(macnae_asd([18,24],scale=2)**2,4*macnae_asd([18,24])**2)
 for f in (0,-1,np.nan):
  with pytest.raises(ValueError):macnae_asd(f)
