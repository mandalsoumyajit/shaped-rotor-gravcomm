import numpy as np
from gravcomm.magnetic_wall import *

def test_free_space_finite_loop():
 for d in (.5,1,2):
  a=.24;H,loss,e=wall_channel(np.array([1.,100.]),0,d,coil_radius=a)
  np.testing.assert_allclose(H,MU0/(2*np.pi*(a*a+d*d)**1.5),rtol=1e-9)
  np.testing.assert_array_equal(e,0)

def test_zero_thickness():
 f=np.array([1.,100.]);a=wall_channel(f,0,1)[0];b=wall_channel(f,1e7,1,thickness=0)[0]
 np.testing.assert_allclose(a,b,rtol=1e-12)
 assert np.all(wall_channel(f,1e7,1,thickness=0)[2]==0)

def test_spectral_limits():
 k=np.array([.1,1,10]);T,R=slab_factors(k,[1,100],0,.1)
 np.testing.assert_allclose(T,1);np.testing.assert_allclose(R,0)
 T,R=slab_factors(k,[1],1e20,.1);assert np.max(abs(T))<1e-10
 np.testing.assert_allclose(R,-1,rtol=1e-5)

def test_passivity_and_quadrature_convergence():
 f=np.geomspace(.01,1000,32)
 for sigma in (4.,1e4,1e7):
  H,L,E=wall_channel(f,sigma,1,order=128)
  H2,L2,E2=wall_channel(f,sigma,1,order=256)
  assert np.all(E>=0)
  np.testing.assert_allclose(H,H2,rtol=.002,atol=1e-25)
  np.testing.assert_allclose(L,L2,rtol=.002)

def test_loading_placement():
 f=np.array([1.,100.]);H,L,_=wall_channel(f,1e6,1,stand_off=.1)
 H2,L2,_=wall_channel(f,1e6,1,stand_off=.8)
 np.testing.assert_allclose(H,H2)
 assert np.all(L>L2)


def test_finite_receiver_reciprocity_and_small_loop_limit():
 f=np.array([1.,100.])
 h=wall_channel(f,1e5,1,coil_radius=.24,receive_radius=.34)[0]
 hr=wall_channel(f,1e5,1,coil_radius=.34,receive_radius=.24)[0]
 np.testing.assert_allclose(h,hr,rtol=1e-12)
 point=wall_channel(f,1e5,1,coil_radius=.24)[0]
 tiny=wall_channel(f,1e5,1,coil_radius=.24,receive_radius=1e-5)[0]
 np.testing.assert_allclose(point,tiny,rtol=1e-8)
