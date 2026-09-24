"""Boundary/energy checks for the dipole-cavity magnetic link model."""
import numpy as np
from scipy.integrate import quad
from gravcomm.magnetic_link import *

def test_vacuum_limit():
 f=np.array([1.,1000.]);g=np.exp(log_field_per_moment(f,0,2,.25))
 np.testing.assert_allclose(g[0],MU0/(2*np.pi*2**3),rtol=1e-12)
 np.testing.assert_allclose(g[1],MU0/(4*np.pi*2**3),rtol=1e-12)

def test_cavity_boundary_matching():
 for t in [.001,.1,1.,10.,100.]:
  x=(1+1j)*t;den=x*x+3*x+3;Dscaled=3/den;Cscaled=-x*x/den
  np.testing.assert_allclose(Dscaled*(1+x),1+Cscaled,atol=1e-13)
  np.testing.assert_allclose(Dscaled*(-2*(1+x)-x*x),-2+Cscaled,atol=1e-13)

def test_reaction_loss_matches_volume_joule_integral():
 for sigma in [.1,100.,1e6]:
  f=10.;a=.25;w=2*np.pi*f;g=(1+1j)*np.sqrt(w*MU0*sigma/2);x=g*a
  # Integral sigma |E_rms|^2 dV, with A_phi=mu m/(4pi)*T*(1+gr)e^-gr/r² sin(theta).
  def integrand(r):return abs(3*(1+g*r)*np.exp(-g*(r-a))/(x*x+3*x+3))**2/r**2
  integral=quad(integrand,a,np.inf,epsabs=1e-10)[0]
  direct=sigma*w*w*(MU0/(4*np.pi))**2*(8*np.pi/3)*integral
  reaction=moment_loss(f,sigma,a,.125,1)[2]
  np.testing.assert_allclose(direct,reaction,rtol=1e-7)

def test_weak_conductor_loss():
 f=1.;a=.25;sigma=1e-8
 expected=MU0**2*sigma*(2*np.pi*f)**2/(6*np.pi*a)
 np.testing.assert_allclose(moment_loss(f,sigma,a,.125,1)[2],expected,rtol=1e-5)

def test_flat_channel_capacity():
 g=np.full(100,np.log(10.));df=np.full(100,.02)
 np.testing.assert_allclose(np.exp(minimum_capacity_logpower(g,df)),2/10*(2**.5-1),rtol=1e-12)

def test_power_monotonic_and_orientation():
 f=np.geomspace(.01,1000,100);g=log_gain(f,4,2)
 assert np.all(np.isfinite(g));assert np.all(moment_loss(f,4,.25,.125,1)[2]>=0)
 assert minimum_capacity_logpower(g.max(axis=0),np.ones(100),2)>minimum_capacity_logpower(g.max(axis=0),np.ones(100),1)
