"""Finite circular coil, air/conducting slab/air, magnetoquasistatic Hankel model."""
import numpy as np
from functools import lru_cache
from scipy.special import roots_laguerre,j1
from numpy.polynomial.legendre import leggauss
from .magnetic_link import MU0,RHO_CU,DENSITY_CU,noise_asd
@lru_cache(None)
def nodes(n):return roots_laguerre(n)

def slab_factors(k,f,sigma,thickness):
 omega=2*np.pi*np.asarray(f)[...,None];k=np.asarray(k)
 q=np.sqrt(k*k+1j*omega*MU0*sigma)
 reflection=-1j*omega*MU0*sigma/(k+q)**2
 e=np.exp(-2*q*thickness);den=1-reflection**2*e
 transmission_scaled=(1-reflection**2)*np.exp(-(q-k)*thickness)/den
 reflected=reflection*(-np.expm1(-2*q*thickness))/den
 return transmission_scaled,reflected

def wall_channel(f,sigma,distance,thickness=.1,coil_radius=.25,copper_mass=1.,stand_off=None,order=128,loss_refinement=1,receive_radius=None):
 """Returns B_rms / m_rms and dissipative watts / m_rms².
 m=N I_rms pi a². Sensor on coil axis beyond wall; normal coil orientation.
 """
 if distance<=thickness:raise ValueError('Distance must exceed wall thickness')
 if stand_off is None:stand_off=(distance-thickness)/2
 if not 0<stand_off<distance-thickness:raise ValueError('Both air stand-offs must be positive')
 f=np.asarray(f);u,w=nodes(order);k=u/distance
 T,_=slab_factors(k,f,sigma,thickness)
 if receive_radius is None:
  H=MU0/(2*np.pi*coil_radius*distance**2)*np.sum(w*u*j1(u*coil_radius/distance)*T,axis=-1)
 else:
  H=MU0/(np.pi*coil_radius*receive_radius*distance)*np.sum(w*j1(u*coil_radius/distance)*j1(u*receive_radius/distance)*T,axis=-1)
 # Composite Gauss quadrature resolves Bessel oscillations even very close to a wall.
 xmax=max(50.,15*coil_radius/stand_off);segments=int(np.ceil(xmax/(np.pi/2)))
 edges=np.linspace(0,xmax,segments+1);gn,gw=leggauss(6*loss_refinement)
 x=((edges[1:]+edges[:-1])[:,None]/2+(edges[1:]-edges[:-1])[:,None]/2*gn).ravel()
 weights=((edges[1:]-edges[:-1])[:,None]/2*gw).ravel()
 k=x/coil_radius;factor=weights*j1(x)**2*np.exp(-2*x*stand_off/coil_radius)/coil_radius
 flat=f.ravel();eddy=np.empty_like(flat)
 for begin in range(0,len(flat),32):
  stop=min(begin+32,len(flat));_,R=slab_factors(k,flat[begin:stop],sigma,thickness)
  eddy[begin:stop]=-(2*np.pi*flat[begin:stop])*MU0/(np.pi*coil_radius**2)*np.sum(factor*R.imag,axis=-1)
 eddy=eddy.reshape(f.shape)
 copper=4*RHO_CU/(coil_radius**2*(copper_mass/DENSITY_CU))
 return H,copper+np.maximum(eddy,0),eddy

def wall_log_gain(f,sigma,distance,ambient=30e-15,sensor='pcb',ambient_knee=0.,**kw):
 H,loss,_=wall_channel(f,sigma,distance,**kw)
 return 2*np.log(np.maximum(abs(H),1e-300))-np.log(loss)-2*np.log(noise_asd(f,ambient,sensor,ambient_knee))
