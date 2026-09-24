"""Magnetoquasistatic conducting-medium link with a centered dipole cavity.
All spectra are one-sided; source moment is RMS. SI units throughout.
"""
import numpy as np
from scipy.special import erfcinv
MU0=4*np.pi*1e-7
RHO_CU=1.724e-8
DENSITY_CU=8960.

def moment_loss(f,sigma,cavity_radius,coil_radius,copper_mass):
 f=np.asarray(f);w=2*np.pi*f;t=cavity_radius*np.sqrt(w*MU0*sigma/2)
 # Imaginary part of x^2/(x^2+3x+3), x=(1+i)t, evaluated without cancellation.
 imaginary=2*t*t*(3+3*t)/((3+3*t)**2+(3*t+2*t*t)**2)
 medium=MU0*w/(2*np.pi*cavity_radius**3)*imaginary
 copper=4*RHO_CU/(coil_radius**2*(copper_mass/DENSITY_CU))
 return copper+medium,copper+np.zeros_like(f),medium

def log_field_per_moment(f,sigma,distance,cavity_radius):
 f=np.asarray(f);g=(1+1j)*np.sqrt(np.pi*f*MU0*sigma);x=g*cavity_radius;u=g*distance
 common=np.log(3*MU0/(4*np.pi))-3*np.log(distance)-np.log(abs(x*x+3*x+3))-g.real*(distance-cavity_radius)
 return np.stack([common+np.log(abs(2*(1+u))),common+np.log(abs(1+u+u*u))])

def noise_asd(f,ambient=30e-15,sensor='pcb',ambient_knee=0.):
 f=np.asarray(f)
 if sensor=='pcb':inst=1.626e-9*np.sqrt((1+20/f)/1.02)/(2*np.pi*f*36.388)
 elif sensor=='ideal_flat':inst=np.zeros_like(f)
 else:raise ValueError(sensor)
 env=ambient*np.sqrt(1+(ambient_knee/f)**2)
 return np.sqrt(inst*inst+env*env)

def log_gain(f,sigma,distance,cavity_radius=.25,coil_fraction=.5,copper_mass=1.,ambient=30e-15,sensor='pcb',ambient_knee=0.):
 if distance<=cavity_radius:raise ValueError('Receiver must be outside the source cavity')
 loss=moment_loss(f,sigma,cavity_radius,cavity_radius*coil_fraction,copper_mass)[0]
 return 2*log_field_per_moment(f,sigma,distance,cavity_radius)-np.log(loss)-2*np.log(noise_asd(f,ambient,sensor,ambient_knee))

def minimum_capacity_logpower(logg,widths,rate=1.):
 """Minimum log average watts for asymptotic Gaussian capacity rate.
 Permits frequency-dependent best orientation; optimistic vector-source bound.
 """
 order=np.argsort(logg)[::-1];g=logg[order];df=widths[order]
 D=np.cumsum(df);S=np.cumsum(df*g);logwater=(rate*np.log(2)-S)/D
 valid=np.flatnonzero(logwater>=-g);k=int(valid[-1]);L=logwater[k]
 factors=-np.expm1(np.minimum(0,-g[:k+1]-L))
 return float(L+np.log(np.sum(df[:k+1]*factors)))

def bpsk_logpower(sigma,distance,rate=1.,ber=.001,**kwargs):
 """Ideal whitened/equalized raised-cosine BPSK benchmark, rolloff .25.
 Unit payload/channel bit rate; framing/acquisition overhead omitted.
 """
 alpha=.25;half=(1+alpha)*rate/2
 carriers=np.geomspace(max(1.,half+.01),1000-half,401)
 # Midpoint pulse-spectrum quadrature, normalized to unit integral.
 z=(np.arange(80)+.5)/80*2*half-half
 pulse=np.ones_like(z)/rate;edge=abs(z)>(1-alpha)*rate/2
 pulse[edge]=.5/rate*(1+np.cos(np.pi/(alpha*rate)*(abs(z[edge])-(1-alpha)*rate/2)))
 weights=pulse/pulse.sum();frequencies=carriers[:,None]+z
 gains=log_gain(frequencies,sigma,distance,**kwargs)
 from scipy.special import logsumexp
 logP=np.log(erfcinv(2*ber)**2*rate)+logsumexp(np.log(weights)-gains,axis=-1)
 axis,index=np.unravel_index(np.argmin(logP),logP.shape)
 return float(logP[axis,index]),float(carriers[index]),int(axis)
