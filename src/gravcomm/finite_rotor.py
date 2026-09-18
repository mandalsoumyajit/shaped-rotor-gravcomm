"""Finite-volume Newtonian rotor using positive spatial quadrature weights.

Unlike the point-mass approximation, this includes the full carrier and inserts.
Rotation is about z; the receiver is at (distance,0,0). Thus x is radial,
y is in-plane transverse and z is axial. Units are SI. Geometry is undeformed.
"""
from dataclasses import dataclass
import math
import numpy as np
from .constants import G


@dataclass(frozen=True)
class FiniteRotor:
    positions_m: np.ndarray
    masses_kg: np.ndarray
    bounding_radius_m: float

    def __post_init__(self):
        p=np.array(self.positions_m,dtype=float,copy=True)
        m=np.array(self.masses_kg,dtype=float,copy=True)
        if p.ndim!=2 or p.shape[1]!=3 or m.shape!=(len(p),) or not len(p):
            raise ValueError('positions must be N by 3, masses length N')
        if not np.all(np.isfinite(p)) or not np.all(np.isfinite(m)) or np.any(m<=0):
            raise ValueError('positions must be finite; quadrature masses finite and positive')
        radius=self.bounding_radius_m
        if not math.isfinite(radius) or radius<=0 or np.max(np.hypot(p[:,0],p[:,1]))>radius*(1+1e-10):
            raise ValueError('bounding radius must enclose every quadrature point')
        p.setflags(write=False);m.setflags(write=False)
        object.__setattr__(self,'positions_m',p);object.__setattr__(self,'masses_kg',m)

    @classmethod
    def from_npz(cls,path):
        with np.load(path,allow_pickle=False) as f:
            return cls(f['positions_m'],f['masses_kg'],float(f['bounding_radius_m']))

    @property
    def total_mass(self):return float(np.sum(self.masses_kg))

    @property
    def max_arm(self):return self.bounding_radius_m

    @property
    def polar_inertia(self):
        return float(self.masses_kg@(self.positions_m[:,0]**2+self.positions_m[:,1]**2))

    @property
    def quadrupole_anisotropy(self):
        """Complex integral m*(x+i*y)^2, not total polar inertia."""
        return complex(self.masses_kg@(self.positions_m[:,0]+1j*self.positions_m[:,1])**2)

    def acceleration_component(self,distance,phase,component='x'):
        if not math.isfinite(distance) or distance<=self.bounding_radius_m:
            raise ValueError('receiver must lie outside the source bounding cylinder')
        if component not in ('x','y','z'):raise ValueError('component must be x, y or z')
        phi=np.asarray(phase,dtype=float)
        if not np.all(np.isfinite(phi)):raise ValueError('phase must be finite')
        result=np.empty(phi.size);flat=phi.ravel();p=self.positions_m
        batch=max(1,min(32,2_000_000//len(self.masses_kg)))
        for start in range(0,len(flat),batch):
            angle=flat[start:start+batch,None];co=np.cos(angle);si=np.sin(angle)
            dx=co*p[:,0]-si*p[:,1]-distance
            y=si*p[:,0]+co*p[:,1]
            numerator={'x':dx,'y':y,'z':p[:,2]}[component]
            result[start:start+batch]=G*np.sum(self.masses_kg*numerator/(dx*dx+y*y+p[:,2]**2)**1.5,axis=1)
        return result.reshape(phi.shape)

    def ac_signal(self,distance,samples=256,component='x'):
        if not isinstance(samples,int) or samples<8:raise ValueError('samples must be an integer >=8')
        phase=np.arange(samples)*2*np.pi/samples
        values=self.acceleration_component(distance,phase,component)
        return phase,values-float(np.mean(values))

    def harmonic_coefficient(self,distance,harmonic=2,samples=256,component='x'):
        if not isinstance(harmonic,int) or harmonic<1 or 2*harmonic>=samples:
            raise ValueError('harmonic must be positive and below the sampling Nyquist limit')
        phase,ac=self.ac_signal(distance,samples,component)
        return complex(2*np.mean(ac*np.exp(-1j*harmonic*phase)))

    def harmonic_amplitude(self,distance,harmonic=2,samples=256,component='x'):
        return abs(self.harmonic_coefficient(distance,harmonic,samples,component))

    def rms_ac_acceleration(self,distance,samples=256,component='x'):
        return float(np.sqrt(np.mean(self.ac_signal(distance,samples,component)[1]**2)))
