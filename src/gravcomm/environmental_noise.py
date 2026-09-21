"""Explicit seismic-background scenarios; not a measured deployment spectrum.

Peterson NLNM/NHNM are ground acceleration spectra. Rayleigh conversion follows
Harms (2015), single horizontal test mass in an isotropic field. Support motion
uses an assumed two-stage base-excited isolation transfer function. A coherent
upper envelope accounts conservatively for correlation of the two seismic paths.
"""
from dataclasses import dataclass
import numpy as np
from .constants import G
from .noise import _PETERSON_NLNM, _PETERSON_NHNM

@dataclass(frozen=True)
class SeismicBackground:
    seismic_model: str = 'NHNM'
    density_kg_m3: float = 2000.
    rayleigh_gamma: float = .8
    sensor_height_m: float = 0.
    rayleigh_speed_m_s: float = 200.
    isolation_frequency_hz: float = .3
    isolation_damping: float = .05
    isolation_stages: int = 2
    high_frequency_extension: str = 'flat_acceleration_asd'

    def components(self, frequency_hz):
        f=np.asarray(frequency_hz,dtype=float)
        if np.any(~np.isfinite(f)) or np.any(f<.1) or np.any(f>100):
            raise ValueError('Scenario domain is 0.1 through 100 Hz')
        if self.seismic_model not in ('NLNM','NHNM'):raise ValueError('Unknown seismic model')
        if self.high_frequency_extension not in ('flat_acceleration_asd','reject'):
            raise ValueError('Unknown extension policy')
        if np.any(f>10) and self.high_frequency_extension=='reject':
            raise ValueError('Peterson does not cover frequencies above 10 Hz')
        if min(self.density_kg_m3,self.rayleigh_gamma,self.rayleigh_speed_m_s,self.isolation_frequency_hz,self.isolation_damping)<=0 or self.sensor_height_m<0:
            raise ValueError('Invalid physical scenario')
        if self.isolation_stages<1 or int(self.isolation_stages)!=self.isolation_stages:
            raise ValueError('Positive integer stage count required')
        table={'NLNM':_PETERSON_NLNM,'NHNM':_PETERSON_NHNM}[self.seismic_model]
        period=1/np.minimum(f,10)
        coefficients=table[np.searchsorted(table[:,0],period,side='right')-1]
        ground=np.sqrt(10**((coefficients[...,1]+coefficients[...,2]*np.log10(period))/10))
        displacement=ground/(2*np.pi*f)**2
        gravity=2*np.pi*G*self.density_kg_m3*self.rayleigh_gamma/np.sqrt(2)*np.exp(-2*np.pi*f*self.sensor_height_m/self.rayleigh_speed_m_s)*displacement
        ratio=f/self.isolation_frequency_hz
        stage=(1+2j*self.isolation_damping*ratio)/(1-ratio**2+2j*self.isolation_damping*ratio)
        transfer=stage**self.isolation_stages
        support=abs(transfer)*ground
        # Same envelope is assumed for horizontal and vertical ground motion.
        # Unknown relative phase/correlation: |a+b| <= |a|+|b|.
        upper=gravity+support
        return dict(ground_acceleration_asd=ground,ground_displacement_asd=displacement,
                    newtonian_acceleration_asd=gravity,support_acceleration_asd=support,
                    support_transfer=transfer,background_upper_psd=upper**2,
                    extrapolated=f>10)

    def acceleration_psd(self,frequency_hz):
        return self.components(frequency_hz)['background_upper_psd']
