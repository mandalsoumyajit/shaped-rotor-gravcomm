"""Ground-motion mitigation scenarios for the revised link study.

The scenarios modify only the support-motion path of the released
SeismicBackground (NHNM ground acceleration, flat above 10 Hz, Rayleigh
Newtonian term). Instrumental noise is still added by CausalReceiver; a
gradiometer adds one further independent instrument through this background.
None of these scenarios is a measured site or isolation specification.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
# Inside the public repository (paper_studies/ground_motion_study) the package sits at the repository root;
# in the authors' working tree it sits in the public_release_v0_5_0 checkout next to the paper folder.
RELEASE = _ROOT if (_ROOT/'paper_studies').is_dir() else _ROOT/'public_release_v0_5_0'
sys.path[:0] = [str(RELEASE/'src'), str(RELEASE/'scripts')]
from gravcomm.constants import G
from gravcomm.environmental_noise import SeismicBackground
from gravcomm.receivers import StructuralOscillator
from gravcomm.finite_rotor import FiniteRotor

ROTOR_NPZ = RELEASE/'fea/results/shaped_waveforms/medium_rounded.npz'
A0 = 5.180259430854157e-10          # s=1 steel rotor at 0.5 m (released constant)
R0 = 0.25                            # normalization radius (m)
TONE_EDGE_HZ = 3*0.14464             # outermost tone offset from carrier
# Released qualified designs: distance -> (carrier, amplitude, scale, mass, inertia, peak torque, peak power)
DESIGNS = {0.5: (24., .8125*A0, 0.9611854, 7.03, 0.259, 1.60, 120.49),
           1.0: (24., .8125*A0, 1.6825, 37.69, 4.262, 26.26, 1980.24),
           2.0: (18., 1.0*A0, 3.0626, 227.32, 85.156, 524.73, 29680.14)}
RELEASED_DESIGNS, RELEASED_TONE_EDGE_HZ = DESIGNS, TONE_EDGE_HZ
# Superseded joint design (hop-limited CF-FSK, K5 binary code, designed labeling, iterative receiver), qualified 2026-11:
# tones k*h/T, k=-3..3, h=3/14, T=224/190 s; torque and power from the exact quintic trajectory
# (drive_hop.py); scale, mass and inertia from sizing.table(0.70, 0.84).
K5_DESIGNS = {0.5: (24., .70*A0, 0.9341619991264621, 6.453575846503782, 0.2245819869821951, 0.611859494654652, 46.65977840075601),
               1.0: (24., .70*A0, 1.6346739293067527, 34.56641404165896, 3.6897229916087366, 10.05241816317332, 766.5871135171993),
               2.0: (18., .84*A0, 2.960317884294902, 205.29664849354262, 71.85477053540292, 195.76380180394054, 11238.957921825311)}
# Joint design (hop-limited CF-FSK, NB-LDPC (43,29) over GF(512), one code symbol per 4-tone word,
# symbol-level iterative receiver), qualified 2026-11-05 at 0.58 / 0.697 A0 (qualify_nb.py):
# tones k*h/T, k=-3..3, h=3/14, T=224/190 s; scale, mass, inertia, torque and power from
# sizing.table(0.58, 0.697, 1.0, 0.5451) (results/sizing_nb.json; torque/power checked with drive_hop.py).
HOP_TONE_EDGE_HZ = 3*(3/14)/(224/190)
HOP_DESIGNS = {0.5: (24., .58*A0, 0.9010437566657981, 5.791238810337648, 0.18749657800956693, 0.5108226310220744, 38.954810661284434),
               1.0: (24., .58*A0, 1.5761135310793883, 30.982996024377382, 3.0745075824072816, 8.376302485171058, 638.7682485770378),
               2.0: (18., .697*A0, 2.854447700523453, 184.04883954571554, 59.89275815922221, 163.1740516379004, 9367.954051734187)}
if __import__('os').environ.get('GC_DESIGN') == 'hop':
    DESIGNS, TONE_EDGE_HZ = HOP_DESIGNS, HOP_TONE_EDGE_HZ


@dataclass(frozen=True)
class MitigatedBackground:
    carrier_hz: float
    isolation: str = 'viscous'          # 'viscous' (released) or 'structural' (low-loss)
    stages: int = 2
    isolation_frequency_hz: float = .3
    damping: float = .05                # viscous damping ratio
    loss_angle: float = .1              # structural loss angle (stage Q = 10)
    feedforward_residual: float = 1.    # amplitude factor left after witness subtraction
    witness_asd: float = 1e-9           # ground witness self-noise (m s^-2 Hz^-1/2)
    common_mode_residual: float = 1.    # gradiometer rejection of shared support motion
    extra_instruments: int = 0          # independent instrument PSDs added (gradiometer: 1)

    def components(self, f):
        f = np.asarray(f, dtype=float)
        c = SeismicBackground().components(f)
        r = f/self.isolation_frequency_hz
        if self.isolation == 'viscous':
            stage = (1 + 2j*self.damping*r)/(1 - r**2 + 2j*self.damping*r)
        elif self.isolation == 'structural':
            stage = (1 + 1j*self.loss_angle)/(1 - r**2 + 1j*self.loss_angle)
        else:
            raise ValueError('unknown isolation model')
        transfer = abs(stage**self.stages)
        support = transfer*c['ground_acceleration_asd']
        if self.feedforward_residual < 1:
            # Imperfect transfer-function model plus witness self-noise projected through H.
            support = np.hypot(self.feedforward_residual*support, transfer*self.witness_asd)
        support = self.common_mode_residual*support
        newtonian = self.common_mode_residual*c['newtonian_acceleration_asd']
        instrument = self.extra_instruments*StructuralOscillator().acceleration_noise_asd(
            50.3 + (f - self.carrier_hz))**2
        return dict(support=support, newtonian=newtonian, instrument_extra=instrument,
                    total_background_psd=instrument + (support + newtonian)**2)

    def acceleration_psd(self, f):
        return self.components(f)['total_background_psd']


SCENARIOS = {
    'baseline':     dict(label='Released two-stage isolation', kw={}),
    'feedforward':  dict(label='Witness feed-forward, 20 dB', kw=dict(feedforward_residual=.1)),
    'structural':   dict(label='Low-loss two-stage isolation', kw=dict(isolation='structural')),
    'three_stage':  dict(label='Three-stage isolation', kw=dict(stages=3)),
    'gradiometer':  dict(label='Gradiometer, 40 dB rejection', kw=dict(common_mode_residual=.01, extra_instruments=1)),
    'combined':     dict(label='Low-loss isolation + feed-forward', kw=dict(isolation='structural', feedforward_residual=.1)),
}


def background(name, carrier_hz):
    return MitigatedBackground(carrier_hz=carrier_hz, **SCENARIOS[name]['kw'])


def total_input_psd(bg, offset_hz, instrument_factor=1.):
    """One-sided physical acceleration PSD at offset from the carrier.

    instrument_factor scales every instrument, including a gradiometer's second sensor.
    """
    from dataclasses import replace
    inst = StructuralOscillator().acceleration_noise_asd(50.3 + offset_hz)**2
    environment = replace(bg, extra_instruments=0).acceleration_psd(bg.carrier_hz + offset_hz)
    return (1 + bg.extra_instruments)*instrument_factor**2*inst + environment


_HOP_PSD = None


def signal_weight(offset_hz):
    """Occupied-band weight.

    Released design: CF-FSK tone span plus 0.15 Hz Gaussian transition shoulders (as published).
    Joint (hop) design: the Welch PSD of the actual baseband waveform with random data
    (results/hop_psd.npz from check_weight.py), i.e. the exact matched-filter SNR weight.
    """
    global _HOP_PSD
    if DESIGNS is HOP_DESIGNS:
        if _HOP_PSD is None:
            z = np.load(Path(__file__).resolve().parent/'results/hop_psd.npz'); _HOP_PSD = (z['f'], z['p'])
        return np.interp(offset_hz, *_HOP_PSD)
    return np.exp(-0.5*(np.maximum(np.abs(offset_hz) - TONE_EDGE_HZ, 0)/0.15)**2)


def equivalent_amplitude_factor(bg, bg_ref, instrument_factor=1.):
    """Amplitude ratio giving equal whitened (matched-filter) SNR over the occupied band."""
    nu = np.linspace(-1, 1, 4001); w = signal_weight(nu)
    ref = np.trapezoid(w/total_input_psd(bg_ref, nu), nu)
    new = np.trapezoid(w/total_input_psd(bg, nu, instrument_factor), nu)
    return float(np.sqrt(ref/new))


class ScaledRotor:
    def __init__(self, path=ROTOR_NPZ):
        self.rotor = FiniteRotor.from_npz(path)
        self.p = np.asarray(self.rotor.positions_m); self.m = np.asarray(self.rotor.masses_kg)

    def coefficient(self, s, point, component=0, samples=128):
        """Complex 2nd-harmonic coefficient of acceleration at an arbitrary receiver point."""
        point = np.asarray(point, dtype=float)
        p = self.p*s; m = self.m*s**3
        phase = np.arange(samples)*2*np.pi/samples; values = np.empty(samples)
        for k, ph in enumerate(phase):
            co, si = np.cos(ph), np.sin(ph)
            rel = np.column_stack((co*p[:, 0] - si*p[:, 1], si*p[:, 0] + co*p[:, 1], p[:, 2])) - point
            values[k] = G*np.sum(m*rel[:, component]/np.sum(rel**2, axis=1)**1.5)
        values -= values.mean()
        return complex(2*np.mean(values*np.exp(-2j*phase)))

    def differential(self, s, d, baseline, layout):
        if layout == 'single':
            return self.coefficient(s, (d, 0, 0))
        if layout == 'radial':
            return self.coefficient(s, (d, 0, 0)) - self.coefficient(s, (d + baseline, 0, 0))
        if layout == 'transverse':
            return self.coefficient(s, (d, baseline/2, 0)) - self.coefficient(s, (d, -baseline/2, 0))
        raise ValueError(layout)

    def invert(self, amplitude, d, baseline=0., layout='single', clearance=.0127):
        """Smallest uniform scale whose (differential) A2 reaches the amplitude."""
        from scipy.optimize import brentq
        f = lambda s: abs(self.differential(s, d, baseline, layout)) - amplitude
        hi = (d - clearance)/R0*0.999
        lo = 0.05
        if f(hi) < 0:
            return None
        return brentq(f, lo, hi, xtol=1e-6)


def sized(d, s, carrier_hz):
    """Mass, inertia, peak torque and power from the released design by geometric scaling."""
    _, _, s0, m0, i0, t0, p0 = DESIGNS[d]
    k = s/s0
    return dict(scale=s, mass_kg=m0*k**3, inertia_kg_m2=i0*k**5, peak_torque_nm=t0*k**5,
                peak_power_w=p0*k**5*carrier_hz/DESIGNS[d][0], diameter_m=2*R0*s,
                rim_ok=bool(s*(carrier_hz + TONE_EDGE_HZ) <= 60))
