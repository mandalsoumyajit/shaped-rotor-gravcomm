"""Paper CF-FSK coefficients plus a separate, non-paper binary illustration.

Binary symbols; tone separation kappa*R at the gravitational second harmonic;
raised-cosine frequency transitions lasting beta/R. Defaults kappa=1,beta=.5.
Torque is rotor-side and excludes drive inertia, losses and elastic dynamics.
"""
import csv
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def transition(inertia,rate,kappa=1.,beta=.5,high_hz=60.):
    separation=kappa*rate
    if not (inertia>0 and rate>0 and kappa>0 and 0<beta<=1 and high_hz>separation):
        raise ValueError('Invalid transition or nonpositive lower tone')
    duration=beta/rate
    peak=inertia*np.pi**2*separation/(2*duration)
    # Integrate an independent sampled frequency trajectory and verify both
    # angular impulse and work against their exact endpoint differences.
    t=np.linspace(0,duration,10001);theta=np.pi*t/duration
    f=high_hz-separation+.5*separation*(1-np.cos(theta))
    omega=np.pi*f
    torque=peak*np.sin(theta)
    energy=.5*inertia*(omega[-1]**2-omega[0]**2)
    integrate=getattr(np,'trapezoid',np.trapz)
    np.testing.assert_allclose(integrate(torque,t),inertia*(omega[-1]-omega[0]),rtol=1e-7)
    np.testing.assert_allclose(integrate(torque*omega,t),energy,rtol=1e-7)
    np.testing.assert_allclose(integrate(torque**2,t)/duration,peak**2/2,rtol=1e-7)
    return dict(raw_binary_rate_bit_s=rate,tone_separation_hz=separation,
        transition_s=duration,low_rpm=30*(high_hz-separation),high_rpm=30*high_hz,
        peak_torque_nm=peak,rms_torque_iid_nm=peak*np.sqrt(beta/4),
        rms_torque_alternating_nm=peak*np.sqrt(beta/2),
        peak_mechanical_power_w=float(np.max(torque*omega)),energy_per_transition_j=energy,
        nominal_solid_12mm_shaft_shear_mpa=16*peak/(np.pi*.012**3)/1e6)


def main():
    out=ROOT/'results'/'operating_envelope';out.mkdir(exist_ok=True)
    rows=[];specs={}
    for name in ('baseline_rounded','medium_rounded','wide_rounded'):
        data=json.loads((ROOT/'results'/'design_sweep'/name/'h0.013_rpm1800'/'metadata.json').read_text())
        s=data['summary'];I=s['aluminum_polar_inertia_kg_m2']+s['steel_polar_inertia_kg_m2']
        specs[name]=dict(inertia_kg_m2=I,torque_coefficient_nm_per_rate_squared=I*np.pi**2,
            energy_at_1800rpm_j=.5*I*(60*np.pi)**2,
            force_per_insert_at_1800rpm_n=s['steel_mass_kg']/2*.2*(60*np.pi)**2)
        for rate in (.1,.5,1.,2.,5.,10.):
            rows.append(dict(case=name,**transition(I,rate)))
    with (out/'torque_vs_rate.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'assumptions.json').write_text(json.dumps(dict(
        status='Conditional actuator sizing; no established safe speed or usable link rate',
        modulation='Binary, uncoded; second-harmonic tone spacing R Hz; raised-cosine frequency transition',
        transition_fraction=.5,upper_signal_frequency_hz=60.,
        iid_transition_probability=.5,shaft_diameter_m=.012,
        exclusions=['motor and shaft inertia','drag','servo tracking error','elastic dynamics',
                    'shaft bending and notches','hub torque transfer','insert retention','receiver BER'],
        rotor_specs=specs),indent=2))
    table=[]
    for r in rows:
        if r['case']!='wide_rounded':continue
        table.append(f"| {r['raw_binary_rate_bit_s']:g} | {r['low_rpm']:g}-1800 | {r['peak_torque_nm']:.3f} | "
            f"{r['rms_torque_iid_nm']:.3f} | {r['peak_mechanical_power_w']/1000:.4f} | "
            f"{r['nominal_solid_12mm_shaft_shear_mpa']:.2f} |")
    paper_rows=[]
    for name,s in specs.items():
        paper_rows.append(f"| {name} | {s['inertia_kg_m2']:.6f} | {np.pi*s['inertia_kg_m2']:.6f} |")
    correction='''# Paper CF-FSK torque sizing: correction to the binary illustration

A subsequent study fixes a numerical seven-tone desktop operating point:
see `../results/desktop_cf_fsk/README.md` and `../sections/desktop_cf_fsk.tex`.
It selects the 30 mm spoke rounded rotor, 8 s symbols, 0.15625 Hz adjacent
spacing, 1.6 s quintic transitions, 0.333 mapped bit/s at 0.5 m, and 1.091 N m
peak inertial torque. That conditional result replaces the unspecified
operating point discussed below; the binary appendix remains illustrative only.

The manuscript defines M-ary CP/CF-FSK in grav_comm.tex, lines 214-250.
The previously reported binary raised-cosine table was an additional assumed
example, not a calculation for the paper's selected modulation framework.
Its 3.032*R^2 coefficient and rate-specific shaft/power conclusions must not
be presented as requirements of the paper's CF-FSK scheme.

The paper's relation for a received-signal frequency change is

    characteristic torque = pi*I*Delta_f_g/T_r.

Using the full shaped-rotor CAD inertia gives:

| Geometry | I, kg m^2 | pi*I, Nm per (Hz/s) |
|---|---:|---:|
'''+ '\n'.join(paper_rows)+'''

For equally spaced M-ary received tones with adjacent spacing delta_f_g,
a hop of j tone intervals has Delta_f_g=j*delta_f_g. The largest hop is M-1
intervals, not one interval. For the discussed seven-tone case the largest
hop is six intervals. If a normalized transition shape s(t/T_r) is specified,
peak torque is pi*I*abs(j)*delta_f_g*max(abs(s'))/T_r. The manuscript's
characteristic torque does not by itself specify that peak shape factor.

Using the manuscript's ideal information-rate convention, R_s=R/log2(M).
If T_r=beta/R_s, then

    characteristic torque = pi*I*abs(j)*delta_f_g*R/(beta*log2(M)).

Only if adjacent tone spacing also scales as delta_f_g=kappa*R_s does this
become proportional to R squared, with coefficient
pi*I*abs(j)*kappa/(beta*log2(M)^2). The general discussion leaves beta and kappa
as design parameters; the new desktop subsection fixes its own transition,
spacing and mapping. Seven-ary CF-FSK was reported as an optimum in the cited
electromechanical study; it is not asserted to be universally optimal here.
Practical seven-ary bit mapping/coding must also be specified to convert the
ideal log2(7) bits/symbol to a net payload rate.

For symmetric seven-tone spacing, the upper tone lies 3*delta_f_g above the
signal center. To keep the modeled upper speed at or below 1800 rpm, require
f_center <= 60-3*delta_f_g Hz, equivalently center rpm <= 1800-90*delta_f_g.
This remains an analysis-envelope condition, not a safe-speed certification.

The existing constrained actuator simulator evaluates a prescribed command;
its demonstration script is not an M-ary CF-FSK encoder or receiver BER test.

---

# Separate binary illustration (not the paper's CF-FSK requirement)

The following retained calculation is conditional on extra assumptions and
must not be used as the paper-specific torque-versus-rate result.

'''
    report='''## Conditional operating envelope and torque sizing

No actual safe maximum speed or modulation torque has yet been established.
1800 rpm is the studied design ceiling, not a certified safe operating limit.
The existing centrifugal FEA fixes the bore and bonds the inserts. It omits
the real shaft/bearings, hub connection, insert retention/contact, acceleration
loads, fatigue and forced response. Averaged stresses are useful comparisons,
but cannot be compared directly to a local yield allowable to certify speed.
The model's 140 MPa aluminium allowable is an assumption, not a material
certificate or validated regional-average acceptance criterion.

For a fixed linear-elastic model, centrifugal stress scales as rpm squared.
That scaling neither supplies missing failure criteria nor identifies the
rotor-bearing critical speeds. The fixed-bore carrier frequencies are not
critical-speed limits or maximum bit rates. See the prior design report.

## Explicit modulation assumptions

R is the raw binary symbol/bit rate before coding. Gravitational signal tone
separation is delta_f=kappa*R, and a changed symbol uses a raised-cosine frequency
transition of duration T=beta/R. These smooth waveforms are illustrative;
tone spacing alone does not prove orthogonality or adequate receiver BER.
At the second harmonic, omega=pi*f_g, so:

    torque = I*pi*df_g/dt
    peak torque = I*pi^2*delta_f/(2*T)
                = I*pi^2*kappa*R^2/(2*beta)

Defaults kappa=1 and beta=0.5 give peak torque = pi^2*I*R^2.
The transition has continuous torque but discontinuous torque slope at its
endpoints. Further smoothing changes the trajectory and must be included in
tracking/receiver calculations. Low-pass filtering cannot remove the required
angular impulse while preserving the same frequency change and transition time.
The absolute lower bound for any transition is I*pi*delta_f/T; the raised-cosine
peak is pi/2 times that bound. This is not a bound on every possible modulation.

The upper tone is held at 60 Hz (1800 rpm), the lower tone is 60-R Hz,
and the nominal center speed is 1800-15*R rpm. Centering the signal at 60 Hz
would exceed the 1800 rpm study ceiling on upper-tone symbols.

## 40 mm spokes / 20 mm rim / 3 mm fillets

Rotor inertia is 0.307167282 kg m^2. Peak torque is 3.03161956*R^2 N m.
Stored energy at 1800 rpm is about 5.457 kJ, and each insert requires about
15.73 kN radial retention force. Neither insert retention nor containment is
validated by the bonded carrier model.

| Raw bit/s | Rotor rpm range | Peak Nm | IID-data RMS Nm | Peak mechanical kW | Nominal shaft shear MPa |
|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(table)+'''

The power magnitude applies to both accelerating and braking transitions.
It is transient mechanical power, not average electrical consumption. Braking
requires absorption/regeneration capacity. RMS assumes independent equally
likely bits (transition probability 1/2); a continuously alternating stream
has RMS torque equal to half the peak under these assumptions. Add losses,
attached inertia and motor/controller margins for a drive specification.

The 12 mm solid-shaft nominal shear is 16*torque/(pi*d^3), excluding bending,
keyways, stress concentrations and hub/setscrew transfer. At 5 bit/s it is
about 223 MPa and at 10 bit/s about 894 MPa; no shaft safety margin is claimed.
This exposes a missing load case even when centrifugal carrier averages are low.

For these assumptions a rotor-side torque limit Tlim implies R <= sqrt(Tlim/3.03162).
For example 10 Nm implies about 1.82 raw bit/s from peak torque alone, provided
the continuous-duty, power, servo, shaft/hub, modal and receiver limits also pass.
This is not a demonstrated safe or decodable data rate.

Baseline and 30 mm spoke inertias and torque coefficients are in assumptions.json;
all three designs are tabulated in torque_vs_rate.csv. Numerical angular-impulse,
work-energy and torque-squared integrals independently check the formulas.

The next structural calculation required for a torque limit is a unit angular-
acceleration load case with the actual shaft/hub connection, followed by combined
centrifugal and alternating modulation loads. Speed qualification additionally
needs a rotor-bearing model, retention design, specified material/fatigue criteria
and experimental vibration/overspeed validation with suitable containment.
Steady centrifugal stress is stationary in the rotating material frame; it
should not automatically be counted as fully reversed fatigue once per turn.
Modulation, starts/stops, gravity, imbalance and support loads need their own
stress histories and cycle counts.

References: [SKF rotor/bearing resonance](https://evolution.skf.com/us/damping-in-a-rolling-bearing/),
[FSK parameter definitions](https://www.mathworks.com/help/comm/ref/comm.fskmodulator-system-object.html).
'''
    (ROOT/'OPERATING_ENVELOPE.md').write_text(correction+report)
    print(json.dumps(specs,indent=2));print('\n'.join(table))


if __name__=='__main__':main()
