# Paper CF-FSK torque sizing: correction to the binary illustration

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
| baseline_rounded | 0.328813 | 1.032996 |
| medium_rounded | 0.316071 | 0.992967 |
| wide_rounded | 0.307167 | 0.964994 |

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

## Conditional operating envelope and torque sizing

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
| 0.1 | 1797-1800 | 0.030 | 0.011 | 0.0057 | 0.09 |
| 0.5 | 1785-1800 | 0.758 | 0.268 | 0.1423 | 2.23 |
| 1 | 1770-1800 | 3.032 | 1.072 | 0.5667 | 8.94 |
| 2 | 1740-1800 | 12.126 | 4.287 | 2.2480 | 35.74 |
| 5 | 1650-1800 | 75.790 | 26.796 | 13.7038 | 223.38 |
| 10 | 1500-1800 | 303.162 | 107.184 | 52.5969 | 893.51 |

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
