"""Rotor stored-energy and centrifugal-load estimates."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .designs import Deployment


def energy_loss_power(kinetic_energy_j: float, fraction_per_hour: float) -> float:
    """Instantaneous P=fraction*E/3600; fraction=0.01 means 1% per hour.

    This is an explicit illustrative loss-rate input, not a vacuum-drag model
    or a prediction that a real machine has this loss rate.
    """
    if any(not math.isfinite(x) or x < 0 for x in (kinetic_energy_j, fraction_per_hour)):
        raise ValueError("energy and fractional loss rate must be finite and nonnegative")
    return kinetic_energy_j * fraction_per_hour / 3600.0


def modulation_step(inertia: float, carrier_hz: float, delta_carrier_hz: float, ramp_seconds: float,
                    harmonic: int = 2) -> dict[str, float]:
    """Linear received-frequency ramp; torque, signed energy and peak drive power.

    Assumes rigid point-mass inertia and ideal reversible actuation. It does
    not certify slew feasibility, motor efficiency, or a bit rate.
    """
    if any(not math.isfinite(x) or x <= 0 for x in (inertia, carrier_hz, ramp_seconds)):
        raise ValueError("inertia, carrier and ramp time must be finite and positive")
    if not isinstance(harmonic, int) or harmonic < 1:
        raise ValueError("harmonic must be a positive integer")
    if not math.isfinite(delta_carrier_hz) or carrier_hz + delta_carrier_hz <= 0:
        raise ValueError("end carrier must be finite and positive")
    start = 2*math.pi*carrier_hz/harmonic
    finish = 2*math.pi*(carrier_hz+delta_carrier_hz)/harmonic
    torque = inertia*(finish-start)/ramp_seconds
    return {"torque_nm": torque,
            "initial_mechanical_power_w": torque*start,
            "peak_absolute_mechanical_power_w": abs(torque)*max(start, finish),
            "energy_change_j": 0.5*inertia*(finish**2-start**2)}


@dataclass(frozen=True)
class RotorMaterial:
    name: str
    allowable_stress: float
    density: float = 0.0
    youngs_modulus: float = 0.0


STEEL = RotorMaterial("steel", 0.5e9, 7850.0, 200e9)
CFRP = RotorMaterial("CFRP", 2.0e9, 1800.0, 150e9)
ALUMINUM = RotorMaterial("aluminum (6061-T6)", 0.14e9, 2700.0, 69e9)


@dataclass(frozen=True)
class RotorMechanicalState:
    material: RotorMaterial
    angular_speed: float
    moment_of_inertia: float
    kinetic_energy: float
    peak_retention_force: float
    required_tie_area: float
    equivalent_tie_diameter: float
    specific_energy: float


def material_for_deployment(deployment: Deployment) -> RotorMaterial:
    """Return the assumed load-bearing arm material for a deployment class."""

    if deployment.key == "heavy":
        return CFRP
    return STEEL


def rotor_mechanical_state(deployment: Deployment) -> RotorMechanicalState:
    """Compute simple point-mass rotor mechanical quantities.

    The peak retention force is the largest centrifugal load carried by one arm
    or tie member. The required tie area assumes the allowable stress already
    includes the selected safety factor.
    """

    rotor = deployment.rotor()
    material = material_for_deployment(deployment)
    angular_speed = deployment.tip_speed / rotor.max_arm
    inertia = rotor.quadrupole_moment
    kinetic_energy = 0.5 * inertia * angular_speed * angular_speed
    forces = [mass * angular_speed * angular_speed * abs(arm) for mass, arm in zip(rotor.masses, rotor.arms)]
    peak_retention_force = max(forces)
    required_tie_area = peak_retention_force / material.allowable_stress
    equivalent_tie_diameter = math.sqrt(4.0 * required_tie_area / math.pi)
    specific_energy = kinetic_energy / rotor.total_mass
    return RotorMechanicalState(
        material=material,
        angular_speed=angular_speed,
        moment_of_inertia=inertia,
        kinetic_energy=kinetic_energy,
        peak_retention_force=peak_retention_force,
        required_tie_area=required_tie_area,
        equivalent_tie_diameter=equivalent_tie_diameter,
        specific_energy=specific_energy,
    )


@dataclass(frozen=True)
class CarrierPlate:
    """Spoked carrier disk that holds the two signal masses as inserts.

    Models a machined disk as an outer rim (annulus) joined to a central hub by
    ``n_spokes`` radial spokes, with the material between spokes cut away.  The
    plate adds rotating inertia (stored energy and modulation torque) and
    provides a distributed load path that retains the inserts. A spoked plate
    is not axisymmetric: it can contribute higher gravitational harmonics, and
    two insert bosses/pockets can also contribute at the second harmonic.
    Use FiniteRotor for the complete shaped-source signal. Inertias
    are polar (about the spin axis); spokes are treated as thin radial bars.
    """

    outer_radius: float
    thickness: float
    rim_width: float
    n_spokes: int
    spoke_width: float
    hub_radius: float
    material: RotorMaterial = ALUMINUM

    def __post_init__(self) -> None:
        if self.rim_width >= self.outer_radius:
            raise ValueError("rim_width must be smaller than outer_radius")
        if self.hub_radius >= self.outer_radius - self.rim_width:
            raise ValueError("hub_radius must be inside the rim")
        if self.n_spokes < 1:
            raise ValueError("n_spokes must be >= 1")

    @property
    def inner_radius(self) -> float:
        """Inner edge of the rim (= outer end of the spokes)."""

        return self.outer_radius - self.rim_width

    @property
    def mass(self) -> float:
        rho_t = self.material.density * self.thickness
        rim = math.pi * (self.outer_radius ** 2 - self.inner_radius ** 2)
        spoke = self.n_spokes * self.spoke_width * (self.inner_radius - self.hub_radius)
        hub = math.pi * self.hub_radius ** 2
        return rho_t * (rim + spoke + hub)

    @property
    def polar_inertia(self) -> float:
        rho_t = self.material.density * self.thickness
        rim_mass = rho_t * math.pi * (self.outer_radius ** 2 - self.inner_radius ** 2)
        rim_inertia = 0.5 * rim_mass * (self.outer_radius ** 2 + self.inner_radius ** 2)
        a, b = self.hub_radius, self.inner_radius
        spoke_mass = self.n_spokes * rho_t * self.spoke_width * (b - a)
        spoke_inertia = spoke_mass * (b ** 3 - a ** 3) / (3.0 * (b - a))
        hub_mass = rho_t * math.pi * self.hub_radius ** 2
        hub_inertia = 0.5 * hub_mass * self.hub_radius ** 2
        return rim_inertia + spoke_inertia + hub_inertia

    @property
    def fill_factor(self) -> float:
        """Mass fraction relative to a solid disk of the same radius/thickness."""

        solid = self.material.density * self.thickness * math.pi * self.outer_radius ** 2
        return self.mass / solid

    def rim_hoop_stress(self, angular_speed: float) -> float:
        """Hoop stress in the spinning rim, rho*(omega*R)^2 (thin-ring approx)."""

        return self.material.density * (angular_speed * self.outer_radius) ** 2


@dataclass(frozen=True)
class DiskRotorMechanicalState:
    angular_speed: float
    insert_inertia: float
    plate_mass: float
    plate_inertia: float
    total_inertia: float
    inertia_overhead: float
    kinetic_energy: float
    modulation_torque_factor: float
    rim_hoop_stress: float
    rim_hoop_margin: float
    insert_retention_force: float
    insert_bearing_stress: float
    insert_bearing_margin: float


def disk_rotor_mechanical_state(
    deployment: Deployment,
    plate: CarrierPlate,
) -> DiskRotorMechanicalState:
    """Mechanical state of a disk-carrier rotor for a deployment class.

    The point-model inserts are the deployment's two masses; the plate adds
    inertia and retains them.  ``inertia_overhead`` and ``modulation_torque_factor``
    quantify the cost of the dead plate (both scale with total/insert inertia).
    The rim-hoop and insert-bearing margins are allowable_stress/stress (>1 is
    nominal margins under the material's built-in safety factor). The bearing
    area rim_width*thickness is a legacy gross-area proxy, not the projected
    diameter*engagement area of the CAD insert or a bonded-interface stress.
    """

    rotor = deployment.rotor()
    angular_speed = deployment.tip_speed / rotor.max_arm
    insert_inertia = rotor.quadrupole_moment
    plate_inertia = plate.polar_inertia
    total_inertia = insert_inertia + plate_inertia
    kinetic_energy = 0.5 * total_inertia * angular_speed * angular_speed
    retention = max(mass * angular_speed * angular_speed * abs(arm) for mass, arm in zip(rotor.masses, rotor.arms))
    bearing_area = plate.rim_width * plate.thickness
    bearing_stress = retention / bearing_area
    hoop_stress = plate.rim_hoop_stress(angular_speed)
    return DiskRotorMechanicalState(
        angular_speed=angular_speed,
        insert_inertia=insert_inertia,
        plate_mass=plate.mass,
        plate_inertia=plate_inertia,
        total_inertia=total_inertia,
        inertia_overhead=plate_inertia / insert_inertia,
        kinetic_energy=kinetic_energy,
        modulation_torque_factor=total_inertia / insert_inertia,
        rim_hoop_stress=hoop_stress,
        rim_hoop_margin=plate.material.allowable_stress / hoop_stress,
        insert_retention_force=retention,
        insert_bearing_stress=bearing_stress,
        insert_bearing_margin=plate.material.allowable_stress / bearing_stress,
    )


def cantilever_bending_stiffness(youngs_modulus: float, radius: float, length: float) -> float:
    """Lateral tip stiffness of a solid circular cantilever, 3 E I / L^3 (N/m)."""

    if length <= 0 or radius <= 0:
        raise ValueError("radius and length must be positive")
    area_moment = math.pi * radius ** 4 / 4.0
    return 3.0 * youngs_modulus * area_moment / length ** 3


def whirl_critical_hz(stiffness: float, mass: float) -> float:
    """First bending (whirl) critical frequency sqrt(k/m)/2pi (Hz)."""

    if stiffness <= 0 or mass <= 0:
        raise ValueError("stiffness and mass must be positive")
    return math.sqrt(stiffness / mass) / (2.0 * math.pi)
