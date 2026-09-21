# Denser-insert option (geometry estimate, not a new qualified design)

User asked whether tungsten could reduce windage by reducing insert volume at fixed mass. Yes. An illustrative tungsten heavy alloy density of 18,000 kg/m3 versus the existing steel 7,850 kg/m3 gives volume/length ratio 0.4361 at fixed cylindrical diameter. Preserve insert mass, radial center and diameter; shorten the axial length. The leading rotating quadrupole and insert polar inertia are preserved for this ideal cylindrical substitution. Finite-volume field corrections, elastic modes, local contact/retention and actual alloy properties must be recomputed.

The insert length reduces from 79.2404 to 34.5576 mm in the reference geometry, while the carrier is 20 mm thick. Total protruding height falls from 59.2404 to 14.5576 mm (about 75.4% reduction). This is the reduction of projected side area above/below the plate for fixed diameter, not a demonstrated reduction in whole-rotor drag or total drive power. Rim/spokes/skin friction remain, and flow interaction and Cd change. For an isolated contribution at fixed speed, P_drag=0.5*rho_air*Cd*A*v^3 motivates the estimate. Shorter overhang may also reduce insert bending demands; radial centrifugal load at fixed mass/radius/speed does not decrease.

Manufacturer density reference: https://www.globaltungsten.com/tungsten-molybdenum-parts/ (Class 3 heavy alloy, nominal 18 g/cm3). Machinable heavy-alloy alternatives: https://tungstenparts.com/products/tungsten-alloys/ . Drag relation: https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/ . Procurement/large-billet availability and alloy-specific strength need assessment; no pricing assumed.

The current steel FEA baseline continues. This is a promising subsequent materials/shape comparison, not a substitution into the frozen designs.

[
  {
    "distance_m": 0.5,
    "insert_diameter_mm": 67.28297892848694,
    "steel_length_mm": 76.16475579216194,
    "tungsten_alloy_length_mm": 33.21629627602618,
    "carrier_thickness_mm": 19.223708265281985,
    "protruding_height_ratio": 0.24573815583808437
  },
  {
    "distance_m": 1.0,
    "insert_diameter_mm": 117.77606278198382,
    "steel_length_mm": 133.3232446424003,
    "tungsten_alloy_length_mm": 58.14374835793569,
    "carrier_thickness_mm": 33.65030365199538,
    "protruding_height_ratio": 0.24573815583808434
  },
  {
    "distance_m": 2.0,
    "insert_diameter_mm": 214.38015373022415,
    "steel_length_mm": 242.67968386036264,
    "tungsten_alloy_length_mm": 105.83530657243594,
    "carrier_thickness_mm": 61.25147249434975,
    "protruding_height_ratio": 0.24573815583808437
  }
]
