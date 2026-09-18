"""Physical constants and plotting unit conversions."""

import math

G = 6.6743e-11
MU0 = 4.0 * math.pi * 1e-7
K_B = 1.380649e-23

GAL = 1e-2
MICROGAL = 1e-8
NANOGAL = 1e-11
# Standard gravity (g) and Gal are different acceleration units.
STANDARD_GRAVITY = 9.80665
NANOG = 1e-9 * STANDARD_GRAVITY

DEFAULT_QUADRUPOLE_COEFFICIENT = 2.25
