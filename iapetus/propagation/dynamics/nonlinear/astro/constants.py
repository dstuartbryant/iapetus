"""
Astrodynamics constants module.

References:
    [1] D. A. Vallado. Fundamentals of astrodynamics and applications.
        Springer, 4th edition, 2013.
"""

SOLAR_DAY = 86400  # [s]
SIDEREAL_DAY = 86164.090517  # [s]

SOLAR_YEAR = SOLAR_DAY * 365

# Gravitational constants [m^3/s^2]
MU = {
    "Sun": 1.32712428e11 * (1e3) ** 3,
    "Earth": 398600.4415 * (1e3) ** 3,
    "Moon": 4902.799 * (1e3) ** 3,
    "Mercury": 2.2032e4 * (1e3) ** 3,
    "Venus": 3.257e5 * (1e3) ** 3,
    "Mars": 4.305e4 * (1e3) ** 3,
    "Jupiter": 1.268e8 * (1e3) ** 3,
    "Saturn": 3.794 * (1e3) ** 3,
    "Uranus": 5.794 * (1e3) ** 3,
    "Neptune": 6.809 * (1e3) ** 3,
    "Pluto": 9.00e2 * (1e3) ** 3,
}

# Orbital periods [s]
PERIOD = {
    "Earth": 0.99997862 * SOLAR_YEAR,
    "Moon": 0.0748 * SOLAR_YEAR,
    "Mercury": 0.24084445 * SOLAR_YEAR,
    "Venus": 0.61518257 * SOLAR_YEAR,
    "Mars": 1.88071105 * SOLAR_YEAR,
    "Jupiter": 11.856525 * SOLAR_YEAR,
    "Saturn": 29.423519 * SOLAR_YEAR,
    "Uranus": 83.747406 * SOLAR_YEAR,
    "Neptune": 163.7232045 * SOLAR_YEAR,
    "Pluto": 248.0208 * SOLAR_YEAR,
}

# Semimajor Axis [m]
SEMIMAJOR_AXIS = {
    "Earth": 149598023e3,
    "Moon": 384400e3,
    "Mercury": 57909083e3,
    "Venus": 108208601e3,
    "Mars": 227939186e3,
    "Jupiter": 778298316e3,
    "Saturn": 1429394133e3,
    "Uranus": 2875038615e3,
    "Neptune": 4504449769e3,
    "Pluto": 5915799000e3,
}

# Mean Equatorial Radius [m]
EQUATORIAL_RADIUS = {"Earth": 6378.1363e3}

# Earth oblateness constants [unitless]
EARTH_SHAPE = {"J2": 0.0010826267, "J3": -0.0000025327, "J4": -0.0000016196}

# Mass [kg]
MASS = {"Sun": 1.9891e30, "Earth": 5.9742e24}

# Plantary rotation rate [rad/ solar sec]
ROTATION = {"Earth": 0.0000729211585530}


# Gravitational constant
GRAV_CONST = 6.673e-20 * (1e3) ** 3  # [m^3/(kg*s^2)]

REL_TOL_0_AND_1 = 1e-6
