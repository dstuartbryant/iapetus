"""Celestial body parameter pytest fixture module."""

import pytest

from iapetus.propagation.dynamics.nonlinear.astro.constants import (
    EARTH_SHAPE,
    EQUATORIAL_RADIUS,
    MASS,
    MU,
    PERIOD,
    ROTATION,
    SEMIMAJOR_AXIS,
)


class EarthParams:
    mu = MU["Earth"]
    period = PERIOD["Earth"]
    semi_major_axis = SEMIMAJOR_AXIS["Earth"]
    equatorial_radius = EQUATORIAL_RADIUS["Earth"]
    j2 = EARTH_SHAPE["J2"]
    j3 = EARTH_SHAPE["J3"]
    j4 = EARTH_SHAPE["J4"]
    mass = MASS["Earth"]
    rotation_rate = ROTATION["Earth"]


@pytest.fixture
def earth_params():
    yield EarthParams()
