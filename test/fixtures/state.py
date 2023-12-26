"""State test fixtures."""


import pytest

from iapetus.propagation.dynamics.nonlinear.astro.eom.payloads import (
    TwoBodyDragState,
)


@pytest.fixture
def sat_state_1():
    yield TwoBodyDragState(
        pi=6656356.11057065,  # [m]
        pj=1700859.15707779,  # [m]
        pk=299734.38071253,  # [m]
        vi=-1794.25660717,  # [mps]
        vj=6353.55570765,  # [mps]
        vk=3792.38315729,  # [mps]
        A=(1 / 3) ** 2,  # [m^2]
        m=50,  # [kg]
        Cd=2.0,  # [unitless]
    )
