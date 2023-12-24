"""State test fixtures."""

import numpy as np
import pytest


class SatelliteState1:
    pi = 6656356.11057065  # [m]
    pj = 1700859.15707779  # [m]
    pk = 299734.38071253  # [m]
    vi = -1794.25660717  # [mps]
    vj = 6353.55570765  # [mps]
    vk = 3792.38315729  # [mps]
    area = 1 / 3**2  # [m^2]
    mass = 50  # [kg]
    Cd = 2.0  # [unitless]

    @property
    def p_vector(self):
        return [self.pi, self.pj, self.pk]

    @property
    def v_vector(self):
        return [self.vi, self.vj, self.vk]

    @property
    def p_mag(self):
        return np.linalg.norm(self.p_vector)

    @property
    def v_mag(self):
        return np.linalg.norm(self.v_vector)


@pytest.fixture
def sat_state_1():
    yield SatelliteState1()
