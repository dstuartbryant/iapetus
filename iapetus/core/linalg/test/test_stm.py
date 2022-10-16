import numpy as np
from iapetus.core.linalg import stm


def test_straight_constant_velocity():
    x = np.array([0, 0, 1, 1])
    dt = 1
    phi = stm.straight_constant_velocity(dt)
    for i in range(4):
        x = np.dot(phi, x)
        assert x[0] == i + 1
        assert x[1] == i + 1
        assert x[2] == 1
        assert x[3] == 1
