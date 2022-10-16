"""Linear propagator tests."""

from iapetus.core import state
from iapetus.core.propagation.linear import StraightConstantVelocity


def test_StraightConstantVelocity():
    X = [0, 0, 1, 1]
    dt = 1

    x_init = state.State.from_array(0, X)
    phi = StraightConstantVelocity(dt)

    x = x_init
    for i in range(4):
        x = phi(x)
        assert x.vector[0] == i + 1
        assert x.vector[1] == i + 1
        assert x.vector[2] == 1
        assert x.vector[3] == 1
