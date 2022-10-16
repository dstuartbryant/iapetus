"""Linear state propagators."""

import numpy as np

from ..linalg.stm import straight_constant_velocity
from ..state import State


class StraightConstantVelocity:
    """Straight (direction) and constant velocity linear propagator class."""

    def __init__(self, dt):
        """Args:
        dt (float): time step, i.e., sample rate
        """
        self.dt = dt
        self.mat = straight_constant_velocity(dt)

    def __call__(self, x: State):
        """Returns state `x` propagated to the next time step.

        Args:
            x_0 (State): starting state
        """

        return State.from_array(x.time + self.dt, np.dot(self.mat, x.vector))
