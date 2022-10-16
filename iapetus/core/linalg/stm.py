"""State transition matrix module.
"""
import numpy as np

"""
TODO: Dev a baseline STM that assumes it's getting np.ndarray inputs - likely
    all that is currently necessary
"""


def straight_constant_velocity(dt):
    """Returns a 2D state transition matrix modeling straight movement with
    constant velocity.

    I.e., an object propagated with this model will move in the same
    direction without turning at constant speed.

    Args:
            dt (float): time increment (aka sample rate)
    """

    return np.block(
        [[np.eye(2), np.eye(2) * dt], [np.zeros((2, 2)), np.eye(2)]]
    )
