"""State models, etc.

TODO: Devise a more descriptive doc string for this file.
"""

from typing import Union

from ..epoch import Epoch
from ..linalg import Vector


class StateVector(Vector):
    """State vector base class.

    Attributes:
        pos (Vector): position vector.
        pos_units (Unit): position vector measurement units.
        vel (Vector): velocity vector.
        vel_units (Unit): velocity vector measurement units.
        acc (Vector, optional): acceleration vector.
        acc_units (Unit, optional): acceleration vector measurement units.
        reference_frame (RefFrame): indicates reference frame for all vector
            coordinates.
        vector: concatenation of all used vectors.

    """

    def __init__(self, pos: Vector, vel: Vector):
        pass


class State:
    """State base class.

    Attributes:
        epoch (Epoch): time component of a state.
        vector (StateVector): vector component of a state.
        reference_frame (RefFrame): indicates reference frame for state vector
            coordinates.

    """

    def __call__(self, t):
        """Propagates state to input epoch and returns a new State
        instantiation at the new epoch.

        Args:
            t (Epoch): timestamp (epoch) to which the state will be propagated
        """
        raise NotImplementedError

    def _propagate(self, t: Union[Epoch, float]):
        """Propagates state vector to a point in time.

        Args:
            t (Epoch or float): point in time to which state vector will be
                propagated.
        """
        raise NotImplementedError
