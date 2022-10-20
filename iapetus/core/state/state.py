from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from ..linalg import Vector


@dataclass
class StateVector:
    """State vector base dataclass.

    Attributes:
        pos: position vector.
        vel: velocity vector.
        acc: acceleration vector, optional.
    """

    pos: Vector
    vel: Vector
    acc: Optional[Vector] = None

    def _full_vector(self):
        x = np.concatenate([self.pos._vector, self.vel._vector])
        if self.acc:
            x = np.concatenate([x, self.acc._vector])
        return x

    @property
    def vector(self):
        return self._full_vector()

    @classmethod
    def from_array(cls, arg):
        """Instantiates a StateVector class from an array (or list).

        Assumes 2D vector space.

        Args:
            arg (list or np.ndarray): Array-like object with 4 or 6 elements.
        """
        if len(arg) == 6:
            return StateVector(
                pos=Vector(arg[:2]), vel=Vector(arg[2:4]), acc=Vector(arg[4:])
            )
        return StateVector(pos=Vector(arg[:2]), vel=Vector(arg[2:4]))


@dataclass
class State:
    """State base class.

    Attributes:
        epoch: time component of a state.
        vector: vector component of a state.

    """

    timestamp: float
    vector_obj: StateVector

    @property
    def time(self):
        return self.timestamp

    @property
    def vector(self):
        return self.vector_obj.vector

    @classmethod
    def from_array(cls, t, x):
        return cls(timestamp=t, vector_obj=StateVector.from_array(x))
