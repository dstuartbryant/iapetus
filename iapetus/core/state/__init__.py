"""State models, etc.

TODO: Devise a more descriptive doc string for this file.
"""

from dataclasses import dataclass
from typing import Optional, Union

from ..epoch import Epoch
from ..linalg import Vector
from ..refframes import ReferenceFrame


@dataclass
class StateVector:
    """State vector base dataclass.

    Attributes:
        pos: position vector.
        vel: velocity vector.
        acc: acceleration vector.
        frame: indicates reference frame for all vector coordinates.
    """

    pos: Vector
    vel: Vector
    acc: Optional[Vector]
    frame: ReferenceFrame


@dataclass
class State:
    """State base class.

    Attributes:
        epoch: time component of a state.
        vector: vector component of a state.
        frame: indicates reference frame for state vector coordinates.

    """

    epoch: Epoch
    vector: StateVector

    @property
    def frame(self):
        """Relays reference frame object from vector attribute."""
        return self.vector.frame
