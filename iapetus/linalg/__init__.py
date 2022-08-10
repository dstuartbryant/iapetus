"""Linear algebra tools."""

from typing import Union

import numpy as np


class Vector(np.ndarray):
    def __init__(self, arg: Union[np.ndarray, list]):
        if isinstance(arg, np.ndarray):
            self._vector = arg
        elif isinstance(arg, list):
            self._vector = np.array(arg)

        raise NotImplementedError("Need to handle vector shape constraints.")

    def __str__(self):
        return str(self._vector)

    def __repr__(self):
        return repr(self._vector)
