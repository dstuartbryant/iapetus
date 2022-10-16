"""Linear algebra tools."""

from typing import List, Union

import numpy as np


class VectorException(Exception):
    pass


class Vector:
    def __init__(self, arg: Union[np.ndarray, List[float]]):
        if isinstance(arg, np.ndarray):
            if arg.ndim != 1:
                raise VectorException(
                    f"Expecting number of dimensions length of 1, found "
                    f"{arg.ndim}."
                )
            self._vector = arg
        elif isinstance(arg, list):
            self._vector = np.array(arg)
            if len([x for x in arg if isinstance(x, list)]) != 0:
                raise VectorException("Not expecting nested lists.")

    def __str__(self):
        return str(self._vector)

    def __repr__(self):
        return repr(self._vector)

    def _dim(self):
        return self._vector.shape[0]

    def _mag(self):
        return np.linalg.norm(self._vector)

    def _unit(self):
        return self._vector / self._mag()

    @property
    def dim(self):
        return self._dim()

    @property
    def mag(self):
        return self._mag()

    @property
    def unit(self):
        return self._unit()
