"""Kalman filter process noise module."""

from abc import ABC, abstractmethod
from typing import Union

import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro import snc as astro_snc
from iapetus.propagation.dynamics.nonlinear.astro.xforms import ntw_matrix


class ProcessNoise(ABC):
    """Kalman filter process noise base class."""

    @abstractmethod
    def __call__(self, dt_s: float, state_mean: np.ndarray) -> np.ndarray:
        """Returns process noise matrix in inertial coordintes."""
        pass


class ZeroProcessNoise(ProcessNoise):
    """Zero process noise subclass.

    Always returns matrix of zeros.
    """

    def __call__(self, dt_s: float, state_mean: np.ndarray) -> np.ndarray:
        """Returns process noise matrix, all zeros."""
        return np.zeros((len(state_mean), len(state_mean)))


class AstroProcessNoise(ProcessNoise):
    """Kalman filter process noise subclass for astrodynamics applications."""

    def __init__(
        self, Q: Union[astro_snc.QInertialRic, astro_snc.QInertialCd]
    ):
        self.Q = Q

    def __call__(self, dt_s: float, state_mean: np.ndarray) -> np.ndarray:
        """Returns process noise matrix in inertial coordintes.

        TODO: This method is FRAGILE because it assumes that the first three
        elements of `state_mean` are position vector elements and that the next
        three elements of `state_mean` are velocity vector elements. This may
        need to be revised somehow.
        """

        R_inertial_to_ntw = ntw_matrix(state_mean[:3], state_mean[3:6])
        return self.Q(dt_s, R_inertial_to_ntw)
