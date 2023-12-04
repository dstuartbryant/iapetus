"""Runge-Kutta base module."""

from abc import ABC, abstractmethod
from typing import Optional

import numpy as np


class DerivativeFunction(ABC):
    @abstractmethod
    def derivative_fcn(self, t: float, y: np.ndarray) -> np.ndarray:
        """Derivative function model for use with Runge-Kutta integrators."""
        pass


class RungeKutta:
    """Runge-Kutta integrator base class."""

    def __init__(
        self,
        derivative_fcn: DerivativeFunction.derivative_fcn,
        tspantol: float,
        extras: Optional[dict] = None,
    ):
        """
        Args:
            derivative_fcn
        """
        self.derivative_fcn = derivative_fcn
        self.tspantol = tspantol
        self.extras = extras
