"""Probabilistic state data model module."""

from dataclasses import dataclass

import numpy as np


@dataclass
class Covariance:
    """Probabilistic state covariance component class."""

    matrix: np.ndarray

    def _standard_deviations(self) -> np.ndarray:
        """Returns an array of the covariance matrix's standard deviations."""
        return np.sqrt(np.diag(self.matrix))

    @property
    def standard_deviations(self) -> np.ndarray:
        return self._standard_deviations()


@dataclass
class State:
    """Probabilistic state composite class."""

    timestamp: float
    mean: np.ndarray
    covariance: Covariance

    def __post_init__(self):
        if not isinstance(self.covariance, Covariance):
            self.covariance = Covariance(self.covariance)

    @property
    def dimension(self) -> int:
        return len(self.mean)
