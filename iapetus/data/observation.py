"""Observation data model module."""

from dataclasses import dataclass
from typing import List

import numpy as np

from .state.probabilistic import State


class ProbabilisticObservation(State):
    pass


@dataclass
class ProbabilisticObservationSet:
    observations: List[ProbabilisticObservation]
    H: np.ndarray
