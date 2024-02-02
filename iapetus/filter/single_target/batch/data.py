"""Batch least squares data modeling module."""

from typing import List

import numpy as np
from pydantic import BaseModel


class BatchDataBaseModel(BaseModel):
    class Config:
        arbitrary_types_allowed = True


class BatchData(BatchDataBaseModel):
    timestamp: float
    timesteps: List[float]
    state_error: np.ndarray
    covariance: np.ndarray
    residuals: List[np.ndarray]
    error_state_transition_matrices: List[np.ndarray]
