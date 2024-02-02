"""Batch least squares data modeling module."""

from typing import List

import numpy as np
from pydantic import BaseModel


class BatchDataBaseModel(BaseModel):
    class Config:
        arbitrary_types_allowed = True


class BatchData(BatchDataBaseModel):
    timestamp: float
    state_error: np.ndarray
    covariance: np.ndarray
    timesteps: List[float]
    residuals: List[np.ndarray]
    obs_error_covariances: List[np.ndarray]
