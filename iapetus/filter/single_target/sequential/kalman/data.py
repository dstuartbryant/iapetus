"""Kalman filter data models module."""

import numpy as np
from pydantic import BaseModel


class KalmanFilterDataPointBaseModel(BaseModel):
    class Config:
        arbitrary_types_allowed = True


class LkfDataPoint(KalmanFilterDataPointBaseModel):
    k: int
    timestamp: float
    predicted_state: np.ndarray
    updated_state: np.ndarray
    state_error: np.ndarray
    error_state_transition_matrix: np.ndarray
    covariance: np.ndarray
    process_noise: np.ndarray


class SmoothedLkfDataPoint(KalmanFilterDataPointBaseModel):
    k: int
    state_error: np.ndarray
    covariance: np.ndarray
