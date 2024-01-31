"""Kalman filter data models module."""

import numpy as np
from pydantic import BaseModel


class KalmanFilterDataPointBaseModel(BaseModel):
    class Config:
        arbitrary_types_allowed = True


class LkfDataPoint(KalmanFilterDataPointBaseModel):
    """Linearized Kalman filter data point model class."""

    k: int  # Filter iteration step
    state_k: np.ndarray
    state_error_k: np.ndarray
    covariance_k: np.ndarray
    error_state_transition_matrix_k_plus_1_k: np.ndarray
    state_k_plus_1_given_k: np.ndarray
    state_error_k_plus_1_given_k: np.ndarray
    Q: np.ndarray
    covariance_k_plus_1_given_k: np.ndarray
    prefit_resid: np.ndarray
    state_k_plus_1: np.ndarray
    state_error_k_plus_1: np.ndarray
    covariance_k_plus_1: np.ndarray


class LkfDataPoint2(KalmanFilterDataPointBaseModel):
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
