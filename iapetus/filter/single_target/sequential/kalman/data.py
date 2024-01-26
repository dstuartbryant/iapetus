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
