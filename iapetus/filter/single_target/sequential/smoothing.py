"""Smoothing module for Kalman filtering.

References:
[1] [1] Born, et al., Statistical Orbit Determination
"""

from typing import List

import numpy as np

from .kalman.data import LkfDataPoint, LkfDataPoint2, SmoothedLkfDataPoint


class SmootherError(Exception):
    pass


def lkf_smoother(data: List[LkfDataPoint2]):
    """Smoother for linearize Kalman filter (LKF).

    NOTE: Assumptions:

    [1] The last LKF iteration is indexed by k, and smoother iterations are
        indexed by ell, so to start, assume L = k + 1, i.e., k = L - 1.

    [2] Smoother iterates all the way to k=0, which assumes that the LKF had an
        observation at t_k = t_0, i.e., the initial state X0 provided to the
        LKF was already propagated to the first observation time at t_0.

    [3] In general, L > k
    """

    L = data[-1].k

    # Stash last foward-filter state error and covariance at step L
    smoothed_data = [
        SmoothedLkfDataPoint(
            k=L, state_error=data[L].state_error, covariance=data[L].covariance
        )
    ]

    # First smoother iteration
    # Set k = L - 1
    k = L - 1
    xhat_L_L = data[L].state_error
    xhat_k_k = data[k].state_error
    P_L_L = data[L].covariance
    P_k_k = data[k].covariance
    Phi = data[k + 1].error_state_transition_matrix
    Q = data[k + 1].process_noise
    P_update = Phi @ P_k_k @ Phi.T + Q
    S_k = P_k_k @ Phi.T @ np.linalg.inv(P_update)
    x_k_L = xhat_k_k + S_k @ (xhat_L_L - Phi @ xhat_k_k)
    P_k_L = P_k_k + S_k @ (P_L_L - P_update) @ S_k.T

    smoothed_data.append(
        SmoothedLkfDataPoint(k=k, state_error=x_k_L, covariance=P_k_L)
    )

    while k > 0:
        k -= 1
        P_k_k = data[k].covariance
        Phi = data[k + 1].error_state_transition_matrix
        Q = data[k + 1].process_noise
        P_update = Phi @ P_k_k @ Phi.T + Q
        S_k = P_k_k @ Phi.T @ np.linalg.inv(P_update)
        x_k_L = xhat_k_k + S_k @ (x_k_L - Phi @ xhat_k_k)
        P_k_L = P_k_k + S_k @ (P_k_L - P_update) @ S_k.T

        smoothed_data.append(
            SmoothedLkfDataPoint(k=k, state_error=x_k_L, covariance=P_k_L)
        )

    return smoothed_data
