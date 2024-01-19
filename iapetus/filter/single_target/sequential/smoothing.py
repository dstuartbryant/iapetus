"""Smoothing module for Kalman filtering.

References:
[1] [1] Born, et al., Statistical Orbit Determination
"""

from typing import List

import numpy as np


class SmootherError(Exception):
    pass


class Smoother:
    """Kalman filtering smoother class."""

    def __init__(
        self,
        X: List[np.ndarray],
        P: List[np.ndarray],
        PHI: List[np.ndarray],
        Q: List[np.ndarray],
    ):
        """
        Args:
            X (List[np.ndarray]): List of ell number of state estimates in
                chronological order, i.e., len(X) = ell
            P (List[np.ndarray]): List of ell number of state covariances in
                chronological order, i.e., len(P) = ell
            PHI (List[np.ndarray]): List of ell-1 number of state transition
                matrices in chronological order, i.e., len(PHI) = ell-1
            Q (List[np.ndarray]): List of ell-1 number of process noise
                matrices (that encapsulate their respective process noise
                transition matrices)in chronological order, i.e.,
                len(Q) = ell-1
        """
        ell = len(X)
        if len(P) != ell:
            raise SmootherError(
                "The nuumber of covariance matrices input to smoother must be "
                "equivalent to number of state estimate."
            )
        if len(PHI) != ell - 1:
            raise SmootherError(
                "The number of state transition matrices input to smoother "
                "must be equivalent to number of state estimates minus 1."
            )
        if len(Q) != ell - 1:
            raise SmootherError(
                "The number of process noise matrices input to smoother must "
                "be equivalent to number of state estimates minus 1."
            )
        self.ell = ell
        X.sort(reverse=True)
        P.sort(reverse=True)
        PHI.sort(reverse=True)
        Q.sort(reverse=True)
        self.X = X
        self.P = P
        self.PHI = PHI
        self.Q = Q
        self.X_smoothed = [self.X[-1]]
        self.P_intermediate = []
        self.S = []

    def __call__(self):
        """Runs smoothing algorithm."""

        ell = self.ell
        count = 0
        for x_k_given_k in self.X[1:]:
            count += 1

            # If ell = 3, k = 2 (first iteration)
            k = ell - count

            # If ell = 3, k = 2, we need P_k|k = P_2|2, but P_3|3 is last value
            # in self.P (first iteration). So we don't want self.P[0], we want
            # self.P[1], and k - 1 = 2 - 1 = 1.
            P_k_given_k = self.P[k - 1]

            # If ell = 3, k = 2, we need Phi(t_3, t_2), which is last value in
            # self.PHI (first iteration). So wee need self.PHI[0], and
            # ell - k - 1 = 3 - 2 - 1 = 0
            Phi_k_plus_1_from_k = self.Phi[ell - k - 1]

            # If ell = 3, k = 2, we need Q_2, which is self.Q[1]
            # (first iteration). ell - k = 3 - 2 = 1
            Q_k = self.Q[ell - k]

            P_k_plus_1_given_k = (
                Phi_k_plus_1_from_k @ P_k_given_k @ Phi_k_plus_1_from_k.T + Q_k
            )
            S_k = (
                P_k_given_k
                @ Phi_k_plus_1_from_k.T
                @ np.linalg.inv(P_k_plus_1_given_k)
            )
            x_k_given_ell = x_k_given_k + S_k @ (
                self.X_smoothed[-1] - Phi_k_plus_1_from_k @ x_k_given_k
            )

            self.P_intermediate.append(P_k_plus_1_given_k)
            self.S.append(S_k)
            self.X_smoothed.append(x_k_given_ell)

        P_k_given_ell = 

        return x_k_given_ell
