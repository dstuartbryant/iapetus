"""Smoothing module for Kalman filtering.

References:
[1] [1] Born, et al., Statistical Orbit Determination
"""

from typing import List

import numpy as np

from .kalman.data import LkfDataPoint


class SmootherError(Exception):
    pass


def Reverse(lst):
    new_lst = lst[::-1]
    return new_lst


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
            X (List[np.ndarray]): List state estimates in chronological order
            P (List[np.ndarray]): List state covariances in chronological order
            PHI (List[np.ndarray]): List state transition matrices in
                chronological order
            Q (List[np.ndarray]): List of process noise matrices (that
                encapsulate their respective process noise transition matrices)
                in chronological order

        NOTE: The first elements of input (unsorted) PHI and Q are not used.
        """
        ell = len(X)
        checks = [
            (P, "covariance matrices"),
            (PHI, "state transition matrices"),
            (Q, "process noise matrices"),
        ]
        for check in checks:
            if len(check[0]) != ell:
                msg = f"The number of {check[1]} input to smoother must be "
                msg += "equivalent to the number of input state estiamtes."
                raise SmootherError(msg)
        self.ell = ell
        self.X = Reverse(X)
        self.P = Reverse(P)
        self.PHI = Reverse(PHI)
        self.Q = Reverse(Q)
        self.X_smoothed = [self.X[0]]
        self.P_smoothed = [self.P[0]]
        self.P_intermediate = []
        self.S = []

    def __call__(self):
        """Runs smoothing algorithm."""

        ell = self.ell
        count = 0
        for x_k_given_k in self.X[1:]:
            count += 1

            # If ell = 3, k = 2 (first iteration)
            #
            # If ell = 3, k = 1 (second iteration)
            k = ell - count

            # If ell = 3, k = 2, we need P_k|k = P_2|2, but P_3|3 is first
            # value in self.P (first iteration). So we don't want self.P[0],
            # we want self.P[1], and ell - k = 3 - 2 = 1.
            #
            # If ell = 3, k =  1 (second iteration), we need P_k|k = P_1|1,
            # which is self.P[2], and ell - k = 3 - 1 = 2
            P_k_given_k = self.P[ell - k]

            # If ell = 3, k = 2, we need Phi(t_3, t_2), which is last value in
            # self.PHI (first iteration). So wee need self.PHI[0], and
            # ell - k - 1 = 3 - 2 - 1 = 0.
            #
            # If ell = 3, k = 1 (second iteration), we need Phi(t_2, t_1),
            # which is self.PHI[1], and ell - k - 1 = 3 - 1 - 1 = 1
            Phi_k_plus_1_from_k = self.PHI[ell - k - 1]

            # If ell = 3, k = 2, we need Q_2, which is self.Q[1]
            # (first iteration). ell - k = 3 - 2 = 1
            #
            # If ell = 3, k = 1 (second iteration), we need Q_1, which is
            # self.Q[2], and ell - k = 3 - 1= 2.
            Q_k = self.Q[ell - k]

            P_k_plus_1_given_k = (
                Phi_k_plus_1_from_k @ P_k_given_k @ Phi_k_plus_1_from_k.T + Q_k
            )
            S_k = (
                P_k_given_k
                @ Phi_k_plus_1_from_k.T
                @ np.linalg.inv(P_k_plus_1_given_k)
            )

            # If ell = 3, k = 2 (first iteration), we need x_k|k = x_2|2, i.e.,
            # we need self.X[1]. ell - k = 3 - 2 = 1. But, we don't have to
            # create/use an index for this variable because we're iterating
            # over its enumaration.

            # If ell = 3, k = 2 (first iteration), we need x_k+1|ell = x_3|3,
            # which we stashed in self.X_smoothed, so all we need to do is pull
            # the last element of that (and ensure we continue appending to it)
            x_k_plus_1_given_ell = self.X_smoothed[-1]

            # If ell = 3, k = 2 (first iteration), we need P_k+1|ell = P_3|3,
            # which we stashed in self.P_smoothed, so all we need to do is pull
            # the last element of that (and ensure we continue appending to it)
            P_k_plus_1_given_ell = self.P_smoothed[-1]

            x_k_given_ell = x_k_given_k + S_k @ (
                x_k_plus_1_given_ell - Phi_k_plus_1_from_k @ x_k_given_k
            )

            P_k_given_ell = (
                P_k_given_k
                + S_k @ (P_k_plus_1_given_ell - P_k_plus_1_given_k) @ S_k.T
            )

            self.P_intermediate.append(P_k_plus_1_given_k)
            self.S.append(S_k)
            self.X_smoothed.append(x_k_given_ell)
            self.P_smoothed.append(P_k_given_ell)

        return x_k_given_ell, P_k_given_ell


def lkf_smoother(data: List[LkfDataPoint]):
    """Smoother for linearize Kalman filter (LKF).

    NOTE: Assumptions:

    [1] The last LKF iteration is indexed by k, and smoother iterations are
        indexed by ell, so to start, assume L = k + 1, i.e., k = L - 1.

    [2] Smoother iterates all the way to k=0, which assumes that the LKF had an
        observation at t_k = t_0, i.e., the initial state X0 provided to the
        LKF was already propagated to the first observation time at t_0.

    [3] In general, L > k
    """
    # Resort data points by k steps in descending order
    data.sort(key=lambda x: x.k, reverse=True)

    L = data[0].k + 1
    k = L
    X_smoothed = [data[0].state_k_plus_1]
    P_smoothed = [data[0].covariance_k_plus_1]
    while k >= 0:
        k -= 1
        dp = [x for x in data if x.k == k][0]  # fwd_filter_data_point

        print(f"dp: {dp.dict()}")

        break

        # P_L-1|L-1 = P_k|k, which is covariance_k
        P_k_k = dp.covariance_k

        # Phi(t_L, t_L-1) = Phi(t_k+1, t_k), which is
        # error_state_transition_matrix_k_plus_1_k
        Phi = dp.error_state_transition_matrix_k_plus_1_k

        # Q_L-1 = Q_k which is just the Q attr
        Q = dp.Q

        # P_L|L-1 = P_k+1|k
        P_kp1_k = Phi @ P_k_k @ Phi.T + Q

        # S_L-1 = S_k
        S_k = P_k_k @ Phi.T @ np.linalg.inv(P_kp1_k)

        # X_L-1|L-1 = X_k|k
        X_k_k = dp.state_k

        # X_k+1|l is last element smoothed X stash
        X_kp1_L = X_smoothed[-1]

        # X_L-1|L = X_k|k+1
        X_k_L = X_k_k + S_k @ (X_kp1_L - Phi @ X_k_k)

        # P_k+1|L
        P_kp1_L = P_smoothed[-1]

        # P_k|L
        P_k_L = P_k_k + S_k @ (P_kp1_L - P_kp1_k) @ S_k.T

        X_smoothed.append(X_k_L)
        P_smoothed.append(P_k_L)

    return X_kp1_L, P_k_L
