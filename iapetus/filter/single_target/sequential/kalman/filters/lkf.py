"""Linearized Kalman filter (LKF) module."""

from typing import Callable, List

import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State

from ..data import LkfDataPoint2
from ..process_noise import ProcessNoise


class LinearizedKalmanFilter:
    """Linearized Kalman filter (LKF) class.

    NOTE: Assumes that initial state is at same time as first observation.
    """

    def __init__(
        self,
        propagator: Callable,
        process_noise: ProcessNoise,
        H_fcn: Callable,
    ):
        self.propagator = propagator
        self.Q = process_noise
        self.H_fcn = H_fcn

    def predict(
        self, t_k: float, X_k_minus_1: State, xbar_k_minus_1: np.ndarray
    ):
        _, mean_k, Phi = self.propagator(
            t_k_minus_1=X_k_minus_1.timestamp,
            mean_k_minus_1=X_k_minus_1.mean,
            t_k=t_k,
        )
        xbar = Phi @ xbar_k_minus_1
        dt = t_k - X_k_minus_1.timestamp
        Q = self.Q(dt, mean_k)
        P = Phi @ X_k_minus_1.covariance.matrix @ Phi.T + Q

        return xbar, mean_k, P, Phi, Q

    def update(
        self,
        X: np.ndarray,
        P: np.ndarray,
        xbar: np.ndarray,
        z: ProbabilisticObservation,
    ):
        y, H_tilde = self.H_fcn(z.timestamp, X)
        Y = z.mean - y
        K = (
            P
            @ H_tilde.T
            @ np.linalg.inv(H_tilde @ P @ H_tilde.T + z.covariance.matrix)
        )
        xhat = xbar + K @ (Y - H_tilde @ xbar)
        P = (np.eye(len(X)) - K @ H_tilde) @ P @ (
            np.eye(len(X)) - K @ H_tilde
        ).T + K @ z.covariance.matrix @ K.T

        return xhat, P

    def __call__(
        self, X0: State, Z: ProbabilisticObservationSet
    ) -> List[LkfDataPoint2]:
        """
        Args:
            X0 (State): initial state with timestamp, mean, and covariance
                attributes
            Z (ProbabilisticObservationSet): Collection of observations

        """
        k = 0
        T = [x.timestamp for x in Z.observations]
        xbar = np.zeros((len(X0.mean),))
        Q0 = self.Q(dt_s=0, state_mean=X0.mean)

        # Process first observation
        xbar, P = self.update(
            X0.mean, X0.covariance.matrix, xbar, Z.observations[0]
        )

        X_k_minus_1 = State(timestamp=T[0], mean=X0.mean + xbar, covariance=P)
        data = [
            LkfDataPoint2(
                k=k,
                timestamp=T[0],
                predicted_state=X0.mean,
                updated_state=X0.mean + xbar,
                state_error=xbar,
                error_state_transition_matrix=np.eye(len(X0.mean)),
                covariance=X0.covariance.matrix,
                process_noise=Q0,
            )
        ]

        # Process remaining observations
        for idx in range(1, len(Z.observations)):
            k += 1

            xbar, mean_k, P, Phi, Q = self.predict(
                T[idx], X_k_minus_1=X_k_minus_1, xbar_k_minus_1=xbar
            )
            xbar, P = self.update(mean_k, P, xbar, Z.observations[idx])
            data.append(
                LkfDataPoint2(
                    k=k,
                    timestamp=T[idx],
                    state=mean_k + xbar,
                    state_error=xbar,
                    error_state_transition_matrix=Phi,
                    covariance=P,
                    process_noise=Q,
                )
            )

        return data
