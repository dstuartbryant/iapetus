"""Kalman filter module."""

import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State
from iapetus.propagation.dynamics.linear import LinearStateTransitionMatrix


class DiscrKalmanFilterError(Exception):
    pass


class DiscreteKalmanFilter:
    """Discrete Kalman filter class.

    Use Case: Linear systems, not nonlinear systems that have been linearized.
    """

    def __init__(
        self,
        initial_state: State,
        propagator: LinearStateTransitionMatrix,
        process_noise: np.ndarray,
        observations: ProbabilisticObservationSet,
    ):
        self.X0 = initial_state
        self.propagator = propagator
        self.Q = process_noise
        self.Z = observations
        self.H = self.Z.H
        self.Identity = np.eye(len(self.X0.mean))

    def predict(self, t_k: float, state_k_minus_1: State) -> State:
        dt = t_k - state_k_minus_1.timestamp
        state_transition_matrix = self.propagator(dt)
        mean_k = state_transition_matrix @ state_k_minus_1.mean
        covariance_k = (
            state_transition_matrix
            @ state_k_minus_1.covariance.matrix
            @ state_transition_matrix.T
            + self.Q
        )
        return State(timestamp=t_k, mean=mean_k, covariance=covariance_k)

    def update(
        self, t_k: float, z_k: ProbabilisticObservation, state_k: State
    ) -> State:
        prefit_residual = z_k.mean - self.H @ state_k.mean
        inner = (
            self.H @ state_k.covariance.matrix @ self.H.T
            + z_k.covariance.matrix
        )
        outer = state_k.covariance.matrix @ self.H.T
        if inner.shape[0] == 1:
            K_k = outer / inner
            mean_k_given_k = state_k.mean + K_k * prefit_residual
        else:
            inv_inner = np.linalg.inv(inner)
            K_k = state_k.covariance.matrix @ self.H.T @ inv_inner
            mean_k_given_k = state_k.mean + K_k @ prefit_residual

        covariance_k_given_k = (
            self.Identity - K_k @ self.H
        ) @ state_k.covariance.matrix
        return (
            State(
                timestamp=t_k,
                mean=mean_k_given_k,
                covariance=covariance_k_given_k,
            ),
            prefit_residual,
        )

    def __call__(self):
        """
        TODO: Add iterative call method and mechanisms to store intermediate
        values like residuals, predicted mean and covariance, etc.

        """
        raise NotImplementedError("See TODO in method docstring.")


class LinearizedKalmanFilter:
    """Linearized Kalman filter class.

    Use Case: Non-linear systems.
    """

    def __init__(self):
        raise NotImplementedError


class ExtendedKalmanFilter:
    """Extended Kalman filter class."""

    def __init__(
        self,
        initial_state: State,
        propagator: object,
        process_noise: np.ndarray,
        observations: ProbabilisticObservationSet,
        num_postponing_obs: int = 3,
    ):
        """
        Args:
            initial_state (State): initial mean and covariance
            propagator (object): performs integration to predict state and
                generate state transition matrices
            process_noise (np.ndarray): process noise matrix
            observations (ProbabilisticObservationSet): sequence of
                observations
            num_postponing_obs (int): Number of observations to process before
                updating the reference trajectory. This can be necessary
                especially when the observations contain significant noise.
                After a few observations have been processed, the estimages of
                state error will stabilize and the trajectory update can be
                initiated - see [Born, p. 209]
        """
        self.X0 = initial_state
        self.propagator = propagator
        self.Q = process_noise
        self.Z = observations
        self.H_fcn = self.Z.H_fcn
        self.H = None
        self.postpone_index = num_postponing_obs
        self.iteration_idx = 0

    def predict(self, t_k: float, state_k_minus_1: State) -> State:
        self.iteration_idx += 1
        self.H = self.H_fcn(t_k, state_k_minus_1.mean)
        dt = t_k - state_k_minus_1.timestamp
        mean_k, state_transition_matrix = self.propagator(
            dt, state_k_minus_1.mean
        )
        covariance_k = (
            state_transition_matrix
            @ state_k_minus_1.covariance.matrix
            @ state_transition_matrix.T
            + self.Q
        )
        return State(timestamp=t_k, mean=mean_k, covariance=covariance_k)

    def update(
        self, t_k: float, z_k: ProbabilisticObservation, state_k: State
    ) -> State:
        prefit_residual = z_k.mean - self.H @ state_k.mean

        if self.iteration_idx < self.postpone_index:
            return state_k, prefit_residual

        R = z_k.covariance.matrix
        inner = self.H @ state_k.covariance.matrix @ self.H.T + R
        outer = state_k.covariance.matrix @ self.H.T
        if inner.shape[0] == 1:
            K_k = outer / inner
        else:
            inv_inner = np.linalg.inv(inner)
            K_k = state_k.covariance.matrix @ self.H.T @ inv_inner

        delta_mean_k_given_k = K_k @ prefit_residual
        mean_k_given_k = state_k.mean + delta_mean_k_given_k
        covariance_k_given_k = (
            state_k.covariance.matrix
            - K_k @ self.H.T @ state_k.covariance.matrix
        )

        return (
            State(
                timestamp=t_k,
                mean=mean_k_given_k,
                covariance=covariance_k_given_k,
            ),
            prefit_residual,
        )
