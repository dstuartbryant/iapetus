"""Kalman filter module."""

from typing import Callable, List

import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State
from iapetus.propagation.dynamics.linear import LinearStateTransitionMatrix

from .data import LkfDataPoint
from .process_noise import ProcessNoise
from .propagators import ExtendedKalmanAstroPropagator


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

    def __init__(
        self,
        propagator: object,
        process_noise: ProcessNoise,
        H_fcn: Callable,
    ):
        """
        Args:
            propagator (object): performs integration to predict state and
                generate state transition matrices
            process_noise (np.ndarray): process noise matrix
            H_fcn (Callable): method for taking partial derivatives wrt
                "reference trajectory" in order to form the H matrix which
                maps the state space to the observation space
        """
        self.propagator = propagator
        self.Q = process_noise
        self.H_fcn = H_fcn
        self.dt = 0

    def propagate_reference_trajectory(
        self, t0: float, x0: np.ndarray, T: list[float]
    ):
        """Propagates initial state across all observation timestamps and
        returns trajectory and error state transition matrices.
        """
        _, means, stms = self.propagator(
            t0=t0,
            mean_k_minus_1=x0,
            T=T,
        )
        return means, stms

    def predict(self, idx: int, dx: np.ndarray, P: np.ndarray):
        x = self.predicted_states[idx]
        phi = self.error_state_transition_matrices[idx]
        Q = self.Q(self.dt, x)
        dx = phi @ dx
        P = phi @ P @ phi.T + Q

        return dx, x, P, phi, Q

    def update(
        self,
        dx: np.ndarray,
        x: np.ndarray,
        P: np.ndarray,
        z: ProbabilisticObservation,
    ):
        H = self.H_fcn(z.timestamp, x)
        b_tilde = z.mean - H @ x
        K = P @ H.T @ np.linalg.inv(H @ P @ H.T + z.covariance.matrix)
        dx = dx + K @ (b_tilde - H @ dx)
        P = P - K @ H @ P
        x = x + dx

        return dx, x, P, b_tilde

    def __call__(
        self, initial_state: State, observations: ProbabilisticObservationSet
    ) -> List[LkfDataPoint]:
        T = [x.timestamp for x in observations.observations]
        self.dt = T[1] - T[0]
        predicted_states, error_stms = self.propagate_reference_trajectory(
            initial_state.timestamp, initial_state.mean, T
        )
        self.predicted_states = predicted_states
        self.error_state_transition_matrices = error_stms
        data_points = []
        k = 0
        x_k = initial_state.mean
        dx_k = np.zeros((len(initial_state.mean),))
        P_k = initial_state.covariance.matrix
        for idx, z_k_plus_1 in enumerate(observations.observations):
            print(f"k: {k}")
            (
                dx_k_plus_1_given_k,
                x_k_plus_1_given_k,
                P_k_plus_1_given_k,
                Phi_k_plus_1_k,
                Q,
            ) = self.predict(idx, dx_k, P_k)

            dx_k_plus_1, x_k_plus_1, P_k_plus_1, prefit_resid = self.update(
                dx_k_plus_1_given_k,
                x_k_plus_1_given_k,
                P_k_plus_1_given_k,
                z_k_plus_1,
            )
            data_points.append(
                LkfDataPoint(
                    k=k,
                    state_k=x_k,
                    state_error_k=dx_k,
                    covariance_k=P_k,
                    error_state_transition_matrix_k_plus_1_k=Phi_k_plus_1_k,
                    state_k_plus_1_given_k=x_k_plus_1_given_k,
                    state_error_k_plus_1_given_k=dx_k_plus_1_given_k,
                    Q=Q,
                    covariance_k_plus_1_given_k=P_k_plus_1_given_k,
                    prefit_resid=prefit_resid,
                    state_k_plus_1=x_k_plus_1,
                    state_error_k_plus_1=dx_k_plus_1,
                    covariance_k_plus_1=P_k_plus_1,
                )
            )
            k += 1
            dx_k = dx_k_plus_1_given_k
            P_k = P_k_plus_1_given_k

        return data_points


class ExtendedKalmanFilter:
    """Extended Kalman filter class."""

    def __init__(
        self,
        initial_state: State,
        propagator: ExtendedKalmanAstroPropagator,
        process_noise: ProcessNoise,
        H_fcn: Callable,
        num_postponing_obs: int = 3,
    ):
        """
        Args:
            initial_state (State): initial mean and covariance
            propagator (object): performs integration to predict state and
                generate state transition matrices
            process_noise (np.ndarray): process noise matrix
            H_fcn (Callable): method for taking partial derivatives wrt
                "reference trajectory" in order to form the H matrix which
                maps the state space to the observation space
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
        self.H_fcn = H_fcn
        self.H = None
        self.postpone_index = num_postponing_obs
        self.iteration_idx = 0

        # For storing values to be used in smoothing
        self.X_stash = []
        self.P_stash = []
        self.Phi_stash = []
        self.Q_stash = []

    def predict(self, t_k: float, state_k_minus_1: State) -> State:
        self.iteration_idx += 1
        self.H = self.H_fcn(t_k, state_k_minus_1.mean)
        _, mean_k, state_transition_matrix = self.propagator(
            t_k_minus_1=state_k_minus_1.timestamp,
            mean_k_minus_1=state_k_minus_1.mean,
            t_k=t_k,
        )
        self.Phi_stash.append(state_transition_matrix)
        dt = t_k - state_k_minus_1.timestamp
        Q = self.Q(dt, mean_k)
        self.Q_stash.append(Q)
        covariance_k = (
            state_transition_matrix
            @ state_k_minus_1.covariance.matrix
            @ state_transition_matrix.T
            + Q
        )
        return State(timestamp=t_k, mean=mean_k, covariance=covariance_k)

    def update(
        self, t_k: float, z_k: ProbabilisticObservation, state_k: State
    ) -> State:
        prefit_residual = z_k.mean - self.H @ state_k.mean

        if self.iteration_idx < self.postpone_index:
            self.X_stash.append(state_k.mean)
            self.P_stash.append(state_k.covariance.matrix)
            return state_k, prefit_residual

        R = z_k.covariance.matrix
        inner = self.H @ state_k.covariance.matrix @ self.H.T + R
        outer = state_k.covariance.matrix @ self.H.T
        if inner.shape[0] == 1:
            K_k = outer / inner
        else:
            # print("SEE ME?")
            inv_inner = np.linalg.inv(inner)
            K_k = state_k.covariance.matrix @ self.H.T @ inv_inner

        # print(f"K_k: {K_k}\n")
        # print(f"self.H.T: {self.H.T}\n")

        delta_mean_k_given_k = K_k @ prefit_residual
        mean_k_given_k = state_k.mean + delta_mean_k_given_k
        covariance_k_given_k = (
            state_k.covariance.matrix
            - K_k @ self.H @ state_k.covariance.matrix
        )
        self.X_stash.append(mean_k_given_k)
        self.P_stash.append(covariance_k_given_k)

        return (
            State(
                timestamp=t_k,
                mean=mean_k_given_k,
                covariance=covariance_k_given_k,
            ),
            prefit_residual,
        )
