"""Kalman filter propagators modules."""

from abc import ABC, abstractmethod
from typing import Tuple

import numpy as np

from iapetus.propagation.integrators import IntegrateConfig
from iapetus.propagation.propagators import AstroProp


class KalmanPropagatorError(Exception):
    pass


class KalmanPropagator(ABC):
    """Kalman filter propagator base class."""

    @abstractmethod
    def __call__(
        self, t_k_minus_1: float, mean_k_minus_1: np.ndarray, t_k: float
    ) -> Tuple[float, np.ndarray, np.ndarray]:
        """Calls to run a propagation process.

        Args:
            t_k_minus_1 (float): "Initial" timestamp associated with "initial"
                state
            mean_k_minus_1 (np.ndarray): "Initial" state vector at time
                t_k_minus_1
            t_k (float): Next subsequent timestamp to which the initial state
             is propagated; t_k can be >, =, or < t_k_minus_1

        Returns:
            t_k (float): Next subsequent timestamp to which the initial state
             was propagated
             mean_k (np.ndarray): state vector at time t_k
             Phi (np.ndarray): State transition matrix from t_k_minus_1 to t_k

        Raises:
            KalmanPropagatorError
        """
        pass

    @property
    @abstractmethod
    def iconfig(self) -> IntegrateConfig:
        pass


# class KalmanAstroPropagator(KalmanPropagator):
#     def __init__(self, init_state: dict, aprop: AstroProp, tspantol: float):
#         self.init_state_dict = init_state
#         self.aprop = aprop
#         self.state_context_manager = self.StateContextManager(init_state)
#         self.tspantol = tspantol

#     def __call__(
#         self, t_k_minus_1: float, mean_k_minus_1: np.ndarray, t_k: float
#     ) -> Tuple[float, np.ndarray, np.ndarray]:
#         tspan = [t_k_minus_1, t_k]
#         dt = t_k - t_k_minus_1
#         ui_state = self.scm.derivative_fcn_to_ui_context(mean_k_minus_1)
#         T, X, Phi = self.aprop(
#             tspan, dt, self.tspantol, ui_state=ui_state.dict()
#         )
#         if T[-1] != t_k:
#             raise KalmanPropagatorError(
#                 f"Final timestamp {T[-1]} does not match t_k {t_k}"
#             )
#         return T[-1], X[-1], Phi[-1]

#     @property
#     def StateContextManager(self):
#         return self.aprop.prop_init.state_context_manager_factory()

#     @property
#     def scm(self):
#         return self.state_context_manager


class KalmanAstroPropagator(KalmanPropagator):
    def __init__(
        self,
        init_state: dict,
        aprop: AstroProp,
        iconfig: IntegrateConfig,
    ):
        self.init_state_dict = init_state
        self.aprop = aprop
        self.state_context_manager = self.StateContextManager(init_state)
        self._iconfig = iconfig

    def __call__(
        self, t_k_minus_1: float, mean_k_minus_1: np.ndarray, t_k: float
    ) -> Tuple[float, np.ndarray, np.ndarray]:
        tspan = [t_k_minus_1, t_k]
        ui_state = self.scm.derivative_fcn_to_ui_context(mean_k_minus_1)
        T, X, Phi = self.aprop(
            tspan, self.iconfig.dt, self.iconfig.tol, ui_state=ui_state.dict()
        )
        if T[-1] != t_k:
            raise KalmanPropagatorError(
                f"Final timestamp {T[-1]} does not match t_k {t_k}"
            )
        return T[-1], X[-1], Phi[-1]

    @property
    def StateContextManager(self):
        return self.aprop.prop_init.state_context_manager_factory()

    @property
    def scm(self):
        return self.state_context_manager

    @property
    def iconfig(self) -> IntegrateConfig:
        """Returns integrator configuration model"""
        return self._iconfig
