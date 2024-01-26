"""Kalman filter propagators modules."""

from abc import ABC, abstractmethod
from typing import List, Tuple

import numpy as np

from iapetus.propagation.propagators import AstroProp

# class KalmanPropagator(ABC):
#     """Kalman filter propagator base class."""

#     @abstractmethod
#     def __call__(
#         self, t_k_minus_1: float, mean_k_minus_1: np.ndarray, t_k: float
#     ) -> Tuple[float, np.ndarray, np.ndarray]:
#         pass


class KalmanAstroPropagatorError(Exception):
    pass


def check_dt_consistency(T: List[float]):
    dts = []
    for a, b in zip(T[1:], T[:-1]):
        dts.append(a - b)
    if len(set(dts)) > 1:
        raise KalmanAstroPropagatorError(
            "Input timestamps are not consistently incremented."
        )


class ExtendedKalmanAstroPropagator:
    def __init__(self, init_state: dict, aprop: AstroProp, tspantol: float):
        self.init_state_dict = init_state
        self.aprop = aprop
        self.state_context_manager = self.StateContextManager(init_state)
        self.tspantol = tspantol

    def __call__(
        self, t_k_minus_1: float, mean_k_minus_1: np.ndarray, t_k: float
    ) -> Tuple[float, np.ndarray, np.ndarray]:
        tspan = [t_k_minus_1, t_k]
        dt = t_k - t_k_minus_1
        ui_state = self.scm.derivative_fcn_to_ui_context(mean_k_minus_1)
        T, X, Phi = self.aprop(
            tspan, dt, self.tspantol, ui_state=ui_state.dict()
        )
        if T[-1] != t_k:
            raise KalmanAstroPropagatorError(
                f"Final timestamp {T[-1]} does not match t_k {t_k}"
            )
        return T[-1], X[-1], Phi[-1]

    @property
    def StateContextManager(self):
        return self.aprop.prop_init.state_context_manager_factory()

    @property
    def scm(self):
        return self.state_context_manager


class LinearizedKalmanAstroPropagator:
    def __init__(self, init_state: dict, aprop: AstroProp, tspantol: float):
        self.init_state_dict = init_state
        self.aprop = aprop
        self.state_context_manager = self.StateContextManager(init_state)
        self.tspantol = tspantol

    def __call__(
        self, t0: float, mean_k_minus_1: np.ndarray, T: List[float]
    ) -> Tuple[float, np.ndarray, np.ndarray]:
        if len(T) <= 1:
            raise KalmanAstroPropagatorError(
                "LKF propagation times list must have more than 1 timestamp "
                "in it."
            )
        check_dt_consistency(T)
        dt = T[1] - T[0]
        ui_state = self.scm.derivative_fcn_to_ui_context(mean_k_minus_1)
        if T[0] - t0 != dt:
            raise KalmanAstroPropagatorError(
                "Initial timestamp increment inconsistent."
            )
        T = np.insert(T, 0, t0)
        T, X, Phi = self.aprop(T, dt, self.tspantol, ui_state=ui_state.dict())

        return T[1:], X[1:], Phi[1:]

    @property
    def StateContextManager(self):
        return self.aprop.prop_init.state_context_manager_factory()

    @property
    def scm(self):
        return self.state_context_manager
