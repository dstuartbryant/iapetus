"""Kalman filter propagators modules."""

from abc import ABC, abstractmethod
from typing import Tuple

import numpy as np

from iapetus.propagation.propagators import AstroProp


class KalmanPropagator(ABC):
    """Kalman filter propagator base class."""

    @abstractmethod
    def __call__(
        self, t_k_minus_1: float, mean_k_minus_1: np.ndarray, t_k: float
    ) -> Tuple[float, np.ndarray, np.ndarray]:
        pass


class KalmanAstroPropagatorError(Exception):
    pass


class KalmanAstroPropagator(KalmanPropagator):
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
