"""Propagators module.
"""

from typing import List, Optional

import numpy as np

from .dynamics.nonlinear.astro import (
    Astrodynamics,
    AstroInit,
    AstroPropIniit,
    unpack_stm,
)
from .integrators import rk4, rk45, rk78

INTEGRATORS = {"rk4": rk4.rk4, "rk45": rk45.rk45, "rk78": rk78.rk78}


class PropagatorConfigError(Exception):
    pass


def get_integrator(name: str):
    if name not in INTEGRATORS.keys():
        raise PropagatorConfigError(
            f"Unexpected integrator name found: {name}"
        )

    return INTEGRATORS[name]


class AstroProp:
    """
    User initializes this, say, as "ap".

    Then calls ap.user_defined_numerical_initial_state_model to get container
    in which to input initial state vector.

    Then use self.prop_init.convert_ui_state_to_eom_state to convert to eom
    context state.


    """

    def __init__(self, ui_config: AstroInit):
        self.dynamics = Astrodynamics(
            state_vector=ui_config.state_vector_content,
            celestial_body=ui_config.celestial_body,
            perturbations=ui_config.perturbations,
        )
        self.prop_init = AstroPropIniit(
            dynamics=self.dynamics,
            stm_flag=ui_config.stm_flag,
        )
        self.integrator = get_integrator(ui_config.integrator)
        self.StateContextManager = (
            self.prop_init.state_context_manager_factory()
        )

    def init_stm(self, y: np.ndarray):
        n = len(y)
        phi = np.eye(n)
        return np.concatenate((y, phi.reshape(n**2)), axis=0)

    def state_context_manager(self, ui_state: dict):
        """Returns an instantiated StateContextManager object.

        Args:
            ui_state (dict): dictionary of UiContext state elements
        """
        return self.StateContextManager(ui_state)

    def __call__(
        self,
        tspan: List[float],
        dt: float,
        tspantol: float,
        ui_state: Optional[dict] = None,
        state_context_manager: Optional[object] = None,
    ):
        if not ui_state and not state_context_manager:
            raise PropagatorConfigError(
                "`ui_state` and `state_context_manager` args cannot both be "
                "None"
            )
        if ui_state:
            scm = self.StateContextManager(ui_state)
        else:
            scm = state_context_manager

        # Switch state to DerivativeFunctionContext model
        y0 = scm.ui_input_to_derivative_fcn_context(scm.init_state)

        if self.prop_init.stm_flag:
            y0 = self.init_stm(y0)
        # y0, eom_state = self.redefine_user_input_numerical_state(u_state)
        derivative_fcn = self.prop_init.update_der(scm)
        T, Y = self.integrator(derivative_fcn, y0, tspan, dt, tspantol)

        Phi = []
        if self.prop_init.stm_flag:
            X = []
            for y in Y:
                x, phi = unpack_stm(y)
                X.append(x)
                Phi.append(phi)
        else:
            X = Y
        return T, X, Phi
