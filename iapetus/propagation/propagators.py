"""Propagators module.
"""

import inspect

from .dynamics.nonlinear.astro import Astrodynamics, AstroInit, AstroPropIniit
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

    def redefine_user_input_numerical_state(self, u_state):
        eom_state_model = self.prop_init.select_eom_state_model()
        eom_state = self.prop_init.convert_ui_state_to_eom_state(u_state)

        # Form a list of state elements based on Astrodynamics definitions
        eom_model_args = [
            x
            for x in inspect.getfullargspec(eom_state_model).args
            if "self" not in x
        ]
        init_state = [getattr(eom_state, x) for x in eom_model_args]
        return init_state, eom_state

    def __call__(self, u_state, tspan, dt, tspantol):
        y0, eom_state = self.redefine_user_input_numerical_state(u_state)
        derivative_fcn = self.prop_init.update_der(eom_state)
        T, Y = self.integrator(derivative_fcn, y0, tspan, dt, tspantol)

        return T, Y

    @property
    def user_defined_numerical_initial_state_model(self):
        return self.prop_init.dynamic_ui_state_vector_model()
