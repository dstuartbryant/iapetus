"""Astrodynamics (dynamics) configuration interface module."""

from typing import List, Optional, Tuple, Union

import numpy as np
from pydantic import BaseModel, create_model, root_validator, validator

from .constants import MU
from .eom import Eom
from .eom import configure as eom_configure
from .eom.payloads import extras
from .eom.payloads.call_states import TwoBodyDragState, TwoBodyState
from .eom.payloads.ui import IMPLEMENTED_PERTURBATION_NAMES, PERTURBATION_NAMES

UI_STATE_VECTOR_OPTIONS = ["translational", "rotational", "Bstar", "Cd"]


class AstrodynamicsConfigError(Exception):
    pass


def dynamic_field(key: str):
    if key in ["Bstar", "Cd"]:
        return {"Cd": (float, ...), "A_m2": (float, ...), "m_kg": (float, ...)}
    else:
        return {key: (float, ...)}


class Astrodynamics(BaseModel):
    """Astrodynamics user interface configuration model."""

    state_vector: List[str]
    celestial_body: str
    perturbations: Optional[List[str]] = []

    @validator("state_vector")
    def check_state_vector_input(cls, value):
        if not all([v in UI_STATE_VECTOR_OPTIONS for v in value]):
            raise AstrodynamicsConfigError(
                "Found unexpected state_vector input."
            )

        if "Bstar" in value and "translational" not in value:
            raise AstrodynamicsConfigError(
                "State vector cannot be configured with Bstar "
                "without translational dynamics."
            )

        if "Bstar" in value and "Cd" in value:
            raise AstrodynamicsConfigError(
                "State vector cannote be configured with both Bstar and Cd at "
                "the same time."
            )

        if "Cd" in value and "translational" not in value:
            raise AstrodynamicsConfigError(
                "State vector cannot be configured with drag coefficient Cd "
                "without translational dynamics."
            )

        return value

    @validator("perturbations")
    def check_perturbations_input(cls, value):
        if not all([v in PERTURBATION_NAMES for v in value]):
            raise AstrodynamicsConfigError(
                "Found unexpected perturbations model input."
            )
        return value

    @validator("perturbations")
    def check_if_perturbations_implemented(cls, value):
        bad_perts = []
        for v in value:
            if v not in IMPLEMENTED_PERTURBATION_NAMES:
                bad_perts.append(v)
        if len(bad_perts) > 0:
            raise AstrodynamicsConfigError(
                f"The following perturbations are not implemented yet: "
                f"{', '.join(bad_perts)}."
            )
        return value

    @validator("celestial_body")
    def check_if_celestial_body_in_list(cls, value):
        if value not in list(MU.keys()):
            raise AstrodynamicsConfigError(
                f"Unexpected celestial body found: {value}"
            )

        return value

    @root_validator
    def cross_check_celestial_body_and_perturbations(cls, values):
        if (
            "atmospheric-drag" in values["perturbations"]
            and values["celestial_body"] != "Earth"
        ):
            raise AstrodynamicsConfigError(
                "Non-Earth atmospheric perturbations are unavailable."
            )

        return values

    @root_validator
    def cross_check_state_vector_and_perturbations_input(cls, values):
        if (
            "Cd" in values["state_vector"]
            and "atmospheric-drag" not in values["perturbations"]
        ):
            raise AstrodynamicsConfigError(
                "Cd cannot be in state vector if atmospheric perturbations are"
                " not used."
            )
        return values

    def form_state_vector_list(self):
        sv_list = []
        if "translational" in self.state_vector:
            sv_list += [
                "position_i_m",
                "position_j_m",
                "position_k_m",
                "velocity_i_mps",
                "velocity_j_mps",
                "velocity_k_mps",
            ]
        if "rotational" in self.state_vector:
            sv_list += [
                "angle_1_rad",
                "angle_2_rad",
                "angle_3_rad",
                "angle_rate_1_radps",
                "angle_rate_2_radps",
                "angle_rate_3_radps",
            ]
        if "Cd" in self.state_vector:
            sv_list += ["Cd"]

        return sv_list

    def _state_variable_abbreviations_map(self):
        """Maps abbreviated state variable names to their longer/full length
        versions.
        """
        return {
            "position_i_m": "pi",
            "position_j_m": "pj",
            "position_k_m": "pk",
            "velocity_i_mps": "vi",
            "velocity_j_mps": "vj",
            "velocity_k_mps": "vk",
            "accelration_i_mps2": "ai",
            "accelration_j_mps2": "aj",
            "accelration_k_mps2": "ak",
            "angle_1_rad": "t1",  # t for theta
            "angle_2_rad": "t2",
            "angle_3_rad": "t3",
            "angle_rate_1_radps": "w1",  # w for omega
            "angle_rate_2_radps": "w2",
            "angle_rate_3_radps": "w3",
            "angular_acceleration_1_radps2": "wdot1",
            "angular_acceleration_2_radps2": "wdot2",
            "angular_acceleration_3_radps2": "wdot3",
            "Bstar": "Bstar",
            "Bstar_rate": "Bstar_rate",
            "Cd": "Cd",
            "Cd_rate": "Cd_rate",
            "A_m2": "A",
            "m_kg": "m",
        }

    def form_abbreviated_state_vector_list(self):
        """Forms an abbreviated state vector list that aligns with attribute
        naming in equations of motion (EOM) modules.

        The intention is to reduce the abount of code by reducing the length of
        variable names.
        """
        sv_list = self.form_state_vector_list()
        a_sv_list = []
        for item in sv_list:
            a_sv_list.append(self.state_abbrv_map[item])
        return a_sv_list

    def form_state_vector_derivative_list(self):
        sv_list = []
        if "translational" in self.state_vector:
            sv_list += [
                "velocity_i_mps",
                "velocity_j_mps",
                "velocity_k_mps",
                "accelration_i_mps2",
                "accelration_j_mps2",
                "accelration_k_mps2",
            ]
        if "rotational" in self.state_vector:
            sv_list += [
                "angle_rate_1_radps",
                "angle_rate_2_radps",
                "angle_rate_3_radps",
                "angular_acceleration_1_radps2",
                "angular_acceleration_2_radps2",
                "angular_acceleration_3_radps2",
            ]
        if "Cd" in self.state_vector:
            sv_list += ["Cd_rate"]

        return sv_list

    def form_abbreviated_state_vector_derivative_list(self):
        d_list = self.form_state_vector_derivative_list()
        a_d_list = []
        for item in d_list:
            a_d_list.append(self.state_abbrv_map[item])
        return a_d_list

    def form_partial_derivatives_matrix_map(self):
        pdm_map_list = []
        for idx, dvar in enumerate(self.abbrv_state_vector_derivative_list):
            for jdx, var in enumerate(self.abbrv_state_vector_list):
                pdm_map_list.append(("d{dvar}_d{var}", idx, jdx))

        return pdm_map_list

    @property
    def state_vector_list(self):
        return self.form_state_vector_list()

    @property
    def state_abbrv_map(self):
        return self._state_variable_abbreviations_map()

    @property
    def state_vector_derivative_list(self):
        return self.form_state_vector_derivative_list()

    @property
    def partial_derivatives_map(self):
        return self.form_partial_derivatives_matrix_map()

    @property
    def abbrv_state_vector_list(self):
        return self.form_abbreviated_state_vector_list()

    @property
    def abbrv_state_vector_derivative_list(self):
        return self.form_abbreviated_state_vector_derivative_list()


class DerivativeFunctionError(Exception):
    pass


def generate_state_only_derviative_fcn(
    x_list: List[str],
    xdot_list: List[str],
    eom: Eom,
    eom_state_model: Union[TwoBodyState, TwoBodyDragState],
    extras_model: Optional[extras.AtmosphericDragExtras] = None,
):
    def der(t, y):
        print(f"\nY: {y}\n")
        # Place incoming state y into appropriate container for EOM processing
        state_kwargs = {}
        for idx, item in enumerate(x_list):
            state_kwargs[item] = y[idx]
        if extras_model:
            for k_extra, v_extra in extras_model():
                state_kwargs[k_extra] = v_extra
        print(f"\nstate_kwargs: {state_kwargs}\n")
        y_modeled = eom_state_model(**state_kwargs)

        # EOM computations
        accels, parts = eom(y_modeled)

        # Assemble ydot output
        ydot = []
        for idx, xdot in enumerate(xdot_list):
            if xdot in x_list:
                ydot.append(y[x_list.index(xdot)])
            else:
                try:
                    ydot.append(accels(xdot))
                except Exception as exc:
                    raise DerivativeFunctionError(
                        f"Unexpected error in derivative function: "
                        f"{str(exc.args[0])}"
                    )
        return ydot

    return der


def generate_state_and_stm_derivative_fcn(
    x_list: List[str],
    xdot_list: List[str],
    partials_map: List[Tuple],
    eom: Eom,
    eom_state_model: Union[TwoBodyState, TwoBodyDragState],
    extras_model: Optional[extras.AtmosphericDragExtras] = None,
):
    def der(t, y):
        # Reshape y and stm
        state_length = len(x_list)
        y_state = y[:state_length]
        y_stm = y[state_length:]
        phi = np.zeros((len(y_state), len(y_state)))
        for idx, item in enumerate(partials_map):
            phi[item[1], item[2]] = y_stm[idx]

        # Place incoming state y into appropriate container for EOM processing
        state_kwargs = {}
        for idx, item in enumerate(x_list):
            state_kwargs[item] = y_state[idx]
        if extras_model:
            for k_extra, v_extra in extras_model():
                state_kwargs[k_extra] = v_extra
        y_modeled = eom_state_model(**state_kwargs)

        # EOM computations
        accels, parts = eom(y_modeled)

        # Form A matrix from partials
        A = np.zeros((len(y_state), len(y_state)))
        for item in partials_map:
            A[item[1], item[2]] = parts(item[0])

        # Compute stm derivative
        phi_dot = A @ phi

        # Assemble ydot output
        ydot = []
        for idx, xdot in enumerate(xdot_list):
            if xdot in x_list:
                ydot.append(y[x_list.index(xdot)])
            else:
                try:
                    ydot.append(accels(xdot))
                except Exception as exc:
                    raise DerivativeFunctionError(
                        f"Unexpected error in derivative function: "
                        f"{str(exc.args[0])}"
                    )
        ydot = np.concatenate(
            (
                ydot,
                phi_dot.reshape(
                    len(y_state) ** 2,
                ),
            ),
            axis=0,
        )

        return ydot

    return der


class PropagatorInit:
    """Produces a propagator to later be instantiated with timing
    parameters.

    This includes dynamically generating derivative functions for integrators.
    """

    def __init__(
        self,
        dynamics: Astrodynamics,
        stm_flag: bool,
    ):
        self.dynamics = dynamics
        self.stm_flag = stm_flag
        self._der = None
        self.Cd_flag = False
        self.Bstar_flag = False

    def compile_eom_init_dict(self) -> dict:
        user_config = {}
        user_config["mu"] = MU[self.dynamics.celestial_body]
        if self.stm_flag:
            user_config["partials_flag"] = True
            if "Bstar" in self.dynamics.state_vector:
                user_config["Bstar_flag"] = True
                self.Bstar_flag = True
            else:
                user_config["Bstar_flag"] = False
                self.Bstar_flag = False
            if "Cd" in self.dynamics.state_vector:
                user_config["Cd_flag"] = True
                self.Cd_flag = True
            else:
                user_config["Cd_flag"] = False
                self.Cd_flag = False
        else:
            user_config["partials_flag"] = False
            user_config["Bstar_flag"] = False
            self.Bstar_flag = False
            user_config["Cd_flag"] = False
            self.Cd_flag = False
        return user_config

    def compile_eom(self) -> Eom:
        user_config = self.compile_eom_init_dict()
        return eom_configure(
            user_config=user_config, perturbations=self.dynamics.perturbations
        )

    def select_eom_state_model(self):
        if "atmospheric-drag" in self.dynamics.perturbations:
            return TwoBodyDragState
        else:
            return TwoBodyState

    def get_extras_model(
        self, init_state: Union[TwoBodyState, TwoBodyDragState]
    ):
        if "atmospheric-drag" in self.dynamics.perturbations:
            return extras.AtmosphericDragExtras(
                init_state=init_state,
                Cd_flag=self.Cd_flag,
                Bstar_flag=self.Bstar_flag,
            )

    def compile_derivative_fcn(
        self, init_state: Union[TwoBodyState, TwoBodyDragState]
    ):
        eom = self.compile_eom()
        eom_state_model = self.select_eom_state_model()
        eom_extras_model = self.get_extras_model(init_state=init_state)
        if self.stm_flag:
            return generate_state_and_stm_derivative_fcn(
                x_list=self.dynamics.abbrv_state_vector_list,
                xdot_list=self.dynamics.abbrv_state_vector_derivative_list,
                partials_map=self.dynamics.partial_derivatives_map,
                eom=eom,
                eom_state_model=eom_state_model,
                extras_model=eom_extras_model,
            )
        else:
            return generate_state_only_derviative_fcn(
                x_list=self.dynamics.abbrv_state_vector_list,
                xdot_list=self.dynamics.abbrv_state_vector_derivative_list,
                eom=eom,
                eom_state_model=eom_state_model,
                extras_model=eom_extras_model,
            )

    def dynamic_ui_state_vector_model(self):
        """Dynamically generates a state vector model for users to input
        initial state values.

        Uses un-abbreviated state vector list so users are forced to include
        units of state vector elements to re-enforce units mindfulness.

        This context is different than the state vector from the state space
        representation model used to define state vector derivative lists and
        state transition matrix partials.

        The context here is 'all of what the user needs to input' to process
        the state space representation models.

        E.g., if atmospheric-density is used, but Cd is not in state space
        representation, the user still needs to provide Cd, A, and m for models
        to function.
        """
        config = {}
        state_vector_list = self.dynamics.state_vector_list
        if "atmospheric-drag" in self.dynamics.perturbations:
            if (
                "Cd" not in state_vector_list
                or "Bstar" not in state_vector_list
            ):
                state_vector_list.append("Cd")
        for item in state_vector_list:
            for k, v in dynamic_field(item).items():
                config[k] = v
        return create_model("DynamicUserInputStateVectorModel", **config)

    def convert_ui_state_to_eom_state(self, dynamic_ui_model):
        eom_state_model = self.select_eom_state_model()
        kwargs = {}
        for k in dynamic_ui_model.__fields__.keys():
            new_k = self.dynamics.state_abbrv_map[k]
            kwargs[new_k] = getattr(dynamic_ui_model, k)
        return eom_state_model(**kwargs)

    # @property
    # def der(self):
    #     if isinstance(self._der, type(None)):
    #         self._der = self.compile_derivative_fcn()
    #     return self._der

    def update_der(self, init_state: Union[TwoBodyState, TwoBodyDragState]):
        self._der = self.compile_derivative_fcn(init_state)
        return self._der
