"""Astrodynamics (dynamics) configuration interface module."""

import inspect
from typing import List, Optional, Tuple

import numpy as np
from pydantic import BaseModel, create_model, root_validator, validator

from .constants import MU
from .eom import Eom
from .eom import configure as eom_configure
from .eom import state
from .eom.payloads.ui import IMPLEMENTED_PERTURBATION_NAMES, PERTURBATION_NAMES

UI_STATE_VECTOR_OPTIONS = ["translational", "rotational", "Bstar", "Cd"]


class AstrodynamicsConfigError(Exception):
    pass


def dynamic_field(key: str):
    if key in ["Cd"]:
        return {"Cd": (float, ...), "A_m2": (float, ...), "m_kg": (float, ...)}
    else:
        return {key: (float, ...)}


def dynamic_field_just_copy(key: str):
    return {key: (float, ...)}


def get_class_arg_names(obj):
    return [x for x in inspect.getfullargspec(obj).args if "self" not in x]


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

        if "Bstar" in self.state_vector:
            sv_list += ["Bstar"]

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

        if "Bstar" in self.state_vector:
            sv_list += ["Bstar_rate"]

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
                pdm_map_list.append((f"d{dvar}_d{var}", idx, jdx))

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


class StateContextManagerFactory:
    """Dynamically creates relational methods between the different types of
    state represenation contexts and the input/output types they use.
    """

    def __init__(self, dynamics: Astrodynamics, stm_flag: bool):
        self.dynamics = dynamics
        self.stm_flag = stm_flag
        self._ui_input_init_context_model = None
        self._abbrv_ui_input_init_context_model = None
        self._eom_context_model = None
        self._ui_output_context_model = NotImplemented

    def generate_ui_state_vector_model(self):
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
                and "Bstar" not in state_vector_list
            ):
                state_vector_list.append("Cd")
        for item in state_vector_list:
            for k, v in dynamic_field(item).items():
                config[k] = v
        return create_model("DynamicUserInputStateVectorModel", **config)

    def define_abbrv_ui_state_model(self):
        """Abbreviated version of initial state model is needed for
        `derivative_fcn_to_eom_context` to function properly.
        """
        ui_model = self.ui_input_init_context_model
        args = [x for x in ui_model.__fields__.keys()]
        abbrv_args = [self.dynamics.state_abbrv_map[x] for x in args]
        config = {}
        for item in abbrv_args:
            for k, v in dynamic_field_just_copy(key=item).items():
                config[k] = v
        return create_model("AbbrvDynamicUserInputStateVectorModel", **config)

    def define_populate_abbrv_ui_state_model_method(self):
        state_abbrv_map = self.dynamics.state_abbrv_map
        abbrv_ui_input_context_model = self.abbrv_ui_input_context_model

        def populate_abbrv_ui_state_model(self, init_state):
            args = {
                k: getattr(init_state, k) for k in init_state.__fields__.keys()
            }
            abbrv_args = {}
            for k, v in args.items():
                abbrv_args[state_abbrv_map[k]] = v
            return abbrv_ui_input_context_model(**abbrv_args)

        return populate_abbrv_ui_state_model

    def define_ui_input_to_derivative_fcn_method(self):
        state_vector_list = self.dynamics.state_vector_list
        ui_input_init_context_model = self.ui_input_init_context_model

        def ui_input_to_derivative_fcn_context(
            self, s: ui_input_init_context_model
        ) -> np.ndarray:
            return np.array([getattr(s, x) for x in state_vector_list])

        return ui_input_to_derivative_fcn_context

    def select_eom_context_model(self):
        if "atmospheric-drag" not in self.dynamics.perturbations:
            if not self.stm_flag:
                return state.TwoBodyWithoutStm
            else:
                return state.TwoBodyWithStm
        else:
            if not self.stm_flag:
                if (
                    "Bstar" not in self.dynamics.state_vector
                    and "Cd" not in self.dynamics.state_vector
                ):
                    return state.TwoBodyDragWithoutStm
                elif "Bstar" in self.dynamics.state_vector:
                    return state.TwoBodyDragBstarWithoutStm
                elif "Cd" in self.dynamics.state_vector:
                    return state.TwoBodyDragCdWithoutStm
            else:
                if (
                    "Bstar" not in self.dynamics.state_vector
                    and "Cd" not in self.dynamics.state_vector
                ):
                    return state.TwoBodyDragWithStm
                elif "Bstar" in self.dynamics.state_vector:
                    return state.TwoBodyDragBstarWithStm
                elif "Cd" in self.dynamics.state_vector:
                    return state.TwoBodyDragCdWithStm

    def define_derivative_fcn_to_eom_method(self):
        """Defines a method to convert essentially an np.ndarray (vector,
        DerivativeFunctionContext) into the appropriate EomContext state used
        for EOM processing.

        Some elements of the EomContext model may not be present in the
        DerivativeFunctionContext vector, e.g., Area and mass when Cd is
        included in the state-space representation of satellite state.

        Therefore, this method requires access to the initial state definitions
        of those values.

        How this can be done? Note that the method produced by this method will
        become a method of the class that is dynamically defined by this class.
        So, we have to make sure that said dynamically defined class also
        includes an attribute that can be called to access such values not
        contained in the DerivativeFunctionContext vector.
        """
        abbrv_state_vector_list = self.dynamics.abbrv_state_vector_list
        eom_context_model = self.eom_context_model
        eom_arg_names = get_class_arg_names(eom_context_model)
        # perts = self.dynamics.perturbations
        args_not_in_deriv_fcn_context = list(
            set(eom_arg_names).difference(set(abbrv_state_vector_list))
        )
        # Build EomContextModel kwargs-index map
        # Keys are kwarg names. Values are indices of DerivativeFunctionContext
        # vector elements that map to those names.
        # NOTE: Only builds map of items that exist in both the
        # DerivativeFunctionContext and EomContext models.
        kwargs_idx_map = {
            k: abbrv_state_vector_list.index(k)
            for k in eom_arg_names
            if k in abbrv_state_vector_list
        }

        def derivative_fcn_to_eom_context(
            self, s: np.ndarray
        ) -> self.eom_context_model:
            """
            NOTE: ASSUMPTION: self.init_state in this method's context
                  assumes the class that is parent to this method has an
                  `init_state` attribute that has the values we need to fill in
                  the EomContextModel, if necessary.
            """
            kwargs = {k: s[v] for k, v in kwargs_idx_map.items()}
            for arg in args_not_in_deriv_fcn_context:
                kwargs[arg] = getattr(self.abbrv_init_state, arg)
            return eom_context_model(**kwargs)

        return derivative_fcn_to_eom_context

    def __call__(self):
        ui_init_state_type = self.ui_input_init_context_model
        method_1 = self.define_ui_input_to_derivative_fcn_method()
        method_2 = self.define_derivative_fcn_to_eom_method()
        method_3 = self.define_populate_abbrv_ui_state_model_method()

        def constructor(self, ui_init_state: dict):
            self.init_state = ui_init_state_type(**ui_init_state)
            self.abbrv_init_state = self.populate_abbrv_ui_state_model(
                self.init_state
            )

        StateContextManager = type(
            "StateContextManager",
            (object,),
            {
                "__init__": constructor,
                # Methods
                "populate_abbrv_ui_state_model": method_3,
                "ui_input_to_derivative_fcn_context": method_1,
                "derivative_fcn_to_eom_context": method_2,
            },
        )
        return StateContextManager

    @property
    def ui_input_init_context_model(self):
        if isinstance(self._ui_input_init_context_model, type(None)):
            self._ui_input_init_context_model = (
                self.generate_ui_state_vector_model()
            )
        return self._ui_input_init_context_model

    @property
    def abbrv_ui_input_context_model(self):
        if isinstance(self._abbrv_ui_input_init_context_model, type(None)):
            self._abbrv_ui_input_init_context_model = (
                self.define_abbrv_ui_state_model()
            )
        return self._abbrv_ui_input_init_context_model

    @property
    def eom_context_model(self):
        if isinstance(self._eom_context_model, type(None)):
            self._eom_context_model = self.select_eom_context_model()
        return self._eom_context_model


class DerivativeFunctionError(Exception):
    pass


def derivative_function_factory(
    x_list: List[str],
    xdot_list: List[str],
    eom: Eom,
    state_context_manager: object,
    stm_flag: bool,
    partials_map: Optional[List[Tuple]] = None,
):
    """Generates a configured derivative function for use with Runge-Kutta
    integrators.

    Args:
        x_list (List[str]): Abbreviated list of state vector elements sourced
            from an instantiated Astrodynamics class
        x_dot_list (List[str]): Abbreviated list of state vector derivative
            elements sourced from an instantiated Astrodynamics class
        eom (Eom): An instantiated equations of motion (EOM) object
        state_context_manager (StateContextManager): Instantiated state context
            manager class that handles state context switching
        stm_flag (bool): If True, indicates the need to perform state
            transition matrix (STM) computations
        partials_map (Optional[List[Tuple]]): Provided when stm_flag == True;
            used to map partial derivatives output for A matrix used in
            Phi_dot = A @ Phi computation.

    Returns:
        der (method): derivative function method for use with Runge-Kutta
            integrators

    """
    if stm_flag and not partials_map:
        raise DerivativeFunctionError(
            "Must provide `partials_map` input when `stm_flag` is set to True."
        )

    def parse_input_state_vector(y, n):
        """Used if stm_flag == True.

        Returns:
            y_state (np.ndarray): the actual state vector
            phi (np.ndarray): nxn state transition matrix, n = lenght of state
                vector
        """
        y_state = y[:n]
        y_stm = y[n:]
        phi = np.zeros((n, n))
        for idx, item in enumerate(partials_map):
            phi[item[1], item[2]] = y_stm[idx]
        return y_state, phi

    def switch_input_state_vector_to_eom_context(y):
        return state_context_manager.derivative_fcn_to_eom_context(y)

    def compute_phi_dot(phi, parts, n):
        # Form A matrix from partials
        A = np.zeros((n, n))
        for item in partials_map:
            A[item[1], item[2]] = parts(item[0])

        # Compute stm derivative
        return A @ phi

    def assemble_ydot_output(y, accels):
        ydot = []
        for idx, xdot in enumerate(xdot_list):
            if xdot in x_list:
                ydot.append(y[x_list.index(xdot)])
            else:
                try:
                    if xdot in ["Cd_rate", "Bstar_rate"]:
                        ydot.append(0)
                    else:
                        ydot.append(accels(xdot))
                except Exception as exc:
                    raise DerivativeFunctionError(
                        f"Unexpected error in derivative function: "
                        f"{str(exc.args[0])}"
                    )
        return ydot

    def reshape_and_concat_phi_dot_to_ydot(phi_dot, ydot, n):
        return np.concatenate(
            (
                ydot,
                phi_dot.reshape(
                    n**2,
                ),
            ),
            axis=0,
        )

    def der(t, y):
        """Derivative function used in Runge-Kutta integrator implementation
        that does not include state transition matrix (STM) computations.

        Args:
            t (float): time in seconds
            y (np.ndarray): state vector

        Returns:
            ydot (np.ndarray): state vector derivative
        """
        y_modeled = switch_input_state_vector_to_eom_context(y)

        # EOM computations
        accels, _ = eom(y_modeled)

        return assemble_ydot_output(y, accels)

    def der_stm(t, y):
        """Derivative function used in Runge-Kutta integrator implementation
        that includes state transition matrix (STM) computations.

        Args:
            t (float): time in seconds
            y (np.ndarray): state vector

        Returns:
            ydot (np.ndarray): state vector derivative
        """
        n = len(x_list)
        y, phi = parse_input_state_vector(y, n)

        y_modeled = switch_input_state_vector_to_eom_context(y)

        # EOM computations
        accels, parts = eom(y_modeled)

        phi_dot = compute_phi_dot(phi, parts, n)

        ydot = assemble_ydot_output(y, accels)

        return reshape_and_concat_phi_dot_to_ydot(phi_dot, ydot, n)

    if stm_flag:
        return der_stm
    else:
        return der


class PropagatorInit:
    """Produces a propagator to later be instantiated with timing
    parameters and an initial state.

    This includes dynamically generating derivative functions for integrators.
    """

    def __init__(
        self,
        dynamics: Astrodynamics,
        stm_flag: bool,
    ):
        self.dynamics = dynamics
        self.stm_flag = stm_flag
        self.state_context_manager_factory = StateContextManagerFactory(
            dynamics=dynamics, stm_flag=stm_flag
        )
        self._der = None

    def compile_eom_init_dict(self) -> dict:
        user_config = {}
        user_config["mu"] = MU[self.dynamics.celestial_body]
        if self.stm_flag:
            user_config["partials_flag"] = True
            if "Bstar" in self.dynamics.state_vector:
                user_config["Bstar_flag"] = True
            else:
                user_config["Bstar_flag"] = False
            if "Cd" in self.dynamics.state_vector:
                user_config["Cd_flag"] = True
            else:
                user_config["Cd_flag"] = False
        else:
            user_config["partials_flag"] = False
            user_config["Bstar_flag"] = False
            user_config["Cd_flag"] = False
        return user_config

    def compile_eom(self) -> Eom:
        user_config = self.compile_eom_init_dict()
        return eom_configure(
            user_config=user_config, perturbations=self.dynamics.perturbations
        )

    def compile_derivative_fcn(
        self,
        state_context_manager: object,
    ):
        eom = self.compile_eom()
        der_factory_kwargs = {
            "x_list": self.dynamics.abbrv_state_vector_list,
            "xdot_list": self.dynamics.abbrv_state_vector_derivative_list,
            "eom": eom,
            "state_context_manager": state_context_manager,
            "stm_flag": False,
        }
        if self.stm_flag:
            der_factory_kwargs["stm_flag"] = True
            der_factory_kwargs[
                "partials_map"
            ] = self.dynamics.partial_derivatives_map

        return derivative_function_factory(**der_factory_kwargs)

    def update_der(self, state_context_manager: object):
        self._der = self.compile_derivative_fcn(state_context_manager)
        return self._der
