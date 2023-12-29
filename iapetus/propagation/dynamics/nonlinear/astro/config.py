"""Astrodynamics (dynamics) configuration interface module."""

from typing import List, Optional

from pydantic import BaseModel, create_model, root_validator, validator

from .eom.payloads.ui import PERTURBATION_NAMES

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
                pdm_map_list.append((f"d{dvar}_d{var}", idx, jdx))

        return pdm_map_list

    def dynamic_ui_state_vector_model(self):
        """Dynamically generates a state vector model for users to input
        initial state values.

        Uses un-abbreviated state vector list so users are forced to include
        units of state vector elements to re-enforce units mindfulness.
        """
        config = {}
        for item in self.state_vector_list:
            for k, v in dynamic_field(item).items():
                config[k] = v
        return create_model("DynamicUserInputStateVectorModel", **config)

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
