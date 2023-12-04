"""Astrodynamics (dynamics) configuration interface module."""

from typing import List

from pydantic import BaseModel, root_validator, validator

UI_STATE_VECTOR_OPTIONS = ["translational", "rotational", "Cd"]
DYNAMICS_MODELS = ["two-body", "non-spherical", "atmospheric-drag"]


class AstrodynamicsConfigError(Exception):
    pass


class Astrodynamics(BaseModel):
    """Astrodynamics user interface configuration model."""

    state_vector: List[str]
    dynamics: List[str] = ["two-body"]

    @validator("state_vector")
    def check_state_vector_input(cls, value):
        if not all([v in UI_STATE_VECTOR_OPTIONS for v in value]):
            raise AstrodynamicsConfigError(
                "Found unexpected state_vector input."
            )

        if "Cd" in value and "translational" not in value:
            raise AstrodynamicsConfigError(
                "State vector cannot be configured with drag coefficient Cd "
                "without translational dynamics."
            )

        return value

    @validator("dynamics")
    def check_dynamics_input(cls, value):
        if not all([v in DYNAMICS_MODELS for v in value]):
            raise AstrodynamicsConfigError(
                "Found unexpected dynamics model input."
            )
        if "two-body" not in value:
            value = ["two-body"] + value
        return value

    @root_validator
    def cross_check_state_vector_and_dynamics_input(cls, values):
        if (
            "Cd" in values["state_vector"]
            and "atmospheric-drag" not in values["dynamics"]
        ):
            raise AstrodynamicsConfigError(
                "Cd cannot be in state vector if atmospheric dynamics are not"
                " used."
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

    def form_partial_derivatives_matrix_map(self):
        pdm_map_list = []
        for idx, dvar in enumerate(self.state_vector_derivative_list):
            for jdx, var in enumerate(self.state_vector_list):
                pdm_map_list.append((dvar, var, idx, jdx))

        return pdm_map_list

    @property
    def state_vector_list(self):
        return self.form_state_vector_list()

    @property
    def state_vector_derivative_list(self):
        return self.form_state_vector_derivative_list()

    @property
    def partial_derivatives_map(self):
        return self.form_partial_derivatives_matrix_map()
