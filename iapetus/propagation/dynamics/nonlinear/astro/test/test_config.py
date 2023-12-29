"""Astrodynamics configuration test module."""

import pytest

from iapetus.propagation.dynamics.nonlinear.astro.config import (
    Astrodynamics,
    AstrodynamicsConfigError,
)


def test_state_vector_validation():
    Astrodynamics(**{"state_vector": ["translational"]})
    Astrodynamics(**{"state_vector": ["rotational"]})
    Astrodynamics(
        **{
            "state_vector": ["translational", "Cd"],
            "perturbations": ["atmospheric-drag"],
        }
    )
    Astrodynamics(
        **{
            "state_vector": ["translational", "rotational", "Cd"],
            "perturbations": ["atmospheric-drag"],
        }
    )

    with pytest.raises(AstrodynamicsConfigError) as excinfo:
        Astrodynamics(**{"state_vector": ["warp"]})

    assert str(excinfo.value) == "Found unexpected state_vector input."

    with pytest.raises(AstrodynamicsConfigError) as excinfo:
        Astrodynamics(**{"state_vector": ["translational", "warp"]})

    assert str(excinfo.value) == "Found unexpected state_vector input."

    with pytest.raises(AstrodynamicsConfigError) as excinfo:
        Astrodynamics(**{"state_vector": ["rotational", "Cd"]})

    assert (
        str(excinfo.value)
        == "State vector cannot be configured with drag coefficient Cd without"
        " translational dynamics."
    )


def test_dynamics_validation():
    ac = Astrodynamics(
        **{
            "state_vector": ["translational"],
            "perturbations": ["non-spherical"],
        }
    )

    assert "non-spherical" in ac.perturbations


def test_cross_check_validation():
    with pytest.raises(AstrodynamicsConfigError) as excinfo:
        Astrodynamics(
            **{"state_vector": ["translational", "rotational", "Cd"]}
        )

    assert (
        str(excinfo.value)
        == "Cd cannot be in state vector if atmospheric perturbations are not used."
    )


def test_state_vector_formation():
    ac = Astrodynamics(**{"state_vector": ["translational"]})
    assert ac.state_vector_list == [
        "position_i_m",
        "position_j_m",
        "position_k_m",
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
    ]

    ac = Astrodynamics(**{"state_vector": ["translational", "rotational"]})
    assert ac.state_vector_list == [
        "position_i_m",
        "position_j_m",
        "position_k_m",
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
        "angle_1_rad",
        "angle_2_rad",
        "angle_3_rad",
        "angle_rate_1_radps",
        "angle_rate_2_radps",
        "angle_rate_3_radps",
    ]

    ac = Astrodynamics(
        **{
            "state_vector": ["translational", "Cd"],
            "perturbations": ["atmospheric-drag"],
        }
    )
    assert ac.state_vector_list == [
        "position_i_m",
        "position_j_m",
        "position_k_m",
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
        "Cd",
    ]

    ac = Astrodynamics(
        **{
            "state_vector": ["translational", "rotational", "Cd"],
            "perturbations": ["atmospheric-drag"],
        }
    )
    assert ac.state_vector_list == [
        "position_i_m",
        "position_j_m",
        "position_k_m",
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
        "angle_1_rad",
        "angle_2_rad",
        "angle_3_rad",
        "angle_rate_1_radps",
        "angle_rate_2_radps",
        "angle_rate_3_radps",
        "Cd",
    ]


def test_state_vector_derivative_formation():
    ac = Astrodynamics(**{"state_vector": ["translational"]})
    assert ac.state_vector_derivative_list == [
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
        "accelration_i_mps2",
        "accelration_j_mps2",
        "accelration_k_mps2",
    ]

    ac = Astrodynamics(**{"state_vector": ["translational", "rotational"]})
    assert ac.state_vector_derivative_list == [
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
        "accelration_i_mps2",
        "accelration_j_mps2",
        "accelration_k_mps2",
        "angle_rate_1_radps",
        "angle_rate_2_radps",
        "angle_rate_3_radps",
        "angular_acceleration_1_radps2",
        "angular_acceleration_2_radps2",
        "angular_acceleration_3_radps2",
    ]

    ac = Astrodynamics(
        **{
            "state_vector": ["translational", "Cd"],
            "perturbations": ["atmospheric-drag"],
        }
    )
    assert ac.state_vector_derivative_list == [
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
        "accelration_i_mps2",
        "accelration_j_mps2",
        "accelration_k_mps2",
        "Cd_rate",
    ]

    ac = Astrodynamics(
        **{
            "state_vector": ["translational", "rotational", "Cd"],
            "perturbations": ["atmospheric-drag"],
        }
    )
    assert ac.state_vector_derivative_list == [
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
        "accelration_i_mps2",
        "accelration_j_mps2",
        "accelration_k_mps2",
        "angle_rate_1_radps",
        "angle_rate_2_radps",
        "angle_rate_3_radps",
        "angular_acceleration_1_radps2",
        "angular_acceleration_2_radps2",
        "angular_acceleration_3_radps2",
        "Cd_rate",
    ]
