"""Astrodynamics state context manager and factory test module.

Need to test both the state context manager factory and its output state
context manager class for the following contexts:
1. TwoBody, No Drag, No STM
2. TwoBody, No Drag, STM
3. TwoBody, Drag, No STM
4. TwoBody, Drag, STM
5. TwoBody, Drag, No STM, Cd in state vector
6. TwoBody, Drag, STM, Cd in state vector
7. TwoBody, Drag, No STM, Bstar in state vector
8. TwoBody, Drag, STM, Bstar in state vector

"""

import numpy as np
import pytest
from pydantic.error_wrappers import ValidationError

from iapetus.propagation.dynamics.nonlinear.astro.config import (
    Astrodynamics,
    StateContextManagerFactory,
)

AC_NO_DRAG = Astrodynamics(
    **{"state_vector": ["translational"], "celestial_body": "Earth"}
)

AC_DRAG = Astrodynamics(
    **{
        "state_vector": ["translational"],
        "celestial_body": "Earth",
        "perturbations": ["atmospheric-drag"],
    }
)

AC_DRAG_CD = Astrodynamics(
    **{
        "state_vector": ["translational", "Cd"],
        "celestial_body": "Earth",
        "perturbations": ["atmospheric-drag"],
    }
)

AC_DRAG_BSTAR = Astrodynamics(
    **{
        "state_vector": ["translational", "Bstar"],
        "celestial_body": "Earth",
        "perturbations": ["atmospheric-drag"],
    }
)


def compare_ui_context(tested_object, truth_object, mode: str = None):
    ui_keys = [
        "position_i_m",
        "position_j_m",
        "position_k_m",
        "velocity_i_mps",
        "velocity_j_mps",
        "velocity_k_mps",
    ]
    if mode == "Cd":
        ui_keys.append("Cd")
    elif mode == "Bstar":
        ui_keys.append("Bstar")
    for key in ui_keys:
        assert getattr(tested_object, key) == truth_object[key]


def compare_derivative_fcn_context(
    tested_object, truth_object, mode: str = None
):
    dc_key_pairs = [
        (0, "position_i_m"),
        (1, "position_j_m"),
        (2, "position_k_m"),
        (3, "velocity_i_mps"),
        (4, "velocity_j_mps"),
        (5, "velocity_k_mps"),
    ]
    if mode == "Cd":
        dc_key_pairs.append((6, "Cd"))
    elif mode == "Bstar":
        dc_key_pairs.append((6, "Bstar"))

    for item in dc_key_pairs:
        assert tested_object[item[0]] == truth_object[item[1]]


def compare_eom_context(tested_object, truth_object, mode: str = None):
    eom_key_pairs = [
        ("pi", "position_i_m"),
        ("pj", "position_j_m"),
        ("pk", "position_k_m"),
        ("vi", "velocity_i_mps"),
        ("vj", "velocity_j_mps"),
        ("vk", "velocity_k_mps"),
    ]
    if mode == "Cd":
        eom_key_pairs.append(("Cd", "Cd"))
        eom_key_pairs.append(("A", "A_m2"))
        eom_key_pairs.append(("m", "m_kg"))
    elif mode == "Bstar":
        eom_key_pairs.append(("Bstar", "Bstar"))

    for item in eom_key_pairs:
        assert getattr(tested_object, item[0]) == truth_object[item[1]]


def test_context_1(ui_state_two_body):
    """TwoBody, No Drag, No STM"""
    scmf = StateContextManagerFactory(dynamics=AC_NO_DRAG, stm_flag=False)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyWithoutStm"

    scm_class = scmf()
    scm = scm_class(ui_state_two_body)

    compare_ui_context(scm.init_state, ui_state_two_body)

    assert scm.abbrv_init_state.pi == ui_state_two_body["position_i_m"]
    assert scm.abbrv_init_state.pj == ui_state_two_body["position_j_m"]
    assert scm.abbrv_init_state.pk == ui_state_two_body["position_k_m"]
    assert scm.abbrv_init_state.vi == ui_state_two_body["velocity_i_mps"]
    assert scm.abbrv_init_state.vj == ui_state_two_body["velocity_j_mps"]
    assert scm.abbrv_init_state.vk == ui_state_two_body["velocity_k_mps"]

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    compare_derivative_fcn_context(dc, ui_state_two_body)

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body)

    p = np.linalg.norm(
        [
            ui_state_two_body["position_i_m"],
            ui_state_two_body["position_j_m"],
            ui_state_two_body["position_k_m"],
        ]
    )

    assert ec.p == p
    assert ec.p3 == p**3
    with pytest.raises(AttributeError):
        ec.p5


def test_context_2(ui_state_two_body):
    """TwoBody, No Drag, STM"""
    scmf = StateContextManagerFactory(dynamics=AC_NO_DRAG, stm_flag=True)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyWithStm"

    scm_class = scmf()
    scm = scm_class(ui_state_two_body)

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    compare_derivative_fcn_context(dc, ui_state_two_body)

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body)

    p = np.linalg.norm(
        [
            ui_state_two_body["position_i_m"],
            ui_state_two_body["position_j_m"],
            ui_state_two_body["position_k_m"],
        ]
    )

    assert ec.p5 == p**5


def test_context_3(ui_state_two_body, ui_state_two_body_drag):
    """TwoBody, Drag, No STM"""
    scmf = StateContextManagerFactory(dynamics=AC_DRAG, stm_flag=False)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyDragWithoutStm"

    scm_class = scmf()
    with pytest.raises(ValidationError):
        scm_class(ui_state_two_body)

    scm = scm_class(ui_state_two_body_drag)

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    assert len(dc) == 6

    compare_derivative_fcn_context(dc, ui_state_two_body_drag)

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body_drag, mode="Cd")

    with pytest.raises(AttributeError):
        ec.p5


def test_context_4(ui_state_two_body_drag):
    """TwoBody, Drag, STM"""
    scmf = StateContextManagerFactory(dynamics=AC_DRAG, stm_flag=True)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyDragWithStm"

    scm_class = scmf()
    scm = scm_class(ui_state_two_body_drag)

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    compare_derivative_fcn_context(dc, ui_state_two_body_drag)

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body_drag, mode="Cd")

    p = np.linalg.norm(
        [
            ui_state_two_body_drag["position_i_m"],
            ui_state_two_body_drag["position_j_m"],
            ui_state_two_body_drag["position_k_m"],
        ]
    )

    assert ec.p5 == p**5


def test_context_5(ui_state_two_body_drag):
    """TwoBody, Drag, No STM, Cd in state vector"""
    scmf = StateContextManagerFactory(dynamics=AC_DRAG_CD, stm_flag=False)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyDragCdWithoutStm"

    scm_class = scmf()
    scm = scm_class(ui_state_two_body_drag)

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    compare_derivative_fcn_context(dc, ui_state_two_body_drag, mode="Cd")

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body_drag, mode="Cd")

    with pytest.raises(AttributeError):
        ec.p5


def test_context_6(ui_state_two_body_drag):
    """TwoBody, Drag, STM, Cd in state vector"""
    scmf = StateContextManagerFactory(dynamics=AC_DRAG_CD, stm_flag=True)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyDragCdWithStm"

    scm_class = scmf()
    scm = scm_class(ui_state_two_body_drag)

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    compare_derivative_fcn_context(dc, ui_state_two_body_drag, mode="Cd")

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body_drag, mode="Cd")

    p = np.linalg.norm(
        [
            ui_state_two_body_drag["position_i_m"],
            ui_state_two_body_drag["position_j_m"],
            ui_state_two_body_drag["position_k_m"],
        ]
    )

    assert ec.p5 == p**5


def test_context_7(ui_state_two_body_drag, ui_state_two_body_drag_Bstar):
    """TwoBody, Drag, No STM, Bstar in state vector"""
    scmf = StateContextManagerFactory(dynamics=AC_DRAG_BSTAR, stm_flag=False)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyDragBstarWithoutStm"

    scm_class = scmf()
    with pytest.raises(ValidationError):
        scm_class(ui_state_two_body_drag)

    scm = scm_class(ui_state_two_body_drag_Bstar)

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    compare_derivative_fcn_context(
        dc, ui_state_two_body_drag_Bstar, mode="Bstar"
    )

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body_drag_Bstar, mode="Bstar")

    with pytest.raises(AttributeError):
        ec.Cd

    with pytest.raises(AttributeError):
        ec.A

    with pytest.raises(AttributeError):
        ec.m

    with pytest.raises(AttributeError):
        ec.p5


def test_context_8(ui_state_two_body_drag_Bstar):
    """TwoBody, Drag, STM, Bstar in state vector"""
    scmf = StateContextManagerFactory(dynamics=AC_DRAG_BSTAR, stm_flag=True)
    eom_model = scmf.eom_context_model

    assert eom_model.__name__ == "TwoBodyDragBstarWithStm"

    scm_class = scmf()
    scm = scm_class(ui_state_two_body_drag_Bstar)

    dc = scm.ui_input_to_derivative_fcn_context(scm.init_state)

    compare_derivative_fcn_context(
        dc, ui_state_two_body_drag_Bstar, mode="Bstar"
    )

    ec = scm.derivative_fcn_to_eom_context(dc)

    compare_eom_context(ec, ui_state_two_body_drag_Bstar, mode="Bstar")

    with pytest.raises(AttributeError):
        ec.Cd

    with pytest.raises(AttributeError):
        ec.A

    with pytest.raises(AttributeError):
        ec.m

    p = np.linalg.norm(
        [
            ui_state_two_body_drag_Bstar["position_i_m"],
            ui_state_two_body_drag_Bstar["position_j_m"],
            ui_state_two_body_drag_Bstar["position_k_m"],
        ]
    )

    assert ec.p5 == p**5
