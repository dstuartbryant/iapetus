"""State test fixtures."""


from copy import deepcopy

import pytest

from iapetus.propagation.dynamics.nonlinear.astro.eom import state

# from iapetus.propagation.dynamics.nonlinear.astro.eom.payloads import (
#     TwoBodyDragState,
# )


PI = 6656356.11057065
PJ = 1700859.15707779
PK = 299734.38071253
VI = -1794.25660717
VJ = 6353.55570765
VK = 3792.38315729
A = (1 / 3) ** 2
M = 50
CD = 2.0
BSTAR = 0.5 * CD * A / M

two_body_input_base = {
    "pi": PI,
    "pj": PJ,
    "pk": PK,
    "vi": VI,
    "vj": VJ,
    "vk": VK,
}

two_body_drag_input = deepcopy(two_body_input_base)
two_body_drag_input["A"] = A
two_body_drag_input["Cd"] = CD
two_body_drag_input["m"] = M

two_body_drag_Bstar_input = deepcopy(two_body_input_base)
two_body_drag_Bstar_input["Bstar"] = BSTAR

ui_two_body_input_base = {
    "position_i_m": PI,
    "position_j_m": PJ,
    "position_k_m": PK,
    "velocity_i_mps": VI,
    "velocity_j_mps": VJ,
    "velocity_k_mps": VK,
}

ui_two_body_drag_input = deepcopy(ui_two_body_input_base)
ui_two_body_drag_input["A_m2"] = A
ui_two_body_drag_input["Cd"] = CD
ui_two_body_drag_input["m_kg"] = M

ui_two_body_drag_Bstar_input = deepcopy(ui_two_body_input_base)
ui_two_body_drag_Bstar_input["Bstar"] = BSTAR


@pytest.fixture
def state_two_body_without_stm():
    yield state.TwoBodyWithoutStm(**two_body_input_base)


@pytest.fixture
def state_two_body_with_stm():
    yield state.TwoBodyWithStm(**two_body_input_base)


@pytest.fixture
def state_two_body_drag_without_stm():
    yield state.TwoBodyDragWithoutStm(**two_body_drag_input)


@pytest.fixture
def state_two_body_drag_with_stm():
    yield state.TwoBodyDragWithStm(**two_body_drag_input)


@pytest.fixture
def state_two_body_drag_cd_with_stm():
    yield state.TwoBodyDragCdWithStm(**two_body_drag_input)


@pytest.fixture
def state_two_body_drag_bstar_with_stm():
    yield state.TwoBodyDragBstarWithStm(**two_body_drag_Bstar_input)


@pytest.fixture
def ui_state_two_body():
    yield ui_two_body_input_base


@pytest.fixture
def ui_state_two_body_drag():
    yield ui_two_body_drag_input


@pytest.fixture
def ui_state_two_body_drag_Bstar():
    yield ui_two_body_drag_Bstar_input
