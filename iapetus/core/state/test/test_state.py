from types import NoneType

import numpy as np
from iapetus.core.linalg import Vector
from iapetus.core.state import state

X = [1, 2, 3, 4]
X_acc = X + [5, 6]


def test_state_vector():
    sv = state.StateVector.from_array(X)
    assert isinstance(sv.pos, Vector)
    assert isinstance(sv.vel, Vector)
    assert isinstance(sv.acc, NoneType)
    np.testing.assert_array_equal(sv.vector, np.array(X))

    sv = sv = state.StateVector.from_array(X_acc)
    assert isinstance(sv.acc, Vector)
    np.testing.assert_array_equal(sv.vector, np.array(X_acc))


def test_state():
    sv = state.StateVector.from_array(X)
    y = state.State(timestamp=1, vector_obj=sv)
    assert y.time == 1
    np.testing.assert_array_equal(sv.vector, np.array(X))
