"""Probabilistic state test module."""

import numpy as np

from iapetus.data.state.probabilistic import State


def test_state():
    m = np.array([1, 2, 3, 4])
    P = np.array([[1, 2, 2, 2], [2, 4, 2, 2], [2, 2, 9, 2], [2, 2, 2, 16]])
    s = State(mean=m, covariance=P)

    assert s.dimension == 4
    assert np.array_equal(m, s.mean)
    assert np.array_equal(
        s.covariance.standard_deviations, np.array([1, 2, 3, 4])
    )
