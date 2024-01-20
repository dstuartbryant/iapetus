"""Smoother test module."""

import numpy as np

from iapetus.filter.single_target.sequential.smoothing import Reverse, Smoother


def test_smoother(smooth_input):
    d = smooth_input
    print(f"X: {d['X']}")
    smoother = Smoother(d["X"], d["P"], d["PHI"], d["Q"])

    _x_1_given_3, _P_1_given_3 = smoother()

    X = Reverse(d["X"])
    P = Reverse(d["P"])
    PHI = Reverse(d["PHI"])
    Q = Reverse(d["Q"])

    P_3_given_2 = PHI[0] @ P[1] @ PHI[0].T + Q[1]
    assert np.array_equal(P_3_given_2, smoother.P_intermediate[0])

    S_2 = P[1] @ PHI[0].T @ np.linalg.inv(P_3_given_2)
    assert np.array_equal(S_2, smoother.S[0])

    x_2_given_3 = X[1] + S_2 @ (X[0] - PHI[0] @ X[1])
    assert np.array_equal(x_2_given_3, smoother.X_smoothed[1])

    P_2_given_3 = P[1] + S_2 @ (P[0] - P_3_given_2) @ S_2.T
    assert np.array_equal(P_2_given_3, smoother.P_smoothed[1])

    P_2_given_1 = PHI[1] @ P[2] @ PHI[1].T + Q[2]
    assert np.array_equal(P_2_given_1, smoother.P_intermediate[1])

    S_1 = P[2] @ PHI[1].T @ np.linalg.inv(P_2_given_1)
    assert np.array_equal(S_1, smoother.S[1])

    x_1_given_3 = X[2] + S_1 @ (x_2_given_3 - PHI[1] @ X[2])
    assert np.array_equal(x_1_given_3, _x_1_given_3)

    P_1_given_3 = P[2] + S_1 @ (P_2_given_3 - P_2_given_1) @ S_1.T
    assert np.array_equal(P_1_given_3, _P_1_given_3)
