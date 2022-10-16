import numpy as np
from iapetus.core.linalg import process_noise


def test_continuous_white_noise_acceleration():
    q = 1
    dt = 1
    Q = process_noise.continuous_white_noise_acceleration(q, dt)

    Q_truth = np.array(
        [
            [1 / 3, 0, 1 / 2, 0],
            [0, 1 / 3, 0, 1 / 2],
            [1 / 2, 0, 1, 0],
            [0, 1 / 2, 0, 1],
        ]
    )
    np.testing.assert_array_equal(Q, Q_truth)
