import numpy as np
from iapetus.core.state import state, trajectory

STATE_VECTORS = [[0, 0, 1, 1], [1, 1, 1, 1], [2, 2, 1, 1]]


def test_trajectory():
    states_list = []
    for idx, x in enumerate(STATE_VECTORS):
        states_list.append(state.State(idx, state.StateVector.from_array(x)))
    label = "test-trajectory-label"
    traj = trajectory.Trajectory(label, states_list)
    svs = traj.state_vectors
    sts = traj.state_timestamps

    for idx, sv in enumerate(svs):
        np.testing.assert_array_equal(sv, STATE_VECTORS[idx])
        assert idx == sts[idx]
