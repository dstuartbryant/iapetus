import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State
from iapetus.filter.single_target.sequential.kalman import DiscreteKalmanFilter
from iapetus.propagation.dynamics.linear import StmStraightLine

propagator = StmStraightLine(1)
Q = np.array([[1, 0, 0], [0, 0.01, 0], [0, 0, 0.01]])

observations = [
    ProbabilisticObservation(5, np.array([385]), np.array([25])),
    ProbabilisticObservation(10, np.array([234]), np.array([25])),
    ProbabilisticObservation(15, np.array([85]), np.array([25])),
    ProbabilisticObservation(20, np.array([2]), np.array([25])),
]


H = np.array([[1, 0, 0]])
observation_set = ProbabilisticObservationSet(observations, H)

X0 = np.array([[600], [-45], [1.5]])
P0 = np.array([[225, 0, 0], [0, 25, 0], [0, 0, 1]])
initial_state = State(0, X0, P0)
dkf = DiscreteKalmanFilter(initial_state, propagator, Q, observation_set)


residuals = []
position_stds = []
velocity_stds = []
acceleration_stds = []

for idx in range(0, len(observations)):
    tk = observations[idx].timestamp
    if idx == 0:
        state_k_minus_1 = initial_state
    state_k = dkf.predict(tk, state_k_minus_1)
    state_k_given_k, resid = dkf.update(tk, observations[idx], state_k)

    residuals.append(list(resid)[0][0])
    position_stds.append(np.sqrt(state_k_given_k.covariance.matrix[0][0]))
    velocity_stds.append(np.sqrt(state_k_given_k.covariance.matrix[1][1]))
    acceleration_stds.append(np.sqrt(state_k_given_k.covariance.matrix[2][2]))

    state_k_minus_1 = state_k_given_k


print(f"Residuals: {residuals}")
print(f"Position sigmas: {position_stds}")
print(f"Velocity sigmas: {velocity_stds}")
print(f"Acceleration sigmas: {acceleration_stds}")
