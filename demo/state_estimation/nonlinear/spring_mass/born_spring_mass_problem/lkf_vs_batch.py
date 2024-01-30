"""Compares Linearized Kalman filter (LKF) and batch processing results using
example in Section 4.8.2, pg. 216 of Ref. [1].

References:
[1] Born, Stat OD text.

"""

import json
import math
from copy import deepcopy
from os import path
from typing import List

import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State as Pstate
from iapetus.filter.single_target.batch.batch_processor import batch_processor
from iapetus.filter.single_target.sequential.kalman.filters import (
    LinearizedKalmanFilter,
)
from iapetus.filter.single_target.sequential.kalman.process_noise import (
    ZeroProcessNoise,
)
from iapetus.filter.single_target.sequential.smoothing import lkf_smoother
from iapetus.propagation.integrators.rk45 import rk45

CURR_DIR = path.dirname(path.abspath(__file__))
OBS_PATH = path.join(CURR_DIR, "obs_data.json")


def propagator(t0: float, X0: np.ndarray, t: float):
    """Spring mass propagator based on Eq. 4.8.5 from Ref. [1].

    ASSUMPTION: t0 = 0.
    """
    if t0 != 0:
        raise ValueError("t0 should be zero.")

    k1 = 2.5
    k2 = 3.7
    m = 1.5
    w = np.sqrt((k1 + k2) / m)

    x0 = X0[0]
    v0 = X0[1]
    x = x0 * math.cos(w * t) + v0 / w * math.sin(w * t)
    v = v0 * math.cos(w * t) - x0 * w * math.sin(w * t)

    Phi = np.array(
        [
            [math.cos(w * t), 1 / w * math.sin(w * t)],
            [-w * math.sin(w * t), math.cos(w * t)],
        ]
    )

    return np.array([x, v]), Phi


def add_stm_to_state_vector(y, phi):
    return np.concatenate(
        (
            y,
            phi.reshape(
                2**2,
            ),
        ),
        axis=0,
    )


def remove_stm_from_state_vector(y):
    x = y[:2]
    phi = y[2:].reshape((2, 2))
    return x, phi


def deriv_fcn(t, y):
    """Derivative function for sequential propagator integration."""

    k1 = 2.5
    k2 = 3.7
    m = 1.5
    w = np.sqrt((k1 + k2) / m)

    y, phi = remove_stm_from_state_vector(y)

    A = np.array([[0, 1], [-(w**2), 0]])
    ydot = A @ y

    phi_dot = A @ phi
    ydot = add_stm_to_state_vector(ydot, phi_dot)

    return ydot


def sequential_propagator(t0: float, X0: np.ndarray, t: float):
    phi = np.eye(2)
    tspan = [t0, t]
    dt = t - t0
    X0 = add_stm_to_state_vector(X0, phi)
    _, Y = rk45(deriv_fcn, X0, tspan, dt, 1e-3)
    x, phi = remove_stm_from_state_vector(Y[-1])
    return x, phi


def propagator_wrapper(t0: float, mean_k_minus_1: np.ndarray, T: List[float]):
    """Wrapper for LKF implementation."""
    if t0 != T[0]:
        raise ValueError("First time stamp in T should match t0")

    states = [mean_k_minus_1]
    stms = [np.eye(2)]
    X_k = mean_k_minus_1
    for t_kp1, t_k in zip(T[1:], T[:-1]):
        x, phi = sequential_propagator(t_k, X_k, t_kp1)
        states.append(x)
        stms.append(phi)
        X_k = x

    return T, states, stms


def H_fcn(t: float, X: np.ndarray) -> np.ndarray:
    """Observation matrix based on Eqs. 4.8.4 and H_tilde function (just above
    Eq. 4.8.5) in Ref. [1].
    """
    h = 5.4
    x = X[0]
    v = X[1]
    rho = np.sqrt(x**2 + h**2)
    rho_dot = x * v / rho

    H_tilde = np.array(
        [[x / rho, 0], [(v / rho - x**2 * v / rho**3), x / rho]]
    )
    return np.array([rho, rho_dot]), H_tilde


# ----------------- Read in obs data ---------------------
data_loaded = json.load(open(OBS_PATH, "r"))
data = []
R = np.eye(2)
for x in data_loaded:
    data.append(
        ProbabilisticObservation(
            timestamp=x[0], mean=np.array([x[1], x[2]]), covariance=R
        )
    )
Z = ProbabilisticObservationSet(observations=data)

# ----------------- Form initial state ---------------------
t0 = 0
X0 = np.array([4.0, 0.2])
P0 = np.array([[1000, 0], [0, 100]])
xbar0 = np.array([0, 0])


# ----------------- Run Batch Filter ---------------------
X0_i = deepcopy(X0)
P0_i = deepcopy(P0)
xbar = deepcopy(xbar0)
xhats = []
for i in range(1):
    (
        xhat,
        xhat_fwd,
        P,
        batch_residuals,
        batch_predicted_states,
        batch_stms,
    ) = batch_processor(t0, X0_i, P0_i, xbar, Z, propagator, H_fcn)
    xhats.append(xhat)
    X0_i += xhat
    xbar = xbar - xhat

propagated_batch_state, _ = propagator(t0, X0_i, Z.observations[-1].timestamp)

# ----------------- Initilize LKF ---------------------


# lkf_T, lkf_predicted_states, lkf_stms = propagator_wrapper(
#     t0, X0, [z.timestamp for z in Z.observations]
# )


# def lkf_predict(idx, dx0, P0):
#     phi = lkf_stms[idx]
#     dx = phi @ dx0
#     P = phi @ P0 @ phi.T

#     return dx, P


# def lkf_update(idx, dxi, Pi, zk):
#     xi = lkf_predicted_states[idx]
#     y, H = H_fcn(zk.timestamp, xi)
#     b = zk.mean - y
#     K = Pi @ H.T @ np.linalg.inv(H @ Pi @ H.T + R)
#     dxk = dxi + K @ (b - H @ dxi)
#     Pk = (np.eye(2) - K @ H) @ Pi
#     xk = xi + dxk
#     return dxk, xk, Pk, b


# dxk = np.array([0, 0])
# Pk = deepcopy(P0_i)
# lkf_residuals = []
# for idx in range(len(Z.observations)):
#     dxi, Pi = lkf_predict(idx, dxk, Pk)
#     dxk, xk, Pk, bk = lkf_update(idx, dxi, Pi, Z.observations[idx])

#     lkf_residuals.append(bk)


LKF = LinearizedKalmanFilter(
    propagator=propagator_wrapper,
    process_noise=ZeroProcessNoise(),
    H_fcn=H_fcn,
)

X_init_LKF = Pstate(timestamp=0, mean=X0, covariance=P0)

# ----------------- Run LKF ---------------------

LKF_data = LKF(initial_state=X_init_LKF, observations=Z)

lkf_ref_traj = [x.state_k_plus_1_given_k for x in LKF_data]


# last_lkf_mean = LKF_data[-1].state_k_plus_1


# # ----------------- Run LKF smoother ---------------------
X0_smooth, P0_smooth = lkf_smoother(LKF_data)


# ----------------- Smoother troubleshooting ------------------

"""

Observations at: t_0, t_1, t_2

L = 2

X_L|L is updated estimate at t_L = t_2

If k = L - 1 = 1,
* X_k|L = Phi(t_k, t_L) X_L_L


"""

data = deepcopy(LKF_data)
data.sort(key=lambda x: x.k, reverse=True)
L = data[0].k + 1
k = L - 1
dp = [x for x in data if x.k == k][0]


# P_k_k = dp.covariance_k

# # Phi(t_L, t_L-1) = Phi(t_k+1, t_k), which is
# # error_state_transition_matrix_k_plus_1_k
# Phi = dp.error_state_transition_matrix_k_plus_1_k

# # Q_L-1 = Q_k which is just the Q attr
# Q = dp.Q

# # P_L|L-1 = P_k+1|k

# P_kp1_k = Phi @ P_k_k @ Phi.T + Q

# # S_L-1 = S_k
# S_k = P_k_k @ Phi.T @ np.linalg.inv(P_kp1_k)

# # X_L-1|L-1 = X_k|k
# X_k_k = dp.state_k

# # X_k+1|l is last element smoothed X stash
# X_kp1_L = dp.state_k_plus_1_given_k

# # X_L-1|L = X_k|k+1
# X_k_L = X_k_k + S_k @ (X_kp1_L - Phi @ X_k_k)

# # # P_k+1|L
# # P_kp1_L = P_smoothed[-1]

# # # P_k|L
# # P_k_L = P_k_k + S_k @ (P_kp1_L - P_kp1_k) @ S_k.T
