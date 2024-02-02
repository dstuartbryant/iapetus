"""Compares Linearized Kalman filter (LKF) and batch processing results using
example in Section 4.8.2, pg. 216 of Ref. [1].

References:
[1] Born, Stat OD text.

"""

import json
import math
from copy import deepcopy
from os import path

import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State as Pstate
from iapetus.filter.single_target.batch import batch_processor
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


def sequential_propagator(
    t_k_minus_1: float, mean_k_minus_1: Pstate, t_k: float
):
    phi = np.eye(2)
    tspan = [t_k_minus_1, t_k]
    dt = t_k - t_k_minus_1
    X0 = add_stm_to_state_vector(mean_k_minus_1, phi)
    tout, Y = rk45(deriv_fcn, X0, tspan, dt, 1e-3)
    x, phi = remove_stm_from_state_vector(Y[-1])
    return tout, x, phi


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

LKF = LinearizedKalmanFilter(
    propagator=sequential_propagator,
    process_noise=ZeroProcessNoise(),
    H_fcn=H_fcn,
)

X_init_LKF = Pstate(timestamp=0, mean=X0, covariance=P0)

# ----------------- Run LKF ---------------------

LKF_data = LKF(initial_state=X_init_LKF, observations=Z)

# lkf_ref_traj = [x.state_k_plus_1_given_k for x in LKF_data]
lkf_ref_traj = [x.predicted_state for x in LKF_data]


# last_lkf_mean = LKF_data[-1].state_k_plus_1

last_lkf_mean_1 = LKF_data[-1].predicted_state + LKF_data[-1].state_error
last_lkf_mean_2 = LKF_data[-1].updated_state


# ----------------- Run Smoother  ---------------------

smooth_data = lkf_smoother(LKF_data)

smooth_state_error = smooth_data[-1].state_error
x0_smooth_correction = X0 + smooth_state_error
