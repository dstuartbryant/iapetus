"""Spring-Mass problem demonstration based on example in Section 4.8.2,
pg. 216 of Ref. [1].

This demo serves as an instructive path to developing a variety of propagators
and filter implementations, e.g., determining how to model propagators for
different types of dynamic systems.

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
from iapetus.filter.single_target.batch.batch_processor import batch_processor

# from typing import List


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

# # ----------------- Test dynamics prop ---------------------
# # Test dynamics propagation by attempting to recreate obs data on file
# # (for troubleshooting)

# X0_true = np.array([3.0, 0.0])
# T_true = np.arange(0, 11, 1)
# X_true = [X0_true]
# Z_true = []
# for t in T_true:
#     if t == 0:
#         continue
#     x, phi = propagator(0, X0_true, t)
#     X_true.append(x)

# for idx, x in enumerate(X_true):
#     z, h = H_fcn(T_true[idx], x)
#     Z_true.append(z)

# Z_error = []
# for idx, z in enumerate(Z_true):
#     z_file = np.array(data_loaded[idx][1:])
#     dz = z_file - z
#     Z_error.append(dz)
#     print(dz)


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
for i in range(4):
    xhat, xhat_fwd, P, residuals = batch_processor(
        t0, X0_i, P0_i, xbar, Z, propagator, H_fcn
    )
    xhats.append(xhat)
    X0_i += xhat
    xbar = xbar - xhat


# ----------------- Compare to Ref. [1] results -------------------------------

# After 4 iterations
X0_text = np.array([3.00019, 1.18181e-3])

dX0 = X0_text - X0_i

assert dX0[0] < 1e-6
assert dX0[0] < 1e-9

sigma_x0_text = 0.411
sigma_v0_text = 0.765
rho_text = 0.0406
resid_rho_mean_text = -4.3e-5
resid_rhodot_mean_text = -1.76e-6
resid_rho_rms = 1.16e-4
resid_rhodot_rms = 4.66e-4

sigma_x0 = np.sqrt(P[0, 0])
sigma_v0 = np.sqrt(P[1, 1])

dsigma_x0 = sigma_x0_text - sigma_x0
dsigma_v0 = sigma_v0_text - sigma_v0

assert dsigma_x0 < 1e-3
assert dsigma_v0 < 1e-3

rho = P[0, 1] / (sigma_x0 * sigma_v0)
drho = rho_text - rho

assert drho < 1e-6

rho_resid = [x[0] for x in residuals]
rhodot_resid = [x[1] for x in residuals]

rho_resid_mean = np.mean(rho_resid)
rhodot_resid_mean = np.mean(rhodot_resid)

d_rho_resid_mean = resid_rho_mean_text - rho_resid_mean
d_rhodot_resid_mean = resid_rhodot_mean_text - rhodot_resid_mean

assert d_rho_resid_mean < 1e-9
assert d_rhodot_resid_mean < 1e-9

rho_reside_rms = np.sqrt(np.mean(np.square(rho_resid)))
rhodot_reside_rms = np.sqrt(np.mean(np.square(rhodot_resid)))

d_rho_resid_rms = resid_rho_rms - rho_reside_rms
d_rhodot_resid_rms = resid_rhodot_rms - rhodot_reside_rms

assert d_rho_resid_rms < 1e-7
assert d_rhodot_resid_rms < 1e-7
