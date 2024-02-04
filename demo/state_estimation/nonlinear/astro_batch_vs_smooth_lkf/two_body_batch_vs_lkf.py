"""Demo for two-body state estimation with batch filter."""

import json
from copy import deepcopy
from typing import List

import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State as Pstate
from iapetus.filter.single_target.batch import BatchIterator
from iapetus.filter.single_target.batch.propagators import BatchAstroPropagator
from iapetus.filter.single_target.sequential.kalman.filters import (
    LinearizedKalmanFilter,
)
from iapetus.filter.single_target.sequential.kalman.process_noise import (
    ZeroProcessNoise,
)
from iapetus.filter.single_target.sequential.kalman.propagators import (
    KalmanAstroPropagator,
)
from iapetus.filter.single_target.sequential.smoothing import lkf_smoother
from iapetus.propagation.dynamics.nonlinear.astro import snc
from iapetus.propagation.dynamics.nonlinear.astro.xforms import ntw_matrix
from iapetus.propagation.integrators import IntegrateConfig
from iapetus.propagation.propagators import AstroInit, AstroProp
from iapetus.view.error_plots import (
    ErrorStdInput,
    Plot3dErrorWithStd,
    Plot3dErrorWithStdInput,
)

np.random.seed(0)

X0 = {
    "position_i_m": 6656356.11057065,
    "position_j_m": 1700859.15707779,
    "position_k_m": 299734.38071253,
    "velocity_i_mps": -1794.25660717,
    "velocity_j_mps": 6353.55570765,
    "velocity_k_mps": 3792.38315729,
}

# T = np.arange(0, 1 * 90 * 60)
T = np.arange(0, 100)  # Used for shorter runs when degugging is needed

tspan = T

iconfig = IntegrateConfig(dt=0.1, tol=1e-1)


a_init = AstroInit(
    state_vector_content=["translational"],
    celestial_body="Earth",
    stm_flag=False,
    integrator="rk45",
    perturbations=[],
)
aprop = AstroProp(a_init)

_, Y_truth, _ = aprop(tspan, iconfig.dt, iconfig.tol, ui_state=X0)

# --------------------- Form Observations -----------------
sigma_pos = 10  # [m]
R = np.diag([sigma_pos**2, sigma_pos**2, sigma_pos**2])
pos_errors = np.random.normal(0, sigma_pos, size=(len(T), 3))
obs_data = []
for idx, pe in enumerate(pos_errors):

    obs_data.append(
        ProbabilisticObservation(
            timestamp=T[idx], mean=Y_truth[idx][:3] + pe, covariance=R
        )
    )
Z = ProbabilisticObservationSet(observations=obs_data)

# --------------------- Initialize state for propagation  -----------------
sigma_vel = sigma_pos * 1e-4
variances = [
    sigma_pos**2,
    sigma_pos**2,
    sigma_pos**2,
    sigma_vel**2,
    sigma_vel**2,
    sigma_vel**2,
]
P0 = np.diag(variances)

# Propagator and Initial State
a_init2 = deepcopy(a_init)
a_init2.stm_flag = True
aprop2 = AstroProp(a_init2)
baprop = BatchAstroPropagator(X0, aprop2)
init_state = Pstate(
    timestamp=T[0],
    mean=baprop.scm.ui_input_to_derivative_fcn_context(baprop.scm.init_state),
    covariance=P0,
)

# --------------------- Initialize batch filter  -----------------

# Add some error to initial state
init_state.mean = np.random.multivariate_normal(
    init_state.mean, init_state.covariance.matrix
)


# Inflate initial covariance a little
for i in range(6):
    init_state.covariance.matrix[i, i] = (
        init_state.covariance.matrix[i, i] * 100
    )

init_state_batch = deepcopy(init_state)

# Map state to observation space
H = np.hstack([np.eye(3), np.zeros((3, 3))])


def H_fcn(t, mean):
    # y = H @ mean
    y = mean[:3]
    return y, H


batch = BatchIterator(batch_iter_tol=1e-6, batch_max_iter=1)

# --------------------- Run batch filter  -----------------

init_xbar = np.zeros((len(init_state_batch.mean),))

X_out_batch, batch_data = batch(
    init_state=init_state_batch,
    xbar=init_xbar,
    obs=Z,
    propagator=baprop,
    H_fcn=H_fcn,
    iconfig=iconfig,
)

# --------------------- Initialize LKF ---------------------

# Propagator and Initial State
a_init3 = deepcopy(a_init)
a_init3.stm_flag = True
aprop3 = AstroProp(a_init3)
kaprop = KalmanAstroPropagator(init_state=X0, aprop=aprop3, iconfig=iconfig)
init_state_lkf = deepcopy(init_state)


# Process Noise
Q = ZeroProcessNoise()


# Init LKF
LKF = LinearizedKalmanFilter(propagator=kaprop, process_noise=Q, H_fcn=H_fcn)


# --------------------- Run LKF ---------------------

data = LKF(initial_state=init_state, observations=Z)

# --------------------- Run LKF Smoother ---------------------

smooth_data = lkf_smoother(data)

smooth_state_error = smooth_data[-1].state_error
x0_smooth_correction = init_state_lkf.mean + smooth_state_error


dx = X_out_batch.mean - x0_smooth_correction
for d in dx:
    print(d)

dxe = (batch_data.state_error - smooth_state_error).tolist()


# lkf_resids = [x.residual for x in data]
# for b, L in zip(batch_data.residuals, lkf_resids):
#     print(b - L)
