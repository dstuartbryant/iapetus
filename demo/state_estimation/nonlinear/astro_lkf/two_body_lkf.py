"""Demo for two-body state estimation with EKF.

1. Use propagator to generate ground truth states
2. Apply Gaussian error to ground truth to form representative observations
3. Run filter

"""

import json
from copy import deepcopy
from typing import List

import numpy as np

from iapetus.data.observation import (
    ProbabilisticObservation,
    ProbabilisticObservationSet,
)
from iapetus.data.state.probabilistic import State as Pstate
from iapetus.filter.single_target.sequential.kalman.filters import (
    LinearizedKalmanFilter,
)
from iapetus.filter.single_target.sequential.kalman.process_noise import (
    AstroProcessNoise,
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
T = np.arange(0, 10)  # Used for shorter runs when degugging is needed

# --------------------- Propagate Truth -----------------
tspan = T
# dt = T[1] - T[0]
# tspantol = 1e-1
iconfig = IntegrateConfig(dt=0.1, tol=1e-2)

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
pos_errors = np.random.normal(0, sigma_pos, size=(len(T), 3))
Z = []
for idx, pe in enumerate(pos_errors):
    Z.append(Y_truth[idx][:3] + pe)

# --------------------- Initialize LKF ---------------------
# Covariance
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
kaprop = KalmanAstroPropagator(init_state=X0, aprop=aprop2, iconfig=iconfig)
init_state = Pstate(
    timestamp=T[0],
    mean=kaprop.scm.ui_input_to_derivative_fcn_context(kaprop.scm.init_state),
    covariance=P0,
)

# Add some error to initial state
init_state.mean = np.random.multivariate_normal(
    init_state.mean, init_state.covariance.matrix
)

# Inflate initial covariance a little
for i in range(6):
    init_state.covariance.matrix[i, i] = (
        init_state.covariance.matrix[i, i] * 100
    )


# Process Noise
Q_fcn = snc.QInertialRic(qr_mps2=1.2e-9, qi_mps2=2.5e-10, qc_mps2=3e-9)
Q = AstroProcessNoise(Q=Q_fcn)
Q = ZeroProcessNoise()

# Map state to observation space
H = np.hstack([np.eye(3), np.zeros((3, 3))])

# Measurement Error
R = np.diag([sigma_pos**2, sigma_pos**2, sigma_pos**2])


def H_fcn(t, mean):
    y = H @ mean
    return y, H


# Init LKF
LKF = LinearizedKalmanFilter(propagator=kaprop, process_noise=Q, H_fcn=H_fcn)

# # --------------------- Test LKF Propagator ---------------------

# t0 = 0
# mean_k_minus_1 = init_state.mean
# T_LKF = T[1:4]

# T_LKF_out, X, Phi = kaprop(t0, mean_k_minus_1, T_LKF)


# --------------------- Format observations ---------------------
time_index = -1
x_k_minus_1 = init_state
obs = []
for z in Z:
    time_index += 1
    t_k = T[time_index]
    obs.append(ProbabilisticObservation(timestamp=t_k, mean=z, covariance=R))

Z_lkf = ProbabilisticObservationSet(observations=obs)

# --------------------- Run LKF ---------------------

data = LKF(initial_state=init_state, observations=Z_lkf)

# --------------------- Run LKF Smoother ---------------------

smooth_data = lkf_smoother(data)

smooth_state_error = smooth_data[-1].state_error
x0_smooth_correction = init_state.mean + smooth_state_error

fpath = "/tmp/smooth_est.json"
with open(fpath, "w") as f:
    json.dump(x0_smooth_correction.tolist(), f, indent=4)

# # --------------------- Propagate final state estimate to t0 -----------------
# # For comparison with batch results, if needed
# xf = data[-1].updated_state
# TF = data[-1].timestamp
# XF = {s
#     "position_i_m": xf[0],
#     "position_j_m": xf[1],
#     "position_k_m": xf[2],
#     "velocity_i_mps": xf[3],
#     "velocity_j_mps": xf[4],
#     "velocity_k_mps": xf[5],
# }
# a_init3 = deepcopy(a_init)
# aprop3 = AstroProp(a_init3)
# kaprop2 = KalmanAstroPropagator(init_state=XF, aprop=aprop3, iconfig=iconfig)
# init_state = Pstate(
#     timestamp=TF,
#     mean=kaprop2.scm.ui_input_to_derivative_fcn_context(kaprop.scm.init_state),
#     covariance=data[-1].covariance,
# )
# tspan_f = [TF, T[0]]
# _, Y_f, _ = aprop(tspan_f, iconfig.dt, iconfig.tol, ui_state=XF)

# x0_filtered = Y_f[-1]

# # -------------------- Parse LKF Output ---------------------------------

# x_estimates = []
# P_estimates = []
# for d in data:
#     x_estimates.append(d.state_k_plus_1)
#     P_estimates.append(d.covariance_k_plus_1)


# # ------------------ Stash collection smoother input components ----------

# # Don't need to stash all data, just a few elements for testing purposes


# def shape_serialize_array(A: np.ndarray) -> List[float]:
#     n, m = A.shape
#     B = A.reshape(n * m)
#     return B.tolist()


# X_stash = []
# P_stash = []
# Phi_stash = []
# Q_stash = []

# num_items = len(EKF.X_stash)
# last_idx = num_items - 1
# num_to_stash = 3
# for idx in range(num_to_stash):
#     stash_idx = last_idx - (num_to_stash - 1) + idx
#     X_stash.append(EKF.X_stash[stash_idx].tolist())
#     P_stash.append(shape_serialize_array(EKF.P_stash[stash_idx]))
#     Phi_stash.append(shape_serialize_array(EKF.Phi_stash[stash_idx]))
#     Q_stash.append(shape_serialize_array(EKF.Q_stash[stash_idx]))

# fpath = "/tmp/smoother_input_data.json"
# smooth_input_data = {
#     "X_stash": X_stash,
#     "P_stash": P_stash,
#     "Phi_stash": Phi_stash,
#     "Q_stash": Q_stash,
# }

# with open(fpath, "w") as f:
#     json.dump(smooth_input_data, f, indent=4)


# # ------------------ Compute IJK error from truth ----------

# est_pos_errors = []
# standard_deviations = []
# for idx, x in enumerate(x_estimates):
#     true_pos = Y_truth[idx + 1][:3]
#     est_pos = x[:3]
#     est_pos_errors.append(true_pos - est_pos)
#     variances = np.diag(P_estimates[idx])
#     standard_deviations.append(np.sqrt(variances[:3]))

# # ------------------ Compute RIC (NTW) Error ----------

# est_ric_pos_errors = []
# ric_standard_deviations = []
# for idx, x in enumerate(x_estimates):
#     R_ECI_to_NTW = ntw_matrix(x[:3], x[3:])
#     est_ric_pos_errors.append(R_ECI_to_NTW @ est_pos_errors[idx])
#     P_eci = P_estimates[idx][:3, :3]
#     P_ric = R_ECI_to_NTW.T @ P_eci @ R_ECI_to_NTW
#     ric_standard_deviations.append(np.sqrt(np.diag(P_ric)))


# # ------------------ Save data to facilitate plot troubleshooting -------------
# fpath = "/tmp/two_body_ekf_data.json"

# data = {
#     "timestamps": T.tolist(),
#     "est_pos_errors": [x.tolist() for x in est_pos_errors],
#     "standard_deviations": [x.tolist() for x in standard_deviations],
#     "residuals": [x.tolist() for x in residuals],
#     "est_ric_pos_errors": [x.tolist() for x in est_ric_pos_errors],
#     "ric_standard_deviations": [x.tolist() for x in ric_standard_deviations],
# }

# with open(fpath, "w") as f:
#     json.dump(data, f, indent=4)

# # ---------------- Unpack Saved Data for Plotting ---------------------------
# data = json.load(open(fpath, "r"))

# T = data["timestamps"]
# est_pos_errors = data["est_pos_errors"]
# standard_deviations = data["standard_deviations"]
# residuals = data["residuals"]
# est_ric_pos_errors = data["est_ric_pos_errors"]
# ric_standard_deviations = data["ric_standard_deviations"]


# ------------------ Plot IJK Error -----------------------------


# top = ErrorStdInput(
#     x=list(T[1:]),
#     y_error=[x[0] for x in est_pos_errors],
#     y_std=[x[0] for x in standard_deviations],
#     x_axis_title="Time elapsed past epoch [s]",
#     y_axis_title="I-axis error [m]",
# )

# middle = ErrorStdInput(
#     x=list(T[1:]),
#     y_error=[x[1] for x in est_pos_errors],
#     y_std=[x[1] for x in standard_deviations],
#     x_axis_title="Time elapsed past epoch [s]",
#     y_axis_title="J-axis error [m]",
# )

# bottom = ErrorStdInput(
#     x=list(T[1:]),
#     y_error=[x[2] for x in est_pos_errors],
#     y_std=[x[2] for x in standard_deviations],
#     x_axis_title="Time elapsed past epoch [s]",
#     y_axis_title="K-axis error [m]",
# )

# fig_input = Plot3dErrorWithStdInput(
#     **{"top": top, "middle": middle, "bottom": bottom}
# )

# plotter = Plot3dErrorWithStd(upper=True, lower=True)
# figure = plotter(data=fig_input)

# # --------------- Plot RIC (NTW) Error -----------------------------


# top = ErrorStdInput(
#     x=list(T[1:]),
#     y_error=[x[0] for x in est_ric_pos_errors],
#     y_std=[x[0] for x in ric_standard_deviations],
#     x_axis_title="Time elapsed past epoch [s]",
#     y_axis_title="Radial error [m]",
# )

# middle = ErrorStdInput(
#     x=list(T[1:]),
#     y_error=[x[1] for x in est_ric_pos_errors],
#     y_std=[x[1] for x in ric_standard_deviations],
#     x_axis_title="Time elapsed past epoch [s]",
#     y_axis_title="In-Track error [m]",
# )

# bottom = ErrorStdInput(
#     x=list(T[1:]),
#     y_error=[x[2] for x in est_ric_pos_errors],
#     y_std=[x[2] for x in ric_standard_deviations],
#     x_axis_title="Time elapsed past epoch [s]",
#     y_axis_title="Cross-Track error [m]",
# )

# fig_input = Plot3dErrorWithStdInput(
#     **{"top": top, "middle": middle, "bottom": bottom}
# )

# plotter = Plot3dErrorWithStd(upper=True, lower=True)
# figure_ric = plotter(data=fig_input)
