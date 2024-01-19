"""Demo for two-body state estimation with EKF.

1. Use propagator to generate ground truth states
2. Apply Gaussian error to ground truth to form representative observations
3. Run filter

"""

import json

import numpy as np

from iapetus.data.observation import ProbabilisticObservation
from iapetus.data.state.probabilistic import State as Pstate
from iapetus.filter.single_target.sequential.kalman.filters import (
    ExtendedKalmanFilter,
)
from iapetus.filter.single_target.sequential.kalman.process_noise import (
    AstroProcessNoise,
)
from iapetus.filter.single_target.sequential.kalman.propagators import (
    KalmanAstroPropagator,
)
from iapetus.propagation.dynamics.nonlinear.astro import snc
from iapetus.propagation.dynamics.nonlinear.astro.xforms import ntw_matrix
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

T = np.arange(0, 1 * 90 * 60)
# T = np.arange(0, 10)  # Used for shorter runs when degugging is needed

# --------------------- Propagate Truth -----------------
tspan = T
dt = T[1] - T[0]
tspantol = 1e-1

a_init = AstroInit(
    state_vector_content=["translational"],
    celestial_body="Earth",
    stm_flag=False,
    integrator="rk45",
    perturbations=[],
)
aprop = AstroProp(a_init)

_, Y_truth, _ = aprop(tspan, dt, tspantol, ui_state=X0)

# --------------------- Form Observations -----------------
sigma_pos = 500  # [m]
pos_errors = np.random.normal(0, sigma_pos, size=(len(T), 3))
Z = []
for idx, pe in enumerate(pos_errors):
    if idx == 0:
        # First state does not need an observation
        continue
    Z.append(Y_truth[idx][:3] + pe)

# --------------------- Initialize EKF ---------------------
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
a_init2 = a_init
a_init2.stm_flag = True
aprop2 = AstroProp(a_init2)
kaprop = KalmanAstroPropagator(X0, aprop2, tspantol)
init_state = Pstate(
    timestamp=T[0],
    mean=kaprop.scm.ui_input_to_derivative_fcn_context(kaprop.scm.init_state),
    covariance=P0,
)

# Add some error to initial state
init_state.mean = np.random.multivariate_normal(
    init_state.mean, init_state.covariance.matrix
)


# Process Noise
Q_fcn = snc.QInertialRic(qr_mps2=1.2e-9, qi_mps2=2.5e-10, qc_mps2=3e-9)
Q = AstroProcessNoise(Q=Q_fcn)

# Map state to observation space
H = np.hstack([np.eye(3), np.zeros((3, 3))])

# Measurement Error
R = np.diag([sigma_pos**2, sigma_pos**2, sigma_pos**2])


def H_fcn(t, mean):
    return H


# Init EKF
EKF = ExtendedKalmanFilter(
    initial_state=init_state, propagator=kaprop, process_noise=Q, H_fcn=H_fcn
)

# --------------------- Run EKF ---------------------
x_predicted = []
x_estimates = []
residuals = []

time_index = 0
x_k_minus_1 = init_state
for z in Z:
    time_index += 1
    print(f"time index: {time_index}")
    t_k = T[time_index]
    x_k = EKF.predict(t_k=t_k, state_k_minus_1=x_k_minus_1)
    z_k = ProbabilisticObservation(timestamp=t_k, mean=z, covariance=R)
    x_k_given_k, resid = EKF.update(t_k=t_k, z_k=z_k, state_k=x_k)

    x_predicted.append(x_k)
    residuals.append(resid)
    x_estimates.append(x_k_given_k)

    x_k_minus_1 = x_k_given_k

# ------------------ Compute IJK error from truth ----------

est_pos_errors = []
standard_deviations = []
for idx, x in enumerate(x_estimates):
    true_pos = Y_truth[idx + 1][:3]
    est_pos = x.mean[:3]
    est_pos_errors.append(true_pos - est_pos)
    variances = np.diag(x.covariance.matrix)
    standard_deviations.append(np.sqrt(variances[:3]))

# ------------------ Compute RIC (NTW) Error ----------

est_ric_pos_errors = []
ric_standard_deviations = []
for idx, x in enumerate(x_estimates):
    R_ECI_to_NTW = ntw_matrix(x.mean[:3], x.mean[3:])
    est_ric_pos_errors.append(R_ECI_to_NTW @ est_pos_errors[idx])
    P_eci = x.covariance.matrix[:3, :3]
    P_ric = R_ECI_to_NTW.T @ P_eci @ R_ECI_to_NTW
    ric_standard_deviations.append(np.sqrt(np.diag(P_ric)))


# ------------------ Save data to facilitate plot troubleshooting -------------
fpath = "/tmp/two_body_ekf_data.json"

data = {
    "timestamps": T.tolist(),
    "est_pos_errors": [x.tolist() for x in est_pos_errors],
    "standard_deviations": [x.tolist() for x in standard_deviations],
    "residuals": [x.tolist() for x in residuals],
    "est_ric_pos_errors": [x.tolist() for x in est_ric_pos_errors],
    "ric_standard_deviations": [x.tolist() for x in ric_standard_deviations],
}

with open(fpath, "w") as f:
    json.dump(data, f, indent=4)

# ---------------- Unpack Saved Data for Plotting ---------------------------
data = json.load(open(fpath, "r"))

T = data["timestamps"]
est_pos_errors = data["est_pos_errors"]
standard_deviations = data["standard_deviations"]
residuals = data["residuals"]
est_ric_pos_errors = data["est_ric_pos_errors"]
ric_standard_deviations = data["ric_standard_deviations"]


# ------------------ Plot IJK Error -----------------------------


top = ErrorStdInput(
    x=list(T[1:]),
    y_error=[x[0] for x in est_pos_errors],
    y_std=[x[0] for x in standard_deviations],
    x_axis_title="Time elapsed past epoch [s]",
    y_axis_title="I-axis error [m]",
)

middle = ErrorStdInput(
    x=list(T[1:]),
    y_error=[x[1] for x in est_pos_errors],
    y_std=[x[1] for x in standard_deviations],
    x_axis_title="Time elapsed past epoch [s]",
    y_axis_title="J-axis error [m]",
)

bottom = ErrorStdInput(
    x=list(T[1:]),
    y_error=[x[2] for x in est_pos_errors],
    y_std=[x[2] for x in standard_deviations],
    x_axis_title="Time elapsed past epoch [s]",
    y_axis_title="K-axis error [m]",
)

fig_input = Plot3dErrorWithStdInput(
    **{"top": top, "middle": middle, "bottom": bottom}
)

plotter = Plot3dErrorWithStd(upper=True, lower=True)
figure = plotter(data=fig_input)

# --------------- Plot RIC (NTW) Error -----------------------------


top = ErrorStdInput(
    x=list(T[1:]),
    y_error=[x[0] for x in est_ric_pos_errors],
    y_std=[x[0] for x in ric_standard_deviations],
    x_axis_title="Time elapsed past epoch [s]",
    y_axis_title="Radial error [m]",
)

middle = ErrorStdInput(
    x=list(T[1:]),
    y_error=[x[1] for x in est_ric_pos_errors],
    y_std=[x[1] for x in ric_standard_deviations],
    x_axis_title="Time elapsed past epoch [s]",
    y_axis_title="In-Track error [m]",
)

bottom = ErrorStdInput(
    x=list(T[1:]),
    y_error=[x[2] for x in est_ric_pos_errors],
    y_std=[x[2] for x in ric_standard_deviations],
    x_axis_title="Time elapsed past epoch [s]",
    y_axis_title="Cross-Track error [m]",
)

fig_input = Plot3dErrorWithStdInput(
    **{"top": top, "middle": middle, "bottom": bottom}
)

plotter = Plot3dErrorWithStd(upper=True, lower=True)
figure_ric = plotter(data=fig_input)
