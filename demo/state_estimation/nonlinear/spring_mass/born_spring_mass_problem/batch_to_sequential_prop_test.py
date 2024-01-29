"""Need to create a sequential propagator out of given batch propagator
for this demo and this script is for prototyping and testing that.
"""

import math

import numpy as np

from iapetus.propagation.integrators.rk45 import rk45


def batch_propagator(t0: float, X0: np.ndarray, t: float):
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


def computing_equilibrium_distance(accel: float, pos: float):
    k1 = 2.5
    k2 = 3.7
    m = 1.5

    xbar = -1 * (m * accel / (-k1 + k2) - pos)
    return xbar


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


def sequential_propagator_deriv_fcn(t, y):
    """Spring mass propagator based on Eq. 4.8.5 from Ref. [1]."""

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


def sequential_propagator(X0, T):
    phi_0 = np.eye(2)
    X_k = X0
    for t_kp1, t_k in zip(T[1:], T[:-1]):
        tspan = [t_k, t_kp1]
        dt = t_kp1 - t_k
        X0 = add_stm_to_state_vector(X0, phi)
        T_out, Y_out = rk45(sequential_propagator_deriv_fcn, X0, T, 1, 1e-3)

    states = []
    stms = []
    for idx in range(len(Y)):
        x, phi = remove_stm_from_state_vector(Y[idx])
        states.append(x)
        stms.append(phi)
    return T, states, stms


X0 = np.array([3.0, 0.0])
T = np.arange(0, 11, 1)

# ----------------- Run Batch Prop ------------------------


X_batch = []
Phi_batch = []
for t in T:
    x, phi = batch_propagator(0, X0, t)
    X_batch.append(x)
    Phi_batch.append(phi)


# ----------------- Run Sequential Prop ------------------------

T_out, X_seq, Phi_seq = sequential_propagator(X0, T)


# ---------------- Compare Output without using STM -----------------------

dx = []
dv = []
for idx in range(len(X_batch)):
    dx.append(X_batch[idx][0] - X_seq[idx][0])
    dv.append(X_batch[idx][1] - X_seq[idx][1])

mean_dx = np.mean(dx)
mean_dv = np.mean(dv)

assert mean_dx < 1e-9
assert mean_dv < 1e-9


# ---------------- Compare Output By using STMs -----------------------
X_stm_states = []
X_k = X0
for idx in range(len(Phi_seq)):
    X_kp1 = Phi_seq[idx] @ X_k
    X_stm_states.append(X_kp1)
    X_k = X_kp1
