"""Batch least squares module."""


from typing import Callable

import numpy as np

from iapetus.data.observation import ProbabilisticObservationSet
from iapetus.data.state.probabilistic import State as Pstate

from .data import BatchData


def batch_processor(
    init_state: Pstate,
    xbar: np.ndarray,
    obs: ProbabilisticObservationSet,
    propagator: Callable,
    H_fcn: Callable,
):
    T = [x.timestamp for x in obs.observations]
    P = init_state.covariance.matrix
    Lam = np.linalg.inv(P)
    N = Lam @ xbar
    COV = []
    # Propagate reference trajectory
    X_ref = []
    STMS = []
    for t in T:
        x, phi = propagator(t0, X0, t)
        X_ref.append(x)
        STMS.append(phi)

    residuals = []
    for idx, z in enumerate(obs.observations):
        x_ref = X_ref[idx]
        phi = STMS[idx]
        y, H_tilde = H_fcn(z.timestamp, x_ref)
        resid = z.mean - y
        residuals.append(resid)
        H = H_tilde @ phi
        R = z.covariance.matrix
        R_inv = np.linalg.inv(R)
        Lam += H.T @ R_inv @ H
        N += H.T @ R_inv @ resid
        COV.append(np.linalg.inv(Lam))

    xhat = np.linalg.inv(Lam) @ N
    P = np.linalg.inv(Lam)

    return xhat, P, residuals, X_ref, STMS
