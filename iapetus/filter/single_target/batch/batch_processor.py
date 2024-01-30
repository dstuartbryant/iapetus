"""Batch least squares processor."""

from typing import Callable

import numpy as np

from iapetus.data.observation import ProbabilisticObservationSet


def batch_processor(
    t0: float,
    X0: np.ndarray,
    P0: np.ndarray,
    xbar: np.ndarray,
    obs: ProbabilisticObservationSet,
    propagator: Callable,
    H_fcn: Callable,
):
    T = [x.timestamp for x in obs.observations]
    P = P0
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
    xhat_mapped_fwd = phi @ xhat
    P = np.linalg.inv(Lam)
    # X0 += xhat
    # xbar = xbar - xhat

    return xhat, xhat_mapped_fwd, P, residuals, X_ref, STMS
