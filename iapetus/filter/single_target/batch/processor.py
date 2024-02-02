"""Batch least squares module."""


from typing import Callable

import numpy as np

from iapetus.data.observation import ProbabilisticObservationSet
from iapetus.data.state.probabilistic import State as Pstate
from iapetus.propagation.integrators import IntegrateConfig

from .data import BatchData
from .propagators import BatchPropagator


def batch_processor(
    init_state: Pstate,
    xbar: np.ndarray,
    obs: ProbabilisticObservationSet,
    propagator: BatchPropagator,
    iconfig: IntegrateConfig,
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
    Tout, X_ref, STMS = propagator(
        t0=init_state.timestamp, X0=init_state.mean, tspan=T, iconfig=iconfig
    )

    residuals = []
    obs_error_covs = []
    for idx, z in enumerate(obs.observations):
        x_ref = X_ref[idx]
        phi = STMS[idx]
        y, H_tilde = H_fcn(z.timestamp, x_ref)
        resid = z.mean - y
        residuals.append(resid)
        H = H_tilde @ phi
        R = z.covariance.matrix
        obs_error_covs.append(R)
        R_inv = np.linalg.inv(R)
        Lam += H.T @ R_inv @ H
        N += H.T @ R_inv @ resid
        COV.append(np.linalg.inv(Lam))

    xhat = np.linalg.inv(Lam) @ N
    P = np.linalg.inv(Lam)

    return BatchData(
        timestamp=init_state.timestamp,
        state_error=xhat,
        covariance=P,
        timesteps=T,
        residuals=residuals,
        obs_error_covariances=obs_error_covs,
    )
