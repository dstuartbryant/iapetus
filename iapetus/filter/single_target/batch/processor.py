"""Batch least squares module."""

from typing import Callable, Optional

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
    H_fcn: Callable,
    iconfig: Optional[IntegrateConfig] = None,
):
    """Batch least squares processor.

    Args:
        init_state  (Pstate): initial state with timestamp, mean state vector,
            and covariance matrix attributes
        xbar (np.ndarray): initial state error vector
        obs (ProbabilisticObservationSet): Set of observations
        propagator (BatchPropagator): propagator instance that generate
            reference trajectory
        H_fcn (Callable): encapsulates mapping from state space to observation
            space
        iconfig (IntegrateConfig, optional): integrator configuration, if
            the propagator includes an integrator without default configs

    Returns:
        (BatchData): state estimation results in a BatchData model
    """
    T = [x.timestamp for x in obs.observations]
    P = init_state.covariance.matrix
    Lam = np.linalg.inv(P)
    N = Lam @ xbar
    COV = []
    X_ref = []
    STMS = []

    if not iconfig:
        Tout, X_ref, STMS = propagator(
            t0=init_state.timestamp, X0=init_state.mean, tspan=T
        )
    else:
        Tout, X_ref, STMS = propagator(
            t0=init_state.timestamp,
            X0=init_state.mean,
            tspan=T,
            iconfig=iconfig,
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
