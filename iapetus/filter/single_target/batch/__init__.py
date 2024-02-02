"""Batch least squares module."""

from copy import deepcopy
from typing import List

from .processor import Callable  # BatchData,
from .processor import (
    BatchPropagator,
    IntegrateConfig,
    ProbabilisticObservationSet,
    Pstate,
    batch_processor,
    np,
)


class BatchIterator:
    """Batch processing iterator class."""

    def __init__(
        self,
        batch_iter_tol: float,
        batch_max_iter: int,
        integrator_time_step: float,
        integrator_time_step_tolerance: float,
    ):
        """
        Args:
            batch_iter_tol (float): Iteration error tolerance; determines when
                batch processor iterations have converged
            batch_max_iter (int): Maximum number times to iterate on batch
                processor solution; guards agains endless loops in presence of
                non-convergence
            integrator_time_step (float): time step size for integrator
            integrator_time_step_tolerance (float): used to evenly divide
                integrator time steps
        """
        self.iter_tol = batch_iter_tol
        self.max_iter = batch_max_iter
        self.iconfig = IntegrateConfig(
            dt=integrator_time_step, tol=integrator_time_step_tolerance
        )
        self.num_iter = 0
        self.residuals_rms = []
        self.diff_residuals_rms = []

    def compute_residuals_rms(
        self, residuals: List[np.ndarray], obs_covarances: List[np.ndarray]
    ) -> float:
        """Computes root-mean-square (RMS) of observation residuals."""

        L = len(residuals)
        p = len(residuals[0])
        m = L * p
        summand = 0
        for epsilon, R in zip(residuals, obs_covarances):
            summand += epsilon.T @ np.linalg.inv(R) @ epsilon
        return np.sqrt(summand / m)

    def __call__(
        self,
        init_state: Pstate,
        xbar: np.ndarray,
        obs: ProbabilisticObservationSet,
        propagator: BatchPropagator,
        H_fcn: Callable,
    ):
        state = deepcopy(init_state)
        for i in range(self.max_iter):
            self.num_iter += 1
            bd = batch_processor(
                init_state, xbar, obs, propagator, self.iconfig, H_fcn
            )
            self.residuals_rms.append(
                self.compute_residuals_rms(
                    residuals=bd.residuals,
                    obs_covarances=bd.obs_error_covariances,
                )
            )
            state.mean += bd.state_error
            xbar = xbar - bd.state_error

            if i > 0:
                self.diff_residuals_rms.append(
                    abs(self.residuals_rms[-1] - self.residuals_rms[-2])
                )

            if self.num_iter == self.max_iter:
                state.covariance.matrix = bd.covariance
                return state

            if self.num_iter > 1:
                if self.diff_residuals_rms[-1] < self.iter_tol:
                    print("Stopping due to tolerance.")
                    state.covariance.matrix = bd.covariance
                    return state

        state.covariance.matrix = bd.covariance
        return state
