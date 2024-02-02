"""Batch least squares propagators interfacing module."""

from abc import ABC, abstractmethod
from typing import List, Optional, Tuple

import numpy as np

from iapetus.propagation.integrators import IntegrateConfig
from iapetus.propagation.propagators import AstroProp


class BatchPropagatorError(Exception):
    pass


class BatchPropagator(ABC):
    """Batch processor propagator base class."""

    @abstractmethod
    def __call__(
        self,
        t0: float,
        X0: np.ndarray,
        tspan: List[float],
        iconfig: Optional[IntegrateConfig] = None,
    ) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """Calls to run a propagation process.

        Args:
            t0 (float): an initial reference timestamp
            X0 (np.ndarray): initial state vector at time t0
            tspan (List[float]): list of reference trajectory timestamps that
                match observation epochs; times in this list do not need to be
                consistently incremented
            iconfig (IntegrateConfig): integrator configurations

        Returns:
            T (List[float]): list of reference trajectory timestamps output
                from the integrator
            Xref (List[np.ndarray]): list of reference trajectory state vectors
                ordered in accordance with tspan
            Phis (List[np.ndarray]): list of error state transition matrices
                ordered in accordance with tspan

        Raises:
            BatchPropagatorError if initial time t0 not first element in tspan.
            NOTE: Rationale: If configured without t0 input, the burden on a
            user then becomes to remember everytime they invoke this propagator
            that they'll have to remember to place t0 within tspan, which
            comes with inherit risks due to general forgetfulnees and the
            complexity of the tools with which this feature is used. Rather,
            force the user to explicitly define t0 and leverage this exception
            as a reminder to ensure that t0 is also the first element in tspan.
        """
        if t0 != tspan[0]:
            raise BatchPropagatorError(
                "Expected to find `t0` at beginning of `tspan`"
            )


class BatchAstroPropagator(BatchPropagator):
    """Batch processor astrodynamics propagator class."""

    def __init__(self, init_state: dict, aprop: AstroProp):
        """
        Args:
            init_state (dict): Initial state dictionary that conforms to
                astrodynamics state models
            aprop (AstroProp): astrodynamics propagator instance
        """
        self.init_state_dict = init_state
        self.aprop = aprop
        self.state_context_manager = self.StateContextManager(init_state)

    def __call__(
        self,
        t0: float,
        X0: np.ndarray,
        tspan: List[float],
        iconfig: IntegrateConfig,
    ) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        ui_state = self.scm.derivative_fcn_to_ui_context(X0)
        T, X_ref, Phis = self.aprop(
            tspan=tspan,
            dt=iconfig.dt,
            tspantol=iconfig.tol,
            ui_state=ui_state.dict(),
        )

        return T, X_ref, Phis

    @property
    def StateContextManager(self):
        return self.aprop.prop_init.state_context_manager_factory()

    @property
    def scm(self):
        return self.state_context_manager
