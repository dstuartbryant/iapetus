"""Extras module for astrodynamics propagation.

Due to nature of having an state-space-context-state used for dynamically
configuring partial derivative matrices, and that the Runge-Kutta integrators
require having state vectors as input, i.e., unlabeled numpy ndarray, it
is sometimes necessary to provide eom-context-states "extra" variables for
underlying computations, e.g., Area, mass, and Cd for atmospheric-drag
perturbation modeling.

Hence the use of "extras" for derivative function configuration.
"""

from .call_states import TwoBodyDragState


class AtmosphericDragExtras:
    """Atmospheric drag extras class.

    Scenario Coverage
    -----------------
    1. When Cd and Bstar are NOT in the state-space-context-state,
        - the underlying assumption is that Cd, A, m, and therefore Bstar,
          are known with high level of certainty, i.e., we assume we know them
          perfectly,
        - so, Cd, A, and must be provided as extras,

    2. When Cd IS (and therefore Bstar is not) in the
       state-space-context-state,
        - the underlying assumption is that A and m are known with a high level
          of certainty,
        - we want the state-space-context-state Cd in the incoming y vector to
          the derivative function to be used in computations, and NOT the Cd
          provided in the initial state because the Cd value will change over
          time via state estimation,
        - so we only want A and m to be provided as extras,

    3. When Bstar IS (and therefore Cd is not) in the
       state-space-context-state,
        - the underlying assumption is that we want Bstar to "soak" up
          uncertainties in the state estimation, which means it will change
          over time via state estimation,
        - this is the case when we have a high level of uncertainty in A, m,
          and Cd,
        - therefore, we do NOT want to provide any extras to the underlying EOM


    """

    def __init__(
        self, init_state: TwoBodyDragState, Cd_flag: bool, Bstar_flag: bool
    ):
        if not Cd_flag and not Bstar_flag:
            self.extras_list = ["Cd", "A", "m"]
        elif Cd_flag and not Bstar_flag:
            self.extras_list = ["A", "m"]
        elif Bstar_flag:
            self.extras_list = []
        self.Cd = init_state.Cd
        self.A = init_state.A
        self.m = init_state.m

    def __call__(self):
        return {k: getattr(self, k) for k in self.extras_list}
