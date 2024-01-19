"""State noise compensation (SNC) module.


References:
[1] Born, et al., Statistical Orbit Determination
"""

import numpy as np


class Gamma:
    """Process noise transition matrix class.

    The process noise transition matrix is commonly denoted as the greek
    capital letter Gamma, see Ref. [1], pg. 228.
    """

    def __init__(self, m: int):
        """
        Args:
            m (int): Determines size/configuration of transition
                matrix. If a state vector consists only of a position and
                velocity components, and k denoted the dimension, i.e., number
                of elements in either position or velocity vectors, then the
                dimension of the statevector would be n = 2*k, and m would
                then be m = k, thus resulting in an nxm Gamma matrix.
                Examples:
                    1. k = 1, n = 2, m = 1, Gamma: 2x1 matrix
                    2. k = 2, n = 4, m = 2, Gamma: 4x2 matrix
                    3. k = 3, n = 6, m = 3, Gamma: 6x3 matrix
        """
        self.m = m
        self.I_matrix = np.eye(self.m)

        raise Exception(
            "This model needs granular revision, as it's not general enough. "
            "It's only *a* version based on STAT OD 2 notes. See Shalom's "
            "secion on Motion models, e.g., pg. 43 CWNA and how it compares "
            " to DWNA on pg. 48. DWNA is like this, i.e., like STAT OD 2. "
            " But why? What's the difference? It has something to do with "
            "assumptions about dynamics and sampling frequency, I think."
        )

    def __call__(self, dt: float) -> np.ndarray:
        return dt * np.vstack([dt / 2 * self.I_matrix, self.I_matrix])
