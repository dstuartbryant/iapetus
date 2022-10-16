"""Linear process noise models.

References:
[1] Bar-Shalom, Yaakov, Peter K. Willett, and Xin Tian. Tracking and data
    fusion. Vol. 11. Storrs, CT, USA:: YBS publishing, 2011.

[2] Blackman, Samuel, and Robert Popoli. "Design and analysis of modern 
    tracking systems(Book)." Norwood, MA: Artech House, 1999. (1999).

"""

import numpy as np


def continuous_white_noise_acceleration(q, dt):
    """Returns a 2D white, discrete process noise covariance matrix.

    According to Ref. [1], `q` is the power spectral density of the continous
    time white nosie that models the motion uncertainty, which is a filter
    design parameter - the "noise density".

    See Refs. [1-2].

    Args:
        q (float): noise density
        dt (float): time step, i.e., sample rate
    """

    return q * np.block(
        [
            [1 / 3 * np.eye(2) * dt**3, 1 / 2 * np.eye(2) * dt**2],
            [1 / 2 * np.eye(2) * dt**2, np.eye(2) * dt],
        ]
    )
