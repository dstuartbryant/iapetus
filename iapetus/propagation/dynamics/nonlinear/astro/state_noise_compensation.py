"""Astrodynamics state noise compensation (SNC) module.

References:
[1] Born et al., Statistical Orbit Determination,
[2] Duncan, M. and Long, A., Realistic Covariance Prediction for the Earth
    Science Constellation, AIAA/AAS Astrodynamics Specialist Conference and
    Exhibit, 21-24 August 2006, Keystone, Colorado.
[3] Bar-Shalom et al., Tracking and Data Fusion, YBS Publishing, 2011.
"""

import numpy as np


class QRic:
    """Radial, in-track, cross-track (RIC) proces noise matrix.

    This model follows the approach in Ref. [2], where the process noise
    transition matrix follows a continuous white noise acceleration (CWNA)
    model. See Ref. [3], p. 43 for more info on CWNA model.

    """

    def __init__(self, qr_mps2: float, qi_mps2: float, qc_mps2: float):
        """
        Args:
            qr_mps2 (float): radial component of acceleration error vector
                [m/s^2]
            qi_mps2 (float): in-track component of acceleration error vector
                [m/s^2]
            qc_mps2 (float): cross-track component of acceleration error vector
                [m/s^2]
        """
        self.qr_mps2 = qr_mps2
        self.qi_mps2 = qi_mps2
        self.qc_mps2 = qc_mps2
        self.q = np.array([qr_mps2**2, qi_mps2**2, qc_mps2**2])

    def q_matrix(self, dt_s: float, q: np.ndarray) -> np.ndarray:
        """Returns a Q_RIC matrix.

        Args:
            dt_s (float): change in time in seconds
            q (np.ndarray): process noise accelerations vector

        Returns:
            (np.ndarray): matrix
        """
        b11_diag = dt_s**4 / 3 * q
        b12_diag = dt_s**3 / 2 * q
        b22_diag = dt_s**2 * q
        return np.block(
            [
                [np.diag(b11_diag), np.diag(b12_diag)],
                [np.diag(b12_diag), np.diag(b22_diag)],
            ]
        )

    def __call__(self, dt_s: float) -> np.ndarray:
        return self.q_matrix(dt_s, self.q)


class QInertialRic(QRic):
    """Inertial process noise matrix from RIC process noise matrix."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.q_ric = self.q

    def q_inertial(self, R_inertial_to_RIC: np.ndarray) -> np.ndarray:
        """Returns RIC process noise accelerations tranformed into inertial
        reference frame coordinates.

        Args:
            R_inertial_to_RIC (np.ndarray): coordinate transformation matrix:
                inertial to RIC

        Returns:
            (np.ndarray): process noise accelerations vector in interial
                coordinates
        """

        return R_inertial_to_RIC.T @ self.q_ric

    def __call__(
        self, dt_s: float, R_inertial_to_RIC: np.ndarray
    ) -> np.ndarray:
        q_inertial = self.q_inertial(R_inertial_to_RIC)

        return self.q_matrix(dt_s, q_inertial)


class QInertialCd(QInertialRic):
    """Radial, in-track, cross-track (RIC) and drag coefficient (Cd) proces
    noise matrix in inertial reference frame.

    This model follows the approach in Ref. [2].
    """

    def __init__(
        self, qr_mps2: float, qi_mps2: float, qc_mps2: float, qcd: float
    ):
        """
        Args:
            qr_mps2 (float): radial component of acceleration error vector
                [m/s^2]
            qi_mps2 (float): in-track component of acceleration error vector
                [m/s^2]
            qc_mps2 (float): cross-track component of acceleration error vector
                [m/s^2]
            qcd (float): drag coefficient process noise element
        """
        super().__init__(qr_mps2, qi_mps2, qc_mps2)
        self.qcd = qcd

    def __call__(self, dt_s: float, R_inertial_to_RIC: np.ndarray):
        q_inertial = self.q_inertial(R_inertial_to_RIC)
        q_matrix = self.q_matrix(dt_s, q_inertial)
        return np.block(
            [
                [q_matrix, np.zeros((6, 1))],
                [np.zeros((1, 6)), self.qcd**2],
            ]
        )
