"""State transition matrix (STM) module."""

import numpy as np

from .constants import EARTH_SHAPE, EQUATORIAL_RADIUS, MU, ROTATION
from .dynamics.atmospheric_models import ExponentialAtmosphericModel
from .integrators.rk45 import rk45


class AstrodynamicsModelError(Exception):
    pass


class AstrodynamicsBaseModel:
    """Base class for astrodynamics models."""

    def accelerations(self, *args, **kwargs):
        """Computes acceleration components."""
        raise NotImplementedError(
            "This base class method needs to be overriden by a subclass."
        )

    def derivative_fcn(self, *args, **kwargs):
        """The derivative function method is passed into a Runge-Kutta
        integrator.
        """
        raise NotImplementedError(
            "This base class method needs to be overriden by a subclass."
        )
    
    def __call__(self, tspan: np.ndarray, dt: float, tspan_tolerance: float):


class TwoBody:
    """Two-body problem dynamics class.

    Integrates for propagation of state vector.
    """

    def __init__(self, X0: np.ndarray, mu: float):
        """
        Args:
            X0 (np.ndarray): 6x1 vector [x_m, y_m, z_m, u_mps, v_mps, w_mps]
            mu (float): gravitational constant in m^3/s^2


        """
        if len(X0) != 6:
            raise AstrodynamicsModelError(
                "Initial state vector must have 6 components."
            )
        self.X0 = X0
        self.mu = mu
        self.M = 6

    def accelerations(self, x_m, y_m, z_m, r_m):
        return (
            -self.mu * x_m / r_m**3,
            -self.mu * y_m / r_m**3,
            -self.mu * z_m / r_m**3,
        )

    def derivative_fcn(self, t, y):
        x_m = y[0]
        y_m = y[1]
        z_m = y[2]
        u_mps = y[3]
        v_mps = y[4]
        w_mps = y[5]
        phi = y[6:].reshape(6, 6)

        r_m = np.linalg.norm([x_m, y_m, z_m])

        ax_mps2, ay_mps2, az_mps2 = self.accelerations(x_m, y_m, z_m, r_m)

        A = self.A_matrix(x_m, y_m, z_m, r_m)

        # phi_dot = np.dot(A, phi)
        phi_dot = A @ phi

        ydot = [u_mps, v_mps, w_mps, ax_mps2, ay_mps2, az_mps2]
        ydot = np.concatenate(
            (
                ydot,
                phi_dot.reshape(
                    36,
                ),
            ),
            axis=0,
        )
        return ydot

    def state_only_derivative_fcn(self, t, y):
        x_m = y[0]
        y_m = y[1]
        z_m = y[2]
        u_mps = y[3]
        v_mps = y[4]
        w_mps = y[5]

        r_m = np.linalg.norm([x_m, y_m, z_m])

        ax_mps2, ay_mps2, az_mps2 = self.accelerations(x_m, y_m, z_m, r_m)

        ydot = [u_mps, v_mps, w_mps, ax_mps2, ay_mps2, az_mps2]

        return ydot

    def propagate_state_only(
        self, tspan: np.ndarray, dt: float, tspan_tolerance: float
    ):
        y0 = self.X0
        T, Y = rk45(
            self.state_only_derivative_fcn, y0, tspan, dt, tspan_tolerance
        )

        return T, Y

    def __call__(self, tspan: np.ndarray, dt: float, tspan_tolerance: float):
        phi0 = np.eye(6)
        y0 = np.concatenate(
            (
                self.X0,
                phi0.reshape(
                    36,
                ),
            ),
            axis=0,
        )
        T, Y = rk45(self.derivative_fcn, y0, tspan, dt, tspan_tolerance)
        X = [x[:6] for x in Y]
        Phi = [x[6:].reshape(6, 6) for x in Y]

        return T, X, Phi

    class TwoBodyStm:
        """Two-body problem dynamics class.

        Integrates for propagation of state vector and solving for state transition
        matrix.
        """

        def __init__(self, X0: np.ndarray, mu: float):
            """
            Args:
                X0 (np.ndarray): 6x1 vector [x_m, y_m, z_m, u_mps, v_mps, w_mps]
                mu (float): gravitational constant in m^3/s^2


            """
            if len(X0) != 6:
                raise AstrodynamicsModelError(
                    "Initial state vector must have 6 components."
                )
            self.X0 = X0
            self.mu = mu
            self.M = 6

        def A_matrix(self, x_m, y_m, z_m, r_m):
            block_11 = np.zeros((3, 3))
            block_12 = np.eye(3)

            block_22 = np.zeros((3, 3))

            block_21_11 = (
                3 * x_m**2 * self.mu / r_m**5 - self.mu / r_m**3
            )
            block_21_22 = (
                3 * y_m**2 * self.mu / r_m**5 - self.mu / r_m**3
            )
            block_21_33 = (
                3 * z_m**2 * self.mu / r_m**5 - self.mu / r_m**3
            )

            block_21_12 = block_21_21 = 3 * x_m * y_m * self.mu / r_m**5
            block_21_13 = block_21_31 = 3 * x_m * z_m * self.mu / r_m**5
            block_21_23 = block_21_32 = 3 * y_m * z_m * self.mu / r_m**5

            block_21 = np.array(
                [
                    [block_21_11, block_21_12, block_21_13],
                    [block_21_21, block_21_22, block_21_23],
                    [block_21_31, block_21_32, block_21_33],
                ]
            )
            return np.block([[block_11, block_12], [block_21, block_22]])

        def accelerations(self, x_m, y_m, z_m, r_m):
            return (
                -self.mu * x_m / r_m**3,
                -self.mu * y_m / r_m**3,
                -self.mu * z_m / r_m**3,
            )

        def derivative_fcn(self, t, y):
            x_m = y[0]
            y_m = y[1]
            z_m = y[2]
            u_mps = y[3]
            v_mps = y[4]
            w_mps = y[5]
            phi = y[6:].reshape(6, 6)

            r_m = np.linalg.norm([x_m, y_m, z_m])

            ax_mps2, ay_mps2, az_mps2 = self.accelerations(x_m, y_m, z_m, r_m)

            A = self.A_matrix(x_m, y_m, z_m, r_m)

            # phi_dot = np.dot(A, phi)
            phi_dot = A @ phi

            ydot = [u_mps, v_mps, w_mps, ax_mps2, ay_mps2, az_mps2]
            ydot = np.concatenate(
                (
                    ydot,
                    phi_dot.reshape(
                        36,
                    ),
                ),
                axis=0,
            )
            return ydot

        def state_only_derivative_fcn(self, t, y):
            x_m = y[0]
            y_m = y[1]
            z_m = y[2]
            u_mps = y[3]
            v_mps = y[4]
            w_mps = y[5]

            r_m = np.linalg.norm([x_m, y_m, z_m])

            ax_mps2, ay_mps2, az_mps2 = self.accelerations(x_m, y_m, z_m, r_m)

            ydot = [u_mps, v_mps, w_mps, ax_mps2, ay_mps2, az_mps2]

            return ydot

        def propagate_state_only(
            self, tspan: np.ndarray, dt: float, tspan_tolerance: float
        ):
            y0 = self.X0
            T, Y = rk45(
                self.state_only_derivative_fcn, y0, tspan, dt, tspan_tolerance
            )

            return T, Y

        def __call__(
            self, tspan: np.ndarray, dt: float, tspan_tolerance: float
        ):
            phi0 = np.eye(6)
            y0 = np.concatenate(
                (
                    self.X0,
                    phi0.reshape(
                        36,
                    ),
                ),
                axis=0,
            )
            T, Y = rk45(self.derivative_fcn, y0, tspan, dt, tspan_tolerance)
            X = [x[:6] for x in Y]
            Phi = [x[6:].reshape(6, 6) for x in Y]

            return T, X, Phi
