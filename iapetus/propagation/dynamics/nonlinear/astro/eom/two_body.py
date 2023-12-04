"""Two-body equations of motion module."""

import numpy as np


class TwoBodyAccelerations:
    """Computes accelerations under two-body motion."""

    def __init__(self, mu_m3ps2: float):
        """
        Args:
            mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]
        """
        self.mu = mu_m3ps2

    def ax(self, *args, **kwargs):
        pass

    def ay(self, *args, **kwargs):
        pass

    def az(self, *args, **kwargs):
        pass


class TwoBodyPartialDerivatives:
    """Computes partial derivative elements under two-body motion."""

    def __init__(self, mu_m3ps2: float):
        """
        Args:
            mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]
        """
        self.mu = mu_m3ps2
        self.dvi_dpi = 0.0

    def _dvidpi(self, *args, **kwargs):
        pass

    def _dvjdpj(self, *args, **kwargs):
        pass

    def _dvkdpk(self, *args, **kwargs):
        pass

    def _dvidpj(self, *args, **kwargs):
        pass


class _TwoBodyEom:
    def __init__(self, mu_m3ps2: float):
        """
        Args:
            mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]
        """
        self.mu = mu_m3ps2
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.u = 0.0
        self.v = 0.0
        self.w = 0.0
        self.r = 0.0
        self.r3 = 0.0
        self.r5 = 0.0

    def update_r(self, x_m: float, y_m: float, z_m: float):
        """Updates Earth-centered position vector magnitude and pre-computes
        exponents thereof for more efficient computation.

        Args:
            x_m (float): i-component position in meters [m]
            y_m (float): j-component position in meters [m]
            z_m (float): k-component position in meters [m]
        """
        self.r = np.sqrt(x_m**2 + y_m**2 + z_m**2)
        self.r3 = self.r**3
        self.r5 = self.r**5

    def update_state(
        self,
        x_m: float,
        y_m: float,
        z_m: float,
        u_mps: float,
        v_mps: float,
        w_mps: float,
    ):
        """Updates Earth-centered position vector magnitude and pre-computes
        exponents thereof for more efficient computation.

        Args:
            x_m (float): i-component position in meters [m]
            y_m (float): j-component position in meters [m]
            z_m (float): k-component position in meters [m]
            u_mps (float): i-component velocity in meters per second [mps]
            v_mps (float): j-component velocity in meters per second [mps]
            w_mps (float): k-component velocity in meters per second [mps]
        """
        self.x = x_m
        self.y = y_m
        self.z = z_m
        self.u = u_mps
        self.v = v_mps
        self.w = w_mps

        self.update_r(x_m, y_m, z_m)

    # Accelerations
    def acceleration(self, p_m: float) -> float:
        """Two-body acceleration, variable component.

        Args:
            p_m (float): One of the component position values in meters [m].

        Returns:
            (float): acceleration component in meters per second squared [mps2]
        """
        return -self.mu * p_m / self.r3

    def accelration_i_mps2(self) -> float:
        return self.acceleration(self.x)

    def accelration_j_mps2(self) -> float:
        return self.acceleration(self.y)

    def accelration_k_mps2(self) -> float:
        return self.acceleration(self.z)

    # Partial derivatives matrix components
    def _deriv_accel_wrt_position_same_component(self, p_m) -> float:
        """Partial derivative of acceleration with respect to
        position, where they share the same identity as either the i-th, j-th,
        or k-th component.

        Args:
            p_m (float): One of the component position values in meters [m].

        Returns:
            float
        """
        return 3 * p_m**2 * self.mu / self.r5 - self.mu / self.r3

    def _derivative_accelration_i_mps2_wrt_position_i_m(self) -> float:
        return self._deriv_accel_wrt_position_same_component(self.x)

    def _derivative_accelration_j_mps2_wrt_position_j_m(self) -> float:
        return self._deriv_accel_wrt_position_same_component(self.y)

    def _derivative_accelration_k_mps2_wrt_position_k_m(self) -> float:
        return self._deriv_accel_wrt_position_same_component(self.z)

    def _deriv_accel_wrt_position_different_component(
        self, p1_m: float, p2_m: float
    ) -> float:
        """Partial derivative of acceleration with respect to
        position, where they do not share the same identity as either the i-th,
        j-th, or k-th component.

        Args:
            p1_m (float): One of the component position values in meters [m].
            p2_m (float): One of the other component position values in meters
                [m].

        Returns:
            float
        """
        return 3 * p1_m * p2_m * self.mu / self.r5

    def _derivative_accelration_i_mps2_wrt_position_j_m(self) -> float:
        return self._deriv_accel_wrt_position_different_component(
            self.x, self.y
        )

    def _derivative_accelration_i_mps2_wrt_position_k_m(self) -> float:
        return self._deriv_accel_wrt_position_different_component(
            self.x, self.z
        )


# First row
def derivative_velocity_i_mps_wrt_position_i_m(*args, **kwargs) -> float:
    """Partial derivative of i-component velocity with respect to i-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_i_mps_wrt_position_j_m(*args, **kwargs) -> float:
    """Partial derivative of i-component velocity with respect to j-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_i_mps_wrt_position_k_m(*args, **kwargs) -> float:
    """Partial derivative of i-component velocity with respect to k-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_i_mps_wrt_velocity_i_mps(*args, **kwargs) -> float:
    """Partial derivative of i-component velocity with respect to i-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 1.0


def derivative_velocity_i_mps_wrt_velocity_j_mps(*args, **kwargs) -> float:
    """Partial derivative of i-component velocity with respect to j-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_i_mps_wrt_velocity_k_mps(*args, **kwargs) -> float:
    """Partial derivative of i-component velocity with respect to k-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


# Second row
def derivative_velocity_j_mps_wrt_position_i_m(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to i-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_j_mps_wrt_position_j_m(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to j-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_j_mps_wrt_position_k_m(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to k-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_j_mps_wrt_velocity_i_mps(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to i-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_j_mps_wrt_velocity_j_mps(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to j-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 1.0


def derivative_velocity_j_mps_wrt_velocity_k_mps(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to k-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


# Third row
def derivative_velocity_k_mps_wrt_position_i_m(*args, **kwargs) -> float:
    """Partial derivative of k-component velocity with respect to i-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_k_mps_wrt_position_j_m(*args, **kwargs) -> float:
    """Partial derivative of k-component velocity with respect to j-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_k_mps_wrt_position_k_m(*args, **kwargs) -> float:
    """Partial derivative of k-component velocity with respect to k-component
    position.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_k_mps_wrt_velocity_i_mps(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to i-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_k_mps_wrt_velocity_j_mps(*args, **kwargs) -> float:
    """Partial derivative of k-component velocity with respect to j-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_k_mps_wrt_velocity_k_mps(*args, **kwargs) -> float:
    """Partial derivative of k-component velocity with respect to k-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 1.0


# Fourth row
def derivative_accelration_i_mps2_wrt_position_i_m(
    x_m: float, r_m3: float, r_m5: float, mu_m3ps2: float
) -> float:
    """Partial derivative of i-component acceleration with respect to
    i-component position.

    Args:
        x_m (float): i-component position in meters [m]
        r_m3 (float): cubed position vector magnitude in meters cubed [m3]
        r_m5 (float): position vector magnitude raised to 5th power in
            meters^5 [m5]
        mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]

    Returns:
        float
    """
    return 3 * x_m * mu_m3ps2 / r_m5 - mu_m3ps2 / r_m3


def derivative_accelration_i_mps2_wrt_position_j_m(
    x_m: float, y_m: float, r_m5: float, mu_m3ps2: float
) -> float:
    """Partial derivative of i-component acceleration with respect to
    j-component position.

    Args:
        x_m (float): i-component position in meters [m]
        y_m (float): j-component position in meters [m]
        r_m5 (float): position vector magnitude raised to 5th power in
            meters^5 [m5]
        mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]

    Returns:
        float
    """
    return 3 * x_m * y_m * mu_m3ps2 / r_m5


def derivative_accelration_i_mps2_wrt_position_k_m(
    x_m: float, z_m: float, r_m5: float, mu_m3ps2: float
) -> float:
    """Partial derivative of i-component acceleration with respect to
    k-component position.

    Args:
        x_m (float): i-component position in meters [m]
        z_m (float): k-component position in meters [m]
        r_m5 (float): position vector magnitude raised to 5th power in
            meters^5 [m5]
        mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]

    Returns:
        float
    """
    return 3 * x_m * z_m * mu_m3ps2 / r_m5


def derivative_velocity_k_mps_wrt_velocity_i_mps(*args, **kwargs) -> float:
    """Partial derivative of j-component velocity with respect to i-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_k_mps_wrt_velocity_j_mps(*args, **kwargs) -> float:
    """Partial derivative of k-component velocity with respect to j-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 0.0


def derivative_velocity_k_mps_wrt_velocity_k_mps(*args, **kwargs) -> float:
    """Partial derivative of k-component velocity with respect to k-component
    velocity.

    Args:
        Any

    Returns:
        float
    """
    return 1.0
