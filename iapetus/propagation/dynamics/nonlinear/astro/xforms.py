"""Coordinate transformation module.


References:
[1] Vallado, 5th ed.

"""

import numpy as np


def rsw_matrix(r_m: np.ndarray, v_mps: np.ndarray):
    """Computes the radial, along-track, cross-track (RSW) coordinate
    transformation matrix.

    References:
        See [1], pgs. 159 & 166.

    Args:
        r_m (np.ndarray): 3x1 position vector [meters] in inertial frame
        v_mps (np.ndarray): 3x1 velocity vector [meters per second] in inertial
            frame

    Returns:
        (np.ndarray): 3x3 matrix that transforms INERTIAL coordinates to RSW
    """

    R = r_m / np.linalg.norm(r_m)
    _w = np.cross(r_m, v_mps)
    W = _w / np.linalg.norm(_w)
    S = np.cross(W, R)

    return np.vstack((R, S, W))


def ntw_matrix(r_m: np.ndarray, v_mps: np.ndarray):
    """Computes the radial, in-track, cross-track (NTW) coordinate
    transformation matrix.

    References:
        See [1], pgs. 157 & 164.

    Args:
        r_m (np.ndarray): 3x1 position vector [meters] in inertial frame
        v_mps (np.ndarray): 3x1 velocity vector [meters per second] in inertial
            frame

    Returns:
        (np.ndarray): 3x3 matrix that transforms INERTIAL coordinates to NTW
    """

    T = v_mps / np.linalg.norm(v_mps)
    _w = np.cross(r_m, v_mps)
    W = _w / np.linalg.norm(_w)
    N = np.cross(T, W)

    return np.vstack((N, T, W))
