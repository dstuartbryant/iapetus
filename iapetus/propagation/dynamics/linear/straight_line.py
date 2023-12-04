"""Straight line, constant velocity dynamics module.

"""
import numpy as np


class LinearStateTransitionMatrix:
    """Linear state transition matrix base class."""

    def __init__(self, dimensions: int):
        self.dimensions = dimensions
        self.matrix_fcn = self._generate_stm()

    def _generate_stm(self) -> np.ndarray:
        """Returns an nxn state transition matrix.

        n = self.dimensions
        """
        raise NotImplementedError(
            "This base class method is meant to be overriden by subclass."
        )

    def __call__(self, dt_s: float):
        return self.matrix_fcn(dt_s)


class StmStraightLineConstantVelocity(LinearStateTransitionMatrix):
    """Straight line, constant velocity state transition matrix subclass."""

    def _generate_stm(self) -> np.ndarray:
        """Returns an nxn state transition matrix.

        n = self.dimensions
        """
        Identity = np.eye(self.dimensions)
        Z = np.zeros((self.dimensions, self.dimensions))

        def stm(dt_s: float) -> np.ndarray:
            """Straight line, constant velocity state transition matrix as
            a function of time.

            Args:
                dt_s (float): delta t in seconds
            """
            return np.block([[Identity, Identity * dt_s], [Z, Identity]])

        return stm


class StmStraightLine(LinearStateTransitionMatrix):
    """Straight line (non-constant velocity) state transition matrix class."""

    def _generate_stm(self) -> np.ndarray:
        """Returns an nxn state transition matrix.

        n = self.dimensions
        """
        Identity = np.eye(self.dimensions)
        Z = np.zeros((self.dimensions, self.dimensions))

        def stm(dt_s: float) -> np.ndarray:
            """Straight line, non-constant velocity state transition matrix as
            a function of time.

            Args:
                dt_s (float): delta t in seconds
            """
            return np.block(
                [
                    [Identity, Identity * dt_s, Identity * dt_s**2 / 2],
                    [Z, Identity, Identity * dt_s],
                    [Z, Z, Identity],
                ]
            )

        return stm
