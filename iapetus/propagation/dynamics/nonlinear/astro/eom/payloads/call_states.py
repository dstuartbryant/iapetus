"""Equations of motion (EOM) call state interface models module."""

import numpy as np


def isnone(obj):
    """Checks if object is NoneType."""
    if isinstance(obj, type(None)) or (isinstance(obj, str) and obj == "None"):
        return True
    return False


class TwoBodyState:
    """State base class for Two-Body equations of motion (EOM) processing."""

    def __init__(
        self, pi: float, pj: float, pk: float, vi: float, vj: float, vk: float
    ):
        self.pi = pi
        self.pj = pj
        self.pk = pk
        self.vi = vi
        self.vj = vj
        self.vk = vk
        self.p = np.linalg.norm([self.pi, self.pj, self.pk])
        self.p3 = self.p**3
        self.p5 = self.p**5


class TwoBodyDragState(TwoBodyState):
    def __init__(
        self,
        pi: float,
        pj: float,
        pk: float,
        vi: float,
        vj: float,
        vk: float,
        Cd: float,
        A: float,
        m: float,
    ):
        super().__init__(pi, pj, pk, vi, vj, vk)
        self.Cd = Cd
        self.A = A
        self.m = m
        self._Bstar = None

    def _bstar(self):
        if isnone(self._Bstar):
            self._Bstar = 0.5 * self.A * self.Cd / self.m
        return self._Bstar

    @property
    def Bstar(self):
        return self._bstar()
