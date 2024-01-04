"""Equations of motion (EOM) state contexts module."""

import numpy as np


class TwoBodyWithoutStm:
    """Two body state in EOM context with no STM computations."""

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


class TwoBodyWithStm(TwoBodyWithoutStm):
    """Two body state in EOM context with STM computations."""

    _p5 = None

    @property
    def p5(self):
        if isinstance(self._p5, type(None)):
            self._p5 = self.p**5
        return self._p5


class TwoBodyDragWithoutStm(TwoBodyWithoutStm):
    """Two body state in EOM context without STM computation and without Cd or
    Bstar in the state vector.
    """

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
        self.Bstar = 0.5 * Cd * A / m


class TwoBodyDragWithStm(TwoBodyWithStm):
    """Two body state in EOM context with STM computation and without Cd or
    Bstar in the state vector.
    """

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
        self.Bstar = 0.5 * Cd * A / m


class TwoBodyDragCdWithoutStm(TwoBodyDragWithoutStm):
    """Two body state in EOM context without STM computations and with Cd in
    the state vector.
    """

    pass


class TwoBodyDragCdWithStm(TwoBodyDragWithStm):
    """Two body state in EOM context with STM computation and with Cd in the
    state vector.
    """

    pass


class TwoBodyDragBstarWithoutStm(TwoBodyWithoutStm):
    """Two body state in EOM context without STM computation and with Bstar in
    the state vector.
    """

    def __init__(
        self,
        pi: float,
        pj: float,
        pk: float,
        vi: float,
        vj: float,
        vk: float,
        Bstar: float,
    ):
        super().__init__(pi, pj, pk, vi, vj, vk)
        self.Bstar = Bstar


class TwoBodyDragBstarWithStm(TwoBodyWithStm):
    """Two body state in EOM context with STM computation and with Bstar in
    the state vector.
    """

    def __init__(
        self,
        pi: float,
        pj: float,
        pk: float,
        vi: float,
        vj: float,
        vk: float,
        Bstar: float,
    ):
        super().__init__(pi, pj, pk, vi, vj, vk)
        self.Bstar = Bstar
