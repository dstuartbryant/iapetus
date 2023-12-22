"""Atmospheric drag equations of motion (EOM) module.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro.constants import ROTATION

from .accelerations import AtmosphericDragAccelerations, accelerations
from .atmospheric_density import ExponentialAtmosphericModel
from .derivatives import AtmosphericDragPartials, partials

OMEGA_EARTH = ROTATION["Earth"]


@dataclass
class PreComputed:
    Bstar: float
    vreli: float
    vrelj: float
    vrelk: float
    vrel: float = None

    def __post_init__(self):
        if isinstance(self.vrel, type(None)):
            self.vrel = np.linalg.norm([self.vreli, self.vrelj, self.vrelk])


@dataclass
class PerturbedOutput:
    accelrations: AtmosphericDragAccelerations
    partials: Optional[AtmosphericDragPartials] = None


class Perturbations:
    """Atmospheric drag perturbations class.

    Handles acceleration and partial derivative computations for atmospheric
    drag dynamics.
    """

    def __init__(
        self,
        partials_flag: bool,
        Bstar_flag: bool,
    ):
        self.w = OMEGA_EARTH
        self.partials_flag = partials_flag
        self.Bstar_flag = Bstar_flag
        self.density_model = ExponentialAtmosphericModel(partials_flag)

    def pre_compute(
        self,
        pi: float,
        pj: float,
        vi: float,
        vj: float,
        vk: float,
        A: float,
        Cd: float,
        m: float,
    ):
        return PreComputed(
            Bstar=0.5 * A * Cd / m,
            vreli=vi + self.w * pj,
            vrelj=vj - self.w * pi,
            vrelk=vk,
        )

    def __call__(
        self,
        pi: float,
        pj: float,
        pk: float,
        p: float,
        vi: float,
        vj: float,
        vk: float,
        A: float,
        Cd: float,
        m: float,
    ) -> PerturbedOutput:
        pc = self.pre_compute(pi, pj, vi, vj, vk, A, Cd, m)
        d = self.density_model(pi, pj, pk, p)
        output = PerturbedOutput(
            accelrations=accelerations(
                pc.vreli,
                pc.vrelj,
                pc.vrelk,
                pc.vrel,
                d.rho,
                pc.Bstar,
            )
        )
        if self.partials_flag:
            output.partials = partials(
                pc.vreli,
                pc.vrelj,
                pc.vrelk,
                pc.vrel,
                d.rho,
                pc.Bstar,
                d.partials.drho_dpi,
                d.partials.drho_dpj,
                d.partials.drho_dpk,
                self.Bstar_flag,
            )
        return output
