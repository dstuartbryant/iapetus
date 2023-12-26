"""Atmospheric drag equations of motion (EOM) module.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro.constants import ROTATION

from ...payloads import AtmosphericDragInitConfig, TwoBodyDragState
from .accelerations import AtmosphericDragAccelerations, accelerations
from .atmospheric_density import ExponentialAtmosphericModel
from .derivatives import AtmosphericDragPartials, partials

OMEGA_EARTH = ROTATION["Earth"]


@dataclass
class PreComputed:
    vreli: float
    vrelj: float
    vrelk: float
    vrel: float = None

    def __post_init__(self):
        if isinstance(self.vrel, type(None)):
            self.vrel = np.linalg.norm([self.vreli, self.vrelj, self.vrelk])


@dataclass
class PerturbedOutput:
    accelerations: AtmosphericDragAccelerations
    partials: Optional[AtmosphericDragPartials] = None


class Perturbations:
    """Atmospheric drag perturbations class.

    Handles acceleration and partial derivative computations for atmospheric
    drag dynamics.
    """

    def __init__(self, c: AtmosphericDragInitConfig):
        self.w = OMEGA_EARTH
        self.partials_flag = c.partials_flag
        self.Bstar_flag = c.Bstar_flag
        self.density_model = ExponentialAtmosphericModel(c.partials_flag)

    def pre_compute(self, s: TwoBodyDragState):
        return PreComputed(
            vreli=s.vi + self.w * s.pj,
            vrelj=s.vj - self.w * s.pi,
            vrelk=s.vk,
        )

    def __call__(self, s: TwoBodyDragState) -> PerturbedOutput:
        pc = self.pre_compute(s)
        d = self.density_model(s.pi, s.pj, s.pk, s.p)
        output = PerturbedOutput(
            accelerations=accelerations(
                pc.vreli,
                pc.vrelj,
                pc.vrelk,
                pc.vrel,
                d.rho,
                s.Bstar,
            )
        )
        if self.partials_flag:
            output.partials = partials(
                pc.vreli,
                pc.vrelj,
                pc.vrelk,
                pc.vrel,
                d.rho,
                s.Bstar,
                d.partials.drho_dpi,
                d.partials.drho_dpj,
                d.partials.drho_dpk,
                self.Bstar_flag,
            )
        return output
