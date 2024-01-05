"""Atmospheric drag equations of motion (EOM) module.
"""

from dataclasses import dataclass

import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro.constants import ROTATION

from ...payloads import AtmosphericDragInitConfig
from ...state import STATE_TYPES
from ..models import AtmosphericDragPerturbedOutput, Perturbation
from .accelerations import accelerations
from .atmospheric_density import ExponentialAtmosphericModel
from .derivatives import partials

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


class AtmosphericPerturbation(Perturbation):
    """Atmospheric drag perturbations class.

    Handles acceleration and partial derivative computations for atmospheric
    drag dynamics.
    """

    def __init__(self, c: AtmosphericDragInitConfig):
        self.w = OMEGA_EARTH
        self.partials_flag = c.partials_flag
        self.Bstar_flag = c.Bstar_flag
        self.Cd_flag = c.Cd_flag
        self.density_model = ExponentialAtmosphericModel(c.partials_flag)

    def pre_compute(self, s: STATE_TYPES):
        return PreComputed(
            vreli=s.vi + self.w * s.pj,
            vrelj=s.vj - self.w * s.pi,
            vrelk=s.vk,
        )

    def __call__(self, s: STATE_TYPES) -> AtmosphericDragPerturbedOutput:
        pc = self.pre_compute(s)
        d = self.density_model(s.pi, s.pj, s.pk, s.p)
        output = AtmosphericDragPerturbedOutput(
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
            partials_kwargs = {
                "vreli": pc.vreli,
                "vrelj": pc.vrelj,
                "vrelk": pc.vrelk,
                "vrel": pc.vrel,
                "rho": d.rho,
                "Bstar": s.Bstar,
                "drho_dpi": d.partials.drho_dpi,
                "drho_dpj": d.partials.drho_dpj,
                "drho_dpk": d.partials.drho_dpk,
                "Bstar_flag": self.Bstar_flag,
                "Cd_flag": self.Cd_flag,
            }
            if self.Cd_flag:
                partials_kwargs["A"] = s.A
                partials_kwargs["m"] = s.m
            output.partials = partials(**partials_kwargs)
        return output
