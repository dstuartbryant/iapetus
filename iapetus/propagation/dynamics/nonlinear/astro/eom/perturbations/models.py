"""Perturbations interface models."""

from dataclasses import dataclass
from typing import Optional


@dataclass
class PerturbedAccelerations:
    """Perturbed accelerations output component base class."""

    ai: float
    aj: float
    ak: float


@dataclass
class PerturbedPartials:
    """Perturbed partials output component base class."""

    dai_dpi: float
    dai_dpj: float
    dai_dpk: float
    dai_dvi: float
    dai_dvj: float
    dai_dvk: float
    daj_dpi: float
    daj_dpj: float
    daj_dpk: float
    daj_dvi: float
    daj_dvj: float
    daj_dvk: float
    dak_dpi: float
    dak_dpj: float
    dak_dpk: float
    dak_dvi: float
    dak_dvj: float
    dak_dvk: float


@dataclass
class PerturbedOutput:
    accelerations: PerturbedAccelerations
    partials: Optional[PerturbedPartials] = None


class Perturbation:
    """Astrodynamics perturbations base class."""

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs) -> PerturbedOutput:
        pass


# ------------ Atmospheric Drag --------------
@dataclass
class AtmosphericDragPartials(PerturbedPartials):
    dai_dBstar: Optional[float] = None
    daj_dBstar: Optional[float] = None
    dak_dBstar: Optional[float] = None
    dai_dCd: Optional[float] = None
    daj_dCd: Optional[float] = None
    dak_dCd: Optional[float] = None


@dataclass
class AtmosphericDragPerturbedOutput(PerturbedOutput):
    partials: Optional[AtmosphericDragPartials] = None
