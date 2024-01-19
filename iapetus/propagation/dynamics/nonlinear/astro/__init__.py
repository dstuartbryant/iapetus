"""Astro-dynamics module."""

from typing import List, Optional

from pydantic import BaseModel

from . import state_noise_compensation as snc
from .config import Astrodynamics
from .config import PropagatorInit as AstroPropIniit
from .config import unpack_stm


class AstroInit(BaseModel):
    state_vector_content: List[str]
    celestial_body: str
    stm_flag: bool
    integrator: str
    perturbations: Optional[List[str]] = []


__all__ = ["Astrodynamics", "AstroPropIniit", "AstroInit", "unpack_stm", "scn"]
