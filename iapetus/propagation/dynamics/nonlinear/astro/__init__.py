"""Astro-dynamics module."""

from typing import List, Optional

from pydantic import BaseModel

from .config import Astrodynamics
from .config import PropagatorInit as AstroPropIniit


class AstroInit(BaseModel):
    state_vector_content: List[str]
    celestial_body: str
    stm_flag: bool
    integrator: str
    perturbations: Optional[List[str]] = []
