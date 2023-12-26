"""Equations of motion (EOM) __init__ configuration interface models module."""


from dataclasses import dataclass


@dataclass
class EomPartialsConfig:
    partials_flag: bool  # Indicates whether partial derivatives are calculated


@dataclass
class TwoBodyInitConfig(EomPartialsConfig):
    mu: float  # celestial body gravitational constant [m3ps2]


@dataclass
class PerturbationInitConfig(EomPartialsConfig):
    """Base class for perturbation configuration parameters."""

    pass


@dataclass
class AtmosphericDragInitConfig(PerturbationInitConfig):
    Bstar_flag: bool  # Indicates whether to additionaly compute Bstar partials
