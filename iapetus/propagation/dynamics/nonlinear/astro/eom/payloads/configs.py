"""Equations of motion (EOM) __init__ configuration interface models module."""


from dataclasses import dataclass


@dataclass
class EomPartialsConfig:
    partials_flag: bool  # Indicates whether partial derivatives are calculated


@dataclass
class TwoBodyInitConfig(EomPartialsConfig):
    """Configuration class for initializing Two Body EOM."""

    mu: float  # celestial body gravitational constant [m3ps2]


@dataclass
class PerturbationInitConfig(EomPartialsConfig):
    """Base class for perturbation configuration parameters."""

    pass


@dataclass
class NonSphericalInitConfig(PerturbationInitConfig):
    """Configuration class for initializing non-spherical gravity
    perturbations.
    """

    pass


@dataclass
class ThirdBodyInitConfig(PerturbationInitConfig):
    """Configuration class for initializing third-body gravity
    perturbations.
    """

    pass


@dataclass
class GeneralRelativityInitConfig(PerturbationInitConfig):
    """Configuration class for initializing general relativity
    perturbations.
    """

    pass


@dataclass
class AtmosphericDragInitConfig(PerturbationInitConfig):
    """Configuration class for initializing Atmospheric Drag perturbations."""

    Bstar_flag: bool  # Indicates whether to additionaly compute Bstar partials


@dataclass
class SrpInitConfig(PerturbationInitConfig):
    """Configuration class for initializing solar radiation pressure
    perturbations.
    """

    pass


@dataclass
class ErpInitConfig(PerturbationInitConfig):
    """Configuration class for initializing Earth radiation pressure
    perturbations.
    """

    pass
