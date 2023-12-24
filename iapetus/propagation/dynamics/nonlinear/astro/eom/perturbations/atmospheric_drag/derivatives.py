"""Atmospheric drag partial derivatives for equations of motion (EOM)."""
from dataclasses import dataclass
from typing import Optional

from iapetus.propagation.dynamics.nonlinear.astro.constants import ROTATION

OMEGA_EARTH = ROTATION["Earth"]


@dataclass
class AtmosphericDragPartials:
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
    dai_dBstar: Optional[float] = None
    daj_dBstar: Optional[float] = None
    dak_dBstar: Optional[float] = None


def partials(
    vreli: float,
    vrelj: float,
    vrelk: float,
    vrel: float,
    rho: float,
    Bstar: float,
    drho_dpi: float,
    drho_dpj: float,
    drho_dpk: float,
    Bstar_flag: bool,
) -> AtmosphericDragPartials:
    """Computes partial derivatives for atmospheric drag dynamics.

    Args:
        vreli (float): i-th relative velocity component [mps]
        vrelj (float): j-th relative velocity component [mps]
        vrelk (float): k-th relative velocity component [mps]
        vrel (float): relative velocity magnitude [mps]
        rho (float): atmospheric density [kgpm3]
        Bstar (float): B* = 1/2*Cd*A/m [m2pkg]; Cd is drag coefficient
            [unitless], A is area [m2], and m is mass [kg].
        drho_dpi (float): Partial derivative of rho wrt pi
        drho_dpj (float): Partial derivative of rho wrt pj
        drho_dpk (float): Partial derivative of rho wrt pk
        Bstar_flat (bool): if True, additionaly computes partial derivatives
            wrt to Bstar
    """
    output = AtmosphericDragPartials(
        dai_dpi=-Bstar
        * vreli
        * (drho_dpi * vrel - (OMEGA_EARTH * rho * vrelj) / vrel),
        dai_dpj=-Bstar
        * (
            drho_dpj * vrel * vreli
            + OMEGA_EARTH * rho * (vreli**2 / vrel + vrel)
        ),
        dai_dpk=-Bstar * drho_dpk * vrel * vreli,
        dai_dvi=-Bstar * rho * (vreli**2 / vrel + vrel),
        dai_dvj=-Bstar * rho * vreli * vrelj / vrel,
        dai_dvk=-Bstar * rho * vreli * vrelk / vrel,
        daj_dpi=-Bstar
        * (
            drho_dpi * vrel * vrelj
            - OMEGA_EARTH * rho * (vrelj**2 / vrel + vrel)
        ),
        daj_dpj=-Bstar
        * vrelj
        * (drho_dpj * vrel + OMEGA_EARTH * rho * vreli / vrel),
        daj_dpk=-Bstar * drho_dpk * vrel * vrelj,
        daj_dvi=-Bstar * rho * vreli * vrelj / vrel,
        daj_dvj=-Bstar * rho * (vrelj**2 / vrel + vrel),
        daj_dvk=-Bstar * rho * vrelj * vrelk / vrel,
        dak_dpi=-Bstar
        * vrelk
        * (drho_dpi * vrel - OMEGA_EARTH * rho * vrelj / vrel),
        dak_dpj=-Bstar
        * vrelk
        * (drho_dpj * vrel + OMEGA_EARTH * rho * vreli / vrel),
        dak_dpk=-Bstar * drho_dpk * vrel * vrelk,
        dak_dvi=-Bstar * rho * vreli * vrelk / vrel,
        dak_dvj=-Bstar * rho * vrelj * vrelk / vrel,
        dak_dvk=-Bstar * rho * (vrelk**2 / vrel + vrel),
    )
    if Bstar_flag:
        coeff = rho * vrel
        output.dai_dBstar = coeff * vreli
        output.daj_dBstar = coeff * vrelj
        output.dak_dBstar = coeff * vrelk
    return output
