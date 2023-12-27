"""Atmospheric drag accelerations for equations of motion (EOM)."""


from ..models import PerturbedAccelerations


def accelerations(
    vreli: float,
    vrelj: float,
    vrelk: float,
    vrel: float,
    rho: float,
    Bstar: float,
) -> PerturbedAccelerations:
    """Computes component-specified acceleration following the drag force
    model.

    Args:
        vreli (float): i-th relative velocity component [mps]
        vrelj (float): j-th relative velocity component [mps]
        vrelk (float): k-th relative velocity component [mps]
        vrel (float): relative velocity magnitude [mps]
        rho (float): atmospheric density [kgpm3]
        Bstar (float): B* = 1/2*Cd*A/m [m2pkg]; Cd is drag coefficient
            [unitless], A is area [m2], and m is mass [kg].
    """
    prefix = -Bstar * rho * vrel
    ai = prefix * vreli
    aj = prefix * vrelj
    ak = prefix * vrelk

    return PerturbedAccelerations(ai, aj, ak)
