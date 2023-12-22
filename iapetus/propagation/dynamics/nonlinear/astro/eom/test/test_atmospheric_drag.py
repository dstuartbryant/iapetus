"""Atmospheric drag EOM test module."""


import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro.constants import ROTATION
from iapetus.propagation.dynamics.nonlinear.astro.eom.atmospheric_drag import (
    Perturbations,
)
from iapetus.propagation.dynamics.nonlinear.astro.eom.atmospheric_drag.atmospheric_density import (
    ExponentialAtmosphericModel,
)

SAT_POSITION = [6656356.11057065, 1700859.15707779, 299734.38071253]
SAT_VELOCITY = [-1794.25660717, 6353.55570765, 3792.38315729]
AREA = 1 / 3**2  # m^2
MASS = 50  # kg
CD = 2.0
OMEGA = ROTATION["Earth"]


def test_accelerations():
    p = SAT_POSITION
    v = SAT_VELOCITY
    A = AREA
    m = MASS
    Cd = CD
    w = OMEGA
    pmag = np.linalg.norm(p)
    Bstar = 0.5 * Cd * A / m
    density_model = ExponentialAtmosphericModel(partials_flag=False)
    rho = density_model(p[0], p[1], p[2], pmag).rho

    pert = Perturbations(partials_flag=False, Bstar_flag=False)
    pout = pert(p[0], p[1], p[2], pmag, v[0], v[1], v[2], A, Cd, m)

    vreli = v[0] + w * p[1]
    vrelj = v[1] - w * p[0]
    vrelk = v[2]
    vrel = np.linalg.norm([vreli, vrelj, vrelk])

    ai = -Bstar * rho * vrel * vreli
    aj = -Bstar * rho * vrel * vrelj
    ak = -Bstar * rho * vrel * vrelk

    assert ai == pout.accelrations.ai
    assert aj == pout.accelrations.aj
    assert ak == pout.accelrations.ak


def test_partials_no_Bstar():
    p = SAT_POSITION
    v = SAT_VELOCITY
    A = AREA
    m = MASS
    Cd = CD
    w = OMEGA
    pmag = np.linalg.norm(p)
    Bstar = 0.5 * Cd * A / m
    density_model = ExponentialAtmosphericModel(partials_flag=True)
    d = density_model(p[0], p[1], p[2], pmag)
    rho = d.rho
    drho_dpi = d.partials.drho_dpi
    drho_dpj = d.partials.drho_dpj
    drho_dpk = d.partials.drho_dpk

    pert = Perturbations(partials_flag=True, Bstar_flag=False)
    pout = pert(p[0], p[1], p[2], pmag, v[0], v[1], v[2], A, Cd, m)

    vreli = v[0] + w * p[1]
    vrelj = v[1] - w * p[0]
    vrelk = v[2]
    vrel = np.linalg.norm([vreli, vrelj, vrelk])

    assert pout.partials.dai_dpi == -Bstar * vreli * (
        drho_dpi * vrel - w * rho * vrelj / vrel
    )
    assert pout.partials.dai_dpj == -Bstar * (
        drho_dpj * vrel * vreli + w * rho * (vreli**2 / vrel + vrel)
    )
    assert pout.partials.dai_dpk == -Bstar * drho_dpk * vrel * vreli
    assert pout.partials.dai_dvi == -Bstar * rho * (vreli**2 / vrel + vrel)
    assert pout.partials.dai_dvj == -Bstar * rho * vreli * vrelj / vrel
    assert pout.partials.dai_dvk == -Bstar * rho * vreli * vrelk / vrel
    assert pout.partials.daj_dpi == -Bstar * (
        drho_dpi * vrel * vrelj - w * rho * (vrelj**2 / vrel + vrel)
    )
    assert pout.partials.daj_dpj == -Bstar * vrelj * (
        drho_dpj * vrel + w * rho * vreli / vrel
    )
    assert pout.partials.daj_dpk == -Bstar * drho_dpk * vrel * vrelj
    assert pout.partials.daj_dvi == -Bstar * rho * vreli * vrelj / vrel
    assert pout.partials.daj_dvj == -Bstar * rho * (vrelj**2 / vrel + vrel)
    assert pout.partials.daj_dvk == -Bstar * rho * vrelj * vrelk / vrel
    assert pout.partials.dak_dpi == -Bstar * vrelk * (
        drho_dpi * vrel - w * rho * vrelj / vrel
    )
    assert pout.partials.dak_dpj == -Bstar * vrelk * (
        drho_dpj * vrel + w * rho * vreli / vrel
    )
    assert pout.partials.dak_dpk == -Bstar * drho_dpk * vrel * vrelk
    assert pout.partials.dak_dvi == -Bstar * rho * vreli * vrelk / vrel
    assert pout.partials.dak_dvj == -Bstar * rho * vrelj * vrelk / vrel
    assert pout.partials.dak_dvk == -Bstar * rho * (vrelk**2 / vrel + vrel)
