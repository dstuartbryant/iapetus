"""Atmospheric drag EOM test module."""


import numpy as np

from . import ExponentialAtmosphericModel, Perturbations


def test_accelerations(sat_state_1, earth_params):
    ss1 = sat_state_1
    ep = earth_params
    w = ep.rotation_rate
    Bstar = 0.5 * ss1.Cd * ss1.area / ss1.mass
    density_model = ExponentialAtmosphericModel(partials_flag=False)
    rho = density_model(ss1.pi, ss1.pj, ss1.pk, ss1.p_mag).rho

    pert = Perturbations(partials_flag=False, Bstar_flag=False)
    pout = pert(
        ss1.pi,
        ss1.pj,
        ss1.pk,
        ss1.p_mag,
        ss1.vi,
        ss1.vj,
        ss1.vk,
        ss1.area,
        ss1.Cd,
        ss1.mass,
    )

    vreli = ss1.vi + w * ss1.pj
    vrelj = ss1.vj - w * ss1.pi
    vrelk = ss1.vk
    vrel = np.linalg.norm([vreli, vrelj, vrelk])

    ai = -Bstar * rho * vrel * vreli
    aj = -Bstar * rho * vrel * vrelj
    ak = -Bstar * rho * vrel * vrelk

    assert ai == pout.accelerations.ai
    assert aj == pout.accelerations.aj
    assert ak == pout.accelerations.ak


def test_partials_no_Bstar(sat_state_1, earth_params):
    ss1 = sat_state_1
    ep = earth_params
    w = ep.rotation_rate
    Bstar = 0.5 * ss1.Cd * ss1.area / ss1.mass
    density_model = ExponentialAtmosphericModel(partials_flag=True)
    d = density_model(ss1.pi, ss1.pj, ss1.pk, ss1.p_mag)
    rho = d.rho
    drho_dpi = d.partials.drho_dpi
    drho_dpj = d.partials.drho_dpj
    drho_dpk = d.partials.drho_dpk

    pert = Perturbations(partials_flag=True, Bstar_flag=False)
    pout = pert(
        ss1.pi,
        ss1.pj,
        ss1.pk,
        ss1.p_mag,
        ss1.vi,
        ss1.vj,
        ss1.vk,
        ss1.area,
        ss1.Cd,
        ss1.mass,
    )

    vreli = ss1.vi + w * ss1.pj
    vrelj = ss1.vj - w * ss1.pi
    vrelk = ss1.vk
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
