"""Atmospheric drag EOM test module."""


import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro.eom.payloads import (
    AtmosphericDragInitConfig,
)

from . import AtmosphericPerturbation, ExponentialAtmosphericModel


def test_accelerations(state_two_body_drag_without_stm, earth_params):
    ss1 = state_two_body_drag_without_stm
    ep = earth_params
    w = ep.rotation_rate
    density_model = ExponentialAtmosphericModel(partials_flag=False)
    rho = density_model(ss1.pi, ss1.pj, ss1.pk, ss1.p).rho

    adic = AtmosphericDragInitConfig(
        partials_flag=False, Bstar_flag=False, Cd_flag=False
    )
    pert = AtmosphericPerturbation(adic)
    pout = pert(ss1)

    vreli = ss1.vi + w * ss1.pj
    vrelj = ss1.vj - w * ss1.pi
    vrelk = ss1.vk
    vrel = np.linalg.norm([vreli, vrelj, vrelk])

    ai = -ss1.Bstar * rho * vrel * vreli
    aj = -ss1.Bstar * rho * vrel * vrelj
    ak = -ss1.Bstar * rho * vrel * vrelk

    assert ai == pout.accelerations.ai
    assert aj == pout.accelerations.aj
    assert ak == pout.accelerations.ak


def test_partials_no_Bstar_no_Cd(state_two_body_drag_with_stm, earth_params):
    ss1 = state_two_body_drag_with_stm
    ep = earth_params
    w = ep.rotation_rate
    density_model = ExponentialAtmosphericModel(partials_flag=True)
    d = density_model(ss1.pi, ss1.pj, ss1.pk, ss1.p)
    rho = d.rho
    drho_dpi = d.partials.drho_dpi
    drho_dpj = d.partials.drho_dpj
    drho_dpk = d.partials.drho_dpk
    adic = AtmosphericDragInitConfig(
        partials_flag=True, Bstar_flag=False, Cd_flag=False
    )
    pert = AtmosphericPerturbation(adic)
    pout = pert(ss1)

    vreli = ss1.vi + w * ss1.pj
    vrelj = ss1.vj - w * ss1.pi
    vrelk = ss1.vk
    vrel = np.linalg.norm([vreli, vrelj, vrelk])

    assert pout.partials.dai_dpi == -ss1.Bstar * vreli * (
        drho_dpi * vrel - w * rho * vrelj / vrel
    )
    assert pout.partials.dai_dpj == -ss1.Bstar * (
        drho_dpj * vrel * vreli + w * rho * (vreli**2 / vrel + vrel)
    )
    assert pout.partials.dai_dpk == -ss1.Bstar * drho_dpk * vrel * vreli
    assert pout.partials.dai_dvi == -ss1.Bstar * rho * (
        vreli**2 / vrel + vrel
    )
    assert pout.partials.dai_dvj == -ss1.Bstar * rho * vreli * vrelj / vrel
    assert pout.partials.dai_dvk == -ss1.Bstar * rho * vreli * vrelk / vrel
    assert pout.partials.daj_dpi == -ss1.Bstar * (
        drho_dpi * vrel * vrelj - w * rho * (vrelj**2 / vrel + vrel)
    )
    assert pout.partials.daj_dpj == -ss1.Bstar * vrelj * (
        drho_dpj * vrel + w * rho * vreli / vrel
    )
    assert pout.partials.daj_dpk == -ss1.Bstar * drho_dpk * vrel * vrelj
    assert pout.partials.daj_dvi == -ss1.Bstar * rho * vreli * vrelj / vrel
    assert pout.partials.daj_dvj == -ss1.Bstar * rho * (
        vrelj**2 / vrel + vrel
    )
    assert pout.partials.daj_dvk == -ss1.Bstar * rho * vrelj * vrelk / vrel
    assert pout.partials.dak_dpi == -ss1.Bstar * vrelk * (
        drho_dpi * vrel - w * rho * vrelj / vrel
    )
    assert pout.partials.dak_dpj == -ss1.Bstar * vrelk * (
        drho_dpj * vrel + w * rho * vreli / vrel
    )
    assert pout.partials.dak_dpk == -ss1.Bstar * drho_dpk * vrel * vrelk
    assert pout.partials.dak_dvi == -ss1.Bstar * rho * vreli * vrelk / vrel
    assert pout.partials.dak_dvj == -ss1.Bstar * rho * vrelj * vrelk / vrel
    assert pout.partials.dak_dvk == -ss1.Bstar * rho * (
        vrelk**2 / vrel + vrel
    )

    assert isinstance(pout.partials.dai_dCd, type(None))
    assert isinstance(pout.partials.daj_dCd, type(None))
    assert isinstance(pout.partials.dak_dCd, type(None))

    assert isinstance(pout.partials.dai_dBstar, type(None))
    assert isinstance(pout.partials.daj_dBstar, type(None))
    assert isinstance(pout.partials.dak_dBstar, type(None))


def test_partials_Cd(state_two_body_drag_cd_with_stm, earth_params):
    """Tests atmospheric-drag partials when Cd in state vector."""

    ss1 = state_two_body_drag_cd_with_stm
    ep = earth_params
    w = ep.rotation_rate
    density_model = ExponentialAtmosphericModel(partials_flag=True)
    d = density_model(ss1.pi, ss1.pj, ss1.pk, ss1.p)
    rho = d.rho
    adic = AtmosphericDragInitConfig(
        partials_flag=True, Bstar_flag=False, Cd_flag=True
    )
    pert = AtmosphericPerturbation(adic)
    pout = pert(ss1)

    vreli = ss1.vi + w * ss1.pj
    vrelj = ss1.vj - w * ss1.pi
    vrelk = ss1.vk
    vrel = np.linalg.norm([vreli, vrelj, vrelk])

    assert pout.partials.dai_dCd == -0.5 * ss1.A / ss1.m * rho * vrel * vreli
    assert pout.partials.daj_dCd == -0.5 * ss1.A / ss1.m * rho * vrel * vrelj
    assert pout.partials.dak_dCd == -0.5 * ss1.A / ss1.m * rho * vrel * vrelk

    assert isinstance(pout.partials.dai_dBstar, type(None))
    assert isinstance(pout.partials.daj_dBstar, type(None))
    assert isinstance(pout.partials.dak_dBstar, type(None))


def test_partials_Bstar(state_two_body_drag_bstar_with_stm, earth_params):
    """Tests atmospheric-drag partials when Cd in state vector."""

    ss1 = state_two_body_drag_bstar_with_stm
    ep = earth_params
    w = ep.rotation_rate
    density_model = ExponentialAtmosphericModel(partials_flag=True)
    d = density_model(ss1.pi, ss1.pj, ss1.pk, ss1.p)
    rho = d.rho
    adic = AtmosphericDragInitConfig(
        partials_flag=True, Bstar_flag=True, Cd_flag=False
    )
    pert = AtmosphericPerturbation(adic)
    pout = pert(ss1)

    vreli = ss1.vi + w * ss1.pj
    vrelj = ss1.vj - w * ss1.pi
    vrelk = ss1.vk
    vrel = np.linalg.norm([vreli, vrelj, vrelk])

    assert pout.partials.dai_dBstar == -rho * vrel * vreli
    assert pout.partials.daj_dBstar == -rho * vrel * vrelj
    assert pout.partials.dak_dBstar == -rho * vrel * vrelk

    assert isinstance(pout.partials.dai_dCd, type(None))
    assert isinstance(pout.partials.daj_dCd, type(None))
    assert isinstance(pout.partials.dak_dCd, type(None))
