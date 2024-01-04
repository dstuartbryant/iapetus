"""Two body EOM test module."""

import pytest

from iapetus.propagation.dynamics.nonlinear.astro.eom.payloads import (
    TwoBodyInitConfig,
)
from iapetus.propagation.dynamics.nonlinear.astro.eom.two_body import TwoBody


def test_accelerations(state_two_body_without_stm, earth_params):
    ss1 = state_two_body_without_stm
    tbic = TwoBodyInitConfig(mu=earth_params.mu, partials_flag=False)
    mu = earth_params.mu
    tba = TwoBody(tbic)
    output = tba(ss1)
    assert output.accelerations.ai == -mu * ss1.pi / ss1.p3
    assert output.accelerations.aj == -mu * ss1.pj / ss1.p3
    assert output.accelerations.ak == -mu * ss1.pk / ss1.p3


def test_partial_derivatives(state_two_body_with_stm, earth_params):
    ss1 = state_two_body_with_stm
    mu = earth_params.mu
    tbic = TwoBodyInitConfig(mu=earth_params.mu, partials_flag=True)
    tbpd = TwoBody(tbic)
    output = tbpd(ss1)
    assert output.partials.dvi_dpi == 0.0
    assert output.partials.dvi_dpj == 0.0
    assert output.partials.dvi_dpk == 0.0
    assert output.partials.dvi_dvi == 1.0
    assert output.partials.dvi_dvj == 0.0
    assert output.partials.dvi_dvk == 0.0
    assert output.partials.dvj_dpi == 0.0
    assert output.partials.dvj_dpj == 0.0
    assert output.partials.dvj_dpk == 0.0
    assert output.partials.dvj_dvi == 0.0
    assert output.partials.dvj_dvj == 1.0
    assert output.partials.dvj_dvk == 0.0
    assert output.partials.dvk_dpi == 0.0
    assert output.partials.dvk_dpj == 0.0
    assert output.partials.dvk_dpk == 0.0
    assert output.partials.dvk_dvi == 0.0
    assert output.partials.dvk_dvj == 0.0
    assert output.partials.dvk_dvk == 1.0
    assert output.partials.dai_dpi == pytest.approx(
        -mu / ss1.p3 + 3 * mu * ss1.pi**2 / ss1.p5, abs=1e-21
    )
    assert output.partials.dai_dpj == pytest.approx(
        3 * mu * ss1.pi * ss1.pj / ss1.p5, abs=1e-21
    )
    assert output.partials.dai_dpk == pytest.approx(
        3 * mu * ss1.pi * ss1.pk / ss1.p5, abs=1e-21
    )
    assert output.partials.dai_dvi == 0.0
    assert output.partials.dai_dvj == 0.0
    assert output.partials.dai_dvk == 0.0
    assert output.partials.daj_dpi == pytest.approx(
        3 * mu * ss1.pi * ss1.pj / ss1.p5, abs=1e-21
    )
    assert output.partials.daj_dpj == pytest.approx(
        -mu / ss1.p3 + 3 * mu * ss1.pj**2 / ss1.p5, abs=1e-21
    )
    assert output.partials.daj_dpk == pytest.approx(
        3 * mu * ss1.pj * ss1.pk / ss1.p5, abs=1e-21
    )
    assert output.partials.daj_dvi == 0.0
    assert output.partials.daj_dvj == 0.0
    assert output.partials.daj_dvk == 0.0
    assert output.partials.dak_dpi == pytest.approx(
        3 * mu * ss1.pi * ss1.pk / ss1.p5, abs=1e-21
    )
    assert output.partials.dak_dpj == pytest.approx(
        3 * mu * ss1.pj * ss1.pk / ss1.p5, abs=1e-21
    )
    assert output.partials.dak_dpk == pytest.approx(
        -mu / ss1.p3 + 3 * mu * ss1.pk**2 / ss1.p5, abs=1e-21
    )
    assert output.partials.dak_dvi == 0.0
    assert output.partials.dak_dvj == 0.0
    assert output.partials.dak_dvk == 0.0
