"""EOM test module."""

import pytest

from iapetus.propagation.dynamics.nonlinear.astro.eom import (
    EomError,
    configure,
)
from iapetus.propagation.dynamics.nonlinear.astro.eom.payloads import (
    AtmosphericDragInitConfig,
    TwoBodyInitConfig,
)
from iapetus.propagation.dynamics.nonlinear.astro.eom.payloads.call_states import (
    TwoBodyDragState,
    TwoBodyState,
)
from iapetus.propagation.dynamics.nonlinear.astro.eom.perturbations.atmospheric_drag import (
    AtmosphericPerturbation,
)
from iapetus.propagation.dynamics.nonlinear.astro.eom.two_body import TwoBody


class TestTwoBody:
    s = TwoBodyState(
        pi=6656356.11057065,  # [m]
        pj=1700859.15707779,  # [m]
        pk=299734.38071253,  # [m]
        vi=-1794.25660717,  # [mps]
        vj=6353.55570765,  # [mps]
        vk=3792.38315729,  # [mps]
    )

    def test_accelerations_only(self, earth_params):
        mu = earth_params.mu
        ui_config = {"mu": mu, "partials_flag": False}
        eom = configure(ui_config, [])
        tbic = TwoBodyInitConfig(mu=earth_params.mu, partials_flag=False)
        tba = TwoBody(tbic)

        accels, parts = eom(self.s)
        output = tba(self.s)

        assert output.accelerations.ai == accels("ai")
        assert output.accelerations.aj == accels("aj")
        assert output.accelerations.ak == accels("ak")

        assert parts("dvi_dpi") == 0
        assert parts("dvi_dpj") == 0
        assert parts("dvi_dpk") == 0
        assert parts("dvj_dpi") == 0
        assert parts("dvj_dpj") == 0
        assert parts("dvj_dpk") == 0
        assert parts("dvk_dpi") == 0
        assert parts("dvk_dpj") == 0
        assert parts("dvk_dpk") == 0
        assert parts("dai_dpi") == 0
        assert parts("dai_dpj") == 0
        assert parts("dai_dpk") == 0
        assert parts("daj_dpi") == 0
        assert parts("daj_dpj") == 0
        assert parts("daj_dpk") == 0
        assert parts("dak_dpi") == 0
        assert parts("dak_dpj") == 0
        assert parts("dak_dpk") == 0

    def test_partials(self, earth_params):
        mu = earth_params.mu
        ui_config = {"mu": mu, "partials_flag": True}
        eom = configure(ui_config, [])
        tbic = TwoBodyInitConfig(mu=earth_params.mu, partials_flag=True)
        tba = TwoBody(tbic)

        accels, parts = eom(self.s)
        output = tba(self.s)

        assert parts("dvi_dpi") == output.partials.dvi_dpi
        assert parts("dvi_dpj") == output.partials.dvi_dpj
        assert parts("dvi_dpk") == output.partials.dvi_dpk
        assert parts("dvj_dpi") == output.partials.dvj_dpi
        assert parts("dvj_dpj") == output.partials.dvj_dpj
        assert parts("dvj_dpk") == output.partials.dvj_dpk
        assert parts("dvk_dpi") == output.partials.dvk_dpi
        assert parts("dvk_dpj") == output.partials.dvk_dpj
        assert parts("dvk_dpk") == output.partials.dvk_dpk
        assert parts("dai_dpi") == output.partials.dai_dpi
        assert parts("dai_dpj") == output.partials.dai_dpj
        assert parts("dai_dpk") == output.partials.dai_dpk
        assert parts("daj_dpi") == output.partials.daj_dpi
        assert parts("daj_dpj") == output.partials.daj_dpj
        assert parts("daj_dpk") == output.partials.daj_dpk
        assert parts("dak_dpi") == output.partials.dak_dpi
        assert parts("dak_dpj") == output.partials.dak_dpj
        assert parts("dak_dpk") == output.partials.dak_dpk


class TestTwoBodyDrag:
    s_2body = TwoBodyState(
        pi=6656356.11057065,  # [m]
        pj=1700859.15707779,  # [m]
        pk=299734.38071253,  # [m]
        vi=-1794.25660717,  # [mps]
        vj=6353.55570765,  # [mps]
        vk=3792.38315729,  # [mps]
    )

    s_drag = TwoBodyDragState(
        pi=6656356.11057065,  # [m]
        pj=1700859.15707779,  # [m]
        pk=299734.38071253,  # [m]
        vi=-1794.25660717,  # [mps]
        vj=6353.55570765,  # [mps]
        vk=3792.38315729,  # [mps]
        A=(1 / 3) ** 2,  # [m^2]
        m=50,  # [kg]
        Cd=2.0,  # [unitless]
    )

    def test_accelerations_only(self, earth_params):
        mu = earth_params.mu
        ui_config = {"mu": mu, "partials_flag": False, "Bstar_flag": False}
        eom = configure(ui_config, ["atmospheric-drag"])
        tbic = TwoBodyInitConfig(mu=earth_params.mu, partials_flag=False)
        tba = TwoBody(tbic)
        adic = AtmosphericDragInitConfig(partials_flag=False, Bstar_flag=False)
        ap = AtmosphericPerturbation(adic)

        accels, parts = eom(self.s_drag)
        output1 = tba(self.s_2body)
        output2 = ap(self.s_drag)

        assert output1.accelerations.ai + output2.accelerations.ai == accels(
            "ai"
        )
        assert output1.accelerations.aj + output2.accelerations.aj == accels(
            "aj"
        )
        assert output1.accelerations.ak + output2.accelerations.ak == accels(
            "ak"
        )

        assert parts("dvi_dpi") == 0
        assert parts("dvi_dpj") == 0
        assert parts("dvi_dpk") == 0
        assert parts("dvj_dpi") == 0
        assert parts("dvj_dpj") == 0
        assert parts("dvj_dpk") == 0
        assert parts("dvk_dpi") == 0
        assert parts("dvk_dpj") == 0
        assert parts("dvk_dpk") == 0
        assert parts("dai_dpi") == 0
        assert parts("dai_dpj") == 0
        assert parts("dai_dpk") == 0
        assert parts("daj_dpi") == 0
        assert parts("daj_dpj") == 0
        assert parts("daj_dpk") == 0
        assert parts("dak_dpi") == 0
        assert parts("dak_dpj") == 0
        assert parts("dak_dpk") == 0

    def test_partials_no_bstar(self, earth_params):
        mu = earth_params.mu
        ui_config = {"mu": mu, "partials_flag": True, "Bstar_flag": False}
        eom = configure(ui_config, ["atmospheric-drag"])
        tbic = TwoBodyInitConfig(mu=earth_params.mu, partials_flag=True)
        tba = TwoBody(tbic)
        adic = AtmosphericDragInitConfig(partials_flag=True, Bstar_flag=False)
        ap = AtmosphericPerturbation(adic)

        accels, parts = eom(self.s_drag)
        output1 = tba(self.s_2body)
        output2 = ap(self.s_drag)

        partials1 = output1.partials
        partials2 = output2.partials

        assert parts("dvi_dpi") == partials1.dvi_dpi
        assert parts("dvi_dpj") == partials1.dvi_dpj
        assert parts("dvi_dpk") == partials1.dvi_dpk
        assert parts("dvj_dpi") == partials1.dvj_dpi
        assert parts("dvj_dpj") == partials1.dvj_dpj
        assert parts("dvj_dpk") == partials1.dvj_dpk
        assert parts("dvk_dpi") == partials1.dvk_dpi
        assert parts("dvk_dpj") == partials1.dvk_dpj
        assert parts("dvk_dpk") == partials1.dvk_dpk
        assert parts("dai_dpi") == partials1.dai_dpi + partials2.dai_dpi
        assert parts("dai_dpj") == partials1.dai_dpj + partials2.dai_dpj
        assert parts("dai_dpk") == partials1.dai_dpk + partials2.dai_dpk
        assert parts("daj_dpi") == partials1.daj_dpi + partials2.daj_dpi
        assert parts("daj_dpj") == partials1.daj_dpj + partials2.daj_dpj
        assert parts("daj_dpk") == partials1.daj_dpk + partials2.daj_dpk
        assert parts("dak_dpi") == partials1.dak_dpi + partials2.dak_dpi
        assert parts("dak_dpj") == partials1.dak_dpj + partials2.dak_dpj
        assert parts("dak_dpk") == partials1.dak_dpk + partials2.dak_dpk

        assert isinstance(partials2.dai_dBstar, type(None))
        assert isinstance(partials2.daj_dBstar, type(None))
        assert isinstance(partials2.dak_dBstar, type(None))

        with pytest.raises(EomError) as e:
            parts("dai_dBstar")
        assert (
            str(e.value.args[0])
            == "Partials component dai_dBstar: Attempting to add to NoneType"
        )

        with pytest.raises(EomError) as e:
            parts("daj_dBstar")
        assert (
            str(e.value.args[0])
            == "Partials component daj_dBstar: Attempting to add to NoneType"
        )

        with pytest.raises(EomError) as e:
            parts("dak_dBstar")
        assert (
            str(e.value.args[0])
            == "Partials component dak_dBstar: Attempting to add to NoneType"
        )
