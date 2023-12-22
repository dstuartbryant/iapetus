"""Two body EOM test module."""

import numpy as np
import pytest

from iapetus.propagation.dynamics.nonlinear.astro.constants import (
    EQUATORIAL_RADIUS,
    MU,
)
from iapetus.propagation.dynamics.nonlinear.astro.eom.two_body import (
    TwoBodyAccelerations,
    TwoBodyPartialDerivatives,
)

MU_EARTH = MU["Earth"]
R_EARTH = EQUATORIAL_RADIUS


def test_accelerations():
    # Position vector
    p = [6656356.11057065, 1700859.15707779, 299734.38071253]
    r = np.linalg.norm(p)
    r3 = r**3
    tba = TwoBodyAccelerations(MU_EARTH)
    tba.compute(p[0], p[1], p[2], r3)
    assert tba.ai == -MU_EARTH * p[0] / r3
    assert tba.aj == -MU_EARTH * p[1] / r3
    assert tba.ak == -MU_EARTH * p[2] / r3


def test_partial_derivatives():
    p = [6656356.11057065, 1700859.15707779, 299734.38071253]
    r = np.linalg.norm(p)
    r3 = r**3
    r5 = r**5
    tbpd = TwoBodyPartialDerivatives(MU_EARTH)
    tbpd.compute(p[0], p[1], p[2], r3, r5)
    assert tbpd.dvi_dpi == 0.0
    assert tbpd.dvi_dpj == 0.0
    assert tbpd.dvi_dpk == 0.0
    assert tbpd.dvi_dvi == 1.0
    assert tbpd.dvi_dvj == 0.0
    assert tbpd.dvi_dvk == 0.0
    assert tbpd.dvj_dpi == 0.0
    assert tbpd.dvj_dpj == 0.0
    assert tbpd.dvj_dpk == 0.0
    assert tbpd.dvj_dvi == 0.0
    assert tbpd.dvj_dvj == 1.0
    assert tbpd.dvj_dvk == 0.0
    assert tbpd.dvk_dpi == 0.0
    assert tbpd.dvk_dpj == 0.0
    assert tbpd.dvk_dpk == 0.0
    assert tbpd.dvk_dvi == 0.0
    assert tbpd.dvk_dvj == 0.0
    assert tbpd.dvk_dvk == 1.0
    assert tbpd.dai_dpi == pytest.approx(
        -MU_EARTH / r3 + 3 * MU_EARTH * p[0] ** 2 / r5, abs=1e-21
    )
    assert tbpd.dai_dpj == pytest.approx(
        3 * MU_EARTH * p[0] * p[1] / r5, abs=1e-21
    )
    assert tbpd.dai_dpk == pytest.approx(
        3 * MU_EARTH * p[0] * p[2] / r5, abs=1e-21
    )
    assert tbpd.dai_dvi == 0.0
    assert tbpd.dai_dvj == 0.0
    assert tbpd.dai_dvk == 0.0
    assert tbpd.daj_dpi == pytest.approx(
        3 * MU_EARTH * p[0] * p[1] / r5, abs=1e-21
    )
    assert tbpd.daj_dpj == pytest.approx(
        -MU_EARTH / r3 + 3 * MU_EARTH * p[1] ** 2 / r5, abs=1e-21
    )
    assert tbpd.daj_dpk == pytest.approx(
        3 * MU_EARTH * p[1] * p[2] / r5, abs=1e-21
    )
    assert tbpd.daj_dvi == 0.0
    assert tbpd.daj_dvj == 0.0
    assert tbpd.daj_dvk == 0.0
    assert tbpd.dak_dpi == pytest.approx(
        3 * MU_EARTH * p[0] * p[2] / r5, abs=1e-21
    )
    assert tbpd.dak_dpj == pytest.approx(
        3 * MU_EARTH * p[1] * p[2] / r5, abs=1e-21
    )
    assert tbpd.dak_dpk == pytest.approx(
        -MU_EARTH / r3 + 3 * MU_EARTH * p[2] ** 2 / r5, abs=1e-21
    )
    assert tbpd.dak_dvi == 0.0
    assert tbpd.dak_dvj == 0.0
    assert tbpd.dak_dvk == 0.0
