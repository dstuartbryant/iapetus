"""Two-body equations of motion module."""

from dataclasses import dataclass
from typing import Optional

from .payloads.call_states import TwoBodyState
from .payloads.configs import TwoBodyInitConfig


@dataclass
class TwoBodyAccelerations:
    ai: float
    aj: float
    ak: float


@dataclass
class TwoBodyPartials:
    dvi_dpi: float = 0.0
    dvi_dpj: float = 0.0
    dvi_dpk: float = 0.0
    dvi_dvi: float = 1.0
    dvi_dvj: float = 0.0
    dvi_dvk: float = 0.0
    dvj_dpi: float = 0.0
    dvj_dpj: float = 0.0
    dvj_dpk: float = 0.0
    dvj_dvi: float = 0.0
    dvj_dvj: float = 1.0
    dvj_dvk: float = 0.0
    dvk_dpi: float = 0.0
    dvk_dpj: float = 0.0
    dvk_dpk: float = 0.0
    dvk_dvi: float = 0.0
    dvk_dvj: float = 0.0
    dvk_dvk: float = 1.0
    dai_dpi: float = 0.0
    dai_dpj: float = 0.0
    dai_dpk: float = 0.0
    dai_dvi: float = 0.0
    dai_dvj: float = 0.0
    dai_dvk: float = 0.0
    daj_dpi: float = 0.0
    daj_dpj: float = 0.0
    daj_dpk: float = 0.0
    daj_dvi: float = 0.0
    daj_dvj: float = 0.0
    daj_dvk: float = 0.0
    dak_dpi: float = 0.0
    dak_dpj: float = 0.0
    dak_dpk: float = 0.0
    dak_dvi: float = 0.0
    dak_dvj: float = 0.0
    dak_dvk: float = 0.0


@dataclass
class TwoBodyOutput:
    accelerations: TwoBodyAccelerations
    partials: Optional[TwoBodyPartials] = None


class TwoBodyEomError(Exception):
    pass


class TwoBody:
    """Two body equations of motion (EOM) class."""

    def __init__(self, c: TwoBodyInitConfig):
        """
        Args:
            mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]
            partials_flag (bool): If True, computes and returns partial
                derivatives
        """
        self.mu = c.mu
        self.partials_flag = c.partials_flag

    def _acceleration(self, p_comp: float, p3: float) -> float:
        return -self.mu * p_comp / p3

    def _da_dp_same_components(
        self, p_comp: float, p3: float, p5: float
    ) -> float:
        return 3 * p_comp**2 * self.mu / p5 - self.mu / p3

    def _da_dp_different_components(
        self, p_comp_1: float, p_comp_2: float, p5: float
    ) -> float:
        return 3 * p_comp_1 * p_comp_2 * self.mu / p5

    def __call__(self, s: TwoBodyState) -> TwoBodyOutput:
        output = TwoBodyOutput(
            accelerations=TwoBodyAccelerations(
                self._acceleration(s.pi, s.p3),
                self._acceleration(s.pj, s.p3),
                self._acceleration(s.pk, s.p3),
            )
        )
        if self.partials_flag:
            if not s.p5:
                raise TwoBodyEomError(
                    "`p5` arg (satellite radius to fifth power) is required "
                    "for partial derivative calculations."
                )
            dai_dpj = self._da_dp_different_components(s.pi, s.pj, s.p5)
            daj_dpk = self._da_dp_different_components(s.pj, s.pk, s.p5)
            dai_dpk = self._da_dp_different_components(s.pi, s.pk, s.p5)
            parts = TwoBodyPartials(
                dai_dpi=self._da_dp_same_components(s.pi, s.p3, s.p5),
                daj_dpj=self._da_dp_same_components(s.pj, s.p3, s.p5),
                dak_dpk=self._da_dp_same_components(s.pk, s.p3, s.p5),
                dai_dpj=dai_dpj,
                daj_dpi=dai_dpj,
                daj_dpk=daj_dpk,
                dak_dpj=daj_dpk,
                dai_dpk=dai_dpk,
                dak_dpi=dai_dpk,
            )
            output.partials = parts

        return output
