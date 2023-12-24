"""Two-body equations of motion module."""

from dataclasses import dataclass
from typing import Optional


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
    partials: TwoBodyPartials


class TwoBodyEomError(Exception):
    pass


class TwoBody:
    """Two body equations of motion (EOM) class."""

    def __init__(self, mu_m3ps2: float, partials_flag: bool):
        """
        Args:
            mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]
            partials_flag (bool): If True, computes and returns partial
                derivatives
        """
        self.mu = mu_m3ps2
        self.partials_flag = partials_flag

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

    def __call__(
        self, pi: float, pj: float, pk: float, p3: float, p5: Optional[float]
    ) -> TwoBodyOutput:
        output = TwoBodyOutput(
            accelerations=TwoBodyAccelerations(
                self._acceleration(pi, p3),
                self._acceleration(pj, p3),
                self._acceleration(pk, p3),
            )
        )
        if self.partials_flag:
            dai_dpj = self._da_dp_different_components(pi, pj, p5)
            daj_dpi = dai_dpj
            daj_dpk = self._da_dp_different_components(pj, pk, p5)
            dak_dpj = daj_dpk
            dai_dpk = self._da_dp_different_components(pi, pk, p5)
            dak_dpi = dai_dpk

            if not p5:
                raise TwoBodyEomError(
                    "`p5` arg (satellite radius to fifth power) is required "
                    "for partial derivative calculations."
                )
            output.partials.dai_dpi = self._da_dp_same_components(pi, p3, p5)
            output.partials.daj_dpj = self._da_dp_same_components(pj, p3, p5)
            output.partials.dak_dpk = self._da_dp_same_components(pk, p3, p5)
            output.partials.dai_dpj = dai_dpj
            output.partials.daj_dpi = daj_dpi
            output.partials.daj_dpk = daj_dpk
            output.partials.dak_dpj = dak_dpj
            output.partials.dai_dpk = dai_dpk
            output.partials.dak_dpi = dak_dpi
