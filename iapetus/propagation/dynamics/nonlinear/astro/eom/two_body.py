"""Two-body equations of motion module."""


class TwoBodyAccelerations:
    """Computes accelerations under two-body motion."""

    def __init__(self, mu_m3ps2: float):
        """
        Args:
            mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]
        """
        self.mu = mu_m3ps2
        self.ai = 0
        self.aj = 0
        self.ak = 0

    def _acceleration(self, p: float, r3: float) -> float:
        return -self.mu * p / r3

    def compute(self, x: float, y: float, z: float, r3: float):
        self.ai = self._acceleration(x, r3)
        self.aj = self._acceleration(y, r3)
        self.ak = self._acceleration(z, r3)


class TwoBodyPartialDerivatives:
    """Computes partial derivative elements under two-body motion."""

    def __init__(self, mu_m3ps2: float):
        """
        Args:
            mu_m3ps2 (float): graviational constant in meters cubed per seconds
            squared [m3ps2]
        """
        self.mu = mu_m3ps2
        self.dvi_dpi = 0.0
        self.dvi_dpj = 0.0
        self.dvi_dpk = 0.0
        self.dvi_dvi = 1.0
        self.dvi_dvj = 0.0
        self.dvi_dvk = 0.0
        self.dvj_dpi = 0.0
        self.dvj_dpj = 0.0
        self.dvj_dpk = 0.0
        self.dvj_dvi = 0.0
        self.dvj_dvj = 1.0
        self.dvj_dvk = 0.0
        self.dvk_dpi = 0.0
        self.dvk_dpj = 0.0
        self.dvk_dpk = 0.0
        self.dvk_dvi = 0.0
        self.dvk_dvj = 0.0
        self.dvk_dvk = 1.0
        self.dai_dpi = 0.0
        self.dai_dpj = 0.0
        self.dai_dpk = 0.0
        self.dai_dvi = 0.0
        self.dai_dvj = 0.0
        self.dai_dvk = 0.0
        self.daj_dpi = 0.0
        self.daj_dpj = 0.0
        self.daj_dpk = 0.0
        self.daj_dvi = 0.0
        self.daj_dvj = 0.0
        self.daj_dvk = 0.0
        self.dak_dpi = 0.0
        self.dak_dpj = 0.0
        self.dak_dpk = 0.0
        self.dak_dvi = 0.0
        self.dak_dvj = 0.0
        self.dak_dvk = 0.0

    def _da_dp_same_components(self, p: float, r3: float, r5: float) -> float:
        return 3 * p**2 * self.mu / r5 - self.mu / r3

    def _da_dp_different_components(
        self, p1: float, p2: float, r5: float
    ) -> float:
        return 3 * p1 * p2 * self.mu / r5

    def compute(self, x: float, y: float, z: float, r3: float, r5: float):
        self.dai_dpi = self._da_dp_same_components(x, r3, r5)
        self.daj_dpj = self._da_dp_same_components(y, r3, r5)
        self.dak_dpk = self._da_dp_same_components(z, r3, r5)

        self.dai_dpj = self._da_dp_different_components(x, y, r5)
        self.daj_dpi = self.dai_dpj

        self.daj_dpk = self._da_dp_different_components(y, z, r5)
        self.dak_dpj = self.daj_dpk

        self.dai_dpk = self._da_dp_different_components(x, z, r5)
        self.dak_dpi = self.dai_dpk
