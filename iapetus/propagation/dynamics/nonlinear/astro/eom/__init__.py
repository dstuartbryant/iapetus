"""Astrodynamics equations of motion (EOM) module."""

from typing import List, Optional

import numpy as np
from pydantic import BaseModel

from .two_body import TwoBodyAccelerations, TwoBodyPartialDerivatives


class EomAtmosDragParam(BaseModel):
    """Atmospheric drag EOM parameters class."""

    C_d: float  # Satellite coefficient of drag
    area_m2: float  # Satellite surface area [m^2]
    mass_kg: float  # Satellite mass [kg]
    omega_radps: float  # Celestial body rotational rate [radians/s]
    radius_m: float  # Celestial body equatorial radius [m]


class EomParam(BaseModel):
    """Astrodynamics equations of motion parameters class.

    This class exposes callable methods for computing acceleration components
    and partial derivative components for various combinations of dynamics
    modeling options.
    """

    mu_m3ps2: float
    atmos: Optional[EomAtmosDragParam] = None


class Eom:
    """Astrodynamics equations of motion (EOM) class."""

    def __init__(
        self,
        param: EomParam,
        eom_list: List[str] = None,
    ):
        """
        Args:
            param (dict): Parameters dictionary
            eom_list (List[str]): List of EOM dynamics to use, defaults to None
        """
        self.param = param
        if eom_list:
            self.eom_list = eom_list
        else:
            self.eom_list = []
        self._validate_params()
        self.dynamics = self._configure()

    def _validate_params(self):
        raise NotImplementedError()

    def _configure(self):
        dynamics = []
        for item in self.eom_list:
            if item == "two-body":
                dynamics.append(self._configure_two_body())
        return dynamics

    def _configure_two_body(self):
        return TwoBody(self.param.mu_m3ps2, partials_flag=self.partials_flag)


class TwoBody:
    def __init__(self, mu_m3ps2: float, partials_flag: bool):
        self.accelerations = TwoBodyAccelerations(mu_m3ps2)
        self.partials_flag = partials_flag

        if partials_flag:
            self.partials = TwoBodyPartialDerivatives(mu_m3ps2)

    def compute(self, x: float, y: float, z: float):
        r = np.linalg.norm([x, y, z])
        r3 = r**3
        self.accelerations.compute(x, y, z, r3)
        if self.partials_flag:
            r5 = r**5
            self.partials.compute(x, y, z, r3, r5)
