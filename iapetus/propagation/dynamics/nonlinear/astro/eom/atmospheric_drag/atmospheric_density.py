"""Atmospheric models for orbital dynamics."""

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro.constants import (
    EQUATORIAL_RADIUS,
)

RE = EQUATORIAL_RADIUS["Earth"]


@dataclass
class ExponentialModelTableEntry:
    """Component class for a single entry in the Exponential Atmospheric Model
    table.
    """

    min_alt_m: float
    max_alt_m: float
    nominal_density_kgpm3: float
    scale_height_m: float
    base_alt_m: float = None

    def __post_init__(self):
        if isinstance(self.base_alt_m, type(None)):
            self.base_alt_m = self.min_alt_m


class ExponentialAtmosphericModelError(Exception):
    pass


@dataclass
class ExponentialAtmosphericModelPartialsOutput:
    drho_dpi: float
    drho_dpj: float
    drho_dpk: float


@dataclass
class ExponentialAtmosphericModelOutput:
    rho: float  # atmospheric density
    partials: Optional[ExponentialAtmosphericModelPartialsOutput] = None


class ExponentialAtmosphericModel:
    """Exponential atmospheric model.

    TODO: Add remainder of tabulated data (only portions for initial motivating
    purposes were included at time of initial writing 2023-11-09).

    Reference:
        Vallado, 4th ed., pg. 567.
    """

    def __init__(
        self, partials_flag: bool
    ) -> List[ExponentialModelTableEntry]:
        """
        Args:
            partials_flag (bool): If True, computes and returns density partial
                derivates with respect to position vector components.
        """
        self.Re = RE
        self.partials_flag = partials_flag
        self.table = [
            ExponentialModelTableEntry(350e3, 400e3, 9.518e-12, 53298),
            ExponentialModelTableEntry(400e3, 450e3, 3.725e-12, 58515),
            ExponentialModelTableEntry(450e3, 500e3, 1.585e-12, 60828),
            ExponentialModelTableEntry(500e3, 600e3, 6.967e-13, 63822),
            ExponentialModelTableEntry(600e3, 700e3, 1.454e-13, 71835),
            ExponentialModelTableEntry(700e3, 800e3, 3.614e-14, 88667),
        ]

    def altitude_lookup(self, altitude_m: float) -> ExponentialModelTableEntry:
        for entry in self.table:
            if altitude_m > entry.min_alt_m and altitude_m < entry.max_alt_m:
                return entry
        raise ExponentialAtmosphericModelError(
            f"Table entry not found for altitude of {altitude_m} meters."
        )

    def density(self, p: float):
        """Computes density at altitude.

        Args:
            p (float): satellite's position magnitude
        """
        altitude_m = p - self.Re
        data = self.altitude_lookup(altitude_m)

        rho = data.nominal_density_kgpm3 * np.exp(
            -(altitude_m - data.base_alt_m) / data.scale_height_m
        )
        return rho, data.scale_height_m

    def partials(
        self,
        pi: float,
        pj: float,
        pk: float,
        p: float,
        rho: float,
        scale_height: float,
    ) -> ExponentialAtmosphericModelPartialsOutput:
        """Computes partial derivatives of density wrt position components.

        Args:
            pi (float): i-th position component of satellite [m]
            pj (float): j-th position component of satellite [m]
            pk (float): k-th position component of satellite [m]
            p (float): satellite's position vector magnitude [m]
            rho (float): atmospheric density [kgpm3]
            scale_height (float): scale height [m]
        """
        coeff = -rho / scale_height
        return ExponentialAtmosphericModelPartialsOutput(
            drho_dpi=coeff * pi / p,
            drho_dpj=coeff * pj / p,
            drho_dpk=coeff * pk / p,
        )

    def __call__(self, pi: float, pj: float, pk: float, p: float):
        rho, H = self.density(p)
        output = ExponentialAtmosphericModelOutput(rho)
        if self.partials_flag:
            output.partials = self.partials(pi, pj, pk, p, rho, H)

        return output
