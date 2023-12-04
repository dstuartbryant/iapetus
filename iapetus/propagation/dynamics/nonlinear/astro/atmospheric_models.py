"""Atmospheric models for orbital dynamics."""

from dataclasses import dataclass
from typing import List


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


class ExponentialAtmosphericModel:
    """Exponential atmospheric model.

    TODO: Add remainder of tabulated data (only portions for initial motivating
    purposes were included at time of initial writing 2023-11-09).

    Reference:
        Vallado, 4th ed., pg. 567.
    """

    def __init__(self) -> List[ExponentialModelTableEntry]:
        self.table = [
            ExponentialModelTableEntry(350e3, 400e3, 9.518e-12, 53298),
            ExponentialModelTableEntry(400e3, 450e3, 3.725e-12, 58515),
            ExponentialModelTableEntry(450e3, 500e3, 1.585e-12, 60828),
            ExponentialModelTableEntry(500e3, 600e3, 6.967e-13, 63822),
            ExponentialModelTableEntry(600e3, 700e3, 1.454e-13, 71835),
            ExponentialModelTableEntry(700e3, 800e3, 3.614e-14, 88667),
        ]

    def altitude_lookup(self, altitude_m: float):
        for entry in self.table:
            if altitude_m > entry.min_alt_m and altitude_m < entry.max_alt_m:
                return entry
        raise ExponentialAtmosphericModelError(
            f"Table entry not found for altitude of {altitude_m} meters."
        )

    def __call__(self, altitude_m):
        return self.altitude_lookup(altitude_m)
