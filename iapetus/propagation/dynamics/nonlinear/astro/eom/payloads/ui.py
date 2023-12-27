"""User interface module for equations of motion (EOM) configuration.

Perturbations Options:
* Non-Spherical (Earth)
* ThirdBody
* GeneralRelativity
* AtmosphericDrag
* Srp
* Erp (Earth Radiation Pressure)
* ...

"""
from typing import List

from pydantic import BaseModel

from .configs import AtmosphericDragInitConfig, TwoBodyInitConfig

PERTURBATION_NAMES = [
    "non-spherical",
    "third-body",
    "general-relativity",
    "atmospheric-drag",
    "solar-radiation-pressure",
    "earth-radiation-pressure",
]


class EomUiError(Exception):
    pass


class TwoBody(BaseModel):
    mu: float
    partials_flag: bool

    @property
    def two_body_init_config(self):
        return TwoBodyInitConfig(mu=self.mu, partials_flag=self.partials_flag)


class TwoBodyWithPerts(TwoBody):
    def _get_perturbations_init_config_names(self):
        names = []
        for a in dir(self):
            if a != "two_body_init_config" and "_init_config" in a:
                names.append(a)
        return names

    @property
    def perturbations_init_configs(self):
        pert_names = self._get_perturbations_init_config_names()
        # No, I don't think this is best approach for this, can do better.


class TwoBodyNonSphere(TwoBody):
    pass


class TwoBodyThirdBody(TwoBody):
    pass


class TwoBodyGenRel(TwoBody):
    pass


class TwoBodyAtmosDrag(TwoBody):
    Bstar_flag: bool

    @property
    def atmos_drag_init_config(self):
        return AtmosphericDragInitConfig(
            partials_flag=self.partials_flag, Bstar_flag=self.Bstar_flag
        )


class TwoBodySrp(TwoBody):
    pass


class TwoBodyErp(TwoBody):
    pass


class TwoBodyNonSphereThirdBody(BaseModel):
    pass


class TwoBodyNonSphereThirdBodyGenRel(BaseModel):
    pass


class TwoBodyNonSphereThirdBodyGenRelAtmosDrag(BaseModel):
    pass


class TwoBodyNonSphereThirdBodyGenRelAtmosDragSrp(BaseModel):
    pass


class TwoBodyNonSphereThirdBodyGenRelAtmosDragSrpErp(BaseModel):
    pass


ui_map = [
    {"perturbations": [], "ui": TwoBody},
    {"perturbations": ["non-spherical"], "ui": TwoBodyNonSphere},
    {"perturbations": ["third-body"], "ui": TwoBodyThirdBody},
    {"perturbations": ["general-relativity"], "ui": TwoBodyGenRel},
    {"perturbations": ["atmospheric-drag"], "ui": TwoBodyAtmosDrag},
    {"perturbations": ["solar-radiation-pressure"], "ui": TwoBodySrp},
    {"perturbations": ["earth-radiation-pressure"], "ui": TwoBodyErp},
    {
        "perturbations": ["non-spherical", "third-body"],
        "ui": TwoBodyNonSphereThirdBody,
    },
    {
        "perturbations": ["non-spherical", "third-body", "general-relativity"],
        "ui": TwoBodyNonSphereThirdBodyGenRel,
    },
    {
        "perturbations": [
            "non-spherical",
            "third-body",
            "general-relativity",
            "atmospheric-drag",
        ],
        "ui": TwoBodyNonSphereThirdBodyGenRelAtmosDrag,
    },
    {
        "perturbations": [
            "non-spherical",
            "third-body",
            "general-relativity",
            "atmospheric-drag",
            "solar-radiation-pressure",
        ],
        "ui": TwoBodyNonSphereThirdBodyGenRelAtmosDragSrp,
    },
    {
        "perturbations": [
            "non-spherical",
            "third-body",
            "general-relativity",
            "atmospheric-drag",
            "solar-radiation-pressure",
            "earth-radiation-pressure",
        ],
        "ui": TwoBodyNonSphereThirdBodyGenRelAtmosDragSrpErp,
    },
]


def ui_seclector(perturbations: List[str]):
    """Returns the appropriate ui object depending on items in list of
    perturbations.
    """
    for p in perturbations:
        if p not in PERTURBATION_NAMES:
            raise EomUiError(f"Unexpected perturbation name found: {p}")

    for pair in ui_map:
        if set(perturbations) == set(pair["perturbations"]):
            return pair["ui"]
    raise EomUiError(
        f"No UI Object defined for perturbations: {perturbations}."
    )
