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
from dataclasses import fields
from typing import List

from pydantic import BaseModel

from ..perturbations import AtmosphericPerturbation
from ..two_body import TwoBody as TwoBodyEom
from .configs import (
    AtmosphericDragInitConfig,
    ErpInitConfig,
    GeneralRelativityInitConfig,
    NonSphericalInitConfig,
    SrpInitConfig,
    ThirdBodyInitConfig,
    TwoBodyInitConfig,
)

PERTURBATION_NAMES = [
    "non-spherical",
    "third-body",
    "general-relativity",
    "atmospheric-drag",
    "solar-radiation-pressure",
    "earth-radiation-pressure",
]


def classFromArgs(className, argDict):
    fieldSet = {f.name for f in fields(className) if f.init}
    filteredArgDict = {k: v for k, v in argDict.items() if k in fieldSet}
    return className(**filteredArgDict)


class EomUiError(Exception):
    pass


class TwoBody(BaseModel):
    mu: float
    partials_flag: bool


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

individual_init_config_map = {
    "non-spherical": NonSphericalInitConfig,
    "third-body": ThirdBodyInitConfig,
    "general-relativity": GeneralRelativityInitConfig,
    "atmospheric-drag": AtmosphericDragInitConfig,
    "solar-radiation-pressure": SrpInitConfig,
    "earth-radiation-pressure": ErpInitConfig,
}

individual_eom_map = {
    "non-spherical": NotImplemented,
    "third-body": NotImplemented,
    "general-relativity": NotImplemented,
    "atmospheric-drag": AtmosphericPerturbation,
    "solar-radiation-pressure": NotImplemented,
    "earth-radiation-pressure": NotImplemented,
}


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


def match_ui_config_to_individual_init_configs(ui_config, individual_config):
    required_fields = list(individual_config.__dataclass_fields__.keys())
    arg_dict = {}
    for k in required_fields:
        arg_dict[k] = getattr(ui_config, k)
    return classFromArgs(individual_config, arg_dict)


def configure_eom_from_user_config(
    user_config: dict, perturbations: List[str]
):
    ui_model = ui_seclector(perturbations)
    ui_config = ui_model(**user_config)
    two_body_ic = match_ui_config_to_individual_init_configs(
        ui_config, TwoBodyInitConfig
    )
    two_body_eom = TwoBodyEom(two_body_ic)
    pert_eom = []
    for p in perturbations:
        eom_obj = individual_eom_map[p]
        pert_eom.append(
            eom_obj(
                match_ui_config_to_individual_init_configs(
                    ui_config, individual_init_config_map[p]
                )
            )
        )
    return two_body_eom, pert_eom
