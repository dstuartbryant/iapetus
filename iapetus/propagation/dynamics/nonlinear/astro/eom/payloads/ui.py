"""User interface module for equations of motion (EOM) configuration.

Perturbations Options:
* Non-Spherical (Earth)
* ThirdBody
* GeneralRelativity
* AtmosphericDrag
* Srp
* Erp (Earth Radiation Pressure)
* 

"""

from pydantic import BaseModel


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
