"""Astrodynamics equations of motion (EOM) module."""

from typing import List

from .payloads.ui import configure_eom_from_user_config
from .perturbations.models import Perturbation, PerturbedOutput
from .state import STATE_TYPES
from .two_body import TwoBody


class EomError(Exception):
    pass


def partials_wrapper(name: str, pert_output: PerturbedOutput):
    """Attempts to get a named partial derivative value, returns zero if
    perturbation object does not have an attribute of that name.

    Not all perturbations compute the same partial derivatives, so, by
    convention, it is assumed that if a given partial derivative does not exist
    for a given pertubation model, then the value that is returned is zero,
    which will not effect the summation of partial derivative components across
    all perturbation models.

    """
    try:
        val = getattr(pert_output, name)
    except AttributeError:
        val = 0
    return val


class Eom:
    """Equations of motion (EOM) class."""

    def __init__(
        self, eom_two_body: TwoBody, eom_perturbations: List[Perturbation]
    ):
        self.eom_two_body = eom_two_body
        self.eom_perturbations = eom_perturbations

    def __call__(self, s: STATE_TYPES):
        two_body_output = self.eom_two_body(s)
        perturbations_output = [x(s) for x in self.eom_perturbations]

        def accelerations_getter(component: str):
            val = getattr(two_body_output.accelerations, component)
            for p in perturbations_output:
                val += getattr(p.accelerations, component)
            return val

        def partials_getter(name: str):
            try:
                val = partials_wrapper(name, two_body_output.partials)
                for p in perturbations_output:
                    val += partials_wrapper(name, p.partials)
                return val
            except TypeError as e:
                if (
                    str(e.args[0])
                    == "unsupported operand type(s) for +=: 'int' and "
                    "'NoneType'"
                ):
                    raise EomError(
                        f"Partials component {name}: Attempting to add to "
                        "NoneType"
                    )

        return accelerations_getter, partials_getter


def configure(user_config: dict, perturbations: List[str]) -> Eom:
    """Configures and returns a callable Eom object based on user
    configurations and a list of perturbations.

    Args:
        user_config (dict): Dictionary of initialization paramers defined by
            user input
        perturbations (List[str]): List of perturbations to use, if any. Can be
            an empty list.
    """
    two_body_eom, pert_eom = configure_eom_from_user_config(
        user_config, perturbations
    )
    return Eom(eom_two_body=two_body_eom, eom_perturbations=pert_eom)
