"""EOM UI test module."""

from iapetus.propagation.dynamics.nonlinear.astro.eom.payloads.ui import (
    PERTURBATION_NAMES,
    individual_eom_map,
    individual_init_config_map,
    ui_map,
)


def test_ui_map_pert_names():
    for item in ui_map:
        for p in item["perturbations"]:
            assert p in PERTURBATION_NAMES


def test_init_config_map_pert_names():
    for k in individual_init_config_map.keys():
        assert k in PERTURBATION_NAMES


def test_eom_map_pert_names():
    for k in individual_eom_map.keys():
        assert k in PERTURBATION_NAMES
