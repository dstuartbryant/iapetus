# Iapetus

## Disambiguation: Objects = Targets
The term `target` is used throughout this package in place of `object` to disambiguate the notions of *object tracking* and *Python objects*.

In all software documentation terms like *multi-object filtering* will be written as *multi-target filtering*, and `object` will always refer to a Python `object`.

## Modules

Tentative package architecture.

- time
- state
- observation
- dynamics
- propagation
- scenario
- analysis
- view
- prometheus
- epimetheus


### Time

*Likely rename/rebrand this folder name to disambiguate*.

Some scenarios will be relatively simple and won't require special considerations for `time`. Envisoned, however, is the future development of advanced models and scenarios for which sophisticated `time` tools will be required. Those tools will live in this folder.


### State

All conditions of a target's existence with respect to scenario circumstances are encapsulated in a objects contained in the `state` module.

### Observation

Sensor and sensor platform modeling tools are contained in the `observation` module.


### Dynamics
The `dynamics` module is for modeling and any fiducial data needed to define state change.

### Propagation
The `propagation` module contains tools for enacting state tranformation.


### Scenario
The `scenario` module provides tools for defining and interfacing with objects that establish conditions for filter test cases.


### Analysis
Tools for filter evaluation, assessment, analysis, etc. live in the `analysis` module. E.g., OSPA tools live here.


### View
The `view` module contains tools for plotting.





# For contributing
## Bib citation styles
See https://github.com/citation-style-language/styles