# Iapetus

## Disambiguation: Objects = Targets
The term `target` is used throughout this package in place of `object` to disambiguate the notions of *object tracking* and *Python objects*.

In all software documentation terms like *multi-object filtering* will be written as *multi-target filtering*, and `object` will always refer to a Python `object`.

## Modules

Tentative package architecture.

- dynamics
- propagation
- scenario
- state
- time
- measurement
- view
- evaluation
- prometheus
- epimetheus



### Dynamics
The `dynamics` folder is for modeling and any fiducial data needed to define state change.

### Propagation
The `propagation` folder contains tools for enacting state tranformation.

### Time

*Likely rename/rebrand this folder name to disambiguate*.

Some scenarios will be relatively simple and won't require special considerations for `time`. Envisoned, however, is the future development of advanced models and scenarios for which sophisticated `time` tools will be required. Those tools will live in this folder.

### State



### Scenario

