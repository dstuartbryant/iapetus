"""Sensors module.


SEE Obsidian note: Elements of a Generic Sensor Tracking System
- continue developing notes out there, then leverge them to support
  organization here.





Sensors Can Be
--------------

* abstract - e.g., for academic, demonstrative purposes
    * e.g., an omniscient sensor that can see a 2x2 grid at all times
* concrete - e.g., pertaining to or concerned with actual sensors
    - though, a modeled sensor here would, most likely, in fact not be
      "concrete" per se, but would be more detailed than an abstract sensor
    - PERHAPS a better term would be PARTICULARIZED



Abstract Sensors Can Be
-----------------------
* turned on/off
* omniscient; all-seeing, all-knowing
* variably to perfectly (or nearly perfectly) precise and accurate
    - nearly perfectly most likely; to avoid numerical issues
* fixed at arbitrary locations
    * i.e., defining it's location as "the/a origin" thus establishing a frame
      of reference for coordinates use
    * e.g., origin (or translated from) of a Cartesian coordinate system
* in motion
    - translational
    - rotational
    - both







Bases for Concrete (Particularized) Sensors Include
---------------------------------------------------

* sonar
* radar
* laser
* lidar
* optical
* lidar
* audio


Concrete (Particularized) Sensors Can Be
----------------------------------------
* turned on/off
* variably to nearly-perfectly precise and accurate
    - likely, nearly-perfectly (to perfectly) would be used for sanity checking
      and debugging purposes.

* fixed at arbitrary locations
* in motion; translational, rotational


Sensors Have
------------
* on/off switch
* state (static or in motion)
* observation models
* aperture
* SNR
* Figure(s) of merit
* fields of view
* fields of regard




"""

from pydantic import BaseModel, Field

from ..core.state.state import State


class Sensor(BaseModel):
    """Sensor model base class."""

    name: str
    on_switch: bool = True
    state: State
    static: bool = True
    obs_model: object  # FOR NOW FIGURE THIS OUT FOR SIMPLE DEMO CASES, I.E. as simple as H = [I_2, 0_2]
    aperture: object
    snr: float
    figures_of_merit: object

    def field_of_regard(self):
        raise NotImplementedError

    def field_of_view(self):
        raise NotImplementedError


class OmniSensor(Sensor):
    pass
