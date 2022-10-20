"""Trajectory module.

A trajectory is a sequence of states, indexed in time, that represents an
object's open-ended path that does not exactly - or necessarily - repeat.

"""

from dataclasses import dataclass
from typing import List

from .state import State

# @dataclass
# class Trajectory:
#     """Trajectory base class.

#     How to initialize?

#     * If assume is truth trajectory, then it can be completely populated at
#       once.

#     * If assume a track that is being built by a filter, then it is populated
#       one state at a time.


#     What does/can it do?

#     * Essentially acts as a container to hold a sequence of states for a given object for later use

#     * Should retain some unique marker/tag to identify the object it represents


#     So, really, what would one want to do with a complete trajectory once they have it?

#     * Plot
#         - Vector component vs time: So you'd want to strip out a single component and it's associated time

#         - Vector component vs Vector component: 2D position plot really. Need
#           to strip out the x,y components in order

#     * Compute differences

#     """

#     label: str  # Uniquely identifies the trajectory
#     states: List[State]


class Trajectory:
    def __init__(self, label: str, states: List[State]):
        self.label = label
        states.sort(key=lambda x: x.timestamp)
        self.states = states

    def _state_vectors(self):
        return [x.vector for x in self.states]

    @property
    def state_vectors(self):
        return self._state_vectors()
