"""Scenario management module.


Developing a single-target tracking scenario is fairly straightforward given
that only one trajectory is considered. At each timestep, this one target
can be born (appear in the surveillance region), survive (propagated from 
previous time step), or die (leave the surveillance region). So that's 1 target
and 3 behaviors in need of consideration.

Developing multi-target tracking scenarios can be a far more tedious task. Each
target can be born, die, survive and spawn new targets. And so, that's N
targets and 4 behaviors in need of consideration, where N is the total number
of targets considered in a given scenario.

At times we'll want to analyse a tracking filter's performance under specific
circumstances, e.g., when two targets physically move past each other, or even
move along side each other, at close range. Tuning scenario parameters to
simulate such circumstances adds further complexity. 


Scenarios Have:
* trajectory(ies)
* sensor(s)
* surveillance region(s)
* illustrations




Scenario Creation ConOps
------------------------

Assumptions

    A1: Trajectories can be designed with some "play" at their entry and exit
        points (of a survellience region). I.e., if we need to lop off states
        leading into the surveillance region box b/c we need to tighten in the
        region for some reasonable constraints (even aesthetic), then we need
        to be able to flag states such that the target was not yet "born".
        In the case of "death", we flag states near the end of the trajectory
        as appropriate.
            - Still thinking on this. In simple 2-D cases, we may also be able
              to play with timestamps, but that seems risky and would NOT
              translate well to astrodynamic applications.

"""
