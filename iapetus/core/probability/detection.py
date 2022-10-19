"""Probability of detection module.

For some time after inception, this module may only contain a method commonly
used to simulate sensor detection behavior. 

But it may stand to reason that more complex detection probability models are
developed here.
"""


def uniform_percentage(heur, unismpl):
    """Returns True (detected) or False (missed).

    Detection based on comparison between the input heuristic probability of
    detection and a random sample from the standard uniform distribution.

    Args:
        heur (float): heuristic probability of detection, such that
            0 <= heur <= 1.

        unismpl (float): random sample from standard uniform distribution

    Returns:
        (bool): Detection signal True or False.
    """

    return True if unismpl <= heur else False
