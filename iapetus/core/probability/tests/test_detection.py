import numpy as np
from iapetus.core.probability import detection


def test_uniform_percentage():
    heur = 0.8
    false = heur + 1e-3
    true = heur - 1e-3

    assert detection.uniform_percentage(heur, false) == False
    assert detection.uniform_percentage(heur, true) == True
