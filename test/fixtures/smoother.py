"""Smoother test fixture(s)."""

import json
from os import path
from typing import List

import numpy as np
import pytest

CURR_DIR = path.dirname(path.abspath(__file__))
DATA_DIR = path.join(CURR_DIR, "..", "data")


def reshape_array(A: List[float], n: int) -> np.ndarray:
    B = np.array(A)
    return B.reshape((n, n))


@pytest.fixture
def smooth_input():
    fpath = path.join(DATA_DIR, "smoother_input_data.json")
    loaded_data = json.load(open(fpath, "r"))
    data = {"X": [], "P": [], "PHI": [], "Q": []}
    for idx in range(len(loaded_data["X_stash"])):
        n = len(loaded_data["X_stash"][idx])
        data["X"].append(np.array(loaded_data["X_stash"][idx]))
        data["P"].append(reshape_array(loaded_data["P_stash"][idx], n))
        data["PHI"].append(reshape_array(loaded_data["Phi_stash"][idx], n))
        data["Q"].append(reshape_array(loaded_data["Q_stash"][idx], n))

    yield data
