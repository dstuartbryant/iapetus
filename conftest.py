"""Iapetus conftest.py"""

import sys
from os import path

CURR_DIR = path.dirname(path.abspath(__file__))
TEST_DIR = path.join(CURR_DIR, "test")
sys.path.append(TEST_DIR)

from fixtures import *
