"""Probability models and tools."""

import numpy as np


class Probability(Exception):
    pass


class Random:
    """Random number generating class.

    Incorporates new best practice for seeding Numpy random number generators.
    See:
        - https://builtin.com/data-science/numpy-random-seed
            [Accessed 17-Oct-2022]
        - https://numpy.org/doc/stable/reference/random/generated/numpy.random.seed.html
            [Accessed 17-Oct-2022]
    """

    def __init__(self, seed=0):
        """
        Args:
            seed (int): Any non-negative integer.
        """
        self.seed = seed
        self.r = np.random.default_rng(seed)

    def gauss_multivar(self, mu, cov):
        """Returns sample drawn from a multivariate normal distribution.

        Args:
            mu (array_like): 1-D, of length N
            cov (array_like): 2-D, of shape NxN

        Returns 1-D numpy.ndarray.
        """
        return self.r.multivariate_normal(mu, cov)

    def uniform(self):
        """Returns sample drawn from a standard uniform distribution.

        **Standard** implies the distribution is bounded by 0 and 1.

        Returns:
            (float): random number x, such that 0 <= x <= 1
        """

        return self.r.uniform(low=0, high=1)
