import numpy as np
from iapetus.core import probability


def test_gauss_multivar_sample():
    """Mostly verifies that implementation of random seed is consistent as
    expected, but also serves as check on consistent use of wrapper in this
    repo.
    """

    rng_1 = probability.Random(0)
    rng_2 = probability.Random(0)

    mu = [1, 2]
    cov = 0.5**2 * np.eye(2)

    X1 = rng_1.gauss_multivar(mu, cov)
    X2 = rng_2.gauss_multivar(mu, cov)

    np.testing.assert_array_equal(X1, X2)
    assert X1.ndim == 1


def test_uniform():
    """Simple consistent use of wrapper check."""

    rng = probability.Random(0)
    assert isinstance(rng.uniform(), float)
