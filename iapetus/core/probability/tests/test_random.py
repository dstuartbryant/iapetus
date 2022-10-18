import numpy as np
from iapetus.core import probability


def test_gauss_multivar_sample():
    rand_1 = probability.Random(0)
    rand_2 = probability.Random(0)

    mu = [1, 2]
    cov = 0.5**2 * np.eye(2)

    X1 = rand_1.gauss_multivar_sample(mu, cov)
    X2 = rand_2.gauss_multivar_sample(mu, cov)

    np.testing.assert_array_equal(X1, X2)
    assert X1.ndim == 1
