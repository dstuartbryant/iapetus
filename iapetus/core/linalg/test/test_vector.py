import numpy as np
import pytest
from iapetus.core.linalg import Vector, VectorException

A = [1, 2, 3]
B = np.array([4, 5, 6])


def test_list_init():
    x = Vector(A)
    assert str(x) == "[1 2 3]"
    assert repr(x) == "array([1, 2, 3])"
    assert x.dim == 3
    assert x.mag == np.linalg.norm(A)
    np.testing.assert_array_equal(x.unit, np.array(A) / x.mag)


def test_array_init():
    x = Vector(B)
    assert str(x) == "[4 5 6]"
    assert repr(x) == "array([4, 5, 6])"
    assert x.dim == 3
    assert x.mag == np.linalg.norm(B)
    np.testing.assert_array_equal(x.unit, np.array(B) / x.mag)


def test_dimension_exception():
    m = np.eye(2)
    with pytest.raises(VectorException):
        Vector(m)


def test_nested_lists_exception():
    m = [[1, 2, 3], [4, 5, 6]]
    with pytest.raises(VectorException):
        Vector(m)
