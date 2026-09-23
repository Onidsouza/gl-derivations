"""
TODO: docstrings for this Lie test module
"""

from gl_derivations import lie
from sympy import Rational, ImmutableMatrix
import pytest

@pytest.mark.parametrize("value_n, dim", [
    (1,1),
    (2,4),
    (3,9),
    (4,16),
    (5,25),
])

def test_construction_gl(value_n,dim):
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(value_n)
    assert g.n == value_n
    assert g.dimension == dim

@pytest.mark.parametrize("value_n", [
    (True),
    (1.0),
    (Rational(1)),
    (8/2.0)
])

def test_type_validation_dimension_gl(value_n):
    """
    TODO: docstring for this test function.
    """
    with pytest.raises(TypeError):
        lie.GeneralLinear(value_n)

@pytest.mark.parametrize("value_n", [
    (-1),
    (0)
])

def test_range_validation_dimension_gl(value_n):
    """
    TODO: docstring for this test function.
    """
    with pytest.raises(ValueError):
        lie.GeneralLinear(value_n)

def test_equality_gl():
    """
    TODO: docstring for this test function.
    """
    assert lie.GeneralLinear(2) == lie.GeneralLinear(1+1)
    assert lie.GeneralLinear(3) != lie.GeneralLinear(2)

@pytest.mark.parametrize("value_n", [
    (1),
    (2),
    (5),
    (7)
])

def test_zero_matrix_gl(value_n):
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(value_n)
    z = g.zero()
    assert z.parent == g
    assert z.matrix.shape == (g.n, g.n)
    assert z.matrix == ImmutableMatrix(g.n,g.n, lambda i,j: Rational(0))

@pytest.mark.parametrize("value_n,i,j,matrix", [
    (1,0,0,ImmutableMatrix([[Rational(1)]])),
    (2,0,0,ImmutableMatrix([[Rational(1), Rational(0)], [Rational(0),Rational(0)]])),
    (3,1,2,ImmutableMatrix([
        [Rational(0), Rational(0), Rational(0)],
        [Rational(0), Rational(0), Rational(1)],
        [Rational(0), Rational(0), Rational(0)]
    ]))
])

def test_elementary_matrix_gl(value_n,i,j,matrix):
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(value_n)
    z = g.E(i,j)
    assert z.parent == g
    assert z.matrix.shape == (g.n, g.n)
    assert z.matrix == matrix

def test_gl_basis():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    basis = g.basis()
    expected = (g.E(0,0), g.E(0,1), g.E(1,0), g.E(1,1))
    assert len(basis) == g.dimension
    assert basis == expected