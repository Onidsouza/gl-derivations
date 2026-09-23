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

def test_lie_element_arithmetic():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    expected = ImmutableMatrix([
        [Rational(1), Rational(1)],
        [Rational(0), Rational(0)]
    ])
    z = g.E(0,0) + g.E(0,1)
    assert z.parent == g
    assert z.matrix == expected
    expected = ImmutableMatrix([
        [Rational(-1), Rational(-1)],
        [Rational(0), Rational(0)]
    ])
    z = -z
    assert z.parent == g
    assert z.matrix == expected
    z = z-z 
    assert z == g.zero()
    z = g.E(0,1)*2
    expected = ImmutableMatrix([
        [Rational(0), Rational(2)],
        [Rational(0), Rational(0)]
    ])
    assert z.parent == g
    assert z.matrix == expected
    z = g.E(1,0)*z
    expected = ImmutableMatrix([
        [Rational(0), Rational(0)],
        [Rational(0), Rational(2)]
    ])
    assert z.parent == g
    assert z.matrix == expected
    z = Rational(-3)*z
    expected = ImmutableMatrix([
        [Rational(0), Rational(0)],
        [Rational(0), Rational(-6)]
    ])
    assert z.parent == g
    assert z.matrix == expected

def test_gl_coordinate_conversion():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    z = g.from_coordinates(ImmutableMatrix([0,1,0,0]))
    assert z == g.E(0,1)
    z = g.E(0,0) - g.E(1,0)
    assert g.coordinates(z) == ImmutableMatrix([1,0,-1,0])
    z = g.E(1,1)
    coord = ImmutableMatrix([1,2,Rational(2,3),4])
    assert z == g.from_coordinates(g.coordinates(z))
    assert coord == g.coordinates(g.from_coordinates(coord))

def test_lie_bracket():
    """
    TODO: docstring for this test function.
    """

    g = lie.GeneralLinear(2)
    assert g.E(0,1).bracket(g.E(1,0)) == g.E(0,0) - g.E(1,1)
    assert g.E(0,1).bracket(g.E(0,0)) == -g.E(0,1)
    assert g.E(0,1).bracket(g.E(1,1)) == g.E(0,1)
    z = g.E(0,0) + g.E(1,1)
    for x in g.basis():
        assert z.bracket(x) == g.zero()
    for x in g.basis():
        for y in g.basis():
            assert x.bracket(y) == -y.bracket(x)
            for z in g.basis():
                assert x.bracket(y.bracket(z)) + z.bracket(x.bracket(y)) + y.bracket(z.bracket(x)) == g.zero()

@pytest.mark.parametrize("value_n", [
    (1),
    (2),
    (3),
    (6)
])

def test_lie_shape_factories(value_n):
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(value_n)
    diags = g.diagonals()
    tr_diags = g.traceless_diagonals()
    triang = g.strictly_upper_triangulars()
    assert len(diags) == g.n
    assert len(tr_diags) == g.n - 1
    assert len(triang) == (g.dimension - g.n) // 2
    for x in diags:
        assert x.matrix.trace() == 1 and x.matrix.is_diagonal()
    for x in tr_diags:
        assert x.matrix.trace() == 0 and x.matrix.is_diagonal()