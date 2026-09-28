"""
TODO: docstrings for this Lie test module
"""

from gl_derivations import lie, symmetric, subspaces
from sympy import Rational, Poly, ImmutableMatrix
import pytest

def test_subspace_construction_from_basis():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    d = S.basis()[3]
    V = subspaces.Subspace.from_basis(S,2*u+2*d,b,u+b+d)
    assert V.ambient == S
    assert V.dimension == 2
    assert V.basis_matrix == ImmutableMatrix([[1,0,0,1],[0,1,0,0]]).transpose()
    V = subspaces.Subspace.from_basis(S)
    assert V.ambient == S
    assert V.dimension == 0
    V = subspaces.Subspace.from_basis(S,S.zero())
    assert V.dimension == 0

def test_subspace_construction_from_equations():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    V = subspaces.Subspace.from_equations(S,ImmutableMatrix([
        [1, 0, 0, -1],
        [0, 0, 1, 0]
    ]))
    assert V.ambient == S
    assert V.dimension == 2
    V = subspaces.Subspace.from_equations(S,ImmutableMatrix([
        [7,0,0,-7],
        [0,0,1,0],
        [2,0,3,-2],
        [0,0,0,0]
    ]))
    assert V.dimension == 2

def test_both_subspace_constructions_coincide():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    V = subspaces.Subspace.from_basis(S,u-b,b-c,c-d)
    assert V.dimension == 3
    assert V.equation_matrix() == ImmutableMatrix([[1,1,1,1]])
    W = subspaces.Subspace.from_equations(S,ImmutableMatrix([[1,1,1,1]]))
    assert V == W
    Z = subspaces.Subspace.from_equations(S,ImmutableMatrix([[1,0,1,0]]))
    assert V != Z

def test_check_element_containment():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    V = subspaces.Subspace.from_basis(S,u-b,b-c,c-d)
    assert V.contains(u+b+c-3*d)
    assert all(not V.contains(x) for x in S.basis())

def test_check_cartan_subspace_factories():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    D = subspaces.Subspace.D(g)
    assert D.contains(g.E(0,0)-g.E(1,1))
    assert not D.contains(g.E(1,0))
    D0 = subspaces.Subspace.D0(g)
    assert D0.contains(g.E(0,0)-g.E(1,1))
    assert not D0.contains(g.E(0,0))
    T = subspaces.Subspace.T(g)
    assert T.contains(g.E(0,1))
    assert not T.contains(g.E(1,0))
    assert not T.contains(g.E(0,0)-g.E(1,1))
    assert T.contains(0*g.E(0,0))