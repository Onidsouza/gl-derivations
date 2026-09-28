"""
TODO: docstrings for this Weights test module
"""

from gl_derivations import lie, symmetric, subspaces, weights, actions
from sympy import Rational, Poly, ImmutableMatrix
import pytest

@pytest.mark.parametrize("value_n, value_k, monomial, weight", [
    (2,1,(1,0,0,0),(0,0)),
    (2,1,(0,0,0,1),(0,0)),
    (2,1,(0,1,0,0),(1,-1)),
    (2,1,(0,0,1,0),(-1,1)),
    (2,2,(0,1,1,0),(0,0)),
    (2,2,(0,2,0,0),(2,-2)),
    (2,2,(1,1,0,0),(1,-1)),
    (2,2,(0,1,0,1),(1,-1)),
    (3,3,(0,2,0,0,0,1,0,0,0),(2,-1,-1))
])

def test_weight_from_monomial_computation(value_n,value_k,monomial,weight):
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(value_n)
    S = symmetric.SymmetricPower(g,value_k)
    assert S.weight_from_exponent_tuple(*monomial) == weight

def test_constructor_of_weight_elements():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,2)
    assert len(S.weight_elements(0,0)) == 4
    assert S.basis()[0] in S.weight_elements(0,0)
    assert S.basis()[3] in S.weight_elements(0,0)
    assert S.basis()[1] in S.weight_elements(1,-1)
    assert S.basis()[2] in S.weight_elements(-1,1)
    assert len(S.weight_elements(-2,2)) == 1
    assert len(S.weight_elements(3,0)) == 0

def test_stability_under_diagonals():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    assert subspaces.Subspace.from_basis(S,u+d).is_diagonal_stable() == True
    assert subspaces.Subspace.from_basis(S,u+b).is_diagonal_stable() == False
    assert subspaces.Subspace.from_basis(S,u+d,b).is_diagonal_stable() == True
    assert subspaces.Subspace.from_basis(S).is_diagonal_stable() ==  True # the trivial subspace
    assert subspaces.Subspace.from_equations(S,ImmutableMatrix(0,S.dimension,lambda i,j: Rational(0))).is_diagonal_stable() == True # the whole space

def test_weight_space_decomposition():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    V = subspaces.Subspace.from_basis(S,u+d,b)
    decomposition = V.weight_decomposition()
    assert len(decomposition.keys()) == 2
    assert (0,0) in decomposition.keys()
    assert (1,-1) in decomposition.keys()
    assert (-1,1) not in decomposition.keys()
    tally = 0
    for key in decomposition.keys():
        tally += decomposition[key].dimension
    assert tally == V.dimension
    assert decomposition[(0,0)].contains(u+d)
    assert decomposition[(1,-1)].contains(b)
    assert actions.action_matrix(g.E(0,0),S)*decomposition[(1,-1)].basis_matrix == decomposition[(1,-1)].basis_matrix
    assert actions.action_matrix(g.E(1,1),S)*decomposition[(1,-1)].basis_matrix == -decomposition[(1,-1)].basis_matrix