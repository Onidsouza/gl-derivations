"""
TODO: docstrings for this Evaluation test module
"""

from gl_derivations import lie, symmetric, subspaces, actions
from sympy import Rational, Poly, ImmutableMatrix
import pytest

def test_evaluation():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,2)
    v = 2*S.monomial(2,0,0,0) + Rational(-1,3) * S.monomial(0,1,1,0)
    x = g.E(0,0) + 2*g.E(0,1) + 3*g.E(1,0) + 4*g.E(1,1)
    assert actions.evaluate(v,x) == 0
    v = S.monomial(2,0,0,0) + 2*S.monomial(1,0,0,1) + S.monomial(0,0,0,2)
    assert actions.evaluate(v,x) == 25
    v = S.monomial(1,0,1,0)
    assert actions.evaluate(v,x) == 2

def test_evaluation_kernel():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    x = g.from_coordinates(1,2,3,4)
    assert subspaces.Subspace.whole(S).evaluation_kernel(g.zero()) == subspaces.Subspace.whole(S)
    assert subspaces.Subspace.from_basis(S,u+d,b).evaluation_kernel(x) == subspaces.Subspace.from_basis(S,3*u-5*b+3*d)

def test_complete_quadractic_calculation():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,2)
    X = subspaces.Subspace.from_basis(g,*(g.basis()))
    Y = subspaces.Subspace.whole(S).invariants(X) # central quadractic elements.
    q1 = S.monomial(2,0,0,0) + 2*S.monomial(1,0,0,1) + S.monomial(0,0,0,2)
    q2 = S.monomial(1,0,0,1) - S.monomial(0,1,1,0)
    assert Y == subspaces.Subspace.from_basis(S,q1,q2)
    assert Y == subspaces.Subspace.from_equations(S,Y.equation_matrix())
    assert Y.is_diagonal_stable()
    assert Y == Y.weight_decomposition()[(0,0)]
    x = g.from_coordinates(1,2,3,4)
    assert Y.evaluation_kernel(x) == subspaces.Subspace.from_basis(S,2*q1 + 25*q2)