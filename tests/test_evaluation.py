from sympy import Rational

from glder import *


def test_evaluation():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    v = 2*S.monomial(2,0,0,0) + Rational(-1,3) * S.monomial(0,1,1,0)
    x = g.E(0,0) + 2*g.E(0,1) + 3*g.E(1,0) + 4*g.E(1,1)
    assert evaluate(v,x) == 0
    v = S.monomial(2,0,0,0) + 2*S.monomial(1,0,0,1) + S.monomial(0,0,0,2)
    assert evaluate(v,x) == 25
    v = S.monomial(1,0,1,0)
    assert evaluate(v,x) == 2

def test_evaluation_kernel():
    g = GeneralLinear(2)
    S = SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    x = g.from_coordinates(1,2,3,4)
    assert Subspace.whole(S).evaluation_kernel(g.zero()) == Subspace.whole(S)
    assert Subspace.from_basis(S,u+d,b).evaluation_kernel(x) == Subspace.from_basis(S,3*u-5*b+3*d)

def test_complete_quadractic_calculation():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    X = Subspace.from_basis(g,*(g.basis()))
    Y = Subspace.whole(S).invariants(X) # central quadractic elements.
    q1 = S.monomial(2,0,0,0) + 2*S.monomial(1,0,0,1) + S.monomial(0,0,0,2)
    q2 = S.monomial(1,0,0,1) - S.monomial(0,1,1,0)
    assert Y == Subspace.from_basis(S,q1,q2)
    assert Y == Subspace.from_equations(S,Y.equation_matrix())
    assert Y.is_diagonal_stable()
    assert Y == Y.weight_decomposition()[(0,0)]
    x = g.from_coordinates(1,2,3,4)
    assert Y.evaluation_kernel(x) == Subspace.from_basis(S,2*q1 + 25*q2)