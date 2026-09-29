
from glder import *


def test_invariant_subalgebras():
    g = GeneralLinear(2)
    S = SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    assert Subspace.whole(S).invariants(Subspace.trivial(g)) == Subspace.whole(S)
    assert Subspace.whole(S).invariants(Subspace.from_basis(g,g.E(0,0)+g.E(1,1))) == Subspace.whole(S)
    assert Subspace.whole(S).invariants(Subspace.whole(g)) == Subspace.from_basis(S,u+d)
    assert Subspace.whole(S).invariants(Subspace.T(g)) == Subspace.from_basis(S,u+d,b)

def test_quadractic_solution_gl2():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    C1 = S.monomial(2,0,0,0) + 2* S.monomial(0,1,1,0) + S.monomial(0,0,0,2)
    C2 = S.monomial(1,0,0,1) - S.monomial(0,1,1,0)
    assert Subspace.whole(S).invariants(Subspace.whole(g)) == Subspace.from_basis(S,C1,C2)