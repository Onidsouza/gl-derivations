"""
TODO: docstrings for this Invariants test module
"""


from gl_derivations import lie, subspaces, symmetric


def test_invariant_subalgebras():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,1)
    u = S.basis()[0]
    b = S.basis()[1]
    c = S.basis()[2]
    d = S.basis()[3]
    assert subspaces.Subspace.whole(S).invariants(subspaces.Subspace.trivial(g)) == subspaces.Subspace.whole(S)
    assert subspaces.Subspace.whole(S).invariants(subspaces.Subspace.from_basis(g,g.E(0,0)+g.E(1,1))) == subspaces.Subspace.whole(S)
    assert subspaces.Subspace.whole(S).invariants(subspaces.Subspace.whole(g)) == subspaces.Subspace.from_basis(S,u+d)
    assert subspaces.Subspace.whole(S).invariants(subspaces.Subspace.T(g)) == subspaces.Subspace.from_basis(S,u+d,b)

def test_quadractic_solution_gl2():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,2)
    C1 = S.monomial(2,0,0,0) + 2* S.monomial(0,1,1,0) + S.monomial(0,0,0,2)
    C2 = S.monomial(1,0,0,1) - S.monomial(0,1,1,0)
    assert subspaces.Subspace.whole(S).invariants(subspaces.Subspace.whole(g)) == subspaces.Subspace.from_basis(S,C1,C2)