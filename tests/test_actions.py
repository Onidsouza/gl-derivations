"""
TODO: docstrings for this Lie test module
"""

from gl_derivations import lie,symmetric,actions
from sympy import Rational, Poly, ImmutableMatrix, zeros
import pytest

def test_adjoint_action():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,2)
    x = g.E(0,1)
    v = S.from_terms({(1,0,1,0): 1})
    assert actions.act(x,v) == S.from_terms({(2,0,0,0) : 1, (1,0,0,1) : -1, (0,1,1,0) : -1})
    v = S.from_terms({(2,0,0,0): 1})
    assert actions.act(x,v) == S.from_terms({(1,1,0,0): -2})

def test_adjoint_action_matrix():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,2)
    x = g.E(0,0) + g.E(1,1)
    assert actions.action_matrix(x,S) == zeros(10,10)
    for x in g.basis():
        A = actions.action_matrix(x,S)
        for v in S.basis():
            assert actions.act(x,v) == S.from_coordinates(A*S.coordinates(v))
    x = g.E(0,1)
    A = ImmutableMatrix([
        [0, 0, 1, 0],
        [-1, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, -1, 0]
    ])
    S = symmetric.SymmetricPower(g,1)
    assert actions.action_matrix(x,S) == A

@pytest.mark.parametrize('value_n', [
    (2,),
    (3,)
])

def test_degree_one_agrees_with_bracket(value_n):
    g = lie.GeneralLinear(*value_n)
    S = symmetric.SymmetricPower(g,1)
    for x in g.basis():
        for y in g.basis():
            assert x.bracket(y) == g.from_coordinates(S.coordinates(actions.act(x,S.from_coordinates(g.coordinates(y)))))

@pytest.mark.parametrize('value_n,value_k', [
    (2,1),
    (2,2),
    (3,1),
])

def test_representation_preserves_bracket(value_n,value_k):
    g = lie.GeneralLinear(value_n)
    S = symmetric.SymmetricPower(g,value_k)
    for x in g.basis():
        for y in g.basis():
            for v in S.basis():
                assert actions.act(x,actions.act(y,v)) - actions.act(y,actions.act(x,v)) == actions.act(x.bracket(y),v)