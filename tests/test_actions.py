import pytest
from sympy import ImmutableMatrix, zeros

from glder import *


def test_adjoint_action():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    x = g.E(0,1)
    v = S.from_terms({(1,0,1,0): 1})
    assert act(x,v) == S.from_terms({(2,0,0,0) : 1, (1,0,0,1) : -1, (0,1,1,0) : -1})
    v = S.from_terms({(2,0,0,0): 1})
    assert act(x,v) == S.from_terms({(1,1,0,0): -2})

def test_adjoint_action_matrix():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    x = g.E(0,0) + g.E(1,1)
    assert action_matrix(x,S) == zeros(10,10)
    for x in g.basis():
        A = action_matrix(x,S)
        for v in S.basis():
            assert act(x,v) == S.from_coordinates(A*S.coordinates(v))
    x = g.E(0,1)
    A = ImmutableMatrix([
        [0, 0, 1, 0],
        [-1, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, -1, 0]
    ])
    S = SymmetricPower(g,1)
    assert action_matrix(x,S) == A

@pytest.mark.parametrize('value_n', [
    (2,),
    (3,)
])

def test_degree_one_agrees_with_bracket(value_n):
    g = GeneralLinear(*value_n)
    S = SymmetricPower(g,1)
    for x in g.basis():
        for y in g.basis():
            assert x.bracket(y) == g.from_coordinates(S.coordinates(act(x,S.from_coordinates(g.coordinates(y)))))

@pytest.mark.parametrize('value_n,value_k', [
    (2,1),
    (2,2),
    (3,1),
])

def test_representation_preserves_bracket(value_n,value_k):
    g = GeneralLinear(value_n)
    S = SymmetricPower(g,value_k)
    for x in g.basis():
        for y in g.basis():
            for v in S.basis():
                assert act(x,act(y,v)) - act(y,act(x,v)) == act(x.bracket(y),v)