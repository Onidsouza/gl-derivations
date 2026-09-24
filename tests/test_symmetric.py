"""
TODO: docstrings for this Lie test module
"""

from gl_derivations import lie, symmetric
import pytest

@pytest.mark.parametrize("value_n, value_k, expected_dimension", [
    (1,4,1),
    (2,3,20),
    (3,2,45),
    (5,6,593775)
])

def test_construction_symmetric_power(value_n,value_k,expected_dimension):
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(value_n)
    S = symmetric.SymmetricPower(g,value_k)
    assert S.algebra == g
    assert S.degree == value_k 
    assert S.dimension == expected_dimension
    assert len(S.generators) == S.algebra.dimension

@pytest.mark.parametrize("value_n,value_k, expected_basis", [
    (1,4,((4,),)),
    (2,1,((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1))),
    (2,2,((2,0,0,0),(1,1,0,0),(1,0,1,0),(1,0,0,1),(0,2,0,0),(0,1,1,0),(0,1,0,1),(0,0,2,0),(0,0,1,1),(0,0,0,2))),
    (4,0,((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),))
])

def test_basis_monomial_construction_symmetric_power(value_n,value_k,expected_basis):
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(value_n)
    S = symmetric.SymmetricPower(g,value_k)
    assert S.basis_labels() == expected_basis
    assert len(S.basis_labels()) == S.dimension

def test_monomial_construction_from_tuple():
    """
    TODO: docstring for this test function.
    """
    g = lie.GeneralLinear(2)
    S = symmetric.SymmetricPower(g,2)
    assert S.monomial(2,0,0,0) == S.basis()[0]
    with pytest.raises(TypeError):
        S.monomial(1,4,7,0,3)
    with pytest.raises(TypeError):
        S.monomial(0,0.5,0,0)
    with pytest.raises(TypeError):
        S.monomial(True,True,0,0)
    with pytest.raises(ValueError):
        S.monomial(3,0,0,0)