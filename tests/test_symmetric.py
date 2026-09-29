import pytest
from sympy import ImmutableMatrix, Poly, Rational

from glder import *


@pytest.mark.parametrize("value_n, value_k, expected_dimension", [
    (1,4,1),
    (2,3,20),
    (3,2,45),
    (5,6,593775)
])

def test_construction_symmetric_power(value_n,value_k,expected_dimension):
    g = GeneralLinear(value_n)
    S = SymmetricPower(g,value_k)
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
    g = GeneralLinear(value_n)
    S = SymmetricPower(g,value_k)
    assert S.basis_labels() == expected_basis
    assert len(S.basis_labels()) == S.dimension

def test_monomial_construction_from_tuple():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    assert S.monomial(2,0,0,0) == S.basis()[0]
    with pytest.raises(TypeError):
        S.monomial(1,4,7,0,3)
    with pytest.raises(TypeError):
        S.monomial(0,0.5,0,0)
    with pytest.raises(TypeError):
        S.monomial(True,True,0,0)
    with pytest.raises(ValueError):
        S.monomial(3,0,0,0)

def test_symmetric_element_construction_from_terms():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    expected = S.monomial(2,0,0,0).poly - Rational(1,3)*S.monomial(0,1,1,0).poly
    result = S.from_terms({(2,0,0,0): 1, (0,1,1,0): Rational(-1,3)})
    assert result.poly == expected
    assert S.from_terms(dict()) == S.zero()
    with pytest.raises(TypeError):
        S.from_terms({(2,0,0,0): 0.5})

def test_symmetric_element_construction_from_poly():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    expected = SymmetricElement(S,2*S.monomial(2,0,0,0).poly - Rational(1,3)*S.monomial(0,1,1,0).poly)
    gens = S.generators
    assert S.from_poly(Poly(2*(gens[0]**2) - Rational(1,3)*gens[1]*gens[2], gens[0],gens[1],gens[2])) == expected
    with pytest.raises(ValueError):
        S.from_poly(Poly(0.5*(gens[0]**2),gens[0]))
    with pytest.raises(ValueError):
        S.from_poly(Poly(gens[0]**2 - gens[1],gens[0],gens[1]))
    with pytest.raises(ValueError):
        S.from_poly(Poly(1,gens))

def test_symmetric_element_coordinates():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    element = S.from_terms({(2,0,0,0) : 2, (0,1,1,0): Rational(-1,3)})
    expected = ImmutableMatrix([Rational(2),Rational(0),Rational(0),Rational(0),Rational(0),Rational(-1,3),Rational(0),Rational(0),Rational(0),Rational(0)])
    assert S.coordinates(element) == expected
    assert S.from_coordinates(S.coordinates(element)) == element
    assert S.coordinates(S.from_coordinates(expected)) == expected

def test_symmetric_element_arithmetic():
    g = GeneralLinear(2)
    S = SymmetricPower(g,2)
    element = S.from_terms({(2,0,0,0) : 2, (0,1,1,0) : Rational(-1,3)})
    expected = S.from_terms({(2,0,0,0) : 6, (0,1,1,0) : -1})
    assert element + (-element) == S.zero()
    assert element *3  == expected
    assert 3*element == expected
    assert element + element + element == expected