"""
TODO: docstrings for this validation test module
"""

from gl_derivations._validation import exact_scalar, validate_index
from sympy import Rational
import pytest

@pytest.mark.parametrize("integer_scalar_input,expected_integer_scalar_output", [
    (7,Rational(7)),
    (-43, Rational(-43)),
    (0, Rational(0)),
    (53,Rational(53))
])

def test_integer_becomes_exact_scalar(integer_scalar_input, expected_integer_scalar_output):
    """
    TODO: docstring for this test function.
    """
    value = exact_scalar(integer_scalar_input)
    assert (value == expected_integer_scalar_output) and isinstance(value,Rational)

def test_bool_is_rejected():
    """
    TODO: docstring for this test function.
    """
    with pytest.raises(TypeError):
        exact_scalar(True)

@pytest.mark.parametrize("float_scalar_input", [
    (0.5),
    (-1/3),
    (0.0),
    (2.0),
    (57/1.0)
])

def test_float_is_rejected(float_scalar_input):
    """
    TODO: docstring for this test function.
    """
    with pytest.raises(TypeError):
        exact_scalar(float_scalar_input)

@pytest.mark.parametrize("index,max_range", [
    (0,1),
    (2,3),
    (0,6),
    (8,9)
])

def test_index_is_valid(index,max_range):
    """
    TODO: docstring for this test function.
    """
    assert validate_index(index,max_range)

@pytest.mark.parametrize("index,max_range", [
    (True,7),
    (Rational(4),5),
    (2.0,3)
])

def test_index_type_is_invalid(index,max_range):
    """
    TODO: docstring for this test function.
    """
    with pytest.raises(TypeError):
        validate_index(index,max_range)

@pytest.mark.parametrize("index,max_range", [
    (1,1),
    (-1,7),
    (13,5)
])

def test_index_value_is_invalid(index,max_range):
    """
    TODO: docstring for this test function.
    """
    with pytest.raises(IndexError):
        validate_index(index,max_range)