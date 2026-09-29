"""Internal validation helpers for scalar coefficients and matrix indices.

Functions:
    exact_scalar:
        Convert a supported exact numeric value to a SymPy rational.
        Accepts Python integers, fractions.Fraction, and sympy.Rational
        instances, including SymPy integers. Rejects booleans and
        unsupported types, including floats and strings.
    validate_index:
        Check that an index is a Python integer, excluding booleans,
        within a specified zero-based range. Return True on success
        and raise an exception for an invalid index.

Notes:
    Scalar coefficients and indices have different type requirements:
    SymPy integers are accepted as coefficients but not as indices.

    The caller of validate_index is responsible for supplying a valid
    exclusive upper bound; the bound itself is not validated.
"""
from fractions import Fraction

from sympy import Rational


def exact_scalar(value):
    """
    This function converts its argument to the scalar data structure used for this package, which is sympy's Rational. It perfoms some validation to ensure that the input is valid, and returns it in the proper format.
    
    Args:
        value (int, fractions.Fraction, sympy.Rational):  the value to be converted to a valid sympy.Rational.
    
    Returns:
        sympy.Rational: the converted value.
    """
    if (not isinstance(value, (int,Fraction,Rational))) or (isinstance(value, bool)):
        raise TypeError(f"Expected int, fractions.Fraction or sympy.Rational, got {type(value).__name__}")
    return Rational(value)

def validate_index(i,max):
    """
    This function validates that the entry is a valid integer between 0 and max-1, so it can be used as index.
    """
    if (not isinstance(i, int)) or (isinstance(i, bool)):
        raise TypeError(f"Expected int, got {type(i).__name__}")
    if (i < 0) or (i >= max):
        raise IndexError(f"Invalid index in the given range 0..{max-1}, got {i}")
    return True