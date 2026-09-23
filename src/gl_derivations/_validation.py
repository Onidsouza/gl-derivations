"""
TODO: Docstring for the validation module.
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