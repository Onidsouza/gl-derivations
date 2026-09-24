"""
TODO: docstrings for this Lie test module
"""

from gl_derivations import lie,symmetric,actions
from sympy import Rational, Poly, ImmutableMatrix
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