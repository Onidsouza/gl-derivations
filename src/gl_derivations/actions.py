"""
TODO: docstring for the Lie module
"""

from gl_derivations.lie import GeneralLinear, LieElement
from gl_derivations.symmetric import SymmetricPower, SymmetricElement
from sympy import symbols, Poly, QQ, ImmutableMatrix

def __adjoint_action_on_generator(lie_element,gens,pos):
    """
    Given a generator of the Symmetric algebra and a Lie element, compute the adjoint action of the latter on the former.

    Arguments:
        lie_element (LieElement): the element of gl(n) that will be acting.
        gens: the ordered basis of generators that will be used for the polynomial
        pos: the position of the generator in the list.

    Return:
        (sympy.Poly) on gens and domain QQ that represents the result of the adjoint action.
    """
    i = pos // lie_element.parent.n
    j = pos % lie_element.parent.n
    result = Poly(0,gens,domain=QQ)
    for a in range(0,lie_element.parent.n):
        result = result + lie_element.matrix[a,i]*gens[a*lie_element.parent.n + j] # adds the result of multiplying the given generator on the left by lie_element
    for b in range(0,lie_element.parent.n):
        result = result - lie_element.matrix[j,b]*gens[i*lie_element.parent.n + b] # subtracts the result of multiplying the given generator on the right by lie_element
    return result

def act(x,v):
    """
    Given a Lie element x and a Symmetric element v, returns the result of the adjoint action ad(x)(v). Does validation to check that the action is well-defined.

    Arguments:
        x: (LieElement) the element that is acting.
        v: (SymmetricElement) the element being acted upon.

    Returns:
        (SymmetricElement): the result of the adjoint action.
    """
    if (not isinstance(x,LieElement)) or (not isinstance(v,SymmetricElement)):
        raise TypeError(f"Expected (LieElement,SymmetricElement) argument pair, got ({type(x).__name__},{type(v).__name__})")
    if x.parent != v.parent.algebra:
        raise ValueError(f"Adjoint action of {x.parent} is not defined on symmetric powers of {v.parent.algebra}")
    gens = v.parent.generators
    result = Poly(0,gens,domain=QQ)
    for letter in gens:
        derivative = v.poly.diff(letter)
        if derivative == 0:
            continue # no need to compute action if dv/dz_i_j = 0. We are summing ad(x,z_i_j) * dv/dz_i_j over all i,j.
        result = result + derivative*__adjoint_action_on_generator(x,gens,gens.index(letter))
    return v.parent.from_poly(result)

def action_matrix(x,V):
    """
    Given a Lie element x and a SymmetricPower V, returns the matrix representing the adjoint action of x on V. 

    Arguments:
        x: (LieElement) the element that is acting.
        V: (SymmetricPower) the space where the element is acting.

    Returns:
        (ImmutableMatrix): square of size V.dimension with rational entries, whose columns are the coordinates of x acting on the canonical ordered basis.
    """
    if not isinstance(x,LieElement):
        raise TypeError(f"Expected LieElement, got {type(x).__name__}")
    if not isinstance(V,SymmetricPower):
        raise TypeError(f"Expected SymmetricPower, got {type(V).__name__}")
    if x.parent != V.algebra:
        raise ValueError(f"Adjoint action of {x.parent} is not defined on symmetric powers of {V.algebra}")
    base = V.basis()
    result = ImmutableMatrix([0] * V.dimension)
    for element in base:
        result = result.col_insert(base.index(element)+1,V.coordinates(act(x,element)))
    return ImmutableMatrix(result.col_del(0))