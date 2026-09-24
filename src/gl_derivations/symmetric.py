"""
TODO: docstring for the Symmetric module
"""

from gl_derivations.lie import GeneralLinear
from sympy import symbols, Poly
from math import comb
from itertools import combinations_with_replacement

class SymmetricPower:
    """
    Each instance of this class represents a space S(k) of the k-th symmetric power of gl(n) for some n.

    Instance attributes:

    Constructor:
        SymmetricPower(algebra,degree) (GeneralLinear, int): returns an instance representing the degree-th symmetric power of algebra.

    Methods:
        basis(): returns a tuple of SymmetricElement elements representing the canonical monomial basis of this space in descending lexographic order.
        basis_labels(): returns a tuple of integer tuples representing the monomial degrees of the canonical ordered basis of this space, in descending lexographic order.
        zero(): returns the SymmetricElement representing the zero polynomial.
        monomial(alpha): (tuple of ints) returns the SymmetricElement representing the monomial with the given tuple of exponents. Performs validation.
    """

    def __init__(self,algebra,degree):
        """
        Constructor for the SymmetricPower class

        Arguments:
            algebra (GeneralLinear): the Lie algebra gl(n) of which this instance is a symmetric power of.
            degree (int): which degree polynomials are we considering.
        """
        if not isinstance(algebra,GeneralLinear):
            raise TypeError(f"Expected GeneralLinear object, got {type(algebra).__name__}.")
        if (not isinstance(degree,int)) or (isinstance(degree,bool)):
            raise TypeError(f"Expected int, got {type(degree).__name__}")
        if (degree < 0):
            raise ValueError(f"Degree must be a non-negative integer, got {degree}")
        self.algebra = algebra
        self.degree = degree
        self.dimension = comb(self.algebra.dimension + self.degree -1, self.degree)
        self.generators = tuple()
        for i in range(0,self.algebra.n):
            for j in range(0,self.algebra.n):
                self.generators = self.generators + (symbols(f"z_{i}_{j}"),)
        self.__basis_labels_tuple = None
        self.__basis_tuple = None

    def basis(self):
        """
        Returns a tuple of Polys for the canonical monomial basis of this space, ordered by descending lexicographic order. Computes lazily and stores in cache.
        """
        if not self.__basis_tuple is None:
            return self.__basis_tuple
        self.__basis_tuple = tuple()
        monomial_tuples = combinations_with_replacement(self.generators,self.degree)
        for tup in monomial_tuples:
            f = SymmetricElement(self,Poly(1,self.generators))
            for monom in tup:
                f.poly = f.poly*monom
            self.__basis_tuple = self.__basis_tuple + (f,)
        return self.__basis_tuple

    def basis_labels(self):
        """
        Returns a tuple of exponent tuples for the canonical monomial basis of this space, ordered by descending lexcographic order.

        Returns:
            A tuple of length self.dimension. Each element of this tuple is a tuple of of self.algebra.dimension integers representing the exponents of a monomial in the given ordered basis. Each such tuple sums to self.degree. Computes it lazily once, and stores in cache.
        """
        if not self.__basis_labels_tuple is None:
            return self.__basis_labels_tuple
        self.__basis_labels_tuple = tuple()
        for f in self.basis():
            self.__basis_labels_tuple = self.__basis_labels_tuple + (f.poly.monoms()[0],)
        return self.__basis_labels_tuple

    def zero(self):
        """
        Returns the SymmetricElement representing the zero polynomial in this space.
        """
        return SymmetricElement(self, Poly(0,self.generators))

    def __eq__(self,other):
        """
        Two symmetric powers are the same if they have the same underlying algebra and the same degree.
        """
        if not isinstance(other,SymmetricPower):
            return False
        return (self.algebra == other.algebra) and (self.degree == other.degree)

    def monomial(self,*args):
        """
        Returns the monomial in this space whose exponents are the given tuple of integers.

        Arguments:
            A tuple of self.algebra.dimension non-negative integers. Are checked against self.basis_labels() for validation.

        Returns:
            (SymmetricElement) representing the given monomial.
        """
        if len(args) != self.algebra.dimension:
            raise TypeError(f"Expected {self.algebra.dimension} arguments, received {len(args)} instead.")
        if (any(isinstance(x,bool) for x in args)) or (not all(isinstance(x,int) for x in args)):
            raise TypeError(f"Expected a list of int, got one entry which is not an int.")
        try:
            return self.basis()[self.basis_labels().index(args)]
        except ValueError:
            raise ValueError(f"Expected a valid tuple of {self.algebra.dimension} non-negative integers summing to {self.degree}, got {args}.")

class SymmetricElement:
    """
    Each instance of this class represents an element of the space represented by a SymmetricPower element.

    Instance attributes:
        parent: (SymmetricPower) the space this element belongs to.
        poly: (sympy.Poly) the underlying polynomial this element represents.

    Constructor:
        SymmetricElement(parent,poly): (SymmetricPower, sympy.Poly) initializes this object with the given parent and polynomial. It does not perform validation, it should not be called by itself.

    Methods:
    """

    def __init__(self,parent,poly):
        """
        Constructor for this class. Does not perform validation, should not be called by itself. Use the SymmetricPower methods to construct SymmetricElements.

        Arguments:
            parent: (SymmetricPower) the ambient space this object lives in.
            poly: (sympy.Poly) the underlying representation
        """
        self.parent = parent
        self.poly = poly

    def __str__(self):
        return f"({self.parent}, {self.poly.as_expr})"

    def __eq__(self,other):
        """
        Two elements are equal if they share the same space and the same underlying polynomial.
        """
        if not isinstance(other,SymmetricElement):
            return False
        return (self.parent == other.parent) and (self.poly == other.poly)