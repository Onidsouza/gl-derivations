"""
TODO: docstring for the Symmetric module
"""

from gl_derivations.lie import GeneralLinear
from gl_derivations._validation import exact_scalar
from sympy import symbols, Poly, Rational, ZZ, QQ, ImmutableMatrix
from math import comb
from itertools import combinations_with_replacement

class SymmetricPower:
    """
    Each instance of this class represents a space S(k) of the k-th symmetric power of gl(n) for some n.

    Instance attributes:

    Constructor:
        SymmetricPower(algebra,degree) (GeneralLinear, int): returns an instance representing the degree-th symmetric power of algebra.

    Methods:
        __eq__(other): spaces are the same if they have the same algebra and same degree
        basis(): returns a tuple of SymmetricElement elements representing the canonical monomial basis of this space in descending lexographic order.
        basis_labels(): returns a tuple of integer tuples representing the monomial degrees of the canonical ordered basis of this space, in descending lexographic order.
        zero(): returns the SymmetricElement representing the zero polynomial.
        monomial(alpha): (tuple of ints) returns the SymmetricElement representing the monomial with the given tuple of exponents. Performs validation.
        from_terms(mapping): (dict from tuples to coefficients) returns the SymmetricElement which has the given coefficient at the given tuple in the supplied dictionary.
        from_poly(poly): (sympy.Poly) validates that the given polynomial makes sense (homogeneous, correct degree, correct variables, valid coefficient domain) and instantiates it as a SymmetricElement.
        coordinates(elem): (SymmetricElement) given a symmetric element, returns the self.dimension-by-1 matrix whose entries are the coordinates of elem in the canonical basis.
        from_coordinates(values): (sympy.ImmutableMatrix) given a matrix of size self.dimension-by-1 whose entries are valid scalars, returns the SymmetricElement with these given coordinates in the canonical ordered basis.
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
            f = SymmetricElement(self,Poly(Rational(1),self.generators,domain=QQ))
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
        return SymmetricElement(self, Poly(0,self.generators,domain=QQ))

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

    def from_terms(self,mapping):
        """
        Returns the SymmetricElement in this space whose coefficients in a given exponent are given by this mapping.

        Arguments:
            mapping: (dict from exponent labels to coefficients). Exponent labels must be a tuple of self.algebra.dimension non-negative integers summing to self.degree. Valides through self.monomial. Coefficients must be int or Rational.
        """
        if not isinstance(mapping,dict):
            raise TypeError(f"Expected a dict, got {type(mapping).__name__}")
        f = Poly(0,self.generators,domain=QQ)
        for key in mapping.keys():
            monom = self.monomial(*key).poly
            coef = exact_scalar(mapping[key])
            f = f+ coef*monom
        return SymmetricElement(self,f)

    def from_poly(self,poly):
        """
        Returns the SymmetricElement in this space whose polynomial is the given one.

        Arguments:
            poly: (sympy.Poly) must have 'QQ' domain and self.generators as its variables. Must be homogeneous of degree self.degree.
        """
        if not isinstance(poly,Poly):
            raise TypeError(f"Expected sympy.Poly, got {type(poly).__name__}.")
        if not (poly.domain == ZZ or poly.domain == QQ):
            raise ValueError(f"Expected sympy.Poly with ZZ or QQ coefficients, got {poly.domain}")
        if not poly.is_homogeneous:
            raise ValueError(f"Polynomial is not homogeneous.")
        if not poly.homogeneous_order() == self.degree:
            raise ValueError(f"Expected polynomial of degree {self.degree}, got {poly.homogeneous_order()}")
        if not (Poly(poly,self.generators).domain in (ZZ,QQ)):
            raise ValueError(f"Polynomial has variables not compatible with this domain: {Poly(poly,self.generators).domain}")
        return SymmetricElement(self,Poly(poly,self.generators,domain=QQ))

    def coordinates(self,elem):
        """
        Returns the matrix of coordinates of the given element in the canonical ordered basis of self.

        Arguments:
            elem (SymmetricElement): the element being converted. Must have self as its parent.

        Returns:
            (sympy.ImmutableMatrix) of shape self.dimension-by-1 with sympy.Rational entries.
        """
        if not isinstance(elem,SymmetricElement):
            raise TypeError(f"Expected SymmetricElement, got {type(elem).__name__}")
        if elem.parent != self:
            raise ValueError(f"This symmetric element does not belong to this symmetric power.")
        return ImmutableMatrix([exact_scalar(elem.poly.coeff_monomial(monom.poly.as_expr())) for monom in self.basis()])

    def from_coordinates(self,*args):
        """
        Returns the matrix of coordinates of the given element in the canonical ordered basis of self.

        Arguments:
            *args is either a single ImmutableMatrix of size self.dimension-by-1 or a tuple of self.dimension scalars.

        Returns:
            (SymmetricElement) represented by the given coefficient matrix in the canonical basis.
        """
        if (len(args) == 1) and (isinstance(args[0],ImmutableMatrix)):
            if not args[0].shape == (self.dimension,1):
                raise IndexError(f"Expected a coordinate matrix of size {self.dimension}-by-1, got {args[0].shape}")
            values = tuple(exact_scalar(args[0][i,0]) for i in range(0,self.dimension))
        elif (len(args) == self.dimension):
            values = tuple(exact_scalar(x) for x in args)
        else:
            raise TypeError(f"Expected 1 ImmutableMatrix or {self.dimension} scalars in the argument, got something else instead.")
        base = self.basis()
        element = Poly(0,self.generators,domain=QQ)
        for i in range(0,self.dimension):
            element = element + values[i]*base[i].poly
        return SymmetricElement(self,element)

class SymmetricElement:
    """
    Each instance of this class represents an element of the space represented by a SymmetricPower element.

    Instance attributes:
        parent: (SymmetricPower) the space this element belongs to.
        poly: (sympy.Poly) the underlying polynomial this element represents.

    Constructor:
        SymmetricElement(parent,poly): (SymmetricPower, sympy.Poly) initializes this object with the given parent and polynomial. It does not perform validation, it should not be called by itself.

    Methods:
        __eq__(other): elements are the same if they have the same parent and poly.
        __str__: returns (self.parent, self.poly) tuple.
        __add__,__sub__,__neg__,__mul__,__rmul__(self,other): arithmetical operations. mul and rul accepts multiplication by a scalar.
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
        return f"({self.parent}, {self.poly.as_expr()})"

    def __eq__(self,other):
        """
        Two elements are equal if they share the same space and the same underlying polynomial.
        """
        if not isinstance(other,SymmetricElement):
            return False
        return (self.parent == other.parent) and (self.poly == other.poly)

    def __add__(self,other):
        """
        Returns new symmetric element with added underlying polynomials, assuming they are compatible.
        """
        if not isinstance(other,SymmetricElement):
            return NotImplemented
        if self.parent != other.parent:
            return ValueError(f"Cannot add elements of different symmetric powers.")
        return SymmetricElement(self.parent,self.poly + other.poly)
    
    def __neg__(self):
        """
        Returns the additive negation of itself.
        """
        return SymmetricElement(self.parent,-self.poly)
    
    def __sub__(self,other):
        """
        Returns the new symmetric element with subtracted underlying polynomials, assuming they are compatible
        """
        return self + (-other)

    def __mul__(self,other):
        """
        Multiplies itself by a scalar.
        """
        try:
            result = SymmetricElement(self.parent, self.poly * exact_scalar(other))
        except TypeError:
            return NotImplemented
        return result

    def __rmul__(self,other):
        """
        Multiplies itself by a scalar.
        """
        return self*other