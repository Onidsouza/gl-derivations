"""Symmetric powers of general linear Lie algebras over the rationals.

Represents S^k(gl(n, Q)) using homogeneous polynomials in commuting
symbols corresponding to the matrix units. Provides monomial bases,
coordinate conversion, element arithmetic, and grouping of basis
monomials by their weights for the diagonal adjoint action.

Classes:
    SymmetricPower:
        Parent representing a fixed symmetric power of a GeneralLinear
        algebra. Provides element constructors, ordered monomial bases,
        coordinate conversion, and weight lookup.
    SymmetricElement:
        An element of a SymmetricPower, represented by a SymPy Poly.
        Supports equality, addition, subtraction, negation, and scalar
        multiplication.

Conventions:
    The symbol z_i_j represents the matrix unit E(i, j). Generators are
    ordered first by row i and then by column j, using zero-based indices.

    An exponent tuple has n**2 nonnegative integer entries in generator
    order, summing to k. The canonical basis consists of ordinary
    monomials with coefficient one, ordered by descending lexicographic
    order of their exponent tuples.

    Coordinate vectors are columns in this ordered monomial basis.

    A weight tuple has n integer entries. For an exponent array alpha,
    its r-th entry is the sum of row r minus the sum of column r.
    This is the eigenvalue of E(r, r) under the induced adjoint action.
    Weight entries may be negative and always sum to zero.

    SymmetricPower parents compare equal when their underlying algebras
    compare equal and their degrees agree. Element addition and
    subtraction require equal parents.

Notes:
    The degree k may be zero, in which case the space consists of
    constants and has basis (1,).

    Construct elements through the SymmetricPower methods. The
    SymmetricElement constructor stores its arguments without validation.

    Basis elements and the weight dictionary are cached and shared.
    Treat returned basis elements and the weight dictionary as read-only.
"""

from itertools import combinations_with_replacement
from math import comb

from sympy import QQ, ZZ, ImmutableMatrix, Poly, Rational, symbols

from ._validation import exact_scalar
from .lie import GeneralLinear


class SymmetricPower:
    """
    Each instance of this class represents a space S(k) of the k-th symmetric power of gl(n) for some n.

    Instance attributes:
        algebra (GeneralLinear): the underlying algebra we are taking a power of.
        degree (int): the degree of the symmetric power we are taking, a non-negative integer.
        dimension (int): the dimension of S(k), given by comb(n**+k-1,k) where n**2 is the size of the matrix algebra gl(n).
        generators (tuple of symbols): a tuple of the sympy symbols used to generate this monomials.

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
        weight_from_exponent_tuple(*args): (tuple of self.algebra.n non-negative integers) returns the tuple of self.algebra.n integers representing the weight of a given monomial represented by the given exponent tuple.
        weight_elements(*args): (tuple of self.algebra.n integers) returns a tuple of SymmetricElements in the canonical basis that are of the given weight.
        weight_element_indices(*args): (tuple of self.algebra.n non-negative integers) returns a tuple of indices in the canonical basis that are of the given weight.
        weight_dictionary(): returns a dictionary from tuples of ints (weights) to tuples of exponents with that given weight.
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
        for i in range(self.algebra.n):
            for j in range(self.algebra.n):
                self.generators = self.generators + (symbols(f"z_{i}_{j}"),)
        self.__basis_labels_tuple = None
        self.__basis_tuple = None
        self.__monomial_weight_dict = None

    def basis(self):
        """
        Returns a tuple of SymmetricElements for the canonical monomial basis of this space, ordered by descending lexicographic order. Computes lazily and stores in cache.
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
        Returns a tuple of exponent tuples for the canonical monomial basis of this space, ordered by descending lexographic order.

        Returns:
            A tuple of length self.dimension. Each element of this tuple is a tuple of self.algebra.dimension integers representing the exponents of a monomial in the given ordered basis. Each such tuple sums to self.degree. Computes it lazily once, and stores in cache.
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
            self.algebra.dimension non-negative integers. Are checked against self.basis_labels() for validation.

        Returns:
            (SymmetricElement) representing the given monomial.
        """
        if len(args) != self.algebra.dimension:
            raise TypeError(f"Expected {self.algebra.dimension} arguments, received {len(args)} instead.")
        if (any(isinstance(x,bool) for x in args)) or (not all(isinstance(x,int) for x in args)):
            raise TypeError("Expected a list of int, got one entry which is not an int.")
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
            raise ValueError("Polynomial is not homogeneous.")
        if (not poly.homogeneous_order() == self.degree) and (poly != 0):
            raise ValueError(f"Expected polynomial of degree {self.degree}, got {poly.homogeneous_order()}")
        if not (Poly(poly,self.generators).domain in (ZZ,QQ)):
            raise ValueError(f"Polynomial has variables not compatible with this domain: {Poly(poly,self.generators).domain}")
        return SymmetricElement(self,Poly(poly,self.generators,domain=QQ))

    def coordinates(self,elem):
        """
        Returns the matrix of coordinates of the given element in the canonical ordered basis of self.

        Arguments:
            elem (SymmetricElement): the element being converted. Must compare == to self as its parent.

        Returns:
            (sympy.ImmutableMatrix) of shape self.dimension-by-1 with sympy.Rational entries.
        """
        if not isinstance(elem,SymmetricElement):
            raise TypeError(f"Expected SymmetricElement, got {type(elem).__name__}")
        if elem.parent != self:
            raise ValueError("This symmetric element does not belong to this symmetric power.")
        return ImmutableMatrix([exact_scalar(elem.poly.coeff_monomial(monom.poly.as_expr())) for monom in self.basis()])

    def from_coordinates(self,*args):
        """
        Returns the element whose coordinates in the canonical basis are the given ones.

        Arguments:
            *args is either a single ImmutableMatrix of size self.dimension-by-1 or separate self.dimension scalars.

        Returns:
            (SymmetricElement) represented by the given coefficient matrix in the canonical basis.
        """
        if (len(args) == 1) and (isinstance(args[0],ImmutableMatrix)):
            if not args[0].shape == (self.dimension,1):
                raise IndexError(f"Expected a coordinate matrix of size {self.dimension}-by-1, got {args[0].shape}")
            values = tuple(exact_scalar(args[0][i,0]) for i in range(self.dimension))
        elif (len(args) == self.dimension):
            values = tuple(exact_scalar(x) for x in args)
        else:
            raise TypeError(f"Expected 1 ImmutableMatrix or {self.dimension} scalars in the argument, got something else instead.")
        base = self.basis()
        element = Poly(0,self.generators,domain=QQ)
        for i in range(self.dimension):
            element = element + values[i]*base[i].poly
        return SymmetricElement(self,element)

    def weight_from_exponent_tuple(self,*args):
        """
        Returns the tuple representing how the diagonals act on a monomial with a given exponent tuple on a SymmetricPower.

        Arguments:
        *args (int): the monomial exponents. Must have self.algebra.dimension separate entries.

        Returns:
        tuple of self.algebra.n integers representing how each diagonal element acts.
        """
        if len(args) != self.algebra.dimension:
            raise ValueError(f"Exponent tuple must have length {self.algebra.dimension}, got {len(args)}")
        if (any(isinstance(x,bool) for x in args)) or (not all(isinstance(x,int) for x in args)):
            raise TypeError("Expected a list of int, got one entry which is not an int.")
        if sum(args) != self.degree:
            raise ValueError(f"Degrees should sum up to {self.degree}, instead they add to {sum(args)}")
        if not all(x >= 0 for x in args):
            raise ValueError("Degrees should be non-negative integers.")
        weights = [0] * self.algebra.n
        for index, exponent in enumerate(args):
            # Each exponent correspond to an element of the canonical ordered basis of g. The index of this exponent tells us which matrix element we are looking at. We can explicitly detect E(i,j) from g.n and this index alone.
            i = index // self.algebra.n
            j = index % self.algebra.n
            weights[i] = weights[i] + exponent
            weights[j] = weights[j] - exponent
        return tuple(weights)

    def weight_elements(self,*args):
        """
        For each tuple of self.algebra.n integers summing to zero, returns a tuple containing the SymmetricElements in the basis whose weight is the given one.

        Arguments:
        args: tuple of integers summing to zero.

        Returns:
        tuple of SymmetricElements in the basis with the given weight.
        """
        if tuple(args) in self.weight_dictionary().keys():
            result = tuple()
            for monom in self.weight_dictionary()[tuple(args)]:
                result = result + (self.monomial(*monom),)
            return result
        return tuple()

    def weight_element_indices(self,*args):
        """
        For each tuple of self.algebra.n integers summing to zero, returns a tuple containing the indices in the ordered basis whose weight is the given one.

        Arguments:
        args: tuple of integers summing to zero.

        Returns:
        tuple of indices in the basis with the given weight. Returns the empty tuple if that weight is absent in this space.
        """
        return tuple(self.basis().index(x) for x in self.weight_elements(*args))

    def weight_dictionary(self):
        """
        Computes lazily, caches and store a dict from weights (tuples of integers) to tuples of exponents with that given weight.
        """
        if self.__monomial_weight_dict is None:
            temporary_monomial_weight_dict = dict()
            for label in self.basis_labels():
                weight = self.weight_from_exponent_tuple(*label)
                if weight not in temporary_monomial_weight_dict:
                    temporary_monomial_weight_dict[weight] = [label]
                else:
                    temporary_monomial_weight_dict[weight].append(label)
            self.__monomial_weight_dict = dict()
            for key in temporary_monomial_weight_dict:
                self.__monomial_weight_dict[key] = tuple(temporary_monomial_weight_dict[key])
        return self.__monomial_weight_dict

    def __str__(self):
        return f"S^{self.degree}(gl({self.algebra.n}))"

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
            raise ValueError("Cannot add elements of different symmetric powers.")
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
        if not isinstance(other,SymmetricElement):
                    return NotImplemented
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