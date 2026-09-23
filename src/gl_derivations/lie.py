"""
TODO: docstring for the Lie module
"""

from gl_derivations._validation import exact_scalar, validate_index
from sympy import ImmutableMatrix, Rational

def kronecker_delta(a,b):
    """
    Implementation of the Rational Kronecker delta. Returns Rational(1) if a == b, otherwise returns Rational(0)
    """
    if a == b:
        return Rational(1)
    return Rational(0)

class GeneralLinear:
    """
    This class holds the information about the Lie algebra gl(n) of n-by-n matrices over the rational numbers.

    Instance attributes:
        n (int): size of the matrices in this Lie algebra.
        dimension (int): dimension of gl(n), returns n**2.
        __basis_tuple: private attribute recording the basis elements as a tuple. Only construct if needed since it allocates n**2 matrices of n-by-n matrices, and then does not compute it again if it is not None.


    Constructor:
        GeneralLinear(val) returns an instance with n = val. Validates that the input val is a valid positive integer, rejects booleans and throws a TypeError at any other type.
    
    Methods:
        __eq__(self,other): gl(n) == gl(m) if and only if n == m.
        __str__(self): returns GeneralLinear(n) as a string.
        zero(self): returns the neutral element.
        E(self,i,j): returns the elementary matrix with 1 in the i-th row, j-th column. Indices run from 0 to self.n.
        basis(self): returns a tuple of self.dimension LieElements representing its canonical ordered basis of elementary matrices.
        coordinates(self,elem): returns a column of size n**2 with the coordinates of the element elem in the given basis.
        from_coordinates(self,values): returns a LieElement with this parent whose matrix has the given coordinates.
    """

    def __init__(self,val):
        """
        Constructor for the GeneralLinear class.

        Arguments:
            val (int): the size of the matrices in this Lie algebra. Must be of type int, and must be a positive integer.
        """
        if (not isinstance(val,int)) or (isinstance(val,bool)):
            raise TypeError(f"Expected int, got {type(val).__name__}")
        if (val < 1):
            raise ValueError(f"Size of matrices in GeneralLinear Lie algebra must be a positive integer, got {val}")
        self.n = val
        self.dimension = val**2
        self.__basis_tuple = None

    def __eq__(self,other):
        """
        Compare two GeneralLinear instances. Returns true if and only if other is also an instance and has the same n.
        """
        if not isinstance(other, GeneralLinear):
            return False
        if self.n != other.n:
            return False
        return True

    def zero(self):
        """
        Creates the zero LieElement in this Lie algebra.

        Return:
            LieElement which has this instance as parent and the zero matrix as its matrix.
        """
        matrix = ImmutableMatrix(self.n, self.n, lambda i,j: Rational(0))
        return LieElement(self,matrix)

    def E(self,i,j):
        """
        Creates the LieElement in this Lie algebra representing the matrix whose only non-zero entry is a one in the i-th row, j-th column.

        Return:
            LieElement which has this instance as parent and the corresponding elementary matrix as underlying matrix.
        """
        validate_index(i,self.n)
        validate_index(j,self.n)
        matrix = ImmutableMatrix(self.n, self.n, lambda a,b: kronecker_delta(i,a)*kronecker_delta(j,b))
        return LieElement(self,matrix)

    def basis(self):
        """
        Creates the tuple of self.n**2 LieElements that form the canonical ordered basis of gl(n). Only does the computation and allocation once, and only when it is first called.

        Return:
            tuple of LieElements which are elementary matrices in the canonical ordered basis having this instance as parent.
        """
        if self.__basis_tuple == None:
            self.__basis_tuple = tuple()
            for i in range(0,self.n):
                for j in range(0, self.n):
                    self.__basis_tuple = self.__basis_tuple + tuple([self.E(i,j)])
        return self.__basis_tuple

    def coordinates(self,elem):
        """
        Returns a matrix of size self.n**2-by-1 with the entries being the coordinates of elem in the canonical ordered basis of self.
        """
        if not isinstance(elem, LieElement):
            raise TypeError(f"Expected LieElement, got {type(elem).__name__}")
        if elem.parent != self:
            raise ValueError(f"Expected matrix of size {self.n}, got {elem.parent.n}")
        return ImmutableMatrix(self.dimension, 1, lambda p,j: elem.matrix[p//self.n, p % self.n])

    def from_coordinates(self,*args):
        """
        Returns a LieElement with this parent and the given values as coordinates.

        Arguments:
            *args (ImmutableMatrix or tuple of self.n**2 scalars): the coordinates in the canonical ordered basis.

        Returns:
            LieElement
        """
        if (len(args) == 1) and (isinstance(args[0], ImmutableMatrix)):
            if not args[0].shape == (self.dimension,1):
                raise IndexError(f"Expected a coordinate matrix of size {self.dimension}-by-1, got {args[0].shape}")
            base = self.basis()
            element = self.zero()
            for i in range(0,self.dimension):
                element = element + args[0][i,0]*base[i]
            return element
        elif (len(args) == self.dimension):
            values = tuple(exact_scalar(x) for x in args)
            base = self.basis()
            element = self.zero()
            for i in range(0, self.dimension):
                element = element + values[i]*base[i]
            return element
        else:
            raise TypeError(f"Expected 1 ImmutableMatrix or {self.dimension} scalars in the arguments, got something else instead.")
        

    def __str__(self):
        return f"GeneralLinear({self.n})"

class LieElement:
    """
    This class holds the information about an element of the Lie algebra gl(n) of n-by-n matrices over the rational numbers.

    Instance attributes:
        parent (GeneralLinear): the ambient Lie algebra this element is a part of.
        matrix (sympy.ImmutableMatrix): the underlying n-by-n matrix representing this element.

    Constructor:
        LieElement(parent,matrix) initializes an instance of this class with the given parent and matrix. It does not perform validation, this constructor is not meant to be called by itself. Instead, use the methods of GeneralLinear to produce Lie elements.
    
    Methods:
        __eq__(self,other): two elements are equal if their parents are the same and their matrices are the same.
        __str__(self): returns its parent and its matrix.
        __add__,__sub__,__neg__,__mul__,__rmul__(self,other): arithmetical operations. Rmul accepts left multiplication by a scalar, mul accepts right multiplication by another Lie element or a scalar.
        bracket(self,other): returns the Lie bracket [self,other].
    """
    def __init__(self,parent,matrix):
        """
        Constructor for the LieElement class. Should not be called by itself, as it does not validate the entries. Instead, create Lie elements through the methods of the GeneralLinear class.

        Arguments:
            parent (GeneralLinear): ambient Lie algebra
            matrix (sympy.ImmutableMatrix): immutable matrix representing this element.
        """
        self.parent = parent
        self.matrix = matrix

    def __eq__(self,other):
        """
        Compares two LieElements. They are the same if they have the same parent and matrix.
        """
        if not isinstance(other,LieElement):
            return False
        return (self.parent == other.parent) and (self.matrix == other.matrix)

    def __add__(self,other):
        """
        Returns new Lie element with added underlying matrices, assuming they are compatible.
        """
        if not isinstance(other, LieElement):
            return NotImplemented
        if other.parent != self.parent:
            raise ValueError(f"Expected matrix of size {self.parent.n}, got {other.parent.n}")
        added_matrix = self.matrix + other.matrix
        return LieElement(self.parent, added_matrix)

    def __neg__(self):
        """
        Returns new Lie element with negated matrix entries.
        """
        return LieElement(self.parent, -self.matrix)

    def __sub__(self,other):
        """
        Returns new Lie element with subtracted underlying matrices, assuming they are compatible.
        """
        if not isinstance(other, LieElement):
            return NotImplemented
        if other.parent != self.parent:
            raise ValueError(f"Expected matrix of size {self.parent.n}, got {other.parent.n}")
        sub_matrix =  self.matrix - other.matrix
        return LieElement(self.parent, sub_matrix)
    
    def __mul__(self,other):
        """
        Returns new Lie element with multiplied underlying matrices or multiplied by a scalar. Accepts only LieElements or scalars on the right hand side.
        """
        if not isinstance(other, LieElement):
            try:
                other = exact_scalar(other)
            except TypeError:
                return NotImplemented
            return LieElement(self.parent, other*self.matrix)
        else:
            if other.parent != self.parent:
                raise ValueError(f"Expected matrix of size {self.parent.n}, got {other.parent.n}")
            mul_matrix = self.matrix * other.matrix
            return LieElement(self.parent, mul_matrix)

    def __rmul__(self,other):
        """
        Called if one requests an operation of the type A*B where B is a LieElement and A has not defined what right multiplication by B means. Only accepts scalars as A.
        """
        try:
            other = exact_scalar(other)
        except TypeError:
            return NotImplemented
        return self*other

    def __str__(self):
        return f"({self.parent}, {self.matrix})"

    def bracket(self,other):
        """
        Returns the Lie bracket [self,other] = self*other - other*self.

        Arguments:
            other (LieElement): should be a Lie element with the same parent.

        Returns:
            LieElement: the result of the bracket.
        """
        if not isinstance(other,LieElement):
            raise TypeError(f"Expected LieElement, got {type(other).__name__}")
        if self.parent != other.parent:
            raise ValueError(f"Expected matrix of size {self.parent.n}, got {other.parent.n}")
        return self*other - other*self