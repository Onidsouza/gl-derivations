"""Rational subspaces of general linear Lie algebras and symmetric powers.

Provides subspace construction from spanning elements or homogeneous
linear equations, membership and equality tests, and standard matrix
subalgebras. For symmetric-power ambients, also provides diagonal
stability tests, weight decomposition, invariant vectors, evaluation
kernels, and intersection with symmetric powers of traceless matrices.

Classes:
    Subspace:
        A linear subspace of a GeneralLinear or SymmetricPower parent,
        represented by an independent basis in ambient coordinates.

Conventions:
    All coordinates refer to the canonical ordered basis of the ambient
    parent, and all computations use exact rational arithmetic.

    If the ambient dimension is N and the subspace dimension is d,
    basis_matrix has shape (N, d). Its columns are the ambient
    coordinate vectors of a basis of the subspace.

    equation_matrix() has shape (N - d, N). Its rows encode independent
    homogeneous linear equations on ambient coordinate columns.

    Multiplication by basis_matrix converts coordinates in the stored
    subspace basis into ambient coordinates.

    Weight tuples record eigenvalues of the diagonal matrix units,
    using the convention of SymmetricPower. Only weight spaces of
    positive dimension are included in a decomposition; the zero
    weight is included when its weight space is nonzero.

    Invariant vectors are annihilated by every element of the acting
    subspace. Evaluation uses the trace-form convention implemented
    in actions.evaluate.

    Returned subspaces and weight-decomposition components retain the
    original ambient parent.

Notes:
    Construct subspaces through from_basis, from_equations, whole,
    trivial, D, D0, or T. The direct constructor performs no validation.

    Spanning elements are reduced to a basis, so the stored basis need
    not be the sequence supplied to from_basis.

    Equation matrices are computed lazily and cached. Treat ambient,
    basis_matrix, and dimension as read-only after construction.
"""

from sympy import ImmutableMatrix, Rational

from ._validation import exact_scalar
from .actions import action_matrix, evaluate
from .lie import GeneralLinear, LieElement, kronecker_delta
from .symmetric import SymmetricElement, SymmetricPower


class Subspace:
    """
    Each instance of this class represents a subspace of either GeneralLinear or SymmetricPower.

    Constructor:

    Attributes:
    ambient (GeneralLinear or SymmetricPower): the ambient space this subspace is contained within
    basis_matrix (ImmutableMatrix): a matrix of size ambient.dimension-by-self.dimension representing a basis of this subspace. Each column is one vector in the basis.
    __equation_matrix (ImmutableMatrix): a matrix of size (ambient.dimension-self.dimension)-by-(ambient.dimension) that represents this space as its nullspace.
    dimension (int): the dimension of this space.

    Methods:
    equation_matrix(self): Computes lazily, caches and returns __equation_matrix
    __eq__(self,other): two subspaces are the same if they have the same ambient, same dimension, and self.equation_matrix()*other.basis_matrix is the zero matrix.
    contains(self,element): (LieElement or SymmetricElement) checks if the given element belongs to this subspace.
    is_diagonal_stable(self): checks if itself is stable under the action of diagonal matrices (the subspace Subspace.D(self.ambient.algebra)). Only valid if its parent is a SymmetricPower.
    weight_decomposition(self): returns a dict from weights to Subspaces corresponding to non-zero weight spaces in this decomposition. Only valid if its parent is a SymmetricPower and self.is_diagonal_stable() returns True.
    basis(): returns a tuple of self.dimension LieElements or SymmetricElements which are a basis of this subspace.
    invariants(self,X): (Subspace with GeneralLinear as parent) returns the subspace of all elements invariant under the LieElements in X.
    evaluation_kernel(self,x): (LieElement) use the trace form to return the subspace of all elements that vanish when evaluated at x.
    intersection_with_symmetric_sl(self): returns the subspace that is the intersection of self with the canonical copy of S^k(sl(n)) inside S^k(gl(n))

    Factory methods:
    from_basis(ambient,elements): (GeneralLinear or SymmetricPower, tuple of LieElements or SymmetricElements) returns a new instance of subspace where it lies in the given ambient space and is spanned by the given elements. Performs validation.
    from_equations(ambient,equation_matrix): (GeneralLinear or SymmetricPower, ImmutableMatrix) returns the space given by the nullspace of the given matrix. Performs validation.
    D(algebra): returns the Subspace of diagonal matrices
    D0(algebra): returns the Subspace of traceless diagonal matrices
    T(algebra): returns the Subspace of strictly upper triangular matrices
    whole(ambient): (GeneralLinear or SymmetricPower) returns an instance of this class representing the whole ambient space as its subspace.
    trivial(ambient): (GeneralLinear or SymmetricPower) returns an instance of this class representing the trivial zero-dimensional subspace.
    """

    def __init__(self,ambient,dimension,basis_matrix):
        """
        Constructor for this class. Requires ambient space, dimension and a basis matrix. Does not provide validation, use class methods to construct elements safely.
        """
        self.ambient = ambient
        self.basis_matrix = basis_matrix
        self.__equation_matrix = None
        self.dimension = dimension

    @staticmethod
    def from_basis(ambient,*args):
        """
        Construct the subspace with the given ambient space followed by a list of elements that span it. Does validation.

        Arguments:
        ambient (GeneralLinear or SymmetricPower): the ambient space this subspace lives in.
        *args (tuple of LieElements or SymmetricElements): the elements that span the subspace.

        Returns:
        Subspace: the one with ambient as parent and spanned by the arguments.
        """
        if not (isinstance(ambient,GeneralLinear) or isinstance(ambient,SymmetricPower)):
            raise TypeError(f"Expected GeneralLinear or Symmetric power as ambient, got {type(ambient).__name__}")
        elements = list(args)
        if len(elements) == 0:
            return Subspace(ambient,0,ImmutableMatrix(ambient.dimension, 0, lambda i,j: Rational(0)))
        if not all(isinstance(elem,LieElement) or isinstance(elem,SymmetricElement) for elem in elements):
            raise TypeError("Not all arguments are either LieElement or SymmetricElement")
        if not all(elem.parent == ambient for elem in elements):
            raise ValueError("Not all elements belong to the given ambient space.")
        basis_matrix = ImmutableMatrix([Rational(0)]* ambient.dimension)
        for elem in elements:
            basis_matrix = basis_matrix.col_insert(elements.index(elem)+1,ambient.coordinates(elem))
        list_of_basis_vectors = basis_matrix.col_del(0).transpose().rref()[0].transpose().columnspace()
        basis_matrix = ImmutableMatrix([Rational(0)]* ambient.dimension)
        for column in list_of_basis_vectors:
            basis_matrix = basis_matrix.col_insert(list_of_basis_vectors.index(column)+1,column)
        return Subspace(ambient,basis_matrix.shape[1]-1,ImmutableMatrix(basis_matrix.col_del(0)))

    @staticmethod
    def from_equations(ambient,equation_matrix):
        """
        Construct the subspace given by the nullspace of the given equation_matrix. Does validation.

        Arguments:
        ambient (GeneralLinear or SymmetricPower): the ambient space this subspace lives in.
        equation_matrix (ImmutableMatrix): a matrix with ambient.dimension columns whose rows encode the homogeneous linear equations this space must satisfy.

        Returns:
        Subspace: the subspace object representing the nullspace of the given equation matrix.
        """
        if not (isinstance(ambient,GeneralLinear) or isinstance(ambient,SymmetricPower)):
            raise TypeError(f"Expected GeneralLinear or Symmetric power as ambient, got {type(ambient).__name__}")
        if not isinstance(equation_matrix, ImmutableMatrix):
            raise TypeError(f"Expected ImmutableMatrix for equations, got {type(equation_matrix).__name__}")
        if not (equation_matrix.shape[1] == ambient.dimension):
            raise ValueError(f"Expected matrix with {ambient.dimension} columns, got {equation_matrix.shape[1]}")
        for entry in equation_matrix:
            try:
                entry = exact_scalar(entry) # Doesn't change anything, just tries to catch type errors.
            except TypeError:
                raise TypeError(f"Expected matrix with valid scalar entries, got {type(entry).__name__}")
        if equation_matrix.shape[0] == 0:
            # if there are no equations, i.e. no rows
            return Subspace(ambient,ambient.dimension,ImmutableMatrix(ambient.dimension,ambient.dimension,lambda i,j: kronecker_delta(i,j)))
        list_of_basis_coordinates = equation_matrix.nullspace()
        tuple_of_basis_vectors = (ambient.from_coordinates(x) for x in list_of_basis_coordinates)
        return Subspace.from_basis(ambient,*tuple_of_basis_vectors)

    @staticmethod
    def whole(ambient):
        """
        returns an instance of this class representing the whole ambient space as its subspace.
        """
        if not (isinstance(ambient,GeneralLinear) or isinstance(ambient,SymmetricPower)):
            raise TypeError(f"Expected GeneralLinear or Symmetric power as ambient, got {type(ambient).__name__}")
        return Subspace.from_equations(ambient,ImmutableMatrix(0,ambient.dimension,lambda i,j: 0))
    
    @staticmethod
    def trivial(ambient):
        """
        returns an instance of this class representing the trivial zero-dimensional subspace.
        """
        return Subspace.from_basis(ambient)

    def equation_matrix(self):
        """
        Computes lazily, caches and returns the (self.ambient.dimension-self.dimension)-by-(self.ambient.dimension) matrix whose nullspace is this subspace.
        """
        if self.__equation_matrix is None:
            list_of_equations = self.basis_matrix.transpose().nullspace()
            self.__equation_matrix = ImmutableMatrix([Rational(0)]*self.ambient.dimension).transpose()
            for row in list_of_equations:
                self.__equation_matrix = self.__equation_matrix.row_insert(list_of_equations.index(row)+1,row.transpose())
            list_of_reduced_equations = self.__equation_matrix.rref()[0].rowspace()
            self.__equation_matrix = ImmutableMatrix([Rational(0)]*self.ambient.dimension).transpose()
            for row in list_of_reduced_equations:
                self.__equation_matrix = self.__equation_matrix.row_insert(list_of_reduced_equations.index(row)+1,row)
            self.__equation_matrix = ImmutableMatrix(self.__equation_matrix.row_del(0))
        return self.__equation_matrix

    def __eq__(self,other):
        """
        Two subspaces are the same if they have the same ambient, same dimension, and self.equation_matrix()*other.basis_matrix is a zero matrix.
        """
        if not isinstance(other,Subspace):
            return False
        if self.ambient != other.ambient:
            return False
        if self.dimension != other.dimension:
            return False
        return (self.equation_matrix()*other.basis_matrix).is_zero_matrix

    def contains(self,element):
        """
        Checks whether the given element belongs to this subspace.

        Arguments:
        element (LieElement or SymmetricElement): the element we are testing.

        Returns
        bool
        """
        if not (isinstance(element, LieElement) or isinstance(element,SymmetricElement)):
            raise TypeError(f"Expected LieElement or SymmetricElement, got {type(element).__name__}")
        if self.ambient != element.parent:
            return False
        return (self.equation_matrix()*self.ambient.coordinates(element)).is_zero_matrix

    
    @staticmethod
    def D(algebra):
        """
        Returns the subspace with basis algebra.diagonals()

        Arguments:
        algebra (GeneralLinear): ambient algebra
        """
        if not isinstance(algebra,GeneralLinear):
            raise TypeError(f"Expected GeneralLinear object, got {type(algebra).__name__}")
        return Subspace.from_basis(algebra,*(algebra.diagonals()))

    @staticmethod
    def D0(algebra):
        """
        Returns the subspace with basis algebra.traceless_diagonals()

        Arguments:
        algebra (GeneralLinear): ambient algebra
        """
        if not isinstance(algebra,GeneralLinear):
            raise TypeError(f"Expected GeneralLinear object, got {type(algebra).__name__}")
        return Subspace.from_basis(algebra,*(algebra.traceless_diagonals()))

    @staticmethod
    def T(algebra):
        """
        Returns the subspace with basis algebra.strictly_upper_triangulars()

        Arguments:
        algebra (GeneralLinear): ambient algebra
        """
        if not isinstance(algebra,GeneralLinear):
            raise TypeError(f"Expected GeneralLinear object, got {type(algebra).__name__}")
        return Subspace.from_basis(algebra,*(algebra.strictly_upper_triangulars()))

    def is_diagonal_stable(self):
        """
        Returns true if this subspace is stable under the action of diagonal matrices in self.ambient.algebra. Only valid if self.ambient is a SymmetricPower
        """
        if not isinstance(self.ambient, SymmetricPower):
            raise TypeError("Diagonal stability is only checked for subspaces of SymmetricPower.")
        return all(
            (self.equation_matrix()*action_matrix(h,self.ambient)*self.basis_matrix).is_zero_matrix for h in self.ambient.algebra.diagonals()
        )

    def basis(self):
        """
        Returns a tuple of LieElements or SymmetricElements which are a basis of this subspace.
        """
        result = tuple()
        for i in range(self.basis_matrix.shape[1]):
            result = result + (self.ambient.from_coordinates(self.basis_matrix.col(i)),)
        return result

    def weight_decomposition(self):
        """
        Returns a dictionary from weights to subspaces of this subspace which are weight spaces of the given weight. Only valid if self.ambient is a SymmetricPower
        """
        if self.dimension == 0:
            return {}
        if not self.is_diagonal_stable():
            raise ValueError("Weight space decomposition only exists for diagonal stable subspaces.")
        all_weights = self.ambient.weight_dictionary().keys()
        decomposition = dict()
        for weight in all_weights:
            support = sorted(self.ambient.weight_element_indices(*weight), reverse=True) # the indices, in the canonical basis, of the monomials in which this weight space is supported. ordered from highest to lowest.
            constrained_matrix = ImmutableMatrix(self.basis_matrix)
            for row_index in support:
                constrained_matrix = constrained_matrix.row_del(row_index)
                # in the result of this operations, for each column, we get the coefficients outside the support of the given weight.
            list_of_coordinates = constrained_matrix.nullspace() # these are the coordinates, in the stored basis of self, of the elements with a given weight.
            nullspace_matrix = ImmutableMatrix([Rational(0)] * self.dimension)
            for column in list_of_coordinates:
                nullspace_matrix = nullspace_matrix.col_insert(list_of_coordinates.index(column)+1,column)
            nullspace_matrix = self.basis_matrix*nullspace_matrix.col_del(0) # the columns in this matrix are coordinates in self.ambient of vectors in self that are weight vectors with the given weight.
            basis_elements = tuple()
            for i in range(nullspace_matrix.shape[1]): # iterate through each column
                basis_elements = basis_elements + (self.ambient.from_coordinates(nullspace_matrix.col(i)),)
            if len(basis_elements) > 0:
                decomposition[weight] = Subspace.from_basis(self.ambient,*basis_elements)
        return decomposition

    def invariants(self,X):
        """
        Returns the subspace of invariant elements under the action of each element of X. Only valid if self.ambient is a SymmetricPower.

        Arguments:
        X: (Subspace) ambient space must be GeneralLinear, equal to self.ambient.algebra

        Returns:
        (Subspace) of self.ambient, contained in this subspace, invariant under every element of X.
        """
        if not isinstance(self.ambient,SymmetricPower):
            raise TypeError("Invariants only defined for subspaces of SymmetricPowers")
        if not isinstance(X,Subspace):
            raise TypeError(f"X must be a Subspace, got {type(X).__name__}")
        if X.ambient != self.ambient.algebra:
            raise ValueError(f"X must be a subspace of {self.ambient.algebra}, got {X.ambient}")
        constraints = ImmutableMatrix(self.equation_matrix())
        for x in X.basis():
            constraints = ImmutableMatrix.vstack(constraints,action_matrix(x,self.ambient))
        return Subspace.from_equations(self.ambient,constraints)

    def evaluation_kernel(self,x):
        """
        Returns the subspace of self consisting of all elements who evaluate to 0 at x, using the trace form to map S(k) to homogeneous polynomials of degree k over g. Only valid if self.ambient is a SymmetricPower

        Arguments:
        x: (LieElement) must have the same parent as self.ambient.algebra

        Returns:
        (Subspace) having the same ambient space as self, contained in this subspace, all of its elements evaluate to 0 at x.
        """
        if not isinstance(self.ambient,SymmetricPower):
            raise TypeError("Evaluation kernels are only defined for subspaces of SymmetricPowers")
        if not isinstance(x,LieElement):
            raise TypeError(f"Expected LieElement, got {type(x).__name__}")
        if x.parent != self.ambient.algebra:
            raise ValueError(f"Cannot evaluate elements of {self.ambient} at elements of {x.parent}")
        if self.dimension == 0:
            return self
        e = ImmutableMatrix([evaluate(v,x) for v in self.basis()]).transpose() # evaluation matrix, a row with each entry being some element of self.basis evaluated at x.
        list_of_basis_coordinates = e.nullspace()
        list_of_basis_vectors = tuple()
        for column in list_of_basis_coordinates:
            list_of_basis_vectors = list_of_basis_vectors + (self.ambient.from_coordinates(self.basis_matrix * column),)
        return Subspace.from_basis(self.ambient,*list_of_basis_vectors)

    def intersection_with_symmetric_sl(self):
        """
        Using the natural embeding of S^k(sl(n)) inside S^k(gl(n)), returns the intersection of self with S^k(sl(n))
        """
        if not isinstance(self.ambient,SymmetricPower):
            raise TypeError("Intersections with sl are only defined for subspaces of SymmetricPowers")
        if self.ambient.degree == 0:
            return self # nothing to compute in degree zero
        # The maths is that this is the kernel of the trace differential operator \sum d/d_{z_i_i}. We write down a matrix representing this operator and compute its nullspace.
        lower_degree_space = SymmetricPower(self.ambient.algebra,self.ambient.degree-1)
        trace_differential_matrix = ImmutableMatrix([0]* lower_degree_space.dimension)
        for x in self.basis():
            trace_image = lower_degree_space.zero()
            for i in range(self.ambient.algebra.n):
                z = lower_degree_space.generators[i*self.ambient.algebra.n+i]
                trace_image.poly = trace_image.poly + x.poly.diff(z)
            trace_differential_matrix = trace_differential_matrix.col_insert(self.basis().index(x)+1,lower_degree_space.coordinates(trace_image))
        trace_differential_matrix = trace_differential_matrix.col_del(0)
        list_of_basis_coordinates = trace_differential_matrix.nullspace()
        list_of_basis_vectors = tuple()
        for column in list_of_basis_coordinates:
            list_of_basis_vectors = list_of_basis_vectors + (self.ambient.from_coordinates(self.basis_matrix * column),)
        return Subspace.from_basis(self.ambient,*list_of_basis_vectors)
