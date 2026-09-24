# GeneralLinear

Each instance of this class represents a Lie algebra gl(n) of n-by-n matrices with rational entries.

## Properties

- n: (Integer) the size of the matrices.
- dimension: (Integer) returns n**2.

## Methods

- zero(): (LieElement) returns the zero matrix in this Lie algebra.
- basis(): (tuple of LieElements) returns a tuple of length self.dimension of LieElements forming the canonical basis of this Lie algebra (i.e. the elementary matrices E(i,j)).
- E(i,j): (LieElement, i: Integer, j: Integer) returns the elementary matrix with entry 1 in the i-th row, j-th column, zero everywhere else. Indices following python convention, and run from 0 to self.n -1.
- coordinates(elem): (sympy.ImmutableMatrix, elem: LieElement) given a LieElement elem whose parant is this Lie algebra, returns a column matrix expressing elem in the canonical basis.
- from_coordinates(col): (LieElement, col: sympy.ImmutableMatrix or tuple of sympy.Rational) given a column matrix or a tuple of sympy.Rational of the right size, returns the LieElement in this given Lie algebra whose coordinates in the canonical basis are the given ones.

# LieElement

Each instance of this class represents a matrix inside gl(n). It is a wrapper for an n-by-n matrix with rational coefficients that remembers its parent Lie algebra, and supports acting on other Lie elements through the Bracket as well as Symmetric elements.

## Properties

- parent: (GeneralLinear) the ambient Lie algebra this LieElement belongs to.
- matrix: (sympy.ImmutableMatrix) the actual underlying matrix associated to this element.

# SymmetricPower

Each instance of this class represents a space S(k) of the k-th symmetric power of gl(n) for some n. 

## Properties

- algebra: (GeneralLinear) the Lie algebra gl(n) used for this symmetric power.
- degree: (Integer) the degree of the homogeneous polynomials in this space.
- dimension: (Integer) the dimension of this space, that is, the binomial coefficient \binom{n^2+k-1}{n^2-1}.
- generators: (tuple of sympy.symbol) a tuple with the symbols 'z_i_j' representing the polynomial generators of our space.

## Methods

- zero(): (SymmetricElement) returns the zero element in this symmetric power.
- basis(): (tuple of SymmetricElement) returns the canonical basis of this symmetric power as a tuple of SymmetricElements, following the monomial ordering according to the total ordering of gl(n) in the design decisions document.
- basis_labels(): (tuple of tuple of ints) returns the exponent labels of the canonical ordered basis.
- coordinates(elem): (sympy.ImmutableMatrix, elem: SymmetricElement) given a SymmetricElement elem in this space, returns the column matrix representing this element in the canonical basis.
- from_coordinates(col): (SymmetricElement, col: sympy.ImmutableMatrix or tuple of sympy.Rational) given a column matrix or a tuple of sympy.Rational of the appropriate size, returns the SymmetricElement in this space with the given coordinates in the canonical basis.

# SymmetricElement

## Properties

- parent: (SymmetricPower) the ambient SymmetricPower this SymmetricElement belongs to.
- poly: (sympy.Poly) the underlying read-only polynomial with rational coefficients representing this element.

# Subspace

This class implements a subspace of either a GeneralLinear object or a SymmetricPower object.

## Properties

- ambient: (GeneralLinear or SymmetricProduct) the ambient space this subspace is contained in.
- dimension: (Integer) the dimension of the subspace
- basis_matrix: (sympy.ImmutableMatrix) an immutable matrix of size self.ambient.dimension-by-self.dimension where each column is the expression of a basis element in the canonical basis of self.ambient.

## Methods

- basis(): (tuple of SymmetricElement or LieElement) returns a tuple of self.dimension SymmetricElements or Lie Elements which are a basis of this subspace.
- equations(): (sympy.ImmutableMatrix) returns the m equations defining this subspace as a matrix A of size m-by-self.ambient.dimension, where this subspace is given as the nullspace of A.
- contains(v): (Boolean, v: SymmetricElement or LieElement) tests whether a given SymmetricElement or LieElement v belongs to this subspace or not.