# glder API reference

This document describes the public interface for version 0.1. Read
[TUTORIAL.md](./TUTORIAL.md) for the mathematical development and run
[gl2_walkthrough.py](../examples/gl2_walkthrough.py) for a complete example.

The Python examples below can be run in order in one session.

## Imports and modules

```python
from fractions import Fraction
from sympy import ImmutableMatrix, Poly, QQ, Rational
from glder import (
    GeneralLinear, LieElement,
    SymmetricPower, SymmetricElement,
    Subspace, act, action_matrix, evaluate,
)
```

| Module | Public objects |
| --- | --- |
| `glder.lie` | `GeneralLinear`, `LieElement`, `kronecker_delta` |
| `glder.symmetric` | `SymmetricPower`, `SymmetricElement` |
| `glder.actions` | `act`, `action_matrix`, `evaluate` |
| `glder.subspaces` | `Subspace` |

The eight principal objects in the import example are re-exported by `glder`.
The elementary helper `kronecker_delta` is documented below under its defining
module. `glder._validation` and names beginning with an underscore are internal.
Imports used to implement a module are not additional package API.

## Shared conventions

All vector spaces and coefficients are over $\mathbb Q$. Accepted scalar
inputs are Python `int`, `fractions.Fraction`, and SymPy `Rational` (including
`Integer`). Booleans, floats, strings, and symbolic coefficients are not rational
scalar inputs. Matrix sizes, degrees, indices, and monomial exponents use Python
integers, excluding booleans.

Indices are zero-based. In $\mathfrak{gl}_n$, the ordered basis is
$(E_{00},E_{01},\ldots,E_{0,n-1},E_{10},\ldots,E_{n-1,n-1})$.
Coordinates are `ImmutableMatrix` columns. A method accepting `*args` takes
separate positional arguments: use `method(*values)` to unpack a tuple.

Parent equality is structural: two `GeneralLinear` instances are equal when
their sizes agree; two `SymmetricPower` instances are equal when their algebras
and degrees agree. They need not be the same Python object.

Use factory methods to construct elements and subspaces. The element and
subspace constructors below are low-level and do not validate their data.
Treat public attributes, returned basis elements, and cached dictionaries as
read-only; rebinding them can invalidate cached computations.

## GeneralLinear

**`GeneralLinear(val)`** represents $\mathfrak{gl}_{\mathrm{val}}(\mathbb Q)$.
`val` must be a positive Python integer. A wrong type raises `TypeError`; a
nonpositive value raises `ValueError`.

| Attribute | Type | Meaning |
| --- | --- | --- |
| `n` | `int` | Matrix size. |
| `dimension` | `int` | $n^2$. |

| Method | Return value and behavior |
| --- | --- |
| `zero()` | The zero `LieElement`. |
| `E(i, j)` | The elementary matrix $E_{ij}$ as a `LieElement`; $0\leq i,j<n$. |
| `basis()` | Cached tuple of the $n^2$ elementary matrices, in row order. |
| `coordinates(elem)` | The $n^2\times1$ coordinate column of a compatible `LieElement`. |
| `from_coordinates(*args)` | A `LieElement` from one $n^2\times1$ `ImmutableMatrix`, or $n^2$ separate exact scalars. |
| `diagonals()` | Tuple $(E_{00},\ldots,E_{n-1,n-1})$. |
| `traceless_diagonals()` | Tuple $(E_{ii}-E_{i+1,i+1})_{i=0}^{n-2}$; empty for $n=1$. |
| `strictly_upper_triangulars()` | Tuple of $E_{ij}$ with $i<j$, in row order. |

`E` raises `TypeError` for a noninteger index and `IndexError` for an index
outside the range. `coordinates` checks the element type and parent.
`from_coordinates` raises `IndexError` for an incorrectly shaped coordinate
matrix and `TypeError` for an unsupported argument arrangement or scalar.
A list is not a coordinate matrix; unpack it or use `ImmutableMatrix`.

Equality uses `==`. `str(g)` has the form `GeneralLinear(n)`.

```python
g = GeneralLinear(2)
e = g.E(0, 1)
f = g.E(1, 0)
h = g.E(0, 0) - g.E(1, 1)

assert g.dimension == 4
assert g.coordinates(e) == ImmutableMatrix([0, 1, 0, 0])
assert g.from_coordinates(1, 0, 0, -1) == h
assert g.from_coordinates(g.coordinates(h)) == h
assert g.diagonals() == (g.E(0, 0), g.E(1, 1))
assert g.traceless_diagonals() == (h,)
assert g.strictly_upper_triangulars() == (e,)
```

### LieElement

**`LieElement(parent, matrix)`** stores a `GeneralLinear` parent and an
$n\times n$ rational `ImmutableMatrix`. The constructor is unchecked;
prefer `g.E`, `g.zero`, or `g.from_coordinates`.

| Attribute | Type | Meaning |
| --- | --- | --- |
| `parent` | `GeneralLinear` | Ambient Lie algebra. |
| `matrix` | `ImmutableMatrix` | Underlying matrix. |

For compatible Lie elements `x` and `y`, the public operations are `x == y`,
`str(x)`, `x + y`, `x - y`, `-x`, `a*x`, `x*a`, `x*y`, and
**`x.bracket(y)`**. Here `a` is an exact scalar, `x*y` is ordinary matrix
multiplication, and `x.bracket(y)` is $xy-yx$. These operations return new
elements; equality returns a boolean and `str` returns a string.

Incompatible Lie parents raise `ValueError` in the arithmetic operations.
Unsupported operand types use Python's normal operator handling. There is no
division operator: multiply by a rational reciprocal. When summing elements,
use `sum(elements, g.zero())` to supply the correct additive identity.

```python
assert e.bracket(f) == h
assert h.bracket(e) == 2*e
assert h.bracket(f) == -2*f
assert e*f == g.E(0, 0)
assert Fraction(1, 2)*h == Rational(1, 2)*h
assert sum(g.diagonals(), g.zero()) == g.from_coordinates(1, 0, 0, 1)
```

### kronecker_delta

**`glder.lie.kronecker_delta(a, b)`** returns SymPy `Rational(1)` when
`a == b`, and `Rational(0)` otherwise. It does not perform index validation.

```python
from glder.lie import kronecker_delta

assert kronecker_delta(0, 0) == 1
assert kronecker_delta(0, 1) == 0
```

## SymmetricPower

**`SymmetricPower(algebra, degree)`** represents $S^k(\mathfrak{gl}_n)$,
where `algebra` is a `GeneralLinear` instance and `degree = k` is a nonnegative
Python integer. Wrong argument types raise `TypeError`; a negative degree
raises `ValueError`.

| Attribute | Type | Meaning |
| --- | --- | --- |
| `algebra` | `GeneralLinear` | Underlying Lie algebra. |
| `degree` | `int` | Homogeneous degree $k$. |
| `dimension` | `int` | $\binom{n^2+k-1}{k}$. |
| `generators` | tuple of SymPy symbols | Symbols `z_i_j` representing $E_{ij}$, in row order. |

The basis consists of ordinary commutative monomials, with no factorial
normalization, in descending lexicographic order of their exponent tuples.
Every label has length $n^2$, nonnegative integer entries, and total $k$.

| Method | Return value and behavior |
| --- | --- |
| `basis()` | Cached tuple of `SymmetricElement` basis monomials. |
| `basis_labels()` | Cached tuple of their exponent tuples in the same order. |
| `zero()` | The zero `SymmetricElement` in this degree. |
| `monomial(*args)` | Basis monomial with the $n^2$ specified exponents. |
| `from_terms(mapping)` | Element from a `dict` mapping exponent tuples to exact coefficients; `{}` gives zero. |
| `from_poly(poly)` | Element from a compatible homogeneous SymPy `Poly` over `ZZ` or `QQ`, converted to `QQ`. |
| `coordinates(elem)` | Coordinate column of shape `(dimension, 1)`. |
| `from_coordinates(*args)` | Element from one such `ImmutableMatrix`, or `dimension` separate exact scalars. |

`from_poly` accepts the zero polynomial in every degree. A nonzero polynomial
must be homogeneous of the parent's degree, with compatible variables and
rational coefficients. Pass a `Poly`, not an unwrapped SymPy expression.
`monomial` validates the number, type, and values of its exponents.
Coordinate conversion checks the element's type and parent; an incorrectly
shaped coordinate matrix passed to `from_coordinates` raises `IndexError`.

```python
S2 = SymmetricPower(g, 2)
u, b, c, d = S2.generators
assert S2.dimension == 10
assert S2.basis_labels()[0] == (2, 0, 0, 0)

bc = S2.monomial(0, 1, 1, 0)
q = S2.from_terms({
    (2, 0, 0, 0): 1,
    (1, 0, 0, 1): -2,
    (0, 0, 0, 2): 1,
    (0, 1, 1, 0): 4,
})
assert q == S2.from_poly(Poly((u-d)**2 + 4*b*c, *S2.generators, domain=QQ))
assert S2.from_coordinates(S2.coordinates(q)) == q
assert S2.from_terms({}) == S2.zero()
```

### Weight methods

A weight is a tuple $(\mu_0,\ldots,\mu_{n-1})$ of eigenvalues for
$E_{00},\ldots,E_{n-1,n-1}$. Its entries may be negative. For a monomial
$z^\alpha$,

$$
\mu_i=\sum_j\alpha_{ij}-\sum_j\alpha_{ji},
\qquad \sum_i\mu_i=0.
$$

| Method | Return value and behavior |
| --- | --- |
| `weight_from_exponent_tuple(*args)` | Weight tuple of a validated degree-$k$ exponent tuple of length $n^2$. |
| `weight_elements(*args)` | Tuple of basis monomials of the supplied weight; empty if absent. Supply $n$ separate weight entries. |
| `weight_element_indices(*args)` | Tuple of their zero-based indices in `basis()`. |
| `weight_dictionary()` | Cached `dict` mapping weights to tuples of **exponent labels**, not elements. |

For the traceless diagonal $E_{ii}-E_{i+1,i+1}$, the eigenvalue is
$\mu_i-\mu_{i+1}$. Weight entries and monomial exponents are different data.

```python
assert S2.weight_from_exponent_tuple(0, 2, 0, 0) == (2, -2)
assert S2.weight_elements(2, -2) == (S2.monomial(0, 2, 0, 0),)
assert S2.weight_elements(9, -9) == ()
indices = S2.weight_element_indices(2, -2)
assert tuple(S2.basis()[i] for i in indices) == S2.weight_elements(2, -2)
assert S2.weight_dictionary()[(2, -2)] == ((0, 2, 0, 0),)
```

### SymmetricElement

**`SymmetricElement(parent, poly)`** stores a `SymmetricPower` parent and its
rational homogeneous `Poly`. The constructor is unchecked; prefer the parent
factories above.

| Attribute | Type | Meaning |
| --- | --- | --- |
| `parent` | `SymmetricPower` | Ambient homogeneous component. |
| `poly` | `sympy.Poly` | Polynomial over `QQ` in the parent's generators. |

The public operations are `v == w`, `str(v)`, `v + w`, `v - w`, `-v`,
`a*v`, and `v*a`. Addition and subtraction require equal parents; multiplication
accepts exact scalars. Use `v.poly.as_expr()` for polynomial display.

Multiplication of two `SymmetricElement` objects and exponentiation of such
objects are not provided. Multiply their underlying polynomials and construct
an element in the resulting degree explicitly:

```python
S4 = SymmetricPower(g, 4)
q_squared = S4.from_poly(q.poly * q.poly)
assert q_squared.parent.degree == 4
assert q + (-q) == S2.zero()
assert sum((q, 2*q), S2.zero()) == 3*q
```

## Actions and evaluation

### act(x, v)

Takes a `LieElement` and a `SymmetricElement`, and returns a
`SymmetricElement` in `v.parent`. Their Lie algebras must agree. The action is
determined by

$$
E_{ab}\cdot z_{ij}=\delta_{bi}z_{aj}-\delta_{aj}z_{ib}
$$

and the Leibniz rule. It preserves degree and annihilates constants.

### action_matrix(x, V)

Takes a `LieElement` and a compatible `SymmetricPower`. Returns an
`ImmutableMatrix` of shape `(V.dimension, V.dimension)` whose column $j$ is
`V.coordinates(act(x, V.basis()[j]))`. Consequently,

$$
[x\cdot v]_V=M_x[v]_V.
$$

### evaluate(v, x)

Takes a `SymmetricElement` followed by a compatible `LieElement`. Returns an
exact SymPy rational number. The argument order differs from `act`.

The trace pairing identifies each generator with the linear function

$$
z_{ij}(x)=\operatorname{tr}(E_{ij}x)=x_{ji}.
$$

Thus evaluation substitutes **transposed matrix entries** into the polynomial.
The three functions reject wrong object types with `TypeError` and mismatched
Lie algebras with `ValueError`.

```python
S1 = SymmetricPower(g, 1)
b1 = S1.monomial(0, 1, 0, 0)
assert evaluate(b1, f) == 1
assert evaluate(b1, e) == 0

assert act(e, q) == S2.zero()
M = action_matrix(e, S2)
assert M.shape == (10, 10)
assert M*S2.coordinates(bc) == S2.coordinates(act(e, bc))

x = g.from_coordinates(1, 2, 3, -1)
assert evaluate(q, x) == 28
```

## Subspace

A `Subspace` retains its ambient `GeneralLinear` or `SymmetricPower`.
Write $N$ for the ambient dimension and $d$ for the subspace dimension.

**`Subspace(ambient, dimension, basis_matrix)`** is an unchecked low-level
constructor. Its matrix must have shape $(N,d)$ with independent rational
columns. Prefer the factories below.

| Attribute | Type | Meaning |
| --- | --- | --- |
| `ambient` | `GeneralLinear` or `SymmetricPower` | Original coordinate space. |
| `dimension` | `int` | Subspace dimension $d$. |
| `basis_matrix` | `ImmutableMatrix` | Independent basis columns, expressed in ambient coordinates. |

### Construction

Call these factories on `Subspace` itself:

| Factory | Result |
| --- | --- |
| `Subspace.from_basis(ambient, *args)` | Span of the supplied compatible elements; dependence and zero vectors are allowed. With no elements, gives the zero subspace. |
| `Subspace.from_equations(ambient, equation_matrix)` | Kernel of a rational `ImmutableMatrix` with $N$ columns. Rows may be dependent; zero rows of equations give the whole ambient space. |
| `Subspace.whole(ambient)` | The ambient space, represented as its own subspace. |
| `Subspace.trivial(ambient)` | The zero subspace. |
| `Subspace.D(algebra)` | Diagonal matrices; dimension $n$. |
| `Subspace.D0(algebra)` | Traceless diagonal matrices; dimension $n-1$. |
| `Subspace.T(algebra)` | Strictly upper triangular matrices; dimension $n(n-1)/2$. |

The last three factories take a `GeneralLinear` instance and return subspaces
of that Lie algebra. `from_basis` may replace the supplied spanning list with
a normalized basis. Do not rely on the input basis being returned unchanged.
For `from_equations`, a non-`ImmutableMatrix` argument raises `TypeError`, and
an incorrect number of columns raises `ValueError`.

```python
D = Subspace.D(g)
D0 = Subspace.D0(g)
T = Subspace.T(g)
assert (D.dimension, D0.dimension, T.dimension) == (2, 1, 1)
assert D0 == Subspace.from_basis(g, h, 2*h, g.zero())

trace_zero = Subspace.from_equations(g, ImmutableMatrix([[1, 0, 0, 1]]))
assert trace_zero.dimension == 3
assert trace_zero == Subspace.from_basis(g, e, f, h)
assert Subspace.trivial(g) == Subspace.from_basis(g)
```

### Basis, equations, equality, and membership

**`Y.basis()`** returns a tuple of $d$ elements of `Y.ambient`.

**`Y.equation_matrix()`** returns a cached `ImmutableMatrix` $A$ of shape
$(N-d,N)$ with independent rows such that

$$
Y=\operatorname{im}(B)=\ker(A),\qquad B=Y.\texttt{basis\_matrix}.
$$

In particular, $AB=0$. Whole and trivial spaces use the corresponding empty
matrix dimensions; no artificial basis vector represents a zero-dimensional
space.

**`Y.contains(element)`** checks membership. A `LieElement` or
`SymmetricElement` with an unequal parent returns `False`; an unsupported
object type raises `TypeError`.

**`Y == Z`** compares the ambient spaces and the spanned subspaces, independently
of which bases were used. Equality with another object type is `False`.

```python
assert trace_zero.contains(e)
assert not trace_zero.contains(g.E(0, 0))
assert trace_zero.equation_matrix()*trace_zero.basis_matrix == ImmutableMatrix.zeros(1, 3)
assert len(trace_zero.basis()) == trace_zero.dimension
```

### is_diagonal_stable()

**`Y.is_diagonal_stable()`** returns whether every diagonal matrix sends `Y`
into itself. Requires a `SymmetricPower` ambient; otherwise raises `TypeError`.
For an equation matrix $A$ and basis matrix $B$, this is the condition
$A M_h B=0$ for every $h$ in the diagonal basis.

### weight_decomposition()

**`Y.weight_decomposition()`** returns a dictionary mapping occurring weights
to their nonzero weight subspaces, all retaining `Y.ambient`. Their direct sum
is `Y`. Requires a symmetric-power ambient and diagonal stability; an incorrect
ambient raises `TypeError`, and failure of stability raises `ValueError`.

This method splits a given stable subspace. The methods on `SymmetricPower`
instead enumerate monomials of a weight in the entire ambient space.

```python
Y = Subspace.whole(S2)
assert Y.is_diagonal_stable()
pieces = Y.weight_decomposition()
assert sum(piece.dimension for piece in pieces.values()) == Y.dimension
assert pieces[(2, -2)] == Subspace.from_basis(S2, *S2.weight_elements(2, -2))
```

### invariants(X)

**`Y.invariants(X)`** returns

$$
Y^X=\{v\in Y:x\cdot v=0\text{ for every }x\in X\}.
$$

`Y.ambient` must be a `SymmetricPower`; `X` must be a `Subspace` of its
underlying Lie algebra. A nonsymmetric ambient for `Y` or a non-`Subspace`
argument raises `TypeError`; an incompatible `X.ambient` raises `ValueError`.
The result retains
`Y.ambient`. Neither stability of `Y` nor the condition that `X` be a Lie
subalgebra is required. This is a common kernel, not a test of stability.

Equivalently, with a basis matrix $B$ for `Y`, the desired vectors are $Bc$
where $c$ lies in the common kernel of the matrices $M_xB$.
Taking `X = Subspace.T(g)` extracts vectors killed by positive root operators;
inside a fixed weight space, these are highest-weight vectors.

```python
all_invariants = Y.invariants(Subspace.whole(g))
assert all_invariants.dimension == 2
assert all_invariants.contains(q)

S8 = SymmetricPower(g, 8)
W = Subspace.from_basis(S8, *S8.weight_elements(2, -2))
H = W.invariants(T)
assert (W.dimension, H.dimension) == (16, 4)
```

### evaluation_kernel(x)

**`Y.evaluation_kernel(x)`** returns $\{v\in Y:v(x)=0\}$, using the trace
pairing in `evaluate`. Requires a symmetric-power ambient and a compatible
`LieElement`. Incorrect types raise `TypeError`; an incompatible Lie algebra
raises `ValueError`. The result retains the ambient space. Its dimension is
$d$ or $d-1$, since evaluation is one linear functional on a fixed degree.

```python
K = Y.evaluation_kernel(x)
assert K.dimension == 9
assert all(evaluate(v, x) == 0 for v in K.basis())
```

### intersection_with_symmetric_sl()

**`Y.intersection_with_symmetric_sl()`** returns
$Y\cap S^k(\mathfrak{sl}_n)$ inside the original
$S^k(\mathfrak{gl}_n)$ ambient. It requires a `SymmetricPower` ambient,
raising `TypeError` otherwise. No stability hypothesis is needed.

The embedded symmetric power is the kernel of trace contraction:

$$
\partial_{\mathrm{tr}}=\sum_{i=0}^{n-1}\frac{\partial}{\partial z_{ii}}.
$$

In degree zero the intersection is `Y` itself. For $n=1$ and positive degree,
it is zero. This operation takes an intersection with an embedded subspace;
it does not construct a quotient or simply impose a polynomial equation
$\sum_i z_{ii}=0$.

```python
sl2_square = Y.intersection_with_symmetric_sl()
assert sl2_square.dimension == 6
assert sl2_square.ambient == S2
assert sl2_square.invariants(Subspace.whole(g)).dimension == 1

W_sl = W.intersection_with_symmetric_sl()
assert (W_sl.dimension, W_sl.invariants(T).dimension) == (4, 1)
```

## Degree zero and scope

`SymmetricPower(g, 0)` is the one-dimensional space of constants. Its basis
label is the all-zero exponent tuple, its weight is the zero tuple, the Lie
action is zero, and evaluation returns the constant.

The API provides fixed-degree spaces and the operations listed here. General
tensor products, enveloping-algebra computations, arbitrary coefficient
fields, numerical approximations, and a generic graded multiplication API are
outside version 0.1. For design rationale and release conventions, see
[DECISIONS.md](./DECISIONS.md).
