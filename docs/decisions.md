# Design decisions for glder 0.1

This document records the mathematical conventions and API choices behind
`glder`, together with the procedure for preparing its first GitHub release.
The callable interface is listed in [API.md](./API.md); the worked mathematical
example is in [TUTORIAL.md](./TUTORIAL.md) and
[gl2_walkthrough.py](../examples/gl2_walkthrough.py).

## 1. Package identity and public interface

The Python import package is **`glder`**. The release recipe below also uses
`glder` as the distribution name, version `0.1.0`, and Git tag `v0.1.0`.

The root package exposes the main interface:

```python
from glder import (
    GeneralLinear, LieElement,
    SymmetricPower, SymmetricElement,
    Subspace, act, action_matrix, evaluate,
)
```

The defining modules remain useful for organization and direct imports.
`kronecker_delta` is an elementary helper available from `glder.lie`.
`_validation` and private names are implementation details. A name imported by
a module for its own implementation does not automatically become supported
public API.

The documentation has three complementary roles: the README gives a short
entry point, the tutorial develops the mathematics through runnable examples,
and the API reference records signatures and contracts. Module docstrings
explain each module's purpose and conventions; class and method docstrings
remain the local reference for individual objects.

## 2. Exact rational arithmetic

The coefficient field is $\mathbb Q$. Exact coefficients keep rank, kernel,
membership, equality, and invariant-space computations free of numerical
tolerances. SymPy supplies `Rational`, `ImmutableMatrix`, and `Poly` over `QQ`.

Python integers, `Fraction`, and SymPy rational numbers are accepted scalar
inputs. Floats and booleans are rejected. In particular, Python's `1/2` is a
float; write `Rational(1, 2)` or `Fraction(1, 2)`. Matrix sizes, degrees, indices,
and monomial exponents use Python integers and exclude booleans.

General symbolic coefficient fields and floating-point linear algebra are
outside the scope of 0.1.

## 3. Parents, elements, and validated construction

`GeneralLinear` and `SymmetricPower` describe spaces. `LieElement` and
`SymmetricElement` carry their elements. Each element stores its parent, so
operations can check compatibility before combining data.

Parent equality is structural: matrix size determines equality of general
linear algebras; matrix size and degree determine equality of symmetric
powers. Code should use `==`, not object identity, to express compatibility.

Factory methods validate external data. Direct `LieElement`,
`SymmetricElement`, and `Subspace` constructors are low-level and unchecked.
This keeps internal construction simple, while ordinary users enter through
`E`, `from_coordinates`, `from_terms`, `from_poly`, `from_basis`, or
`from_equations`.

Public attributes expose the underlying exact objects for mathematical work.
They should be treated as read-only, although the Python wrappers do not
enforce complete immutability. Bases, labels, weight dictionaries, and equation
matrices are computed lazily and cached. Mutating returned data or rebinding
parent attributes can make those caches inconsistent.

## 4. Bases, ordering, and coordinates

Matrix indices are zero-based. The canonical basis of $\mathfrak{gl}_n$ lists
$E_{ij}$ by increasing row and then increasing column. Position $in+j$
corresponds to $E_{ij}$.

The symmetric generators `z_i_j` use the same order. A symmetric-power basis
consists of ordinary monomials in descending lexicographic order of their
exponent tuples. There are no divided-power or factorial factors.

Coordinates are columns. A subspace basis matrix has basis vectors as columns;
an equation matrix has linear equations as rows. An action matrix has the
image of basis vector $j$ in column $j$. These conventions give

$$
[x\cdot v]=M_x[v],\qquad
Y=\operatorname{im}(B)=\ker(A),\qquad AB=0.
$$

Subspace construction may normalize a supplied spanning list. Mathematical
comparisons should use subspace equality or membership, rather than comparing
the particular basis vectors returned by two construction paths.

## 5. The adjoint action and polynomial evaluation

The Lie bracket is $[x,y]=xy-yx$. Its extension to the symmetric algebra is
the derivation determined by

$$
E_{ab}\cdot z_{ij}=\delta_{bi}z_{aj}-\delta_{aj}z_{ib}.
$$

It preserves homogeneous degree. The identity matrix acts by zero, so the
adjoint action of $\mathfrak{gl}_n$ factors through its traceless part.

The trace form identifies a Lie element $y$ with the linear function
$x\mapsto\operatorname{tr}(yx)$. Consequently,

$$
z_{ij}(x)=x_{ji}.
$$

`evaluate(v, x)` therefore substitutes transposed matrix entries. The transpose
is part of the convention, not a correction to be applied by the caller.
For instance, the generator representing $E_{01}$ evaluates to one at
$E_{10}$.

`act(x, v)` and `evaluate(v, x)` have different argument orders because they
express different operations: an element acting on a vector, and a polynomial
being evaluated at a point.

## 6. Fixed homogeneous degrees

A `SymmetricPower` represents one $S^k(\mathfrak{gl}_n)$, with dimension

$$
\binom{n^2+k-1}{k}.
$$

Addition and subtraction stay within one degree. `SymmetricElement`
multiplication supports scalars; multiplication of two symmetric elements is
performed explicitly on their `Poly` objects and wrapped in the appropriate
new degree. This keeps the parent of the result explicit.

Degree zero is included: it is the one-dimensional space of constants, with
zero Lie action and zero weight. Its intersection with $S^0(\mathfrak{sl}_n)$
is unchanged.

The products in $S(\mathfrak g)$ are commutative. They are distinct from both
matrix multiplication and products in the universal enveloping algebra.
In particular, the tutorial's $h^2+4ef$ is a symmetric-algebra invariant.

## 7. Weights and highest-weight vectors

Weights are eigenvalue tuples for the full diagonal basis
$(E_{00},\ldots,E_{n-1,n-1})$. A monomial has weight

$$
\mu_i=\sum_j\alpha_{ij}-\sum_j\alpha_{ji}.
$$

Every occurring weight has coordinate sum zero. Its entries may be negative;
they are not monomial exponents. Eigenvalues for the traceless diagonal basis
are the successive differences $\mu_i-\mu_{i+1}$.

`Subspace.D(g)` is the diagonal subalgebra, `Subspace.D0(g)` its traceless
part, and **`Subspace.T(g)` is the strictly upper triangular subalgebra**.
The letter `T` does not denote a diagonal torus in this API. Strictly upper
triangular matrices are the positive root operators for this convention.

`Y.invariants(X)` imposes annihilation by every element of `X`. It does not
merely test that the action preserves `Y`. The method makes sense for any
linear subspace `X` of the Lie algebra, and `Y` need not be stable under `X`.

For a weight space $V_\lambda$, the quantities
$\dim V_\lambda$ and $\dim V_\lambda^T$ differ. The latter counts primitive
vectors; in the finite-dimensional representations considered here, for
dominant $\lambda$ it gives the multiplicity of the irreducible representation
with highest weight $\lambda$.

The degree-eight example makes this distinction concrete:

| Representation | $\dim V_{(2,-2)}$ | $\dim V_{(2,-2)}^T$ |
| --- | ---: | ---: |
| $S^8(\mathfrak{gl}_2)$ | 16 | 4 |
| $S^8(\mathfrak{sl}_2)$ | 4 | 1 |

On $h=E_{00}-E_{11}$, this weight evaluates to $4$. It corresponds to the
five-dimensional irreducible $\mathfrak{sl}_2$ representation.

## 8. Subspaces retain their ambient coordinates

Subspace operations return subspaces of the same ambient parent. Invariant
spaces and evaluation kernels can therefore be compared, evaluated, or acted
on using the original coordinates.

An arbitrary subspace need not split as a direct sum of its intersections
with weight spaces. `weight_decomposition()` requires diagonal stability and
reports a failure when that hypothesis is not met.

The inclusion $\mathfrak{sl}_n\subset\mathfrak{gl}_n$ induces an inclusion of
symmetric powers. In characteristic zero, its image is

$$
S^k(\mathfrak{sl}_n)
=\ker\left(\sum_i\frac{\partial}{\partial z_{ii}}\right)
\subset S^k(\mathfrak{gl}_n).
$$

`intersection_with_symmetric_sl()` computes intersection with this kernel.
It preserves the ambient `SymmetricPower` and does not replace it with a new
parent. Restricting polynomial functions to traceless matrices instead gives
a quotient by the trace polynomial; that is a different construction.

## 9. Scope and performance

Version 0.1 favors transparent exact computations and inspectable SymPy
objects. Bases and action matrices are explicit; their sizes grow quickly
with both $n$ and $k$. Caching avoids some repeated construction, but does not
remove the combinatorial growth in dimension.