# A worked tutorial for General Linear Derivations

This tutorial introduces the public API through exact computations with $\mathfrak{gl}_2(\mathbb Q)$ and its symmetric powers. The complete executable companion is [gl2_walkthrough.py](../examples/gl2_walkthrough.py). It contains assertions for the calculations below and additional short checks covering the public methods.

We will compute the quadratic invariant $$Q=(u-d)^2+4bc\in S^2(\mathfrak{sl}_2)$$ and then calculate the multiplicity of the irreducible representation of highest weight $(2,-2)$ in degree eight. There are two different dimensions to distinguish:

| Ambient representation | Ambient dimension | Weight-space dimension | Dimension annihilated by $T$ |
| --- | ---: | ---: | ---: |
| $S^8(\mathfrak{gl}_2)$ | 165 | 16 | 4 |
| $S^8(\mathfrak{sl}_2)$ | 45 | 4 | 1 |

Here $T=\mathbb Q E_{01}$ is the strictly upper triangular **Lie subalgebra**. Vectors of weight $(2,-2)$ annihilated by $T$ are highest-weight vectors. Their dimension counts copies of the corresponding irreducible representation; the dimension of the full weight space counts all vectors of that weight.

## 1. Install and run

Place this file and `gl2_walkthrough.py` together, for example at the repository root. From the project root, with your development environment active, run:

```bash
python -m pip install -e .
python gl2_walkthrough.py
```

These examples use the names established in earlier sections. The companion script organizes them into `lie_algebra_example`, `quadratic_example`, and `highest_weight_example` so you can read or run the complete sequence.

For expanded computed bases, or a different degree, use:

```bash
python gl2_walkthrough.py --degree 8 --show-basis
python gl2_walkthrough.py --degree 7
```

The walkthrough accepts degrees at least two. Its default is eight. Increasing the degree increases the size of the ambient action matrices; degree eight is already a useful demonstration of the distinction between weights and highest weights. The script has no dependency beyond the package and SymPy, and writes no files.

## 2. Matrices, parents, and exact coefficients

```python
g = GeneralLinear(2)
e00, e, f, e11 = g.basis()
h = e00 - e11
identity = e00 + e11
```

`g` is a parent describing a vector space and its conventions. An element has
attributes `parent` and `matrix`. For example, `e.parent == g` and `e.matrix`
is the immutable SymPy matrix with its only nonzero entry at position `(0, 1)`.
The public parent attributes are `g.n == 2` and `g.dimension == 4`.

All matrix indices are zero-based. The matrix-unit basis is in row-major order:

$$
E_{00},\quad E_{01},\quad E_{10},\quad E_{11}.
$$

```python
assert e == g.E(0, 1)
assert e.bracket(f) == h
assert h.bracket(e) == 2 * e
assert h.bracket(f) == -2 * f
assert e * f == e00
assert e + g.zero() == e
assert e - e == g.zero()
assert -e + e == g.zero()
assert Fraction(1, 2) * h == h * Rational(1, 2)
```

The product `e * f` is matrix multiplication. The Lie bracket is explicitly
`e.bracket(f)`, which computes $ef-fe$.

Scalar arguments accept Python integers, `fractions.Fraction`, and SymPy
`Rational` values, including SymPy integers. Use `Rational(1, 2)` or
`Fraction(1, 2)` for a half: Python's expression `1/2` produces a float, which
the scalar validator rejects. Booleans and strings are also rejected.
Indices, degrees, and exponent entries use Python integers, excluding booleans.

Equality of parents means equal mathematical parameters. In particular,
`GeneralLinear(2) == GeneralLinear(2)` even for separate instances. Element
compatibility uses this equality rather than object identity.

The helper `kronecker_delta` is available from its defining module:

```python
from gl_derivations.lie import kronecker_delta
assert kronecker_delta(0, 0) == Rational(1)
assert kronecker_delta(0, 1) == Rational(0)
```

It is not needed for ordinary user calculations. `_validation` is an internal
module and is not part of the public API used in this tutorial.

## 3. Coordinates and matrix subspaces

```python
x = e00 + 2 * e + 3 * f - e11
column = g.coordinates(x)
assert column == ImmutableMatrix([1, 2, 3, -1])
assert g.from_coordinates(column) == x
assert g.from_coordinates(1, 2, 3, -1) == x
```

Coordinates are columns in the parent's ordered basis. `from_coordinates`
accepts either one immutable column matrix or separate scalar arguments. To
supply a tuple, unpack it: `g.from_coordinates(*values)`.

Standard subalgebras are supplied both as tuples of basis elements and as
`Subspace` objects:

| Basis tuple from `g` | Corresponding `Subspace` | Dimension for $n=2$ |
| --- | --- | ---: |
| `g.diagonals()` | `Subspace.D(g)` | 2 |
| `g.traceless_diagonals()` | `Subspace.D0(g)` | 1 |
| `g.strictly_upper_triangulars()` | `Subspace.T(g)` | 1 |

```python
D = Subspace.D(g)
D0 = Subspace.D0(g)
T = Subspace.T(g)
sl2 = Subspace.from_basis(g, e, f, h)
assert sl2.dimension == 3
assert sl2.contains(h)
assert not sl2.contains(identity)
```

The spanning arguments to `from_basis` need not be independent. They may
include repetitions or zero vectors. Pass them separately, or unpack a
sequence with `*`. The stored basis is normalized by row reduction and may
differ from the sequence you supplied.

The same space can be specified by equations. In our coordinate order,
trace zero means the first and fourth coordinates sum to zero:

```python
sl2_by_equations = Subspace.from_equations(
    g, ImmutableMatrix([[1, 0, 0, 1]])
)
assert sl2_by_equations == sl2
assert Subspace.from_equations(g, sl2.equation_matrix()) == sl2
assert Subspace.from_basis(g, *sl2.basis()) == sl2
```

For a subspace $Y$ of dimension $r$ in an ambient space of dimension $N$,
`Y.basis_matrix` is an $N\times r$ matrix $B$. Its columns are ambient
coordinates. `Y.equation_matrix()` is an $(N-r)\times N$ matrix $A$ satisfying

$$
Y=\mathrm{im}B=\ker A,\qquad AB=0.
$$

If $c$ is a column of coordinates in the stored basis of $Y$, then $Bc$ is its
ambient coordinate column. This is the conversion used by the kernel methods.

```python
assert sl2.ambient == g
assert sl2.basis_matrix.shape == (4, 3)
assert sl2.equation_matrix().shape == (1, 4)
assert (sl2.equation_matrix() * sl2.basis_matrix).is_zero_matrix
assert Subspace.whole(g).dimension == 4
assert Subspace.trivial(g).dimension == 0
assert Subspace.from_basis(g) == Subspace.trivial(g)
```

An equation matrix can have redundant rows. A matrix with shape `(0, N)`
imposes no equations and defines the whole ambient space. Use explicit shapes
when constructing empty matrices. `Subspace.whole` and `Subspace.trivial`
handle these cases for you.

## 4. Symmetric elements and their constructors

```python
V2 = SymmetricPower(g, 2)
u, b, c, d = V2.generators
assert V2.algebra == g
assert V2.degree == 2
assert V2.dimension == 10
```

The symbols represent the vectors $E_{00},E_{01},E_{10},E_{11}$, respectively,
inside the symmetric algebra. They are printed as `z_0_0`, `z_0_1`, `z_1_0`,
and `z_1_1`. The short names in this tutorial are local Python aliases.

The ordinary monomial basis of $S^2(\mathfrak{gl}_2)$ is

$$
u^2,\ ub,\ uc,\ ud,\ b^2,\ bc,\ bd,\ c^2,\ cd,\ d^2.
$$

This is descending lexicographic order on exponent tuples in the generator
order `(u, b, c, d)`. There are no factorial or divided-power normalizations.
`V2.basis()` returns symmetric elements; `V2.basis_labels()` returns their
exponent tuples.

```python
alpha = (1, 0, 1, 0)
uc = V2.monomial(*alpha)
assert uc == V2.from_terms({alpha: 1})
assert uc == V2.from_poly(Poly(u * c, *V2.generators, domain=QQ))
assert uc.parent == V2
assert uc.poly.as_expr() == u * c
assert V2.from_terms({}) == V2.zero()
```

An exponent tuple has $n^2$ nonnegative entries summing to the parent's degree.
`from_terms` accepts a dictionary from such tuples to exact coefficients.
`from_poly` accepts a SymPy `Poly` over `ZZ` or `QQ`, checks homogeneity and
degree, and stores the result over `QQ` using the parent's generators. Supply
the generators explicitly even for constants and zero polynomials.

```python
coords = V2.coordinates(uc)
assert coords.shape == (10, 1)
assert V2.from_coordinates(coords) == uc
assert V2.from_coordinates(*tuple(coords)) == uc
assert uc + V2.zero() == uc
assert uc - uc == V2.zero()
assert -uc + uc == V2.zero()
assert uc * Rational(1, 2) == Rational(1, 2) * uc
```

`SymmetricElement` supports scalar multiplication. To multiply polynomials or
take powers for a new symmetric degree, operate on their SymPy expressions
and pass the result to the appropriate parent's `from_poly`. For example,
degree eight will use `Poly(b**2 * Q_expression**3, ..., domain=QQ)`.

Treat cached basis elements and weight dictionaries as read-only. In
particular, do not reassign `.poly` on an element obtained from `basis()` or
`monomial()`; those elements may be shared by later computations.

## 5. The adjoint action and its matrix

On matrix units, the action is the commutator. On the symmetric algebra it is
extended by the Leibniz rule. For $e=E_{01}$,

$$
e\cdot u=-b,\qquad e\cdot b=0,\qquad
e\cdot c=u-d,\qquad e\cdot d=b.
$$

Consequently $e\cdot(uc)=u^2-ud-bc$:

```python
e = g.E(0, 1)
expected = V2.from_poly(Poly(u**2 - u*d - b*c, *V2.generators, domain=QQ))
assert act(e, uc) == expected
M = action_matrix(e, V2)
assert M.shape == (10, 10)
assert M * V2.coordinates(uc) == V2.coordinates(act(e, uc))
```

Column $j$ of `action_matrix(e, V2)` records the action on basis element $j$.
The matrix multiplies coordinate columns on the left. `act` preserves the
symmetric degree and returns an element with the same parent as its input.

## 6. Diagonal weights and decomposition

The weight of $u^a b^p c^q d^r$ is $(p-q,q-p)$. More generally, the $i$-th
weight entry of a monomial is its exponent row sum minus its exponent column
sum. These entries are the eigenvalues of $E_{00},E_{11}$ in the present case.

```python
assert V2.weight_from_exponent_tuple(1, 0, 1, 0) == (-1, 1)
monomials = V2.weight_elements(1, -1)
indices = V2.weight_element_indices(1, -1)
assert monomials == tuple(V2.basis()[i] for i in indices)
assert len(monomials) == 2  # ub and bd
labels_by_weight = V2.weight_dictionary()
assert (1, 0, 1, 0) in labels_by_weight[(-1, 1)]
assert V2.weight_elements(3, -3) == ()
```

A weight tuple has $n$ signed integer entries summing to zero; an exponent
tuple has $n^2$ nonnegative entries. The weight lookup methods take separate
entries, so use `V2.weight_elements(*weight)` for a tuple variable. A dictionary
value from `weight_dictionary` is a tuple of exponent labels, not a subspace.

```python
whole2 = Subspace.whole(V2)
assert whole2.is_diagonal_stable()
pieces = whole2.weight_decomposition()
assert pieces[(0, 0)].dimension == 4
assert sum(piece.dimension for piece in pieces.values()) == 10
```

The decomposition has dimensions $1,2,4,2,1$ at weights $(-2,2)$, $(-1,1)$,
$(0,0)$, $(1,-1)$, and $(2,-2)$, respectively. Each component retains `V2`
as its ambient parent. Components of dimension zero are omitted; the zero
weight is included when it occurs. `weight_decomposition` requires diagonal
stability and raises `ValueError` otherwise. The identity acts trivially, so
stability under all diagonals and under traceless diagonals is equivalent here.

## 7. Compute the quadratic invariant in the embedded sl(2) symmetric power

Start by computing invariants, rather than supplying the answer:

```python
invariants_gl2 = whole2.invariants(Subspace.whole(g))
assert invariants_gl2.dimension == 2
invariant_line = invariants_gl2.intersection_with_symmetric_sl()
assert invariant_line.dimension == 1
computed = invariant_line.basis()[0]
Q = (Rational(1) / computed.poly.coeff_monomial(u**2)) * computed
Q_expression = (u - d)**2 + 4*b*c
assert Q == V2.from_poly(Poly(Q_expression, *V2.generators, domain=QQ))
assert all(act(y, Q) == V2.zero() for y in g.basis())
```

The basis algorithm can choose any nonzero multiple of the invariant. Dividing
by its $u^2$ coefficient selects the normalization in the displayed formula.

For $Y\subseteq S^k(\mathfrak{gl}_n)$, `Y.invariants(X)` means

$$
Y^X=\{v\in Y:x\cdot v=0\text{ for every }x\in X\}.
$$

The subspace $Y$ need not be stable under $X$, and $X$ need not be a Lie
subalgebra. This computes annihilated vectors, not just a subspace preserved
by the action. In this example, invariance under $\mathfrak{gl}_2$ is equivalent
to invariance under $\mathfrak{sl}_2$, since scalar matrices act as zero.

The two quadratic invariants of $\mathfrak{gl}_2$ can also be written as
$\tau^2=(u+d)^2$ and $\Delta=ud-bc$. The invariant in $S^2(\mathfrak{sl}_2)$ is

$$
Q=\tau^2-4\Delta=(u-d)^2+4bc.
$$

To construct the entire embedded symmetric power, use:

```python
symmetric_sl2 = whole2.intersection_with_symmetric_sl()
assert symmetric_sl2.dimension == 6
assert symmetric_sl2.ambient == V2
assert symmetric_sl2.invariants(Subspace.whole(g)) == invariant_line
```

The intersection retains $S^2(\mathfrak{gl}_2)$ as its ambient. It does not
create a separate symmetric-power parent for $\mathfrak{sl}_2$.
For positive degree, the defining condition is

$$
\left(\frac{\partial}{\partial u}+\frac{\partial}{\partial d}\right)f=0.
$$

Indeed, with $h=u-d$ and $\tau=u+d$, the operator is $2\partial/\partial\tau$;
its kernel is the polynomial algebra in $h,b,c$. It imposes homogeneous linear
equations on polynomial coefficients. Evaluating at the identity alone would
not impose this condition.

## 8. Trace-form evaluation and its kernel

Evaluation uses the identification $y\mapsto(A\mapsto\mathrm{tr}(yA))$.
Thus $z_{ij}(A)=A_{ji}$, with the indices transposed. The quadratic $bc$ hides
this transpose, so a degree-one example makes it visible:

```python
V1 = SymmetricPower(g, 1)
b1 = V1.monomial(0, 1, 0, 0)
assert evaluate(b1, g.E(1, 0)) == 1
assert evaluate(b1, g.E(0, 1)) == 0
assert evaluate(Q, x) == 28
```

For $x=\begin{pmatrix}1&2\\3&-1\end{pmatrix}$, the last value is
$(1-(-1))^2+4\cdot3\cdot2=28$. Notice the argument order: `act(x, v)` but
`evaluate(v, x)`.

```python
kernel = symmetric_sl2.evaluation_kernel(x)
assert kernel.dimension == 5
assert kernel.ambient == V2
assert all(evaluate(v, x) == 0 for v in kernel.basis())
assert invariant_line.evaluation_kernel(x).dimension == 0
assert invariant_line.evaluation_kernel(g.E(0, 1)) == invariant_line
```

Evaluation is a linear map to $\mathbb Q$. Its restricted kernel therefore
has codimension zero or one. It does not require an action-stability assumption.

## 9. Weight (2,-2) in degree eight

```python
V8 = SymmetricPower(g, 8)
weight = (2, -2)
W = Subspace.from_basis(V8, *V8.weight_elements(*weight))
H = W.invariants(T)
W_sl = W.intersection_with_symmetric_sl()
H_sl = W_sl.invariants(T)
assert V8.dimension == 165
assert (W.dimension, H.dimension) == (16, 4)
assert (W_sl.dimension, H_sl.dimension) == (4, 1)
assert H.ambient == V8 and H_sl.ambient == V8
```

The weight $(2,-2)$ has eigenvalue $4$ on $h=E_{00}-E_{11}$. Its highest-weight
irreducible is consequently the **five-dimensional** $\mathfrak{sl}_2$ module
of highest weight $4$. The two-coordinate weight in this package is not the
single fundamental-weight coordinate for $\mathfrak{sl}_2$.

The full weight space is generally not $T$-stable: the raising operator sends
weight $(2,-2)$ to weight $(3,-3)$. The `invariants` method can nevertheless
compute its kernel, which is exactly what we need. It suffices to check the
single basis element $e$ of $T$.

For a hand count, a degree-eight monomial of weight $(2,-2)$ is

$$
u^a b^{j+2}c^j d^r,\qquad a+r=6-2j,\qquad 0\leq j\leq3.
$$

There are $7+5+3+1=16$ possibilities. In $S^8(\mathfrak{sl}_2)$ the corresponding
weight space has the four-element basis

$$
b^2h^6,\quad b^3ch^4,\quad b^4c^2h^2,\quad b^5c^3.
$$

Here $h$ denotes the degree-one symmetric generator $u-d$. Because
$e\cdot h=-2b$, $e\cdot c=h$, and $e\cdot b=0$, a combination with coefficients
$(A_0,A_1,A_2,A_3)$ is killed by $e$ precisely when

$$
-12A_0+A_1=0,\qquad -8A_1+2A_2=0,\qquad -4A_2+3A_3=0.
$$

The solution space is one-dimensional, with coefficients $(1,12,48,64)$.
Its generator is

$$
b^2Q^3=b^2h^6+12b^3ch^4+48b^4c^2h^2+64b^5c^3.
$$

In $S^8(\mathfrak{gl}_2)$ the four highest-weight vectors can be taken as

$$
b^2\tau^6,\qquad b^2\tau^4Q,\qquad b^2\tau^2Q^2,\qquad b^2Q^3.
$$

Both $\tau$ and $Q$ are invariant, and $e$ kills $b$, so every displayed vector
is killed by $e$. They are independent because their powers of $\tau$ differ
and $Q$ depends only on $h,b,c$. For completeness, the same coefficient
recurrence in degree $m$ gives one such $\mathfrak{sl}_2$ highest vector when
$m\geq2$ is even and none when $m$ is odd. Expanding by powers of $\tau$ then
gives precisely these four vectors in total degree eight.

We can compare their span with the computed answer without requiring the
package to return the same normalized basis:

```python
u8, b8, c8, d8 = V8.generators
tau8 = u8 + d8
q8_expression = (u8 - d8)**2 + 4*b8*c8
explicit = tuple(
    V8.from_poly(Poly(b8**2 * tau8**(6 - 2*j) * q8_expression**j,
                     *V8.generators, domain=QQ))
    for j in range(4)
)
assert H == Subspace.from_basis(V8, *explicit)
assert H_sl == Subspace.from_basis(V8, explicit[-1])
```

Finite-dimensional representations of $\mathfrak{sl}_2$ over $\mathbb Q$ are
completely reducible. Each irreducible contributes one line of highest-weight
vectors at its own highest weight. Therefore the dimensions $4$ and $1$ count
copies of the highest-weight-$4$ irreducible in the two degree-eight spaces.
The full weight multiplicities remain $16$ and $4$, respectively.

## 10. Degree zero, unchecked constructors, and errors

Degree zero is supported: $S^0(\mathfrak g)=\mathbb Q$. Its only basis monomial
has every exponent zero, evaluation returns its constant value, and every
Lie element acts as zero.

```python
V0 = SymmetricPower(g, 0)
one = V0.monomial(0, 0, 0, 0)
assert V0.dimension == 1
assert act(e, one) == V0.zero()
assert evaluate(one, x) == 1
assert Subspace.whole(V0).intersection_with_symmetric_sl() == Subspace.whole(V0)
```

The low-level constructors `LieElement(parent, matrix)`,
`SymmetricElement(parent, poly)`, and `Subspace(ambient, dimension, basis_matrix)`
store their arguments without validation. They are documented API points, but
the parent methods and subspace factories are the ordinary construction route.
If you already have valid data, these examples illustrate their contracts:

```python
assert LieElement(g, x.matrix) == x
assert SymmetricElement(V2, uc.poly) == uc
assert Subspace(g, sl2.dimension, sl2.basis_matrix) == sl2
```

The direct element constructors require compatible shapes, coefficient types,
generators, and degree. The direct subspace constructor requires independent
basis columns and a dimension matching their number. Treat parent attributes
and subspace storage as read-only after construction because derived data may
be cached. `str` is available for `GeneralLinear`, `LieElement`, and
`SymmetricElement`; for an unambiguous polynomial display use `.poly.as_expr()`.

Common input rules are:

| Situation | Rule |
| --- | --- |
| Matrix index | Python integer, excluding booleans, between zero and `n - 1`. |
| Exponent tuple | $n^2$ nonnegative Python integers summing to the degree; unpack for `monomial`. |
| Coordinate matrix | `ImmutableMatrix` column with exactly the ambient dimension in rows. |
| Equation matrix | `ImmutableMatrix` with exactly the ambient dimension in columns. |
| Spanning sequence | Pass separate compatible elements, or unpack the sequence with `*`. |
| Scalar coefficient | Exact supported scalar; no floats or booleans. |
| Weight lookup | Supply $n$ signed integer entries summing to zero; an absent weight gives an empty tuple. |
| Weight decomposition | Requires a diagonally stable subspace of a `SymmetricPower`. |
| Membership | A valid element with an incompatible parent gives `False`; an unsupported object type is an error. |

The exact validation exceptions are documented on each method. The walkthrough
uses valid inputs and checks its mathematical outputs; the package's test suite
is the place to test rejected inputs and exception types.

## 11. API coverage map

This table maps every public method in the package modules to a worked
section. It excludes the private action helper and `_validation` helpers.

| Object | Public API | Section |
| --- | --- | --- |
| `GeneralLinear` | Constructor; `n`, `dimension`; equality and string conversion | 2 |
| `GeneralLinear` | `E`, `basis`, `zero` | 2 |
| `GeneralLinear` | `coordinates`, `from_coordinates` | 3 |
| `GeneralLinear` | `diagonals`, `traceless_diagonals`, `strictly_upper_triangulars` | 3 |
| `LieElement` | `parent`, `matrix`; equality, addition, subtraction, negation, scalar and matrix multiplication; `bracket` | 2 |
| `LieElement` | Unchecked constructor; string conversion | 10 |
| `lie.kronecker_delta` | Rational Kronecker delta | 2 |
| `SymmetricPower` | Constructor; `algebra`, `degree`, `dimension`, `generators`; equality | 4; companion script |
| `SymmetricPower` | `basis`, `basis_labels`, `zero`, `monomial`, `from_terms`, `from_poly` | 4 |
| `SymmetricPower` | `coordinates`, `from_coordinates` | 4 |
| `SymmetricPower` | `weight_from_exponent_tuple`, `weight_elements`, `weight_element_indices`, `weight_dictionary` | 6 |
| `SymmetricElement` | `parent`, `poly`; equality, addition, subtraction, negation, scalar multiplication | 4 |
| `SymmetricElement` | Unchecked constructor; string conversion | 10 |
| Actions | `act`, `action_matrix`, `evaluate` | 5, 8 |
| `Subspace` | `from_basis`, `from_equations`, `whole`, `trivial`, `D`, `D0`, `T` | 3 |
| `Subspace` | `ambient`, `dimension`, `basis_matrix`, `basis`, `equation_matrix`, equality, `contains` | 3 |
| `Subspace` | `is_diagonal_stable`, `weight_decomposition` | 6 |
| `Subspace` | `invariants`, `intersection_with_symmetric_sl` | 7, 9 |
| `Subspace` | `evaluation_kernel` | 8 |
| `Subspace` | Unchecked constructor | 10 |

## 12. Expected walkthrough output

The default run ends with the following results. Additional expanded bases
appear only with `--show-basis`; their normalization can differ from the
compact bases displayed in the mathematics.

```text
g = GeneralLinear(2); dimension = 4
dim D, D0, T, sl2 = 2, 1, 1, 3
Quadratic invariant Q = (u - d)**2 + 4*b*c
dim S^2(gl2)^gl2 = 2; dim S^2(sl2)^sl2 = 1
Q([[1, 2], [3, -1]]) = 28
dim evaluation kernel in S^2(sl2) at this matrix = 5

Degree 8; weight (2, -2)
dim S^8(gl2) = 165
gl2: weight-space dimension = 16; T-invariant dimension = 4
sl2: weight-space dimension = 4; T-invariant dimension = 1
The weight has eigenvalue 4 on h = E00 - E11.
Its T-invariant dimension is the multiplicity of the 5-dimensional irreducible with highest weight (2,-2).
Explicit gl2 highest-weight basis:
  b**2 * tau**6 * Q**0
  b**2 * tau**4 * Q**1
  b**2 * tau**2 * Q**2
  b**2 * tau**0 * Q**3
Explicit sl2 highest-weight basis: b**2 * Q**3

All walkthrough checks passed.
```

For SymPy details used here, see its [polynomial reference](https://docs.sympy.org/latest/modules/polys/reference.html)
and [matrix reference](https://docs.sympy.org/latest/modules/matrices/matrices.html).
The invariant and multiplicity calculations above are derived explicitly in
this tutorial rather than relying on a graded-multiplicity formula.
