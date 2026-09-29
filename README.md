# glder

**Exact adjoint derivations on symmetric powers of general linear Lie algebras.**

`glder` is a Python library for computations over the rational numbers in
$\mathfrak{gl}_n$ and $S^k(\mathfrak{gl}_n)$. It constructs Lie elements,
homogeneous symmetric tensors, action matrices, weight spaces, and subspaces
defined by generators or linear equations. SymPy supplies exact matrices,
polynomials, and rational linear algebra.

The action is the adjoint action extended by the Leibniz rule:

$$
x\cdot(y_1\cdots y_k)
=\sum_{r=1}^{k}y_1\cdots[x,y_r]\cdots y_k.
$$

Typical computations include finding invariant polynomials, extracting
highest-weight vectors, imposing evaluation conditions, and intersecting with
the embedded subspace $S^k(\mathfrak{sl}_n)$.

## Installation

Use a Python version supported by the release's `pyproject.toml`. Installing
the package also installs its declared runtime dependencies, including SymPy.

For release **0.1.0**, download the wheel from this repository's GitHub
**Releases** page. From the directory containing the downloaded file, run:

```bash
python -m pip install ./glder-0.1.0-py3-none-any.whl
```

Alternatively, install from the repository root or an extracted source
distribution:

```bash
python -m pip install .
```

For development, use an editable installation:

```bash
python -m pip install -e .
python -m pip install pytest
python -m pytest
```

## A quadratic invariant in sl(2)

Let $e=E_{01}$, $f=E_{10}$, and $h=E_{00}-E_{11}$. Then
$[h,e]=2e$, $[h,f]=-2f$, and $[e,f]=h$. The symmetric tensor

$$
Q=h^2+4ef
$$

spans the invariant line in $S^2(\mathfrak{sl}_2)$. The following code computes
that line and recovers this normalization:

```python
from sympy import Poly, QQ
from glder import GeneralLinear, SymmetricPower, Subspace, act

g = GeneralLinear(2)
S2 = SymmetricPower(g, 2)
u, b, c, d = S2.generators

sl2_symmetric_square = Subspace.whole(S2).intersection_with_symmetric_sl()
invariant_line = sl2_symmetric_square.invariants(Subspace.whole(g))
assert invariant_line.dimension == 1

q = invariant_line.basis()[0]
q = (1 / q.poly.coeff_monomial(u**2)) * q
expected = S2.from_poly(Poly((u - d)**2 + 4*b*c, *S2.generators, domain=QQ))
assert q == expected
assert all(act(x, q) == S2.zero() for x in g.basis())
print(q.poly.as_expr())
```

Here `u`, `b`, `c`, and `d` represent $E_{00}$, $E_{01}$, $E_{10}$, and
$E_{11}$ in the symmetric algebra. These products are commutative products
in $S(\mathfrak{gl}_2)$; matrix multiplication and multiplication in the
universal enveloping algebra are different operations.

## A weight space and its highest-weight vectors

The weight $(2,-2)$ means that $E_{00}$ acts by $2$ and $E_{11}$ acts by $-2$.
To extract its highest-weight vectors, additionally require annihilation by
$T=\mathbb Q E_{01}$, the strictly upper triangular Lie subalgebra:

```python
S8 = SymmetricPower(g, 8)
weight_space = Subspace.from_basis(S8, *S8.weight_elements(2, -2))
highest_vectors = weight_space.invariants(Subspace.T(g))

print(weight_space.dimension)       # 16
print(highest_vectors.dimension)    # 4

sl2_weight_space = weight_space.intersection_with_symmetric_sl()
sl2_highest_vectors = sl2_weight_space.invariants(Subspace.T(g))

print(sl2_weight_space.dimension)       # 4
print(sl2_highest_vectors.dimension)    # 1
```

Thus the weight multiplicity and the dimension of highest-weight vectors are
different quantities:

| Ambient representation | Weight multiplicity of $(2,-2)$ | Dimension annihilated by $T$ in that weight |
| --- | ---: | ---: |
| $S^8(\mathfrak{gl}_2)$ | 16 | 4 |
| $S^8(\mathfrak{sl}_2)$ | 4 | 1 |

In this notation, **`T` means strictly upper triangular matrices**. The diagonal
subalgebra is `Subspace.D(g)`. Since $(2,-2)$ is nonzero, these vectors are not
annihilated by all diagonal matrices.

## Public interface

The main classes and functions are available directly from the package:

```python
from glder import (
    GeneralLinear,
    LieElement,
    SymmetricPower,
    SymmetricElement,
    Subspace,
    act,
    action_matrix,
    evaluate,
)
```

| Object | Purpose |
| --- | --- |
| `GeneralLinear`, `LieElement` | Construct matrices over $\mathbb Q$ and compute Lie brackets. |
| `SymmetricPower`, `SymmetricElement` | Work with a fixed homogeneous degree and its monomial coordinates. |
| `act`, `action_matrix` | Apply an adjoint derivation or represent it by a matrix. |
| `evaluate` | Evaluate a symmetric tensor using the trace pairing. |
| `Subspace` | Construct subspaces, split weights, and impose annihilation or evaluation conditions. |

## Conventions

- Matrix indices start at **zero**. Elementary matrices and polynomial
  generators are ordered by rows.
- Scalars are exact rational numbers. Use `Rational(1, 2)` or `Fraction(1, 2)`
  for one half; the Python expression `1/2` is a float and is not accepted as
  a rational input.
- Coordinates are columns. Column $j$ of an action matrix is the coordinate
  vector of the action on basis vector $j$.
- Evaluation uses $z_{ij}(x)=\mathrm{tr}(E_{ij}x)=x_{ji}$.
- `invariants(X)` means vectors **annihilated by every element of `X`**.
- Construct elements and subspaces through the documented factories, and
  treat their public attributes and cached data as read-only.

The dimension $\binom{n^2+k-1}{k}$ grows rapidly. Version 0.1 focuses on exact,
explicit computations; constructing large bases or dense action matrices can
be expensive.

## Documentation

- [TUTORIAL.md](./docs/TUTORIAL.mdTUTORIAL.md): the mathematical conventions and a worked tour
  of the public API.
- [gl2_walkthrough.py](./examples/gl2_walkthrough.py): executable examples and assertions.
  After installation, run `python ./examples/gl2_walkthrough.py` from the repository root.
- [API.md](./docs/API.md): constructors, public attributes, methods, and return values.
- [DECISIONS.md](./docs/DECISIONS.md): design choices and the
  [0.1.0 release procedure](./docs/DECISIONS.md#releasing-010-on-github).

The tutorial derives the quadratic invariant and explains the degree-eight
dimension counts, including the difference between weight multiplicity and
irreducible multiplicity.
