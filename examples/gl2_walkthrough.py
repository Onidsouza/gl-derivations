#!/usr/bin/env python3
"""Walk through the public API and compute invariants for gl(2) and sl(2).

Read the accompanying TUTORIAL.md for the mathematics and API explanations.
Install the project first, then run:

    python gl2_walkthrough.py
    python gl2_walkthrough.py --degree 8 --show-basis

The assertions check the examples using exact rational arithmetic.
"""

from argparse import ArgumentParser
from fractions import Fraction

from sympy import QQ, ImmutableMatrix, Poly, Rational

from glder import (
    GeneralLinear,
    LieElement,
    Subspace,
    SymmetricElement,
    SymmetricPower,
    act,
    action_matrix,
    evaluate,
)
from glder.lie import kronecker_delta


def lie_algebra_example():
    """Construct gl(2), use its arithmetic, and construct standard subspaces."""
    g = GeneralLinear(2)
    assert (g.n, g.dimension) == (2, 4)
    assert g == GeneralLinear(2)

    e00, e01, e10, e11 = g.basis()
    e, f, h = e01, e10, e00 - e11
    identity = e00 + e11
    assert e == g.E(0, 1)
    assert e.parent == g
    assert e.matrix == ImmutableMatrix([[0, 1], [0, 0]])
    assert e + g.zero() == e
    assert e + (-e) == g.zero()
    assert h.bracket(e) == 2 * e
    assert h.bracket(f) == -2 * f
    assert e.bracket(f) == h
    assert e * f == e00  # This product is matrix multiplication.
    assert Fraction(1, 2) * h == h * Rational(1, 2)
    assert kronecker_delta(0, 0) == Rational(1)
    assert kronecker_delta(0, 1) == Rational(0)

    x = e00 + 2 * e + 3 * f - e11
    coords = g.coordinates(x)
    assert coords == ImmutableMatrix([1, 2, 3, -1])
    assert g.from_coordinates(coords) == x
    assert g.from_coordinates(1, 2, 3, -1) == x

    # Direct element construction is unchecked. Reuse already valid data here.
    assert LieElement(g, x.matrix) == x
    assert isinstance(str(g), str) and isinstance(str(x), str)

    D = Subspace.D(g)
    D0 = Subspace.D0(g)
    T = Subspace.T(g)
    assert (D.dimension, D0.dimension, T.dimension) == (2, 1, 1)
    assert D == Subspace.from_basis(g, *g.diagonals())
    assert D0 == Subspace.from_basis(g, *g.traceless_diagonals())
    assert T == Subspace.from_basis(g, *g.strictly_upper_triangulars())

    # A homogeneous equation on the matrix-entry coordinate column.
    sl2 = Subspace.from_equations(g, ImmutableMatrix([[1, 0, 0, 1]]))
    assert sl2.dimension == 3
    assert sl2 == Subspace.from_basis(g, e, f, h)
    assert sl2.contains(h) and not sl2.contains(identity)
    assert sl2.ambient == g
    assert sl2.basis_matrix.shape == (4, 3)
    assert Subspace(g, sl2.dimension, sl2.basis_matrix) == sl2  # Unchecked.
    assert sl2.equation_matrix().shape == (1, 4)
    assert (sl2.equation_matrix() * sl2.basis_matrix).is_zero_matrix
    assert len(sl2.basis()) == sl2.dimension
    assert Subspace.whole(g).dimension == 4
    assert Subspace.trivial(g).dimension == 0
    assert Subspace.from_basis(g) == Subspace.trivial(g)
    assert Subspace.from_basis(g, e, 2 * e, g.zero()) == T

    print(f"g = {g}; dimension = {g.dimension}")
    print("dim D, D0, T, sl2 = 2, 1, 1, 3")
    return g, T, x


def quadratic_example(g, x):
    """Use all polynomial constructors, actions, weights, and quadratic kernels."""
    V2 = SymmetricPower(g, 2)
    u, b, c, d = V2.generators  # E00, E01, E10, E11, in this order.
    assert V2.algebra == g and V2.degree == 2 and V2.dimension == 10
    assert V2 == SymmetricPower(GeneralLinear(2), 2)
    assert V2.basis_labels() == (
        (2, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1),
        (0, 2, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1), (0, 0, 2, 0),
        (0, 0, 1, 1), (0, 0, 0, 2),
    )
    assert len(V2.basis()) == V2.dimension

    # Three constructors for the same quadratic element.
    alpha = (1, 0, 1, 0)
    uc = V2.monomial(*alpha)
    assert uc == V2.from_terms({alpha: 1})
    assert uc == V2.from_poly(Poly(u * c, *V2.generators, domain=QQ))
    assert uc.parent == V2 and uc.poly.as_expr() == u * c
    assert SymmetricElement(V2, uc.poly) == uc  # Unchecked constructor.
    assert isinstance(str(uc), str)
    assert V2.from_terms({}) == V2.zero()
    assert uc + V2.zero() == uc and uc - uc == V2.zero()
    assert uc + (-uc) == V2.zero()
    assert uc * Rational(1, 2) == Fraction(1, 2) * uc
    assert V2.from_terms({alpha: Fraction(1, 2)}) == Rational(1, 2) * uc

    coords = V2.coordinates(uc)
    assert coords.shape == (10, 1)
    assert V2.from_coordinates(coords) == uc
    assert V2.from_coordinates(*tuple(coords)) == uc

    # The adjoint derivation: e.u = -b and e.c = u-d.
    e = g.E(0, 1)
    expected_action = V2.from_poly(Poly(u**2 - u * d - b * c,
                                      *V2.generators, domain=QQ))
    assert act(e, uc) == expected_action
    M = action_matrix(e, V2)
    assert M.shape == (10, 10)
    assert M * coords == V2.coordinates(act(e, uc))

    # Diagonal weights are signed n-tuples, not exponent tuples.
    assert V2.weight_from_exponent_tuple(*alpha) == (-1, 1)
    labels_by_weight = V2.weight_dictionary()
    assert alpha in labels_by_weight[(-1, 1)]
    weight_monomials = V2.weight_elements(1, -1)
    indices = V2.weight_element_indices(1, -1)
    assert weight_monomials == tuple(V2.basis()[i] for i in indices)
    assert len(weight_monomials) == 2
    assert V2.weight_elements(3, -3) == ()

    whole2 = Subspace.whole(V2)
    assert Subspace.trivial(V2).basis() == ()
    assert whole2.is_diagonal_stable()
    pieces = whole2.weight_decomposition()
    assert {weight: piece.dimension for weight, piece in pieces.items()} == {
        (-2, 2): 1, (-1, 1): 2, (0, 0): 4, (1, -1): 2, (2, -2): 1,
    }
    assert sum(piece.dimension for piece in pieces.values()) == V2.dimension
    for weight, piece in pieces.items():
        assert piece.ambient == V2
        assert piece == Subspace.from_basis(V2, *V2.weight_elements(*weight))
    assert Subspace.from_basis(V2, *whole2.basis()) == whole2
    assert Subspace.from_equations(V2, whole2.equation_matrix()) == whole2

    # Compute the invariant, rather than assuming its formula.
    invariants_gl2 = whole2.invariants(Subspace.whole(g))
    assert invariants_gl2.dimension == 2
    invariant_line = invariants_gl2.intersection_with_symmetric_sl()
    assert invariant_line.dimension == 1
    computed = invariant_line.basis()[0]
    coefficient = computed.poly.coeff_monomial(u**2)
    assert coefficient != 0
    Q = (Rational(1) / coefficient) * computed
    Q_expression = (u - d)**2 + 4 * b * c
    assert Q == V2.from_poly(Poly(Q_expression, *V2.generators, domain=QQ))
    assert all(act(generator, Q) == V2.zero() for generator in g.basis())
    assert invariant_line == Subspace.from_basis(V2, Q)
    assert invariant_line.ambient == V2

    # Recover the same invariant line from the familiar two gl(2) invariants.
    tau_squared = V2.from_poly(Poly((u + d)**2, *V2.generators, domain=QQ))
    determinant = V2.from_poly(Poly(u * d - b * c, *V2.generators, domain=QQ))
    assert invariants_gl2 == Subspace.from_basis(V2, tau_squared, determinant)
    assert Q == tau_squared - 4 * determinant

    symmetric_sl2 = whole2.intersection_with_symmetric_sl()
    assert symmetric_sl2.dimension == 6
    assert symmetric_sl2.ambient == V2
    assert symmetric_sl2.invariants(Subspace.whole(g)) == invariant_line
    assert symmetric_sl2.intersection_with_symmetric_sl() == symmetric_sl2

    # Trace-form evaluation uses z_ij(x) = x_ji.
    V1 = SymmetricPower(g, 1)
    b1 = V1.monomial(0, 1, 0, 0)
    assert evaluate(b1, g.E(1, 0)) == 1
    assert evaluate(b1, g.E(0, 1)) == 0
    assert evaluate(Q, x) == 28
    kernel = symmetric_sl2.evaluation_kernel(x)
    assert kernel.dimension == 5 and kernel.ambient == V2
    assert all(evaluate(v, x) == 0 for v in kernel.basis())
    assert invariant_line.evaluation_kernel(x).dimension == 0
    assert invariant_line.evaluation_kernel(g.E(0, 1)) == invariant_line

    # Degree zero consists of constants, and every Lie element acts as zero.
    V0 = SymmetricPower(g, 0)
    one = V0.monomial(0, 0, 0, 0)
    assert V0.dimension == 1 and one == V0.basis()[0]
    assert act(e, one) == V0.zero()
    assert evaluate(one, x) == 1
    assert Subspace.whole(V0).intersection_with_symmetric_sl() == Subspace.whole(V0)

    print("Quadratic invariant Q = (u - d)**2 + 4*b*c")
    print("dim S^2(gl2)^gl2 = 2; dim S^2(sl2)^sl2 = 1")
    print("Q([[1, 2], [3, -1]]) = 28")
    print("dim evaluation kernel in S^2(sl2) at this matrix = 5")
    return Q


def highest_weight_example(g, T, Q, degree, show_basis=False):
    """Compute weight (2,-2), then impose annihilation by T, in both ambients."""
    V = SymmetricPower(g, degree)
    weight = (2, -2)
    weight_space = Subspace.from_basis(V, *V.weight_elements(*weight))
    highest_gl2 = weight_space.invariants(T)

    # Restrict the weight space first, keeping the computation smaller than
    # constructing the whole embedded S^degree(sl2).
    weight_space_sl2 = weight_space.intersection_with_symmetric_sl()
    highest_sl2 = weight_space_sl2.invariants(T)
    assert highest_gl2.ambient == V and highest_sl2.ambient == V

    J = (degree - 2) // 2
    weight_dimension_gl2 = sum(degree - 2 * j - 1 for j in range(J + 1))
    assert weight_space.dimension == weight_dimension_gl2
    assert highest_gl2.dimension == J + 1
    assert weight_space_sl2.dimension == J + 1
    assert highest_sl2.dimension == (1 if degree % 2 == 0 else 0)

    u, b, c, d = V.generators
    assert V.generators == Q.parent.generators
    tau = u + d
    Q_expr = Q.poly.as_expr()
    explicit_highest = tuple(
        V.from_poly(Poly(b**2 * tau**(degree - 2 - 2 * j) * Q_expr**j,
                         *V.generators, domain=QQ))
        for j in range(J + 1)
    )
    assert highest_gl2 == Subspace.from_basis(V, *explicit_highest)
    assert highest_sl2 == highest_gl2.intersection_with_symmetric_sl()

    if degree % 2 == 0:
        explicit_sl2 = V.from_poly(
            Poly(b**2 * Q_expr**J, *V.generators, domain=QQ)
        )
        assert highest_sl2 == Subspace.from_basis(V, explicit_sl2)
        assert act(g.E(0, 1), explicit_sl2) == V.zero()

    print(f"\nDegree {degree}; weight {weight}")
    print(f"dim S^{degree}(gl2) = {V.dimension}")
    print(f"gl2: weight-space dimension = {weight_space.dimension}; "
          f"T-invariant dimension = {highest_gl2.dimension}")
    print(f"sl2: weight-space dimension = {weight_space_sl2.dimension}; "
          f"T-invariant dimension = {highest_sl2.dimension}")
    print("The weight has eigenvalue 4 on h = E00 - E11.")
    print("Its T-invariant dimension is the multiplicity of the "
          "5-dimensional irreducible with highest weight (2,-2).")
    print("Explicit gl2 highest-weight basis:")
    for j in range(J + 1):
        print(f"  b**2 * tau**{degree - 2 - 2 * j} * Q**{j}")
    if degree % 2 == 0:
        print(f"Explicit sl2 highest-weight basis: b**2 * Q**{J}")
    else:
        print("The sl2 highest-weight space is zero in this odd degree.")
    if show_basis:
        print("Computed gl2 basis, expanded by the package:")
        for vector in highest_gl2.basis():
            print(f"  {vector.poly.as_expr()}")
        print("Computed sl2 basis, expanded by the package:")
        for vector in highest_sl2.basis():
            print(f"  {vector.poly.as_expr()}")


def main():
    """Run the examples; accept an optional degree of at least two."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=8,
                        help="symmetric degree for the weight calculation (default: 8)")
    parser.add_argument("--show-basis", action="store_true",
                        help="also print the package's expanded highest-weight bases")
    args = parser.parse_args()
    if args.degree < 2:
        parser.error("--degree must be at least 2 for this walkthrough")

    g, T, x = lie_algebra_example()
    Q = quadratic_example(g, x)
    highest_weight_example(g, T, Q, args.degree, args.show_basis)
    print("\nAll walkthrough checks passed.")


if __name__ == "__main__":
    main()
