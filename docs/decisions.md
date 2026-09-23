# Record of module design choices.

- All computations have the field Q of rational numbers as their base field.
- All computations use exact arithmetic in Q, no floating point input is accepted.
- The space Q^n denotes the space of *column* vectors with n rows. Matrices act on column vectors by multiplying on the left.
- Variable indices will follow the Python convention: the vector space Q^n has as a basis the elementary vectors e_0,...,e_{n-1}. The elementary matrices E(i,j) has a 1 in the i-th row, j-th column, and 0 everywhere else. Its indices run from 0 to n-1.
- The adjoint action of gl(n) on itself is given by ad(X)(Y) = [X,Y] = XY-YX.
- S(k) denotes the k-th symmetric product of the vector space gl(n). We denote by z(i,j) the element of S(1) corresponding to the matrix E(i,j), since these two objects should not be identified. The vector space S(k) is the span of the commuting products of k elements of the form z(i,j). We extend the adjoint action of gl(n) on S(1) to an action on S(k) by the derivation rule x(z_1z_2) = (xz_1)z_2 + z_1(xz_2) for any x in gl(n), z_1 in S(j) and z_2 in S(k-j).
- The ordering of the elementary matrices follows row order: the lexographic ordering on the pair (i,j), or, in other words, matrices are read left to right, top to bottom. For n = 2, their ordering is E(0,0), E(0,1), E(1,0), E(1,1). This is a total order on the chosen basis of gl(n).
- The monomial order in S(k) is given by labeling the exponents into an n^2-tuple of integers according to the order on the elementary matrices, then following the lexographic ordering on this tuple. For example, using a,b,c,d for the matrices E(0,0), E(0,1), E(1,0) and E(1,1) respectively, the ordering on the monomials in S(2) for n = 2 is: a^2, ab, ac, ad, b^2, bc, bd, c^2, cd, d^2.
- Our convention is that S(0) = Q, and the adjoint action of gl(n) on S(0) is trivial (i.e. zero for every element).