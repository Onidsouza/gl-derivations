# Record of module design choices.

- All computations have the field Q of rational numbers as their base field.
- All computations use exact arithmetic in Q, no floating point input is accepted.
- The space Q^n denotes the space of *column* vectors with n rows. Matrices act on column vectors by multiplying on the left.
- Variable indices will follow the Python convention: the vector space Q^n has as a basis the elementary vectors e_0,...,e_{n-1}. The elementary matrices E(i,j) has a 1 in the i-th row, j-th column, and 0 everywhere else. Its indices run from 0 to n-1.
- The adjoint action of gl(n) on itself is given by ad(X)(Y) = [X,Y] = XY-YX.