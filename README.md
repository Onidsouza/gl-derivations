# General Linear Derivations

The purpose of this package is to implement the action of the rational Lie algebra gl(n) of n-by-n matrices on its k-th symmetric product S(k), using what it can from Python's symbolic computation capabilities. It should be able to:

- Handle addition, subtraction and scalar multiplication for elements of gl(n) and S(k), as well as providing constructors and easy access to canonical bases of each.
- Given x in gl(n) and v in S(k), compute xv in terms of the canonical basis of S(k).
- Implement the diagonal subalgebra D(n), D0(n) of traceless diagonal matrices and the strictly upper triangular subalgebra T(n).
- Implement subspaces of S(k), both in terms of defining equations and given bases. Be able to pass from one to another. Test if a given subspace is stable under the D(n) action or D0(n) action. Decompose this subspace into eigenspaces for both the D(n) and the D0(n) action.
- Implement subspaces of gl(n). Given a subspace X of gl(n), compute the invariant subspace S(k)^X.
- Using the trace form as a map gl(n) -> gl(n)^*, implement the pairing gl(n) x gl(n) -> Q and thus S(k) x gl(n) -> Q by seeing the elements of S(k) as polynomial functions on gl(n) through the trace form.
- Given x in gl(n) and a subspace Y of S(k), compute the kernel of the map Y -> Q given by evaluating Y at x.

# State of development

- [X] Setting up environment, files and configuration.
- [X] Input validation.
- [X] Setting up test environment.
- [X] Lie elements, coordinates and bracket.
- [ ] Symmetric powers and monomial coordinates.
- [ ] The adjoint derivation and action matrices.
- [ ] Subspaces from bases and equations.
- [ ] Weights, stability and eigenspace decomposition.
- [ ] Simultaneous environments.
- [ ] Trace evaluation and evaluation kernels.
- [ ] Teste catalogue and debugging
- [ ] Final README, API and mathematical documentation.
- [ ] Automated checks
- [ ] Build, isolated installation and publication.
