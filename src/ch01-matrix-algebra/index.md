# Matrix Algebra for Linear Models

A linear model is a statement about vectors and matrices, and nearly
every argument in this book is an argument about them. This chapter collects the
matrix algebra the book relies on: rank and the fundamental subspaces, inverses and
generalized inverses, consistent systems of equations, orthogonal matrices, determinants
and partitioned matrices, the spectral theorem and positive definiteness, the singular
value decomposition, idempotent matrices, Kronecker products and the calculus of
quadratic forms. It is a reference chapter, meant to be skimmed first and consulted
later, but every central result is proved.

**What you need.** None beyond the Notation chapter. A first course in linear algebra
helps, but is not assumed.

**A slower road through the same algebra.** This chapter is a reference. It states what
the rest of the book uses, in the order a regression argument needs it, and proves it —
but it moves fast, and it stays with real matrices because that is all a linear model
requires. The companion book [*Linear Algebra*](https://linear-algebra.akizumifox.com/)
builds the same material from the beginning, in more generality (arbitrary fields, linear
maps rather than matrices) and at a gentler pace, with the motivation and the examples
that a reference chapter has to leave out. A reader who wants that road can take the two
side by side:

| Section here | Read there |
|---|---|
| [1.1 Vectors and matrices](01-basics.html) | [Ch. 1, Vector Spaces](https://linear-algebra.akizumifox.com/ch01-vector-spaces/index.html); [§2.1, Three views of a linear system](https://linear-algebra.akizumifox.com/ch02-linear-systems/01-linear-systems.html) |
| [1.2 Rank and the fundamental subspaces](02-rank.html) | [§2.5](https://linear-algebra.akizumifox.com/ch02-linear-systems/05-rank.html); [§3.3, rank–nullity](https://linear-algebra.akizumifox.com/ch03-linear-maps/03-rank-nullity.html); [§3.8](https://linear-algebra.akizumifox.com/ch03-linear-maps/08-rank-factorization.html) |
| [1.3 Inverses and generalized inverses](03-inverses.html) | [§2.4](https://linear-algebra.akizumifox.com/ch02-linear-systems/04-elementary-matrices-and-inverses.html); [§12.11, the Moore–Penrose inverse](https://linear-algebra.akizumifox.com/ch12-psd-and-svd/11-pseudoinverse.html) |
| [1.4 Systems and consistency](04-systems.html) | [§2.2](https://linear-algebra.akizumifox.com/ch02-linear-systems/02-gaussian-elimination.html); [§10.4, least squares and minimum-norm solutions](https://linear-algebra.akizumifox.com/ch10-inner-products/04-least-squares.html) |
| [1.5 Orthogonality](05-orthogonality.html) | [§10.1](https://linear-algebra.akizumifox.com/ch10-inner-products/01-inner-products.html)–[§10.3](https://linear-algebra.akizumifox.com/ch10-inner-products/03-orthogonal-complements-and-projections.html); [§10.8, QR and Householder](https://linear-algebra.akizumifox.com/ch10-inner-products/08-qr-and-householder.html) |
| [1.6 Trace, determinants, partitioned matrices](06-partitioned.html) | [Ch. 6, Determinants](https://linear-algebra.akizumifox.com/ch06-determinants/index.html); [§3.10, trace](https://linear-algebra.akizumifox.com/ch03-linear-maps/10-projections-and-trace.html); [§7.3, Schur complements](https://linear-algebra.akizumifox.com/ch07-block-matrices/03-block-determinants-and-schur-complements.html) |
| [1.7 Eigenvalues and positive definiteness](07-eigen.html) | [Ch. 8](https://linear-algebra.akizumifox.com/ch08-eigenvalues/index.html); [§11.5, the real spectral theorem](https://linear-algebra.akizumifox.com/ch11-spectral-theory/05-spectral-theorem-real.html); [§12.1](https://linear-algebra.akizumifox.com/ch12-psd-and-svd/01-positive-definite-matrices.html)–[§12.2](https://linear-algebra.akizumifox.com/ch12-psd-and-svd/02-square-roots-and-cholesky.html); [§16.1, the Rayleigh quotient](https://linear-algebra.akizumifox.com/ch16-variational/01-the-rayleigh-quotient.html) |
| [1.8 The SVD and conditioning](08-svd.html) | [§12.8](https://linear-algebra.akizumifox.com/ch12-psd-and-svd/08-singular-value-decomposition.html); [§12.10, low-rank approximation](https://linear-algebra.akizumifox.com/ch12-psd-and-svd/10-low-rank-approximation.html); [Ch. 19, perturbation](https://linear-algebra.akizumifox.com/ch19-perturbation/index.html) |
| [1.9 Idempotent matrices](09-idempotent.html) | [§3.10](https://linear-algebra.akizumifox.com/ch03-linear-maps/10-projections-and-trace.html); [§12.12, projections revisited](https://linear-algebra.akizumifox.com/ch12-psd-and-svd/12-projections-revisited.html) |
| [1.10 Kronecker products and vec](10-kronecker.html) | [§7.4](https://linear-algebra.akizumifox.com/ch07-block-matrices/04-kronecker-product.html); [Ch. 14, tensors](https://linear-algebra.akizumifox.com/ch14-tensors/index.html) |

Section 1.11 has no counterpart there: differentiating a quadratic form is calculus
written in matrix notation, and the identities we need are proved where they are used.

## Roadmap

- [Vectors, matrices and the operations we use](01-basics.html)
- [Rank, column space and null space](02-rank.html)
- [Inverses and generalized inverses](03-inverses.html)
- [Systems of equations and consistency](04-systems.html)
- [Orthogonality and orthogonal matrices](05-orthogonality.html)
- [Trace, determinants and partitioned matrices](06-partitioned.html)
- [Eigenvalues, the spectral theorem and positive definiteness](07-eigen.html)
- [The singular value decomposition and conditioning](08-svd.html)
- [Idempotent matrices](09-idempotent.html)
- [Kronecker products and the vec operator](10-kronecker.html)
- [Vector and matrix calculus](11-calculus.html)
- [Summary and notes](12-summary-and-notes.html)
