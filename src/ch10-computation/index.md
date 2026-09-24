# Computation: QR, SVD and Conditioning

Every least squares estimate in this book is written as
\( \hbeta=(\X\T\X)^{-1}\X\T\y \) or as \( \X(\X\T\X)\ginv\X\T\y \). These are the right
formulas to reason with. They are not how a computer should produce the numbers.
This chapter explains what a computer should do instead, and how accurate the result is. It
separates two questions: how sensitive the least squares problem is to small changes in the
data, which depends on \( \kappa(\X) \) and the size of the residual, and whether an algorithm
adds more error than that. The answers give a practical ranking. The Cholesky factorization of \( \X\T\X \) is
fast and is enough when \( \X \) is well conditioned. Householder or Givens QR is the
default. The singular value decomposition is the tool for rank deficiency and near
dependence. Pivoted QR decides the rank cheaply, which is how statistical software
chooses what to do with aliased columns.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html): the Cholesky factorization
(@thm-mat-pd-characterizations), QR (@thm-mat-qr), the singular value decomposition
(@thm-mat-svd), the Eckart–Young inequality (@prp-mat-svd-norms) and the condition number
(@def-mat-condition-number). [Chapter 5](../ch05-model-and-least-squares/index.html): the linear
model and the moments of \( \hbeta \) (@thm-lm-moments). [Chapter 6](../ch06-projections/index.html):
least squares as projection, and the first look at computation in
[Section 6.10](../ch06-projections/10-computation.html). [Section 10.4](04-svd.html) uses the covariance of a linear
transformation from [Chapter 2](../ch02-random-vectors/index.html) (@thm-rv-linear), and that is
the chapter's only call on random vectors. A few remarks use the sampling theory of
[Chapter 7](../ch07-optimality/index.html) and the sums of squares of
[Chapter 9](../ch09-sums-of-squares/index.html).

## Roadmap

- [Why not invert the cross-product matrix](01-why-not-invert.html)
- [The Cholesky factorization](02-cholesky.html)
- [Householder and Givens QR](03-householder-givens.html)
- [The singular value decomposition](04-svd.html)
- [Conditioning and perturbation bounds](05-conditioning.html)
- [Updating and downdating](06-updating.html)
- [Rank-revealing factorizations](07-rank-revealing.html)
- [Summary and notes](08-summary-and-notes.html)
