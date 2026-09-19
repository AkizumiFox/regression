# High-Dimensional Regression: p > n

Every chapter so far has assumed, openly or quietly, that there are more observations than
coefficients. Once \( p>n \) the columns of \( \X \) cannot be linearly independent, least squares
fits the data exactly whatever they are, and no individual coefficient is estimable. This chapter
asks what can still be done. The first answer is to keep least squares and pick the shortest of its
many solutions. Its risk can be computed exactly, and it produces the "double descent" curve. The second answer is to assume that only a few coefficients matter and
to use the lasso. Here a restricted form of the eigenvalue conditions of Part II replaces full rank,
and short non-asymptotic proofs show that the lasso predicts and estimates almost as well as an
oracle that knew which coefficients were zero. Selecting exactly the right variables needs a
stronger condition, which can fail even when prediction is good. The chapter ends with inference.
Intervals computed after selection are wrong, and a one-step correction of the lasso restores
approximately normal estimates of single coefficients.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html) (the singular value decomposition,
@thm-mat-svd, and the Moore–Penrose inverse, @def-mat-moore-penrose),
[Chapter 6](../ch06-projections/index.html) (projections, minimum-norm least squares,
@prp-proj-min-norm, and the Frisch–Waugh–Lovell theorem, @thm-proj-fwl),
[Chapter 8](../ch08-estimability/index.html) (estimable functions, @thm-est-estimable-identifiable),
[Chapter 10](../ch10-computation/index.html) (least squares through the SVD, @thm-cmp-svd-ls),
[Chapter 13](../ch13-multiplicity/index.html) (the union bound, @thm-mc-bonferroni) and
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) (random regressors and the inverse Wishart
mean, @lem-cor-inverse-wishart). From this part: the lasso, its optimality conditions and soft
thresholding ([Chapter 27](../ch27-shrinkage/index.html), @def-shr-lasso and @thm-shr-lasso-orthonormal)
and ridge regression (@thm-shr-ridge). Selection bias and
cross-validation ([Chapter 29](../ch29-model-selection/index.html)) are referred to but not needed.

## Roadmap

- [Minimum-norm interpolation](01-minimum-norm.html)
- [Restricted eigenvalue conditions](02-restricted-eigenvalues.html)
- [Prediction and estimation bounds for the lasso](03-lasso-bounds.html)
- [Support recovery](04-support-recovery.html)
- [Post-selection and debiased inference](05-debiased-inference.html)
- [Summary and notes](06-summary-and-notes.html)
