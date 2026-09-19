# Regularization and Boosting

Ridge regression and the lasso are two members of one family: least squares plus a penalty on the
coefficients. This chapter studies the family as a whole: when a minimizer exists and is unique, how to recognize one,
what prior it corresponds to, and what is gained by giving up convexity. Two penalties then carry structure the lasso
cannot see: the elastic net keeps correlated regressors together, and the group lasso selects whole factors. The second
half turns to boosting, which fits a regression by many small steps along one coordinate at a time. Run to the end, it
reproduces least squares; stopped early, it is a shrinkage estimator whose path, with small steps, is almost the lasso
path. The number of steps plays the part of the penalty.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html) (the singular value decomposition,
@thm-mat-svd, and the Cauchy–Schwarz inequality, @prp-mat-cauchy-schwarz), [Chapter 6](../ch06-projections/index.html)
(least squares as projection, @thm-proj-ls-projection, and the minimum-norm solution, @prp-proj-min-norm),
[Chapter 7](../ch07-optimality/index.html) (the conjugate posterior, @thm-opt-bayes-conjugate and @prp-opt-shrinkage),
and indicator coding of factors ([Chapter 15](../ch15-anova-subspaces/index.html)). From this part:
[Chapter 27](../ch27-shrinkage/index.html) (ridge regression, @thm-shr-ridge, and the lasso, @def-shr-lasso and @thm-shr-lasso-orthonormal) and [Chapter 29](../ch29-model-selection/index.html) (covariance penalties,
@thm-sel-optimism, and cross-validation). [Chapter 28](../ch28-high-dimensional/index.html) is needed only for remarks.

## Roadmap

- [Penalized least squares in general](01-penalized-least-squares.html)
- [The elastic net](02-elastic-net.html)
- [Group penalties](03-group-penalties.html)
- [Componentwise boosting](04-componentwise-boosting.html)
- [Boosting as regularization](05-boosting-as-regularization.html)
- [Summary and notes](06-summary-and-notes.html)
