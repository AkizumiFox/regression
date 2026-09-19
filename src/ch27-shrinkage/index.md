# Shrinkage: Ridge, Principal Components and the Lasso

Least squares is the best linear unbiased estimator, and under normal errors the best unbiased
estimator of any kind. With nearly dependent columns it can still be very bad, because its
variance along poorly determined directions swamps everything else. This chapter gives up unbiasedness on purpose, measures estimators by mean
squared error, and studies three ways of buying a large cut in variance with a little bias. Ridge
regression shrinks every direction of the coefficient space, the weakly determined ones most;
principal component regression discards the weakly determined directions; the lasso shrinks
coefficients with an absolute-value penalty and sets some exactly to zero. The chapter ends with
Stein's paradox: in three or more dimensions a shrinkage estimator beats least squares for every
value of the parameter, for an empirical Bayes reason that also explains the other methods.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html) (the singular value
decomposition, @thm-mat-svd), [Chapter 6](../ch06-projections/index.html) (projections and the
variance inflation identity @eq-proj-vif-preview), [Chapter 7](../ch07-optimality/index.html)
(the Gauss–Markov theorem, @thm-opt-gauss-markov, and the conjugate normal model,
@thm-opt-bayes-conjugate and @prp-opt-shrinkage), [Chapter 10](../ch10-computation/index.html)
(least squares through the SVD, @thm-cmp-svd-ls and @prp-cmp-tsvd),
[Chapter 19](../ch19-theory-of-departures/index.html) (the short regression as a biased
estimator, @thm-dep-mse), [Chapter 20](../ch20-residuals-leverage-influence/index.html) (the
leave-one-out formulas, @thm-res-deletion) and
[Chapter 26](../ch26-collinearity/index.html) (what collinearity does to variances, @thm-col-variance). The James–Stein section uses the noncentral chi-squared distribution of
[Chapter 4](../ch04-quadratic-forms/index.html) (@thm-qf-ncchisq).

## Roadmap

- [The bias–variance trade-off](01-mean-squared-error.html)
- [Ridge regression](02-ridge.html)
- [Principal component regression](03-principal-components.html)
- [The lasso](04-lasso.html)
- [Stein's paradox and James–Stein estimation](05-james-stein.html)
- [Summary and notes](06-summary-and-notes.html)
