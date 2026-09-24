# Optimality: Gauss–Markov, Maximum Likelihood and Minimum Variance

Least squares was introduced as a way of fitting, and [Chapter 6](../ch06-projections/index.html)
showed that it is a projection. Neither fact says it is a *good* way of estimating. This chapter
asks in what precise senses least squares is the best one can do, and where those senses end.
With only means and covariances assumed, least squares has the smallest variance among
linear unbiased estimators (Gauss–Markov), and the proof is one more application of the projection
theorem. With normal errors three more things hold. Least squares is maximum likelihood. It has
the smallest variance among *all* unbiased estimators (Lehmann–Scheffé). And its sampling
distribution is known exactly, which is what Part III builds tests and intervals on. The
same assumptions also show where optimality stops. Nonlinear estimators can win when the errors are
not normal, and biased estimators can win on mean squared error. A Bayesian analysis with a
conjugate prior is one principled biased estimator, and it recovers least squares as the prior
becomes flat.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html): nonnegative definite
matrices with their characterizations and square roots (@def-mat-nnd, @prp-mat-pd-properties,
@thm-mat-pd-characterizations, @thm-mat-square-root), which is what a comparison of two
covariance matrices is a statement about; and, in single proofs, generalized inverses
(@prp-mat-ginverse-props), idempotent matrices (@prp-mat-idempotent-basic) and the Woodbury
identity (@thm-mat-woodbury).
Covariances of linear functions and the mean of a quadratic form
([Chapter 2](../ch02-random-vectors/index.html): @thm-rv-linear and @thm-rv-quadform-mean). The
multivariate normal and its density ([Chapter 3](../ch03-multivariate-normal/index.html):
@thm-mvn-linear and @thm-mvn-density). Chi-squared quadratic forms and their independence
([Chapter 4](../ch04-quadratic-forms/index.html): @thm-qf-chisq and @thm-qf-indep-linear).
The linear model and its first two moments ([Chapter 5](../ch05-model-and-least-squares/index.html):
@def-lm-linear-model and @thm-lm-moments). Least squares as projection, including the rank-deficient
case ([Chapter 6](../ch06-projections/index.html): @thm-proj-ls-projection,
@thm-proj-normal-equations and @thm-proj-invariant-functions).

## Roadmap

- [Linear unbiased estimation](01-linear-unbiased.html)
- [The Gauss–Markov theorem](02-gauss-markov.html)
- [Maximum likelihood under normality](03-maximum-likelihood.html)
- [Sufficiency, completeness and minimum variance](04-minimum-variance.html)
- [Sampling distributions of the estimates](05-sampling-distributions.html)
- [Bayesian thread I: the conjugate linear model](06-bayes-conjugate.html)
- [Summary and notes](07-summary-and-notes.html)
