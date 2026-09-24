# Clustered, Longitudinal and Split-Plot Data

Most data arrive in groups. Plots lie in fields, pupils sit in classrooms,
patients return for repeated visits, firms are followed for twenty years. Two
observations from the same group are more alike than two from different groups,
and the error covariance is no longer \( \sigma^2\I \). This chapter is the
applied half of Part VII: it takes the general Gauss–Markov theory of
[Chapter 31](../ch31-general-gauss-markov/index.html) and the mixed models of
[Chapter 32](../ch32-linear-mixed-models/index.html) and puts them to work on the
three structures that account for most correlated data in practice. A split-plot
experiment randomizes at two levels and therefore has two error terms, so some
comparisons are precise and others are not. A repeated-measures analysis needs a
condition on the covariance that is weaker than equicorrelation and stronger than
nothing, and there are three standard ways of dealing with its failure. A
longitudinal study needs a covariance model chosen from a small catalogue, checked
against the data, and hedged with a standard error that survives its being wrong.
The chapter closes with growth curves, where each subject has its own trajectory
and the whole point is to borrow strength across subjects.

**What you need.** [Chapter 6](../ch06-projections/index.html) (projection in a
general inner product, @def-proj-gls, and when it agrees with least squares,
@thm-proj-kruskal),
[Chapter 9](../ch09-sums-of-squares/index.html) (expected mean
squares, @thm-ss-expected-mean-squares),
[Chapter 16](../ch16-multiway-layouts/index.html) (balanced layouts as Kronecker
products, @lem-tw-kron, @thm-tw-balanced),
[Chapter 18](../ch18-covariance-and-design/index.html) (randomized
blocks, @thm-dsn-rcbd),
[Chapter 19](../ch19-theory-of-departures/index.html) (the effect of a wrong
covariance, @thm-dep-covariance),
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (sandwich
estimators, @thm-het-sandwich, autoregressive errors, @def-het-ar1),
[Chapter 29](../ch29-model-selection/index.html) (@def-sel-aic-bic), and the two
chapters before this one: generalized least squares and its
inference (@def-ggm-gls, @thm-ggm-inference, @thm-ggm-ols-blue) and the mixed
model with its predictors and likelihoods (@def-mix-model,
@thm-mix-blup, @thm-mix-henderson, @def-mix-reml).
The distribution theory is that of Part I, used directly and not only through those chapters: the
noncentral \( F \), the mean and variance of a quadratic form and the independence of orthogonal
projections from [Chapter 4](../ch04-quadratic-forms/index.html) (@def-qf-noncentral-f,
@thm-qf-mean-var, @thm-qf-orthogonal-projections); whitening, the covariance of a linear function
and the mean of a quadratic form from [Chapter 2](../ch02-random-vectors/index.html)
(@prp-rv-whitening, @thm-rv-linear, @thm-rv-quadform-mean); the Woodbury identity and a matrix
square root from [Chapter 1](../ch01-matrix-algebra/index.html) (@thm-mat-woodbury,
@thm-mat-square-root); Aitken's corollary from [Chapter 7](../ch07-optimality/index.html)
(@cor-opt-aitken); and, in [Section 33.3](03-repeated-measures.html), the centred Wishart lemma of
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) (@lem-cor-centred-wishart).

## Roadmap

- [Split-plot designs](01-split-plot.html)
- [Clustered data and compound symmetry](02-clustered-data.html)
- [Repeated measures](03-repeated-measures.html)
- [Covariance models for longitudinal data](04-covariance-models.html)
- [Growth curves](05-growth-curves.html)
- [Summary and notes](06-summary-and-notes.html)
