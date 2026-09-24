# Random Effects and Linear Mixed Models

Part IV treated the levels of a factor as chosen: five nitrogen rates, three
machines, the treatments the experimenter cared about. Many factors are not
like that. The twelve laboratories in a proficiency study, the eleven firms in a
panel, the patients in a trial and the classrooms in a school district are a
*sample* of levels, and the question is not how laboratory 3 compares with
laboratory 7 but how much laboratories vary. Declaring such a factor random
changes three things at once. The parameters of interest become variances. The
individual levels are *predicted* rather than estimated, and the prediction is
shrunk towards the overall mean. And the covariance of the response stops being
\( \sigma^2\I \): observations sharing a level are correlated, which is what
Part VII is about. This chapter develops the model that results, the best linear
unbiased predictor and Henderson's equations that compute it, restricted maximum
likelihood, and what inference survives when the covariance itself has been
estimated from the same data.

**What you need.** [Chapter 6](../ch06-projections/index.html) (projections and
generalized least squares, @def-proj-gls, @thm-proj-kruskal),
[Chapter 7](../ch07-optimality/index.html) (the Gauss–Markov and Aitken theorems,
@thm-opt-gauss-markov and @cor-opt-aitken; maximum likelihood, @thm-opt-mle; the
conjugate posterior, @thm-opt-bayes-conjugate),
[Chapter 4](../ch04-quadratic-forms/index.html) (distributions of quadratic forms,
@thm-qf-chisq and @thm-qf-indep-quadratic),
[Chapter 9](../ch09-sums-of-squares/index.html) (expected mean squares,
@thm-ss-expected-mean-squares and @exr-ss-random-oneway),
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) (best linear
unbiased prediction, @thm-cor-blup, resting on the best linear predictor
@thm-rv-blp of [Chapter 2](../ch02-random-vectors/index.html)),
[Section 15.1](../ch15-anova-subspaces/01-oneway-model.html) (the one-way layout and
the first look at random effects, @def-aov-oneway) and
[Chapter 16](../ch16-multiway-layouts/index.html) (nested and
crossed layouts, @thm-tw-nested). From Part V,
[Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html)
(weighted least squares with estimated weights) and
[Section 23.3](../ch23-resampling-inference/03-confidence-intervals.html) (the
parametric bootstrap, @prp-bs-parametric); from Part VI,
[Chapter 27](../ch27-shrinkage/index.html) (ridge and empirical Bayes,
@thm-shr-ridge and @prp-shr-empirical-bayes) and
[Chapter 29](../ch29-model-selection/index.html) (AIC and BIC, @def-sel-aic-bic).
The generalized least squares theory of
[Chapter 31](../ch31-general-gauss-markov/index.html) is used throughout.
Underneath the whole chapter sits the partitioned-matrix algebra of
[Chapter 1](../ch01-matrix-algebra/index.html) — the Woodbury identity, the block determinant, the
partitioned inverse, the cyclic trace, a matrix square root and the calculus of quadratic forms
with Lagrange multipliers (@thm-mat-woodbury, @eq-mat-woodbury, @thm-mat-block-determinant,
@thm-mat-partitioned-inverse, @thm-mat-trace-cyclic, @thm-mat-square-root, @prp-mat-quadratic-min,
@prp-mat-lagrange) — which is what turns Henderson's equations and the likelihood into something
computable; and, in [Section 32.4](04-likelihood.html), the normal density and its behaviour under
a linear map from [Chapter 3](../ch03-multivariate-normal/index.html) (@thm-mvn-density,
@thm-mvn-linear). [Section 32.3](03-blup.html) points forward to the smoothing spline as a mixed
model ([Section 43.6](../ch43-smoothing/06-mixed-model-and-bayes.html)); that is a promise, not a
prerequisite.

## Roadmap

- [Random effects](01-random-effects.html)
- [The linear mixed model](02-mixed-model.html)
- [Prediction: BLUP and Henderson's equations](03-blup.html)
- [Maximum likelihood and REML](04-likelihood.html)
- [Inference when the covariance is estimated](05-inference.html)
- [Bayesian thread III: mixed models as hierarchical models](06-bayes.html)
- [Summary and notes](07-summary-and-notes.html)
