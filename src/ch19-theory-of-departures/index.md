# What Goes Wrong: A Theory of Departures

Parts II to IV derived everything from a short list of assumptions: the mean is
\( \X\bbeta \) for the regressors actually used, the errors are uncorrelated with a common variance, they are
normal, and the regressors are fixed and measured without error. Real data meet each of these
only approximately. This chapter asks, one assumption at a time, what a specific
departure does to the least squares estimator, to \( s^2 \) and to the tests and
intervals built from them. Omitting a regressor biases the
coefficients, and yet the biased estimator can still have the smaller mean squared error. A wrong
covariance matrix leaves the estimator nearly efficient but can make its standard error
wrong by a factor of two or more. With nonnormal errors the \( t \) and \( F \) tests are
asymptotically valid when no observation has large leverage, while tests and intervals for
\( \sigma^2 \) are not valid in any sample size. A single bad observation can move the
estimate arbitrarily far. With random regressors and a curved mean, the usual standard errors describe
the wrong target. Collinearity breaks no assumption at all, but it inflates variances. Chapters [20](../ch20-residuals-leverage-influence/index.html)–[22](../ch22-transformations/index.html) turn to diagnosis and remedies, and
[Chapter 23](../ch23-resampling-inference/index.html) to inference that relies less on the assumptions.

**What you need.** [Chapter 5](../ch05-model-and-least-squares/index.html) (the moments of
\( \hbeta \) under departures, @prp-lm-misspecified, and omitted regressors, @prp-lm-omitted),
[Chapter 6](../ch06-projections/index.html) (the Frisch–Waugh–Lovell theorem, leverage, Kruskal's
theorem and the population projection), [Chapter 7](../ch07-optimality/index.html) (Gauss–Markov and
Aitken), [Chapter 10](../ch10-computation/index.html) (the deletion formulas, @prp-cmp-loo), and the
tests and intervals of [Chapter 11](../ch11-general-linear-hypothesis/index.html) and
[Chapter 12](../ch12-intervals-and-bands/index.html). [Section 14.1](../ch14-correlation-lack-of-fit-prediction/01-random-regressors.html)
treats random regressors under the normal model, and [Section 19.5](05-random-regressors.html) builds on it.
Two earlier chapters carry more of the weight than the list above suggests.
[Chapter 2](../ch02-random-vectors/index.html) supplies the moments recomputed under every
departure — the covariance of a linear function, the mean and variance of a quadratic form, the
expected residual sum of squares and the law of total covariance (@thm-rv-linear,
@thm-rv-quadform-mean, @thm-rv-quadform-variance, @eq-rv-rss-mean, @prp-rv-total-covariance) — and
[Chapter 1](../ch01-matrix-algebra/index.html) the spectral theorem and extremal Rayleigh quotient
that [Sections 19.2](02-wrong-covariance.html) and [19.6](06-collinearity.html) measure efficiency
with, together with the QR factorization and the Cauchy–Schwarz inequality (@thm-mat-spectral,
@thm-mat-extremal-rayleigh, @thm-mat-qr, @prp-mat-cauchy-schwarz).
[Section 19.2](02-wrong-covariance.html) points forward to feasible generalized least squares
([Section 31.4](../ch31-general-gauss-markov/04-feasible-gls.html)); that is a promise, not a
prerequisite.

## Roadmap

- [Omitting and adding regressors](01-underfitting-overfitting.html)
- [An incorrect covariance matrix](02-wrong-covariance.html)
- [Nonnormal errors](03-non-normal-errors.html)
- [A single bad observation](04-outliers.html)
- [Random regressors and a curved mean](05-random-regressors.html)
- [Collinearity, previewed](06-collinearity.html)
- [Summary and notes](07-summary-and-notes.html)
