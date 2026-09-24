# The General Gauss–Markov Model

Every chapter so far has assumed that the errors are uncorrelated and equally
variable. That assumption is a statement about the *geometry* of observation
space: it says that Euclidean distance is the right way to measure how far the
data are from the model. When observations share a common shock, follow one
another in time, or are measured with different precision, Euclidean distance is
the wrong ruler. This chapter replaces \( \sigma^2\I \) by a general
\( \sigma^2\V \): first with \( \V \) known and invertible, where the theory of
Parts II and III survives a change of coordinates; then in the two places where
it does not survive, when \( \V \) is singular and some functions of \( \bbeta \)
become known exactly, and when \( \V \) has to be estimated from the data that
produced the fit.

**What you need.** [Chapter 6](../ch06-projections/index.html), especially
[Section 6.9](../ch06-projections/09-inner-products.html): projection in the
\( \V^{-1} \) inner product (@def-proj-gls, @thm-proj-A-projection), oblique projections (@prp-proj-oblique)
and Kruskal's theorem (@thm-proj-kruskal).
[Section 7.2](../ch07-optimality/02-gauss-markov.html) for the Gauss–Markov theorem (@thm-opt-gauss-markov)
and Aitken's corollary (@cor-opt-aitken), which this chapter cites
and does not reprove. From [Chapter 1](../ch01-matrix-algebra/index.html): generalized
inverses, square roots (@thm-mat-square-root) and the Woodbury identity (@thm-mat-woodbury).
From [Chapter 2](../ch02-random-vectors/index.html): whitening (@prp-rv-whitening)
and the support of a singular covariance (@thm-rv-cov-nnd).
Estimability ([Chapter 8](../ch08-estimability/index.html)) and the \( F \) test
([Chapter 11](../ch11-general-linear-hypothesis/index.html)) are used throughout, and
the distribution theory of [Section 31.1](01-generalized-least-squares.html) is the quadratic-form
theory of [Chapter 4](../ch04-quadratic-forms/index.html), reached through @thm-opt-sampling applied
to the whitened model rather than cited here directly.
[Section 19.2](../ch19-theory-of-departures/02-wrong-covariance.html) and
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) supply the
diagnostic side: what a wrong covariance does (@thm-dep-covariance), how much
efficiency it costs (@thm-dep-efficiency), and the sandwich estimator (@thm-het-sandwich).
[Section 31.1](01-generalized-least-squares.html) also builds its interval on the \( t \) interval
for an estimable function from [Chapter 12](../ch12-intervals-and-bands/index.html)
(@thm-ci-estimable-interval).

## Roadmap

- [Generalized least squares](01-generalized-least-squares.html)
- [When ordinary least squares is best](02-ols-blue.html)
- [Singular covariance matrices](03-singular-covariance.html)
- [Estimated covariance: feasible generalized least squares](04-feasible-gls.html)
- [Summary and notes](05-summary-and-notes.html)
