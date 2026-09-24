# Non-Normality, Heteroscedasticity and Serial Correlation

The exact theory of Part III rests on three assumptions about the errors: they are normal, they have a common
variance, and they are uncorrelated. [Chapter 19](../ch19-theory-of-departures/index.html) worked out what
least squares does when each of them fails. This chapter is about the data analyst's side of the same
question: how to see a failure in the residuals, how large it must be before it matters, and what to do
about it. The three assumptions do not matter equally. Coefficient inference survives nonnormal errors in
moderately large samples, but prediction intervals do not. Unequal variances leave least squares unbiased but
make its usual standard errors wrong, sometimes by a factor of four or more, and they can be repaired either
by weighting or by a variance estimator that does not trust the model. Serial correlation does the same
damage in time-ordered data, and it often signals a mean function that is missing its dynamics. Every
diagnostic in the chapter is computed from residuals, which are not the errors. How they differ is where the
chapter starts.

**What you need.** Residuals, leverage and generalized least squares from
[Chapter 6](../ch06-projections/index.html) (@prp-proj-fit-residual, @prp-proj-leverage,
@def-proj-gls, @thm-proj-A-projection), the Gauss–Markov and Aitken theorems and the normal sampling theory of
[Chapter 7](../ch07-optimality/index.html) (@thm-opt-gauss-markov, @cor-opt-aitken, @thm-opt-sampling), the
moments of quadratic forms from
[Chapter 2](../ch02-random-vectors/index.html) (@thm-rv-quadform-mean, @thm-rv-quadform-variance) and their normal distribution theory from
[Chapter 4](../ch04-quadratic-forms/index.html) (@thm-qf-canonical, @cor-qf-mahalanobis), the deletion
formulas of [Chapter 10](../ch10-computation/index.html) (@prp-cmp-loo), and prediction intervals from
[Chapter 12](../ch12-intervals-and-bands/index.html) (@thm-ci-prediction-interval) and their large-sample
coverage under nonnormal errors (@prp-ci-pi-limit). From
[Chapter 19](../ch19-theory-of-departures/index.html) we use the covariance of least squares under a general
error covariance (@thm-dep-covariance) and its large-sample behaviour under nonnormal
errors (@thm-dep-nonnormal). From [Chapter 20](../ch20-residuals-leverage-influence/index.html) we use studentized
residuals (@def-res-residuals, @thm-res-external-t). [Section 21.3](03-weighted-least-squares.html) uses the moments of least squares and of \( s^2 \)
from [Chapter 5](../ch05-model-and-least-squares/index.html) (@thm-lm-moments, @thm-lm-sigma2), and
[Section 21.5](05-durbin-watson.html) the extremal Rayleigh quotient of
[Chapter 1](../ch01-matrix-algebra/index.html) (@thm-mat-extremal-rayleigh) together with the
behaviour of a normal vector under a linear map (@thm-mvn-linear,
[Chapter 3](../ch03-multivariate-normal/index.html)). The limit theorems of probability that the chapter needs
are gathered, or pointed to, in [Section 21.2](02-detecting-heteroscedasticity.html).

## Roadmap

- [Checking normality](01-checking-normality.html)
- [Detecting heteroscedasticity](02-detecting-heteroscedasticity.html)
- [Weighted least squares](03-weighted-least-squares.html)
- [Heteroscedasticity-consistent covariance estimators](04-sandwich-estimators.html)
- [Serial correlation and the Durbin–Watson test](05-durbin-watson.html)
- [Autoregressive errors](06-autoregressive-errors.html)
- [Summary and notes](07-summary-and-notes.html)
