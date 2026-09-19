# The Model, Least Squares and the Normal Equations

A regression model is a claim about how the distribution of a response changes
with a set of regressors, and the linear model is the simplest useful claim of
that kind. It says that the mean of the response is a linear combination of known
columns, that the errors around this mean have a common variance and are
uncorrelated, and nothing more. This chapter sets the model up and draws out what
follows from those assumptions alone. It starts with what the assumptions do and do not
say, fits a straight line by calculus, and passes to the matrix form of the model,
where polynomials, harmonics, indicator variables and interactions all turn out to be
linear models. It then derives the least squares estimator and the normal equations,
computes the mean and covariance of the estimator, estimates the error variance
without bias, and asks what a fitted coefficient means. The last section surveys the
models the rest of the book builds on this one. Throughout, the model matrix has full
column rank. [Chapter 6](../ch06-projections/index.html) retells the same story
geometrically and removes that condition.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html): rank
(@prp-mat-rank-product), inverses, nonnegative definite matrices
(@prp-mat-pd-properties), traces (@thm-mat-trace-cyclic) and the derivatives of
linear and quadratic forms (@thm-mat-quadform-derivative).
[Chapter 2](../ch02-random-vectors/index.html): means and covariances of linear
transformations (@thm-rv-linear), the mean and variance of a quadratic form
(@thm-rv-quadform-mean, @thm-rv-quadform-variance) and the laws of total
expectation and covariance (@prp-rv-total-covariance). Normality
([Chapter 3](../ch03-multivariate-normal/index.html)) appears only in remarks.

## Roadmap

- [What a regression model claims](01-what-a-model-claims.html)
- [Simple linear regression as a warm-up](02-simple-regression.html)
- [The multiple linear regression model](03-multiple-regression-model.html)
- [Least squares and the normal equations](04-least-squares.html)
- [Properties of the estimator under second-moment assumptions](05-properties.html)
- [Estimating the error variance](06-error-variance.html)
- [Interpreting coefficients](07-interpreting-coefficients.html)
- [A panorama of regression models to come](08-panorama.html)
- [Summary and notes](09-summary-and-notes.html)
