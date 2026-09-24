# Polynomials, Piecewise Fits and the Road to Splines

For forty-one chapters the mean of the response has been a linear function of
columns chosen in advance. This part keeps least squares, the projection
geometry and every distribution theory built on them, and gives up one thing
only: the insistence that the columns be the regressors themselves. The oldest
way to do that is to let them be powers of a regressor, which makes a polynomial
regression a linear model in a transformed basis, so that Parts II and III apply
word for word — and that is what makes its failures instructive, because they
are failures of the *family of shapes*, not of the fitting.

A polynomial is a global object: raising its degree to capture a feature at one
end of the data changes the curve everywhere, sends the variance soaring at the
edges, and turns any extrapolation into fiction. Worse, the obvious basis
\( 1,x,x^2,\dots \) is one of the worst-conditioned matrices in numerical
analysis. This chapter separates the two complaints. Orthogonal polynomials
repair the arithmetic completely and the statistics not at all; knots repair the
statistics, and lead to splines.

**What you need.** [Chapter 6](../ch06-projections/index.html) (projection,
reparameterization, leverage), [Chapter 9](../ch09-sums-of-squares/index.html)
(sequential sums of squares, @exm-ss-nile-polynomials),
[Section 1.8](../ch01-matrix-algebra/08-svd.html) (the condition
number, @def-mat-condition-number) and
[Chapter 10](../ch10-computation/index.html) (what conditioning does to a fit), [Chapter 11](../ch11-general-linear-hypothesis/index.html)
(restricted least squares), [Chapter 12](../ch12-intervals-and-bands/index.html)
(bands) and [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) (lack
of fit). The later sections borrow the vocabulary of
[Chapter 29](../ch29-model-selection/index.html) without needing its theory, and the
closing discussion uses the bias–variance trade-off of
[Chapter 19](../ch19-theory-of-departures/index.html), the diagnostics of
[Chapter 20](../ch20-residuals-leverage-influence/index.html), the transformations of
[Chapter 22](../ch22-transformations/index.html) and the causal reading of
[Chapter 25](../ch25-causal-interpretation/index.html).
[Section 42.4](04-regression-splines.html) points forward to the penalized spline of
[Section 43.3](../ch43-smoothing/03-penalized-splines.html); that is a promise, not a
prerequisite, and nothing here is proved from it. Of the chapters listed in the closing discussion,
[Chapter 20](../ch20-residuals-leverage-influence/index.html) is a pointer only: it is named for the
diagnostics a fitted spline should be put through, and nothing in this chapter cites it.

## Roadmap

- [Polynomial regression and its instability](01-polynomial-regression.html)
- [Orthogonal polynomials](02-orthogonal-polynomials.html)
- [Piecewise polynomials](03-piecewise-polynomials.html)
- [Regression splines and truncated power bases](04-regression-splines.html)
- [Summary and notes](05-summary-and-notes.html)
