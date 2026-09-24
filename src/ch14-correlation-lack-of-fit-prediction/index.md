# Correlation, Lack of Fit and Prediction Theory

So far the regressors have been fixed numbers. In most observational studies they are random, drawn
along with the response. This chapter shows that the exact tests and intervals of Part III survive
that change, because their conditional laws do not depend on the design, while the laws of the
estimates, and the power of tests, do change. With random rows a regression is a prediction rule, and
the chapter measures how well an estimated predictor predicts new cases and how much its in-sample
error flatters it. It then studies the correlations that measure linear association: the multiple
correlation, whose sample version \( R^2 \) is the overall \( F \) statistic in disguise; the simple
correlation and Fisher's transformation; and the partial correlation, whose test is the \( t \) test of a
coefficient. The last two sections test whether a mean function is adequate: pure error gives an exact
test when rows are replicated, and several exact tests, none dominating the others, replace it when
they are not.

**What you need.** The normal linear model and its sampling distributions from
[Chapter 7](../ch07-optimality/index.html) (@thm-opt-sampling, @cor-opt-t, @cor-opt-quadratic),
nested projections and the Frisch–Waugh–Lovell theorem from
[Chapter 6](../ch06-projections/index.html) (@thm-proj-nested, @thm-proj-fwl), best linear prediction
in the population (@thm-rv-blp, @thm-proj-blp, @prp-proj-consistency) and the multiple correlation (@prp-rv-multiple-correlation) from [Chapter 2](../ch02-random-vectors/index.html) and [Chapter 6](../ch06-projections/index.html),
normal conditional distributions and partial correlations from
[Chapter 3](../ch03-multivariate-normal/index.html) (@thm-mvn-conditional, @def-mvn-partial-correlation), the noncentral distributions of
[Chapter 4](../ch04-quadratic-forms/index.html) (@thm-qf-orthogonal-projections, @thm-qf-f-power),
and \( R^2 \) and partial \( R^2 \) from [Chapter 9](../ch09-sums-of-squares/index.html) (@thm-ss-r2-null, @thm-ss-partial-r2). The \( F \) test of a reduced model
([Chapter 11](../ch11-general-linear-hypothesis/index.html), @thm-glh-f-test) is the template for
every lack-of-fit test. [Sections 14.1](01-random-regressors.html) and
[14.2](02-best-linear-prediction.html) also use the moments of least squares in the simple model
from [Chapter 5](../ch05-model-and-least-squares/index.html) (@cor-lm-simple-moments,
@thm-lm-sigma2), and, from [Chapter 1](../ch01-matrix-algebra/index.html), the Helmert matrix, the
partitioned inverse and a matrix square root (@exm-mat-helmert, @thm-mat-partitioned-inverse,
@thm-mat-square-root).

## Roadmap

- [Random regressors and conditional inference](01-random-regressors.html)
- [Best linear prediction](02-best-linear-prediction.html)
- [Multiple correlation](03-multiple-correlation.html)
- [The correlation coefficient and Fisher's z](04-correlation-coefficient.html)
- [Partial correlation](05-partial-correlation.html)
- [Testing lack of fit with replicates](06-lack-of-fit.html)
- [Near replicates and nonparametric lack-of-fit tests](07-near-replicates.html)
- [Summary and notes](08-summary-and-notes.html)
