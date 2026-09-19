# Estimability, Identifiability and Indicator Variables

When the columns of the model matrix are linearly dependent, the data cannot
tell some coefficient vectors apart, because different vectors produce the same mean. This
chapter asks which questions about \( \bbeta \) still have answers. The answer is
the same from every direction. A linear function \( \blambda\T\bbeta \) can be estimated
without bias exactly when it is determined by the mean vector \( \X\bbeta \), exactly when
\( \blambda \) lies in the row space of \( \X \), and exactly when every least squares
solution gives it the same value. The rest of the chapter uses this criterion.
Reparameterizations and side conditions turn out to be two ways of choosing coordinates,
and neither can change what is estimable. Factors enter a regression through indicator
variables, and each coding of a factor amounts to a choice of which linear functions of the
cell means the coefficients report. Contrasts are the estimable functions of a factor, and
orthogonal sets of contrasts split the factor's sum of squares into separate pieces.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html) (rank, generalized inverses and
the consistency of linear systems, especially @thm-mat-consistency and @cor-mat-invariant),
[Chapter 2](../ch02-random-vectors/index.html) (means and covariances of linear functions,
@thm-rv-linear), [Chapter 4](../ch04-quadratic-forms/index.html) (for the distribution of contrast sums of
squares under normality, @thm-qf-orthogonal-projections), [Chapter 5](../ch05-model-and-least-squares/index.html) (the linear model,
@def-lm-linear-model), and [Chapter 6](../ch06-projections/index.html), above all
[Section 6.2](../ch06-projections/02-subspaces.html) and
[Section 6.4](../ch06-projections/04-least-squares.html), where rank-deficient least squares
and the invariance of \( \blambda\T\hbeta \) for \( \blambda\in\C(\X\T) \) (@thm-proj-invariant-functions)
first appeared. A few remarks use the Gauss–Markov
theorem of [Chapter 7](../ch07-optimality/index.html).

## Roadmap

- [Identifiability](01-identifiability.html)
- [Estimable functions](02-estimable-functions.html)
- [Characterizations of estimability](03-characterizations.html)
- [Reparameterization and side conditions](04-side-conditions.html)
- [Indicator variables and the coding of a factor](05-factor-coding.html)
- [Several factors, interactions and covariates](06-several-factors.html)
- [Contrasts](07-contrasts.html)
- [Summary and notes](08-summary-and-notes.html)
