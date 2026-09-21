# Smoothing: Kernels, Splines and Penalties

[Chapter 42](../ch42-polynomials-piecewise/index.html) reached a natural stopping
point. Piecewise polynomials can follow almost any shape, but only if someone
chooses how many pieces there are and where they join, and that choice is both
awkward and consequential. This chapter removes the choice. It replaces the
question "how many knots?" with a question about a single continuous dial: how
much roughness are we willing to pay for a closer fit? Everything else — the
basis, the number of knots, the algorithm — becomes a matter of convenience once
that dial exists.

Two traditions answer the question, and they meet in the middle. One fits a
simple model in a moving window — a constant, a line, a cubic, weighted by
distance from the point being estimated. Its dial is the width of the window,
and its theory is asymptotic, a bias–variance calculation that ends in the rate
\( n^{-4/5} \). The other keeps a single global fit but adds a penalty on
roughness; its dial is the weight on the penalty, and its theory is the algebra
of [Chapter 27](../ch27-shrinkage/index.html) and
[Chapter 30](../ch30-regularization-boosting/index.html): a penalized least squares
problem, a linear operator, a trace that counts the degrees of freedom spent.
The second tradition is the more useful one here, because a penalty carries over
unchanged to several covariates
([Chapter 44](../ch44-additive-models/index.html)), to non-normal responses and to
parameters other than the mean ([Chapter 45](../ch45-quantile-gamlss/index.html)),
and because it has a second reading as a normal prior, which turns the smoothing
parameter into a variance ratio the data can estimate.

**What you need.** [Chapter 6](../ch06-projections/index.html) (least squares as
projection, @prp-proj-trace-rank, and the leverage @prp-proj-leverage that
returns in [Section 43.5](05-choosing-lambda.html)),
[Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html)
(weighted least squares, @thm-het-wls),
[Chapter 27](../ch27-shrinkage/index.html) (ridge regression @thm-shr-ridge) and
[Section 30.1](../ch30-regularization-boosting/01-penalized-least-squares.html)
(penalized least squares, @def-reg-penalized). From
[Chapter 29](../ch29-model-selection/index.html): the leave-one-out
shortcut @thm-sel-loocv, generalized cross-validation @def-sel-gcv, AIC @prp-sel-aic-bic and
the warning @prp-sel-selection-bias. From
[Chapter 32](../ch32-linear-mixed-models/index.html): the linear mixed
model @def-mix-model, Henderson's equations @thm-mix-henderson, the BLUP @thm-mix-blup
and REML (@def-mix-reml, @thm-mix-reml). [Section 1.8](../ch01-matrix-algebra/08-svd.html)
supplies the condition number @def-mat-condition-number, and
[Chapter 42](../ch42-polynomials-piecewise/index.html) the spline space
itself (@def-ply-piecewise, @thm-ply-spline-space, @prp-ply-knots). [Section 43.6](06-mixed-model-and-bayes.html)
draws on
[Section 39.4](../ch39-glms-in-practice-bayes/04-bayesian-glms-mcmc.html)
for the Bayesian vocabulary (@def-prc-posterior, @prp-prc-priors).

**Notation for this chapter.** There is one covariate \( x \) taking values
\( x_1,\dots,x_n \) in an interval \( [a,b] \), one response, and an unknown
function \( f \) with \( \mathbf{f}=\{f(x_1),\dots,f(x_n)\}\T \). [Section 43.1](01-kernels.html) writes
\( K \) for a kernel function and \( h \) for a bandwidth; from [Section 43.2](02-b-splines.html) on,
\( K \) is instead the *number of interior knots*, written
\( t_1<\dots<t_K \), with \( \kappa_j \) reserved for the extended knot sequence
of the B-spline recurrence; the kernel does not reappear except by name.
\( \B \) is a matrix of basis functions, \( \bgamma \) its coefficients, \( \bP \)
a penalty matrix — a local use, not the projection \( \bP \) of Chapters 4 to 33 —
\( \bD \) a difference matrix and \( \bS_\lambda \) a smoother matrix. In [Section 43.4](04-smoothing-splines.html), \( \Q \) and \( \R \) are
the two band matrices of the Reinsch form, a local use unrelated to the \( \R \)
of Part VII.

## Roadmap

- [Kernel and local polynomial regression](01-kernels.html)
- [B-splines](02-b-splines.html)
- [Penalized splines](03-penalized-splines.html)
- [Smoothing splines](04-smoothing-splines.html)
- [Choosing the smoothing parameter](05-choosing-lambda.html)
- [Mixed-model and Bayesian representations](06-mixed-model-and-bayes.html)
- [Summary and notes](07-summary-and-notes.html)
