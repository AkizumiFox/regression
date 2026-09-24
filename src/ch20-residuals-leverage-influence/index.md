# Residuals, Leverage and Influence

Every result of Parts II–IV was proved for a model that holds for all \( n \) cases. A record is mistyped, one unit is unlike the others, or the model is right
for most of the range and wrong at its edge. Least squares treats every case as equally credible, so a
single bad case can move the fit a long way, and it can do so without leaving a large residual behind.
This chapter builds the tools for finding such cases. Residuals, rescaled so that they can be compared,
show cases that the model fits badly. Leverage shows cases whose position among the regressors gives
them a strong pull on the fit. Deleting one case at a time, cheap by the updating
formulas of [Chapter 10](../ch10-computation/index.html), turns these two ingredients into
exact tests and into measures of influence. Two warnings run through the chapter. Diagnostics that look
at one case at a time can be defeated by a few bad cases acting together. And a flagged case is a
question to be investigated, not a verdict. The last section uses a robust fit as a second opinion:
a fit that the bad cases cannot capture shows them by their residuals.

**What you need.** [Chapter 4](../ch04-quadratic-forms/index.html) (independent projections),
[Chapter 5](../ch05-model-and-least-squares/index.html) (moments of fitted values and
residuals, @prp-lm-fit-moments), [Chapter 6](../ch06-projections/index.html) (the hat matrix and its
properties, @prp-proj-leverage and @prp-proj-leverage-mahalanobis; the Frisch–Waugh–Lovell theorem, @thm-proj-fwl),
[Chapter 10](../ch10-computation/index.html) (the deletion formulas, @prp-cmp-loo),
[Chapter 11](../ch11-general-linear-hypothesis/index.html) (\( t \) and \( F \) tests), and
[Chapter 13](../ch13-multiplicity/index.html) (the Bonferroni inequality). [Chapter 19](../ch19-theory-of-departures/index.html) describes in general
terms what a contaminated case does to least squares (@thm-dep-outlier); this chapter is about finding
such cases in a particular data set.
From [Chapter 1](../ch01-matrix-algebra/index.html), the deletion identity that every
leave-one-out formula is an instance of (@exm-mat-deletion), with the extremal Rayleigh quotient,
a matrix square root and the Cauchy–Schwarz inequality (@thm-mat-extremal-rayleigh,
@thm-mat-square-root, @prp-mat-cauchy-schwarz); from
[Chapter 2](../ch02-random-vectors/index.html), the covariance of a linear transformation
(@thm-rv-linear); and, in [Section 20.3](03-deletion.html), the single-dummy exercise of
[Chapter 8](../ch08-estimability/index.html) (@exr-est-single-dummy), which is what makes deleting
a case the same as adding one indicator.

## Roadmap

- [Kinds of residuals](01-kinds-of-residuals.html)
- [Leverage and the hat matrix](02-leverage.html)
- [Deletion diagnostics](03-deletion.html)
- [Cook's distance, DFFITS and DFBETAS](04-influence.html)
- [Masking and swamping](05-masking.html)
- [Robust fitting as a diagnostic](06-robust.html)
- [Summary and notes](07-summary-and-notes.html)
