# Collinearity: Diagnosis and Consequences

Part VI gives up unbiasedness in exchange for smaller variance, and collinearity is the
situation that first made that trade attractive. When the columns of the model matrix are
nearly linearly dependent, the data answer some questions well and others barely at all, and
the badly answered ones are often those the analyst asked: the separate effects of regressors that
moved together. This chapter says exactly which linear functions lose precision, builds the standard
diagnostics (variance inflation factors and their generalization to groups of columns, condition
indices and the Belsley–Kuh–Welsch variance decomposition), shows why predictions that follow the
pattern of the data survive intact, and ends with the remedies, the last of which, biased estimation,
is where [Chapter 27](../ch27-shrinkage/index.html) begins.

**What you need.** [Chapter 1](../ch01-matrix-algebra/index.html) (the spectral theorem and
Rayleigh quotients, @thm-mat-spectral and @thm-mat-extremal-rayleigh; the singular value
decomposition and the condition number, @thm-mat-svd and @def-mat-condition-number),
[Chapter 6](../ch06-projections/index.html) (the Frisch–Waugh–Lovell theorem and the variance
identity @eq-proj-vif-preview; leverage, @prp-proj-leverage),
[Chapter 8](../ch08-estimability/index.html) (estimability),
[Chapter 10](../ch10-computation/index.html) (conditioning and column scaling,
@thm-cmp-perturbation and @prp-cmp-van-der-sluis) and
[Chapter 12](../ch12-intervals-and-bands/index.html) (confidence ellipsoids and prediction at a
new point). [Section 19.6](../ch19-theory-of-departures/06-collinearity.html) previewed the subject (@prp-dep-collinear-directions);
this chapter builds on that preview rather than repeating it.

## Roadmap

- [What collinearity does to variance](01-variance.html)
- [Variance inflation factors](02-vif.html)
- [Condition indices and variance decomposition](03-condition-indices.html)
- [Collinearity and prediction](04-prediction.html)
- [What to do about collinearity](05-remedies.html)
- [Summary and notes](06-summary-and-notes.html)
