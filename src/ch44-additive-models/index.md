# Additive, Geoadditive and Structured Additive Models

[Chapter 43](../ch43-smoothing/index.html) learned to estimate one unknown curve. Real data
sets have a dozen covariates, and one unknown surface in a dozen dimensions is hopeless: no
neighbourhood of a point holds enough observations to average. The way out is a structural
assumption weaker than linearity but far stronger than nothing — the covariates act
**additively**, each through its own unknown function — which turns a hopeless problem into
a stack of one-dimensional ones held together by a single penalized criterion.

This chapter builds that framework and keeps widening it. A coefficient may vary smoothly
with a second covariate. A term may be a spatial effect, penalized towards its neighbours on
a map or towards a smooth surface over coordinates. A term may be a random intercept, a
penalized term whose penalty is the identity. All of these are *one* model — a structured
additive predictor — with *one* estimation theory: a quadratic penalty per term, effective
degrees of freedom from a trace, smoothing parameters from a restricted likelihood, and a
Bayesian reading in which every penalty is a prior. The last section puts a link and a
response distribution in front of that predictor and closes the circle with
[Part VIII](../ch34-exponential-families-glm/index.html).

**What you need.** [Chapter 6](../ch06-projections/index.html) (projections, the trace as a
rank, @thm-proj-fwl), [Chapter 26](../ch26-collinearity/index.html) (collinearity and the
variance inflation factor, @def-col-vif), [Chapter 27](../ch27-shrinkage/index.html) and
[Chapter 30](../ch30-regularization-boosting/index.html) (ridge, penalized least squares and
effective degrees of freedom, @thm-shr-ridge and @def-reg-penalized),
[Chapter 29](../ch29-model-selection/index.html) (cross-validation and information criteria,
@thm-sel-loocv, @def-sel-gcv, @prp-sel-selection-bias),
[Chapter 32](../ch32-linear-mixed-models/index.html) (mixed models, BLUP and REML,
@def-mix-model, @thm-mix-blup, @thm-mix-reml), [Chapter 34](../ch34-exponential-families-glm/index.html)
(the generalized linear model and IRLS, @def-glm-model, @thm-glm-irls),
[Chapter 38](../ch38-quasi-likelihood/index.html) (quasi-likelihood and dispersion,
@def-ql-quasi, @prp-ql-dispersion), [Chapter 39](../ch39-glms-in-practice-bayes/index.html)
(the Bayesian machinery, @def-prc-posterior, @prp-prc-priors) and, above all,
[Chapter 43](../ch43-smoothing/index.html), whose P-spline blocks and smoothing-parameter
criteria are used here as given.

::: {.remark}
[Notation in this chapter]

A **term** is a pair: a design block \( \Z_j \) of size \( n\times d_j \) and a
nonnegative definite penalty \( \mathbf{K}_j \), with coefficients \( \bgamma_j \) and
smoothing parameter \( \lambda_j \). So \( \Z \) keeps its meaning from
[Part VII](../ch31-general-gauss-markov/index.html), the design of the coefficients that are
shrunk, and [Section 44.5](05-structured-additive.html) shows the two uses are the same.
Unpenalized columns stay in \( \X \) with coefficients \( \bbeta \);
\( \A=\Z\T\W\Z+\mathbf{K} \) is the penalized cross-product, \( \W \) the working
weights of [Chapter 34](../ch34-exponential-families-glm/index.html), \( \bS \) the
smoother matrix and \( \mathbf{C}_j \) the operator producing the \( j \)th
component.
:::

## Roadmap

- [The additive model](01-additive-models.html)
- [Backfitting and concurvity](02-backfitting.html)
- [Varying-coefficient models](03-varying-coefficients.html)
- [Spatial effects](04-spatial-effects.html)
- [Structured additive regression](05-structured-additive.html)
- [Generalized additive models](06-generalized-additive-models.html)
- [Summary and notes](07-summary-and-notes.html)
