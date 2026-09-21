# Quantile Regression and GAMLSS

Forty-four chapters have been spent on one number. Least squares estimates a
conditional mean; the generalized linear model estimates a transformed conditional
mean; splines and additive models let that mean bend. Everything else about the
conditional distribution of the response — its spread, its skewness, its tails —
entered as an assumption to be checked, a nuisance to be corrected for, or a reason
to widen a standard error. This last chapter turns the nuisance into the subject.

There are two ways to do it. One drops the parametric model and estimates the
conditional quantiles directly, each level fitted by its own criterion: that is
quantile regression, whose loss function is as simple as squared error and as easy to
justify. The other keeps a parametric family but lets *every* parameter of it carry
its own structured additive predictor: that is GAMLSS. Between them sits expectile
regression, which buys a differentiable loss at the price of a less transparent
interpretation. The chapter ends where a book about models should end: with how to
tell whether one estimated distribution is better than another.

**What you need.** [Chapter 5](../ch05-model-and-least-squares/index.html) and
[Chapter 6](../ch06-projections/index.html) (the linear model, least squares and the
best linear predictor @thm-proj-blp),
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html)
(heteroscedasticity, weighted least squares and the sandwich @def-het-hc),
[Chapter 22](../ch22-transformations/index.html) (@def-tr-box-cox, the older route to
a well-behaved response), [Chapter 23](../ch23-resampling-inference/index.html) (the
bootstrap @def-bs-bootstrap and its schemes @def-bs-schemes),
[Chapter 25](../ch25-causal-interpretation/index.html) (what a fitted curve does and
does not mean), [Chapter 29](../ch29-model-selection/index.html)
(cross-validation @thm-sel-loocv, AIC @prp-sel-aic-bic, selection bias @prp-sel-selection-bias),
[Chapter 30](../ch30-regularization-boosting/index.html) (penalized
fitting @def-reg-penalized) and [Chapter 34](../ch34-exponential-families-glm/index.html)
(the generalized linear model @def-glm-model and its iteratively reweighted least
squares @thm-glm-irls). The last two sections use the penalized splines of
[Chapter 43](../ch43-smoothing/index.html) and the structured additive predictors of
[Chapter 44](../ch44-additive-models/index.html); the M-estimation machinery of
[Section 20.6](../ch20-residuals-leverage-influence/06-robust.html) (@prp-res-irls)
is the ancestor of two algorithms below.

## Roadmap

- [Why model more than the mean](01-beyond-the-mean.html)
- [Quantile regression](02-quantile-regression.html)
- [Inference for quantile regression](03-inference.html)
- [Expectile regression](04-expectiles.html)
- [Models for location, scale and shape](05-gamlss.html)
- [Comparing distributional models](06-comparing.html)
- [Summary and notes](07-summary-and-notes.html)
