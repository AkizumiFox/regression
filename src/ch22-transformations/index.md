# Transformations

A regression is fitted on a scale, and nothing in least squares says which one. Food expenditure can be
modelled in currency units or in logarithms, a count as it stands or through its square root. The choice
affects all three claims of the normal linear model: a linear mean, a constant variance and normal
errors. When the three failures share a cause, one transformation repairs them together; when they do
not, the chapter says which to fix first and at what price. It treats the choice of scale as estimation:
the Box–Cox profile likelihood for the response, the Box–Tidwell procedure and component-plus-residual
plots for regressors, and the delta method for transformations that stabilize a variance tied to the
mean. It ends with what a transformed model says on the scale the reader cares about: medians, means,
prediction intervals and elasticities.

**What you need.** Least squares as projection and the Frisch–Waugh–Lovell theorem
([Chapter 6](../ch06-projections/index.html), @thm-proj-fwl), maximum likelihood and profile likelihood in the
normal linear model ([Chapter 7](../ch07-optimality/index.html), @thm-opt-mle and @def-opt-profile), the
interpretation of log-scale coefficients ([Section 5.7](../ch05-model-and-least-squares/07-interpreting-coefficients.html)),
prediction intervals (@thm-ci-prediction-interval), the delta method as used in
[Section 14.4](../ch14-correlation-lack-of-fit-prediction/04-correlation-coefficient.html), and the
diagnostics of [Chapter 20](../ch20-residuals-leverage-influence/index.html) and
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html).
[Section 22.5](05-interpreting.html) also uses the moment generating function of the multivariate
normal from [Chapter 3](../ch03-multivariate-normal/index.html) (@thm-mvn-mgf), which is what says
what a fitted value on the log scale estimates.

## Roadmap

- [Transforming the response](01-transforming-the-response.html)
- [The Box–Cox family](02-box-cox.html)
- [Transforming regressors](03-transforming-regressors.html)
- [Variance-stabilizing transformations](04-variance-stabilizing.html)
- [Interpreting a transformed model](05-interpreting.html)
- [Summary and notes](06-summary-and-notes.html)
