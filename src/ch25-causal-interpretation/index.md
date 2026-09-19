# Causal Interpretation of Coefficients

Earlier chapters treated a regression coefficient as a feature of a joint distribution. People who fit regressions
usually want more: what would happen to the response if a regressor were *changed*, if a patient took the drug or a
school hired more teachers. No amount of data on the undisturbed system answers that without further assumptions. This
chapter states the assumptions and proves what they buy: potential outcomes and randomization, graphs and the back-door
criterion for observational data, the bias from adjusting for mediators and colliders, and covariate adjustment in
experiments, which can hurt unless the regression is interacted. This makes precise the
warnings of [Section 5.7](../ch05-model-and-least-squares/07-interpreting-coefficients.html), of the remark after the
Frisch–Waugh–Lovell theorem ([Section 6.6](../ch06-projections/06-fwl.html)) and of
[Section 18.1](../ch18-covariance-and-design/01-ancova.html).

**What you need.** [Chapter 5](../ch05-model-and-least-squares/index.html) (what a coefficient
measures, @def-lm-coefficient-interpretation, and the omitted-variable formula, @prp-lm-omitted), [Chapter 6](../ch06-projections/index.html) (the Frisch–Waugh–Lovell theorem,
@thm-proj-fwl, and projection in the population, @thm-proj-blp and @prp-proj-consistency),
[Chapter 3](../ch03-multivariate-normal/index.html) and [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html)
(partial correlation, @prp-mvn-partial-meaning and @prp-cor-partial-coefficient),
[Chapter 18](../ch18-covariance-and-design/index.html) (analysis of covariance and randomization,
@thm-dsn-ancova, @prp-dsn-precision and @lem-dsn-srs), and, from this part,
[Chapter 19](../ch19-theory-of-departures/index.html) (omitted-variable bias of least squares, @thm-dep-omitted)
and [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (sandwich standard errors, @thm-het-sandwich).

## Roadmap

- [Prediction versus intervention](01-prediction-intervention.html)
- [Potential outcomes and the regression coefficient](02-potential-outcomes.html)
- [Confounding and causal graphs](03-confounding-dags.html)
- [Adjustment sets and the back-door criterion](04-backdoor.html)
- [Bad controls: mediators and colliders](05-bad-controls.html)
- [Regression adjustment in randomized experiments](06-randomized-adjustment.html)
- [Summary and notes](07-summary-and-notes.html)
