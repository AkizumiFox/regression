# Errors in Variables

Every earlier chapter takes the regressors as recorded without error. Often they are not: one clinic reading stands
in for long-run blood pressure, a questionnaire for usual diet. The scientific model concerns the true value \( X \),
but the data contain a noisy \( W \), and least squares on \( W \) estimates the best linear predictor from what was
recorded. This chapter measures the gap and closes it. Classical error flattens a slope by a factor no sample size
removes, makes correctly measured covariates look important, and in the normal model leaves the slope unidentified.
Information from outside the regression rescues the analysis: a known error variance, repeated measurements, a known
variance ratio, an instrument, or validation data. The moment correction, instrumental variables, regression calibration
and SIMEX turn that information into estimates, each with its own conditions and costs.

**What you need.** [Chapter 6](../ch06-projections/index.html) (the population projection, @thm-proj-blp and
@prp-proj-consistency, and @thm-proj-M-formula). [Chapter 5](../ch05-model-and-least-squares/index.html) (moments of least squares
under departures, @prp-lm-misspecified). [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html)
(random regressors, @thm-cor-conditional, and the limit theorems of
[Section 14.4](../ch14-correlation-lack-of-fit-prediction/04-correlation-coefficient.html)). From this part,
[Chapter 19](../ch19-theory-of-departures/index.html) (omitted-variable bias, @thm-dep-omitted),
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (sandwich covariance estimators, @thm-het-sandwich)
and [Chapter 23](../ch23-resampling-inference/index.html) (the pairs bootstrap, @def-bs-bootstrap).

## Roadmap

- [Attenuation bias](01-attenuation.html)
- [The classical measurement error model](02-measurement-error-model.html)
- [Reliability and moment correction](03-moment-correction.html)
- [Instrumental variables](04-instrumental-variables.html)
- [Regression calibration and SIMEX](05-calibration-simex.html)
- [Summary and notes](06-summary-and-notes.html)
