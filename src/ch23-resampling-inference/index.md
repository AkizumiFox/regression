# Bootstrap and Permutation Inference

Every test and interval in Part III took its reference distribution from the normal model. This chapter asks whether the
reference distribution can be taken from the data instead. The bootstrap imitates the process that produced the data by
resampling residuals, whole cases, or residuals with random signs, and each choice assumes a different model. A permutation test rearranges the data in ways
that leave their distribution unchanged under the null hypothesis, and is exact whenever such rearrangements exist. In
regression they exist only for simple hypotheses, and the Freedman–Lane scheme shows how to get close in other cases.

**What you need.** [Chapter 6](../ch06-projections/index.html) (projections, leverage, @prp-proj-leverage,
and the Frisch–Waugh–Lovell theorem, @thm-proj-fwl), [Chapter 7](../ch07-optimality/index.html) (sampling
distributions in the normal model, @thm-opt-sampling and @cor-opt-t),
[Chapter 10](../ch10-computation/index.html) (the deletion formulas, @prp-cmp-loo),
[Chapter 11](../ch11-general-linear-hypothesis/index.html) (the \( F \) test of a linear hypothesis, @thm-glh-f-test),
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) (random regressors and conditional inference) and
[Chapter 18](../ch18-covariance-and-design/index.html) (randomized blocks). From this part:
[Chapter 19](../ch19-theory-of-departures/index.html) (asymptotic normality without normal errors),
[Chapter 20](../ch20-residuals-leverage-influence/index.html) (studentized residuals) and
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (sandwich covariance estimators and serial
correlation). The moments of least squares and the consistency of \( s^2 \) from
[Chapter 5](../ch05-model-and-least-squares/index.html) (@thm-lm-moments, @cor-lm-s2-consistent) are
what [Sections 23.1](01-cases-and-residuals.html) and [23.5](05-freedman-lane.html) compare a
resampling distribution against.

## Roadmap

- [Resampling cases versus resampling residuals](01-cases-and-residuals.html)
- [The wild bootstrap](02-wild-bootstrap.html)
- [Bootstrap confidence intervals](03-confidence-intervals.html)
- [Permutation tests for linear hypotheses](04-permutation-tests.html)
- [Freedman–Lane and related schemes](05-freedman-lane.html)
- [When resampling fails](06-failures.html)
- [Summary and notes](07-summary-and-notes.html)
