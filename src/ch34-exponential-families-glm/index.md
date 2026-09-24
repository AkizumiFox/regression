# Exponential Families and the Generalized Linear Model

Every chapter so far has modelled a response whose scatter about a linear mean is, if not
normal, at least additive and of constant shape. Counts, proportions, waiting times and
insurance claims are none of these: their means are confined to an interval, their spread
grows with their mean, their distributions are skewed by construction. This chapter builds
the framework that handles them — a response from an **exponential dispersion family**, a
**linear predictor** \( \x_{(i)}\T\bbeta \), and a **link** joining the two. Least squares
gives way to maximum likelihood, the normal equations to the likelihood equations, and the
closed-form estimate to an iteration that turns out to be weighted least squares, run until
it stops changing. Exact distribution theory gives way to asymptotics, so every claim here
comes with conditions attached.

**What you need.** [Chapter 5](../ch05-model-and-least-squares/index.html) and
[Chapter 6](../ch06-projections/index.html) (the linear model, projection, leverage),
[Section 7.3](../ch07-optimality/03-maximum-likelihood.html) and
[Section 7.4](../ch07-optimality/04-minimum-variance.html) (maximum
likelihood @thm-opt-mle, the score and information @prp-opt-information, exponential
families @lem-opt-exponential-family), [Chapter 11](../ch11-general-linear-hypothesis/index.html)
(the likelihood ratio test @thm-glh-lrt and the trinity @prp-glh-trinity),
[Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html)
(weighted least squares @thm-het-wls). [Chapter 20](../ch20-residuals-leverage-influence/index.html)
(residuals and influence, and iteratively reweighted least squares @prp-res-irls) and
[Chapter 29](../ch29-model-selection/index.html) (AIC) are used in the last section, which also
needs the chi-squared distribution of a quadratic form from
[Chapter 4](../ch04-quadratic-forms/index.html) (@thm-qf-chisq) for the deviance in the normal case.

## Roadmap

- [Exponential dispersion families](01-exponential-dispersion-families.html)
- [Link functions](02-link-functions.html)
- [The likelihood equations](03-likelihood-equations.html)
- [Iteratively reweighted least squares](04-irls.html)
- [Deviance and asymptotic inference](05-deviance-and-inference.html)
- [Summary and notes](06-summary-and-notes.html)
