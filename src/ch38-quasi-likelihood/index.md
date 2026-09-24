# Quasi-Likelihood and Overdispersion

A generalized linear model asks for a great deal and gives a great deal back. Naming the
family fixes the whole distribution of the response, and with it the exact relation between
its mean and its variance — the part of the bargain the data most often refuse. Counts of
accidents, litters of affected pups, votes in households, durations of strikes are nearly
always more variable than the Poisson or the binomial allows, and no choice of link repairs
it. This chapter takes the least drastic way out. The likelihood equations of
[Chapter 34](../ch34-exponential-families-glm/index.html) use the density only through the
mean and the variance function, so one may write down those two and estimate from the
equations alone. The estimate is unchanged, its consistency survives a wrong variance
function, and its standard errors can be read off the data rather than the model. What is
given up is the likelihood, and with it the likelihood ratio test, the Akaike criterion,
and every statement about a probability.

**What you need.** [Chapter 34](../ch34-exponential-families-glm/index.html) throughout:
the model (@def-glm-model), the variance function (@def-glm-variance), the score and the
information (@thm-glm-score), iteratively reweighted least squares (@thm-glm-irls), the
deviance (@def-glm-deviance, @thm-glm-deviance) and the
asymptotics (@thm-glm-asymptotics, @prp-glm-noncanonical).
[Section 21.4](../ch21-nonnormality-heteroscedasticity-serial/04-sandwich-estimators.html)
for the sandwich estimators (@def-het-hc, @thm-het-sandwich) and
[Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html)
for weighted least squares (@thm-het-wls).
[Chapter 37](../ch37-counts/index.html) for the negative binomial
models (@def-cnt-negbin, @prp-cnt-negbin) and the score test for
overdispersion (@thm-cnt-dispersion-score). The optimality proof of
[Section 38.2](02-estimating-equations.html) uses the projection inequality of
[Chapter 6](../ch06-projections/index.html);
[Chapter 20](../ch20-residuals-leverage-influence/index.html),
[Chapter 23](../ch23-resampling-inference/index.html),
[Chapter 29](../ch29-model-selection/index.html) and
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) are pointed at in single
sections, as context rather than as steps.

## Roadmap

- [Variance functions](01-variance-functions.html)
- [Quasi-likelihood estimating equations](02-estimating-equations.html)
- [Sandwich covariance estimation](03-sandwich.html)
- [Overdispersion](04-overdispersion.html)
- [Summary and notes](05-summary-and-notes.html)
