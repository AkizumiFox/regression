# Binary Responses

A response that can only be zero or one is the smallest departure from the normal
linear model that still changes everything. Its mean is a probability, so it lives
in an interval rather than on a line; its variance is determined by its mean, so
there is no free dispersion parameter; and the exact distribution theory of Parts
III and IV is gone. What remains is the likelihood, and the apparatus of
[Chapter 34](../ch34-exponential-families-glm/index.html) applies word for word.
This chapter is what that apparatus looks like pointed at the simplest possible
response, and what that response does to it: odds ratios that survive sampling
schemes no other model can tolerate, a maximum likelihood estimate that sometimes
fails to exist, and a deviance that, for ungrouped data, carries no information
about fit at all.

**What you need.** [Chapter 34](../ch34-exponential-families-glm/index.html) for the
exponential dispersion family, the likelihood equations, iteratively reweighted
least squares, the deviance and the asymptotics; [Chapter 7](../ch07-optimality/index.html)
for maximum likelihood and Fisher information (@thm-opt-mle, @prp-opt-information);
[Chapter 11](../ch11-general-linear-hypothesis/index.html) for the likelihood ratio
test in the normal model (@thm-glh-lrt, @prp-glh-trinity);
[Chapter 20](../ch20-residuals-leverage-influence/index.html) for residuals, leverage
and influence; [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html)
for weighted least squares; [Chapter 22](../ch22-transformations/index.html) for
empirical logits (@exr-tr-empirical-logit);
[Chapter 25](../ch25-causal-interpretation/index.html) for the causal reading of a
coefficient; [Chapter 27](../ch27-shrinkage/index.html) and
[Chapter 30](../ch30-regularization-boosting/index.html) for penalized estimation;
[Chapter 29](../ch29-model-selection/index.html) for AIC and cross-validation. Only
[Chapter 34](../ch34-exponential-families-glm/index.html) is a prerequisite in the strict sense:
every proof in this chapter that cites an earlier chapter cites that one. The rest are named
because the chapter points at them — for a diagnostic, a penalty, a criterion or a reading of a
coefficient — and a reader may follow the pointer when it arises.

## Roadmap

- [Logistic regression](01-logistic-regression.html)
- [Probit and complementary log–log links](02-probit-cloglog.html)
- [Separation](03-separation.html)
- [Goodness of fit](04-goodness-of-fit.html)
- [Case–control sampling](05-case-control.html)
- [Summary and notes](06-summary-and-notes.html)
