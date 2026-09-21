# GLMs in Practice; the Bayesian Thread

Chapters 34 to 38 built the generalized linear model and its inference. This chapter asks
two practical questions the theory leaves open. The first is how one arrives at a
particular model at all: which family, which link, which scale for the covariates, which
interactions, and how candidates are compared. The second is how the model is criticized
once fitted, since it has three separable pieces — a mean, a link and a variance function —
and a plot that condemns the whole fit is much less useful than one naming the piece at
fault.

The rest of the chapter takes up the Bayesian thread, running through the book since
[Section 7.6](../ch07-optimality/06-bayes-conjugate.html). For the normal linear model the
posterior was available in closed form and the flat prior reproduced least squares
exactly. None of that survives the move to a generalized linear model, for one
identifiable reason: the log-likelihood is no longer a quadratic. What replaces it is
approximation — the Laplace expansion, which also explains where the BIC comes from — and
simulation, in the form of Metropolis–Hastings and Gibbs samplers written out in numpy.
The chapter ends where the Bayesian and penalized-likelihood threads meet: a normal prior
is ridge regression, a Laplace prior the lasso at its mode, and a merely weakly
informative prior is enough to make an estimate exist where maximum likelihood has none.

**What you need.** [Chapter 34](../ch34-exponential-families-glm/index.html)
throughout (@def-glm-model, @def-glm-link, @def-glm-variance, @thm-glm-score, @thm-glm-irls,
@def-glm-deviance, @thm-glm-deviance, @thm-glm-asymptotics, @def-glm-residuals);
[Chapter 35](../ch35-binary-responses/index.html) for separation (@thm-bin-separation),
Firth's penalty (@prp-bin-firth) and the deviance for binary data (@prp-bin-fit);
[Chapter 37](../ch37-counts/index.html) for offsets (@prp-cnt-offsets);
[Chapter 38](../ch38-quasi-likelihood/index.html) for dispersion (@prp-ql-dispersion) and
overdispersion (@prp-ql-overdispersion). From earlier parts:
[Chapter 20](../ch20-residuals-leverage-influence/index.html) (@def-res-cooks,
@thm-res-deletion, @prp-res-partial-residual);
[Chapter 22](../ch22-transformations/index.html) for
@prp-tr-retransformation and @thm-tr-box-tidwell; the parametric bootstrap of
@prp-bs-parametric; [Chapter 29](../ch29-model-selection/index.html) (@def-sel-aic-bic, @prp-sel-aic-bic,
@def-sel-cv, @prp-sel-selection-bias); and penalties as priors
([Chapter 27](../ch27-shrinkage/index.html),
[Chapter 30](../ch30-regularization-boosting/index.html): @def-shr-ridge,
@def-shr-lasso, @prp-reg-map). The Bayesian thread so far:
[Section 7.6](../ch07-optimality/06-bayes-conjugate.html) for
@thm-opt-bayes-conjugate and @cor-opt-flat-prior,
[Section 12.6](../ch12-intervals-and-bands/06-bayes.html) for @def-ci-credible, and
[Section 32.6](../ch32-linear-mixed-models/06-bayes.html) for @prp-mix-bayes.

## Roadmap

- [Building a model](01-model-building.html)
- [Diagnostics](02-diagnostics.html)
- [Bayesian thread IV: what carries over](03-bayesian-linear-revisited.html)
- [Bayesian generalized linear models and MCMC](04-bayesian-glms-mcmc.html)
- [Priors as regularizers](05-priors-as-regularizers.html)
- [Summary and notes](06-summary-and-notes.html)
