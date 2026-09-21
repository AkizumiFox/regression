# GLMMs and Generalized Estimating Equations

Part VII replaced independent errors by a covariance matrix with structure; Part VIII has
so far replaced the normal response by an exponential dispersion family. This chapter does
both at once, and the combination is not the sum of its parts. In a linear model the mean
is unaffected by what happens to the covariance, so a random effect changes the standard
errors and leaves the coefficients alone. Once the link is nonlinear that is false:
averaging a curve over a random intercept does not give the same curve back, so the model
has *two* sets of coefficients — one for a single cluster, one for the population — and
they are not equal.

Two families of methods follow from that fork. If the interest is in what happens within a
cluster, and a full probability model is acceptable, the answer is the generalized linear
mixed model: random effects inside the linear predictor and a likelihood obtained by
integrating them out, numerically, with different approximations giving visibly different
answers. If the interest is in the population average, and no model for the dependence is
wanted, the answer is generalized estimating equations: solve a quasi-score built from a
*working* correlation that is allowed to be wrong, and price the estimate with a sandwich.
The second route buys robustness by giving up the likelihood, and with it the deviance, the
random-effect predictions, and any tolerance for data that go missing for the wrong
reasons.

**What you need.** [Chapter 34](../ch34-exponential-families-glm/index.html) throughout
(the generalized linear model @def-glm-model and @def-glm-link, the likelihood equations
@thm-glm-score, iteratively reweighted least squares @thm-glm-irls, the
asymptotics @thm-glm-asymptotics), [Chapter 35](../ch35-binary-responses/index.html)
(the logistic model @def-bin-logistic, non-collapsibility @exm-bin-noncollapsible
and case-control sampling @thm-bin-casecontrol),
[Chapter 36](../ch36-multinomial-ordinal/index.html) and
[Chapter 37](../ch37-counts/index.html) (Poisson regression @def-cnt-poisson and the
negative binomial @def-cnt-negbin). From Part VII: the linear mixed model @def-mix-model,
best linear unbiased prediction @thm-mix-blup and its shrinkage reading, REML
@def-mix-reml, the boundary problem @thm-mix-boundary, the covariance models
@def-cls-covariance and the working-covariance result @prp-cls-selection. From Part V:
the sandwich covariance @thm-het-sandwich.
[Chapter 29](../ch29-model-selection/index.html) (AIC and BIC, @def-sel-aic-bic)
and [Chapter 25](../ch25-causal-interpretation/index.html) (what a
coefficient means) are used in places. [Chapter 38](../ch38-quasi-likelihood/index.html)
supplies the quasi-likelihood machinery this chapter generalizes, while the
missingness conditions that [Section 40.4](04-gee.html) needs are stated in
[Chapter 41](../ch41-missing-data/index.html), which comes after this one.

## Roadmap

- [Marginal versus conditional models](01-marginal-versus-conditional.html)
- [Generalized linear mixed models](02-glmms.html)
- [Integrated likelihood approximations](03-integrated-likelihood.html)
- [Generalized estimating equations](04-gee.html)
- [Working correlation](05-working-correlation.html)
- [Summary and notes](06-summary-and-notes.html)
