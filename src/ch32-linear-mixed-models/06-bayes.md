# Bayesian thread III: mixed models as hierarchical models

The Bayesian thread began in
[Section 7.6](../ch07-optimality/06-bayes-conjugate.html) with a conjugate prior on
\( (\bbeta,\sigma^2) \) and continued in
[Section 12.6](../ch12-intervals-and-bands/06-bayes.html) with credible sets and
predictive distributions. Mixed models are where the two traditions meet most closely,
because a random effect *is* a prior: the statement
\( \bu\sim\Normal_q(\bzero,\G) \) is formally indistinguishable from a normal prior
on \( q \) unknown constants. What separates the methods of
[Section 32.4](04-likelihood.html) from a fully Bayesian analysis is only that they
estimate \( \G \) and then treat it as known, which is empirical
Bayes (@prp-shr-empirical-bayes), while a Bayesian puts a prior on \( \G \) too and
integrates.

## The hierarchical reading

Write the model in three levels:
\[
\Y\mid\bbeta,\bu,\boldsymbol{\uptheta}
\sim\Normal_n(\X\bbeta+\Z\bu,\R),\quad
\bu\mid\boldsymbol{\uptheta}\sim\Normal_q(\bzero,\G),\quad
(\bbeta,\boldsymbol{\uptheta})\sim\pi .
\]{#eq-mix-hierarchy}

Level one is a linear model with known coefficients; level two is the prior on
\( \bu \) that the experiment's sampling of levels supplies; level three is the part
that is genuinely a matter of choice.

::: {#prp-mix-bayes}
[BLUP as a posterior mean, REML as a marginal posterior]

Consider @eq-mix-hierarchy with \( \G \) and \( \R \) positive definite and
\( \rank(\X)=p \).

::: {.enumerate options="label=(\alph*)"}
1. *(BLUP is a posterior mean.)* Fix \( \boldsymbol{\uptheta} \) and take the
   improper flat prior \( \pi(\bbeta)\propto1 \). Then the joint posterior of
   \( (\bbeta,\bu) \) is
   \[
   \begin{pmatrix}\bbeta\\\bu\end{pmatrix}\Bigm|\y
   \sim\Normal_{p+q}\Bigl(
   \begin{pmatrix}\hbeta\\\hat{\bu}\end{pmatrix},\ \mathbf{C}^{-1}\Bigr),
   \]
   with \( \hbeta,\hat{\bu} \) and \( \mathbf{C} \) as in @thm-mix-henderson. In
   particular \( \E(\bu\mid\y)=\hat{\bu} \) is the BLUP, \( \E(\bbeta\mid\y)=\hbeta \)
   is the generalized least squares estimate, and
   \( \Cov(\bu\mid\y) \) is the prediction error covariance of @eq-mix-prediction-variance.

2. *(REML is a marginal posterior.)* With the flat prior on \( \bbeta \) and any
   prior \( \pi(\boldsymbol{\uptheta}) \), the marginal posterior of the variance
   components is
   \[
   \pi(\boldsymbol{\uptheta}\mid\y)\ \propto\
   \pi(\boldsymbol{\uptheta})\exp\{\ell_R(\boldsymbol{\uptheta})\},
   \]{#eq-mix-marginal-posterior}

   with \( \ell_R \) the restricted log-likelihood of @def-mix-reml. The REML
   estimate is therefore the mode of the marginal posterior under a flat prior on
   \( \boldsymbol{\uptheta} \).

3. *(What the plug-in leaves out.)* Under (a) with \( \boldsymbol{\uptheta} \)
   replaced by an estimate, the reported posterior variance is
   \( \Cov(\bu\mid\y,\hat{\boldsymbol{\uptheta}}) \), whereas the correct marginal
   posterior variance is
   \( \E\{\Cov(\bu\mid\y,\boldsymbol{\uptheta})\mid\y\}
   +\Cov\{\E(\bu\mid\y,\boldsymbol{\uptheta})\mid\y\} \). The plug-in omits the second
   term, so it is too small.
:::

:::

::: {.proof}
(a) With \( \pi(\bbeta)\propto1 \) and \( \boldsymbol{\uptheta} \) fixed, the joint
density of \( (\y,\bu) \) times the prior is proportional to
\( \exp\{-\tfrac12Q(\bbeta,\bu)\} \) with \( Q \) the penalized criterion
@eq-mix-penalized, because
\( -2\log f(\y\mid\bbeta,\bu)=(\y-\X\bbeta-\Z\bu)\T\R^{-1}
(\y-\X\bbeta-\Z\bu) \) and \( -2\log f(\bu)=\bu\T\G^{-1}\bu \), up to terms free of
\( (\bbeta,\bu) \). By @thm-mix-henderson(b), \( Q \) is a quadratic in
\( \bw=(\bbeta\T,\bu\T)\T \) with positive definite matrix \( \mathbf{C} \), minimized
at \( \hat{\bw}=(\hbeta\T,\hat{\bu}\T)\T \), so
\( Q(\bw)=(\bw-\hat{\bw})\T\mathbf{C}(\bw-\hat{\bw})+Q(\hat{\bw}) \). The posterior
is therefore \( \Normal_{p+q}(\hat{\bw},\mathbf{C}^{-1}) \), and
@thm-mix-henderson(c) identifies the blocks.

(b) The joint posterior is
\( \pi(\bbeta,\boldsymbol{\uptheta}\mid\y)\propto
\pi(\boldsymbol{\uptheta})\exp\{\ell(\bbeta,\boldsymbol{\uptheta})\} \) with
\( \ell \) the marginal log-likelihood @eq-mix-loglik, the random effects having
already been integrated out in @eq-mix-marginal. Completing the square as in
@lem-mix-P(b),
\[
(\y-\X\bbeta)\T\V^{-1}(\y-\X\bbeta)
=\y\T\bP\y+(\bbeta-\hbeta)\T\X\T\V^{-1}\X(\bbeta-\hbeta),
\]
so integrating over \( \bbeta\in\Real^p \) contributes the factor
\( (2\pi)^{p/2}\det(\X\T\V^{-1}\X)^{-1/2} \) and
\[
\begin{aligned}
\int e^{\ell(\bbeta,\boldsymbol{\uptheta})}\,d\bbeta
=\exp\bigl[-\tfrac12\{&(n-p)\log2\pi+\log\det\V\\
&+\log\det(\X\T\V^{-1}\X)+\y\T\bP\y\}\bigr].
\end{aligned}
\]
By @eq-mix-reml-determinant the exponent is \( \ell_R(\boldsymbol{\uptheta}) \) up to the
additive constant \( \tfrac12\log\det(\X\T\X) \), which does not involve
\( \boldsymbol{\uptheta} \). That is @eq-mix-marginal-posterior.

(c) This is the decomposition of a variance into the mean of the conditional variance
plus the variance of the conditional mean (@prp-rv-total-covariance), applied to the
posterior distribution with \( \boldsymbol{\uptheta} \) as the conditioning variable.
The second term is nonnegative definite.
:::

Part (b) is Harville's (1974) observation, and it settles a question that the
frequentist derivation leaves slightly mysterious. REML looked like a technical
device for recovering lost degrees of freedom; it is in fact the exact marginal
likelihood of the variance components once the fixed effects have been integrated out
with a flat prior. The "\( n-p \)" is not a correction but a consequence.

::: {.idea}
Three readings of the same computation. Frequentist: \( \hat{\bu} \) minimizes
prediction mean squared error among unbiased linear predictors. Penalized: it
minimizes a sum of squares with a ridge penalty on \( \bu \). Bayesian: it is the
posterior mean of \( \bu \) under a normal prior with covariance \( \G \). Henderson's
equations compute all three at once, and the three readings differ only in what they
say about \( \G \): a covariance to be estimated, a penalty to be tuned, or a prior to
be chosen.
:::

## Priors for variance components

The conditional distributions in @eq-mix-hierarchy make the conjugate choices
obvious. Given \( \bu \) and \( \bbeta \), the sum of squares
\( \norm{\y-\X\bbeta-\Z\bu}^2 \) is an ordinary normal likelihood for
\( \sigma^2 \), so an inverse gamma prior is conjugate; given \( \bu \), the sum
\( \norm{\bu}^2 \) plays the same role for \( \sigma_a^2 \). This is what makes the
Gibbs sampler below so short. What is not obvious is which inverse gamma.

::: {.warning}
The prior \( \sigma_a^2\sim\text{IG}(\epsilon,\epsilon) \) with \( \epsilon \) small
is often described as non-informative. It is not, and the limit does not exist. The
marginal (restricted) likelihood \( \exp\{\ell_R\} \) is continuous and *strictly
positive* at \( \sigma_a^2=0 \), because \( \sigma_a^2=0 \) is a perfectly good model.
Any prior density whose integral diverges at the origin, such as
\( \pi(\sigma_a^2)\propto1/\sigma_a^2 \) — equivalently a flat prior on
\( \log\sigma_a \) — therefore gives an improper posterior:
\( \int\pi(\sigma_a^2)e^{\ell_R}\,d\sigma_a^2=\infty \). For small \( \epsilon \) the
\( \text{IG}(\epsilon,\epsilon) \) posterior is proper but arbitrarily close to
improper, and with few groups the answer depends on \( \epsilon \). Gelman (2006)
made this point and recommended instead a prior that is flat or half-\( t \) on the
*standard deviation* \( \sigma_a \), which has finite mass near zero.
:::

A convenient weakly informative default is the half-Cauchy on \( \sigma_a \) with
scale \( A \), chosen large relative to the plausible size of a group effect. It
keeps the Gibbs sampler conjugate through the scale-mixture representation
\[
\sigma_a^2\mid\xi\sim\text{IG}(\tfrac12,\ 1/\xi),
\qquad
\xi\sim\text{IG}(\tfrac12,\ 1/A^2),
\]
whose marginal for \( \sigma_a \) is half-Cauchy\( (0,A) \). The auxiliary variable
\( \xi \) is sampled alongside everything else.

## A Gibbs sampler in a dozen lines

For the one-way model with a flat prior on \( \mu \), the four full conditionals are
all standard. Given \( \mu,\sigma^2,\sigma_a^2 \), the effects are independent with
\[
a_k\mid\text{rest}\sim
\Normal\Bigl(\frac{(m/\sigma^2)(\bar y_k-\mu)}{m/\sigma^2+1/\sigma_a^2},\
\frac{1}{m/\sigma^2+1/\sigma_a^2}\Bigr),
\]
which is @eq-mix-shrinkage with the current parameter values; given the effects,
\( \mu\sim\Normal(\overline{y-\Z\mathbf{a}},\ \sigma^2/n) \); and the two variances
have inverse gamma conditionals with shapes \( n/2 \) and \( g/2 \) plus their prior
shapes. Nothing needs tuning.

```{.python .run #cell-gibbs-sampler}
import numpy as np

g, m = 12, 4
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

Y = y.reshape(g, m)                             # the REML fit of Section 32.4
ms_between = m * np.sum((Y.mean(axis=1) - Y.mean()) ** 2) / (g - 1)
ms_within = np.sum((Y - Y.mean(axis=1)[:, None]) ** 2) / (n - g)
reml_s2, reml_s2a = ms_within, (ms_between - ms_within) / m

def gibbs(y, g, m, draws=40000, burn=2000, prior="inverse-gamma", A=1.0, seed=0):
    """Posterior draws of (mu, sigma^2, sigma_a^2, a) with a flat prior on mu."""
    rng = np.random.default_rng(seed)
    Y = y.reshape(g, m)
    n = g * m
    mu_, s2, s2a, xi = Y.mean(), Y.var(), 0.1, 1.0
    out = np.empty((draws, 3 + g))
    for t in range(draws + burn):
        prec = m / s2 + 1 / s2a                                  # a_k | rest
        a = rng.normal((m / s2) * (Y.mean(axis=1) - mu_) / prec, np.sqrt(1 / prec))
        mu_ = rng.normal((y - np.repeat(a, m)).mean(), np.sqrt(s2 / n))   # mu | rest
        resid = y - mu_ - np.repeat(a, m)
        s2 = 1 / rng.gamma(0.001 + n / 2, 1 / (0.001 + resid @ resid / 2))
        if prior == "inverse-gamma":                             # sigma_a^2 | rest
            s2a = 1 / rng.gamma(0.001 + g / 2, 1 / (0.001 + a @ a / 2))
        else:                                                    # half-Cauchy(0, A)
            s2a = 1 / rng.gamma((g + 1) / 2, 1 / (1 / xi + a @ a / 2))
            xi = 1 / rng.gamma(1.0, 1 / (1 / A**2 + 1 / s2a))
        if t >= burn:
            out[t - burn] = np.concatenate([[mu_, s2, s2a], a])
    return out
```

```{.python .run #cell-gibbs-demo}
# a shorter run than the 40 000 draws quoted in the text, so the cell finishes quickly
short = gibbs(y, g, m, draws=4000, burn=500, prior="half-Cauchy", seed=12)
print(f"posterior median of sigma_a^2  {np.median(short[:, 2]):.4f}")
print(f"posterior mean of mu           {short[:, 0].mean():.4f}")
print(f"REML estimate of sigma_a^2     {reml_s2a:.4f}")
```

::: {#exm-mix-gibbs}
[The proficiency study, four ways]

Run the sampler for \( 40\,000 \) draws on the data
of @exm-mix-proficiency. Under the \( \text{IG}(0.001,0.001) \) prior the posterior
median of \( \sigma_a^2 \) is \( 0.0356 \); under the
half-Cauchy prior with \( A=1 \) it is \( 0.0438 \). The
REML estimate \( 0.0386 \) sits between them. The two priors
therefore disagree by about a fifth of the estimate, with twelve
laboratories — which is already a reason to report the prior.

The posterior means of the twelve effects track the BLUPs closely: their correlation
exceeds \( 0.99 \) and the largest in absolute value is
\( 0.294 \) against the BLUP's
\( 0.321 \), slightly more shrunk — not because the posterior
sits below the REML estimate, since its mean \( 0.0431 \) is
above it, but because the shrinkage factor
\( B=m\sigma_a^2/(\sigma^2+m\sigma_a^2) \) is concave in \( \sigma_a^2 \): averaging
\( B \) over a posterior with mass near zero weights the data less than
\( \hat\sigma_a^2 \) does. Their posterior
standard deviations average \( 0.1218 \), against the plug-in
prediction standard error \( 0.1156 \) of
@eq-mix-prediction-variance-oneway: \( 5 \) per cent wider, which is the
uncertainty about \( \sigma_a^2 \) that the plug-in throws away, exactly as
@prp-mix-bayes(c) says. The same happens for the fixed effect: the posterior standard
deviation of \( \mu \) is \( 0.0709 \), against the plug-in
\( \sqrt{\hat\lambda_1/n}=0.0672 \).

Panel (b) of [Figure 32.6.1](06-bayes.html#fig-mix-gibbs) repeats the analysis on the
first \( 5 \) laboratories. The posterior medians are now
\( 0.0435 \) and \( 0.0693 \), a
difference of sixty per cent, against a REML estimate of
\( 0.0493 \). With five groups the prior is doing a
large part of the work, and saying so is not optional.
:::

::: {when-format="html"}
![**Figure 32.6.1.** Posterior densities of the between-laboratory standard deviation
\( \sigma_a \) under two priors, from a Gibbs sampler with 40 000 draws. The dotted
line marks the REML estimate. (a) All twelve laboratories. (b) The first five, where
the prior matters a great deal.](gibbs_posterior.svg){#fig-mix-gibbs width=100%}
:::

::: {when-format="pdf"}
![Posterior densities of the between-laboratory standard deviation \( \sigma_a \)
under two priors, from a Gibbs sampler with 40 000 draws. The dotted line marks the
REML estimate. (a) All twelve laboratories. (b) The first five, where the prior
matters a great deal.](gibbs_posterior.pdf){width=100%}
:::

## What the Bayesian version buys, and what it costs

Three gains. The uncertainty about \( \boldsymbol{\uptheta} \) is carried through
automatically, so intervals for fixed effects and prediction intervals for random
effects are honestly wider, with no need for the corrections of
[Section 32.5](05-inference.html). The boundary problem disappears: the posterior of
\( \sigma_a^2 \) is a proper distribution on \( [0,\infty) \), with no point mass and
no chi-bar-squared mixture to calibrate. And the same machinery handles models with no
closed form, which is what
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) needs for richer
covariance structures and Chapter 40 needs once the response is not normal.

Two costs. A prior for the variance components must be chosen, and with few groups it
matters, as panel (b) shows; the only defensible practice is to state it and show the
sensitivity. And the answers are not exactly the frequentist ones: the two frameworks
agree on the point estimates, but they part company on the spread, and the
frequentist coverage of a credible interval for a variance component is not
guaranteed (@prp-ci-average-coverage).

## Exercises

### A. Check your understanding

::: {#exr-mix-prior-is-random-effect}
[A1]

Explain in two sentences why a random effect and a prior are formally the same
object, and name one thing that distinguishes them in practice.
:::

### B. Practice

::: {#exr-mix-improper-posterior}
[B1]

Show that in the one-way model the marginal posterior density of \( \sigma_a^2 \)
under the prior \( \pi(\sigma_a^2)\propto\sigma_a^{-2} \) is not integrable, by
showing that \( \exp\{\ell_R(\sigma_a^2)\} \) tends to a strictly positive limit as
\( \sigma_a^2\downarrow0 \) with \( \sigma^2 \) profiled out. Contrast with the
behaviour of the prior \( \pi(\sigma_a)\propto1 \).
:::

::: {.solution}
By @eq-mix-oneway-reml, \( -2\ell_R \) is continuous in
\( (\lambda_1,\lambda_2) \) on \( \lambda_1\ge\lambda_2>0 \), and at
\( \sigma_a^2=0 \) it equals \( (n-1)\log\lambda_2+(\text{SSB}+\text{SSE})/\lambda_2 \),
finite and minimized at a positive \( \lambda_2 \). So
\( \exp\{\ell_R\} \to c>0 \) as \( \sigma_a^2\downarrow0 \), and
\( \int_0^\delta c\,\sigma_a^{-2}\,d\sigma_a^2=\infty \). With
\( \pi(\sigma_a)\propto1 \) the Jacobian gives
\( \pi(\sigma_a^2)\propto\sigma_a^{-1} \), and
\( \int_0^\delta\sigma_a^{-1}d(\sigma_a^2)=2\int_0^{\sqrt\delta}d\sigma_a<\infty \):
the posterior is proper near the origin.
:::

::: {#exr-mix-conditional-variance}
[B2]

In @prp-mix-bayes(c), suppose \( \boldsymbol{\uptheta} \) takes only two values with
posterior probabilities \( 1/2 \) each, giving conditional posterior means
\( \hat{\bu}_1,\hat{\bu}_2 \) and common conditional variance \( v \) for a
coordinate. Compute the correct marginal posterior variance, and say how large the
disagreement between \( \hat u_1 \) and \( \hat u_2 \) must be before the plug-in
understates the standard deviation by ten per cent.
:::

### C. Going deeper

::: {#exr-mix-reml-prior}
[C1]

Formula @eq-mix-marginal-posterior shows REML is the posterior mode of
\( \boldsymbol{\uptheta} \) under a flat prior on \( \boldsymbol{\uptheta} \). A mode
is not invariant to reparameterization: the mode of the posterior of
\( \sigma_a^2 \) is not the square of the mode of the posterior of \( \sigma_a \).
Reconcile this with @thm-mix-reml(c), which says the REML estimate is well defined.
What is invariant, and what is not?
:::

::: {#exr-mix-empirical-bayes-james-stein}
[C2]

Take the balanced one-way model with \( \sigma^2 \) known, \( m=1 \) so that
\( \bar y_k=y_k \), and \( \mu \) known to be zero. Show that plugging the unbiased
moment estimate of \( \sigma_a^2 \) into @eq-mix-shrinkage gives a shrinkage factor of
the form \( 1-c\,\sigma^2/\sum_ky_k^2 \), and compare with the James–Stein estimator
of @thm-shr-james-stein. Which constant \( c \) does the plug-in give, and which does
James–Stein use?
:::

::: {.solution}
With \( m=1 \) and \( \mu=0 \) the data are \( y_k\sim\Normal(a_k,\sigma^2) \) with
\( a_k\sim\Normal(0,\sigma_a^2) \), so marginally
\( y_k\sim\Normal(0,\sigma^2+\sigma_a^2) \). The moment estimate from
\( \sum_ky_k^2 \) is \( \hat\sigma_a^2=\sum_ky_k^2/g-\sigma^2 \), and
\[
B=\frac{\hat\sigma_a^2}{\hat\sigma_a^2+\sigma^2}
=1-\frac{g\sigma^2}{\sum_ky_k^2},
\]
so \( c=g \) in the stated form. The James–Stein estimator of
@thm-shr-james-stein uses \( 1-(g-2)\sigma^2/\sum_ky_k^2 \). The difference is the
same degrees-of-freedom correction that separates ML from REML: \( g-2 \) is the
constant that makes \( (g-2)\sigma^2/\sum y_k^2 \) unbiased for
\( \sigma^2/(\sigma^2+\sigma_a^2) \), and it is what makes the estimator dominate the
unshrunk one for every \( \boldsymbol{\upalpha} \), not merely on average.
:::
