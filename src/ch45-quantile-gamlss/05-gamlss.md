# Models for location, scale and shape

Quantile and expectile regression estimate one feature of the conditional law at a
time and never write the law down. The opposite strategy writes it down in full:
choose a parametric family with enough parameters to bend as the data require, and
let every one of them depend on the covariates through its own structured additive
predictor. Rigby and Stasinopoulos (2005) named the result generalized additive
models for location, scale and shape.

The gain is a likelihood, and with it penalized estimation, effective degrees of
freedom, AIC and Bayesian posteriors. The cost is a parametric assumption about the
shape of the conditional distribution, which the data must be allowed to
contradict — the subject of [Section 45.6](06-comparing.html).

## The model

::: {#def-qnt-gamlss}
[Generalized additive model for location, scale and shape]

Let \( D(\theta_1,\dots,\theta_K) \) be a family of distributions on the response
space with density or mass function \( d(y;\theta_1,\dots,\theta_K) \); in this
chapter \( \theta_k \) denotes the \( k \)th parameter of that family, not the
natural parameter of @def-glm-model. A **GAMLSS** asserts that the responses
\( Y_1,\dots,Y_n \) are independent with
\[
Y_i\sim D(\theta_{i1},\dots,\theta_{iK}),\qquad
g_k(\theta_{ik})=\eta_{ik},\qquad k=1,\dots,K,
\]
where each \( g_k \) is a strictly increasing **link** mapping the admissible range of
\( \theta_k \) onto \( \Real \) (the identity for an unrestricted parameter, the
logarithm for a positive one, the logit for one in \( (0,1) \)), and each
\( \eta_{ik} \) is a structured additive predictor (@def-add-star)
\[
\eta_{ik}=\beta_{0k}+\sum_{j=1}^{J_k}f_{jk}(v_{ij}),
\qquad f_{jk}(\cdot)\text{ represented by }\Z_{jk}\bgamma_{jk},
\]{#eq-qnt-gamlss-predictor}

with each \( \bgamma_{jk} \) carrying a quadratic penalty
\( \lambda_{jk}\bgamma_{jk}\T\bP_{jk}\bgamma_{jk} \). The predictors for different
parameters may use different covariates and different amounts of structure.
:::

Taking \( K=1 \) with \( D \) an exponential dispersion family and \( \eta \) linear
recovers @def-glm-model; with additive \( \eta \) it recovers @thm-add-gam. The new
content is \( K>1 \): the scale, and the shape, become regression problems.

::: {#exm-qnt-normal-gamlss}
[The heteroscedastic normal model]

The smallest interesting case has \( K=2 \), \( D=\Normal(\mu,\sigma^2) \), identity
link on \( \mu \) and log link on \( \sigma \):
\[
Y_i\sim\Normal(\mu_i,\sigma_i^2),\qquad
\mu_i=\eta_{i1},\qquad \log\sigma_i=\eta_{i2}.
\]{#eq-qnt-locscale-gamlss}

The log link keeps \( \sigma_i>0 \) whatever the predictor does, as it kept a Poisson
mean positive in [Chapter 37](../ch37-counts/index.html). With both predictors linear
this is @prp-qnt-beyond-mean(b) with normal \( \varepsilon \); with penalized splines
it is the fit of @exm-qnt-engel-gamlss.
:::

## Estimation

::: {#thm-qnt-gamlss}
[Penalized likelihood and block-wise Fisher scoring]

Write \( \ell_i(\boldsymbol{\uptheta}_i)=\log d(y_i;\theta_{i1},\dots,\theta_{iK}) \) and
collect the coefficients of predictor \( k \) into \( \bgamma_{k} \) with design
\( \Z_k \), penalty matrix \( \bP_k \) and smoothing parameter \( \lambda_k \),
so that \( \boldsymbol{\upeta}_k=\Z_k\bgamma_k \). The **penalized log-likelihood** is
\[
\ell_p(\bgamma_1,\dots,\bgamma_K)
=\sum_{i=1}^n\ell_i(\boldsymbol{\uptheta}_i)
 -\tfrac12\sum_{k=1}^K\lambda_k\,\bgamma_k\T\bP_k\bgamma_k .
\]{#eq-qnt-penalized-loglik}

::: {.enumerate options="label=(\alph*)"}
1. **Score.** With \( u_{ik}=\partial\ell_i/\partial\eta_{ik} \) and
   \( \bu_k=(u_{1k},\dots,u_{nk})\T \),
   \[
   \frac{\partial\ell_p}{\partial\bgamma_k}=\Z_k\T\bu_k-\lambda_k\bP_k\bgamma_k .
   \]{#eq-qnt-gamlss-score}

2. **Fisher scoring is penalized weighted least squares.** Let
   \( w_{ik}=-\E(\partial^2\ell_i/\partial\eta_{ik}^2) \) and
   \( \W_k=\diag(w_{1k},\dots,w_{nk}) \). Updating block \( k \) with the others held
   fixed by one Fisher scoring step gives
   \[
   \bgamma_k^{\text{new}}
   =\bigl(\Z_k\T\W_k\Z_k+\lambda_k\bP_k\bigr)^{-1}\Z_k\T\W_k\bz_k,
   \qquad \bz_k=\boldsymbol{\upeta}_k+\W_k^{-1}\bu_k,
   \]{#eq-qnt-gamlss-step}

   the penalized version of the iteratively reweighted least squares
   of @thm-glm-irls, with \( \bz_k \) the working response for parameter \( k \). Cycling
   over \( k=1,\dots,K \) until nothing changes is
   backfitting (@thm-add-backfitting) applied to the \( K \) blocks.

3. **The normal location-scale model is orthogonal.**
   For @eq-qnt-locscale-gamlss,
   \[
   u_{i1}=\frac{y_i-\mu_i}{\sigma_i^{2}},\quad w_{i1}=\frac{1}{\sigma_i^{2}},
   \qquad
   u_{i2}=\frac{(y_i-\mu_i)^{2}}{\sigma_i^{2}}-1,\quad w_{i2}=2,
   \]{#eq-qnt-normal-scores}

   and the expected cross-information
   \( -\E(\partial^2\ell_i/\partial\eta_{i1}\partial\eta_{i2}) \) is zero. The two
   blocks are therefore orthogonal: neither Fisher scoring step needs the other's
   derivatives, the working response for \( \mu \) is \( \bz_1=\y \), and the working
   response for \( \log\sigma \) is
   \( \bz_2=\boldsymbol{\upeta}_2+\tfrac12\bu_2 \).
:::

:::

::: {.proof}
(a) By the chain rule \( \partial\ell_i/\partial\bgamma_k
=(\partial\ell_i/\partial\eta_{ik})\,\partial\eta_{ik}/\partial\bgamma_k \), and
\( \partial\eta_{ik}/\partial\bgamma_k \) is the \( i \)th row of \( \Z_k \).
Summing \( u_{ik} \) times that row over \( i \) gives \( \Z_k\T\bu_k \), and
differentiating the quadratic penalty gives \( -\lambda_k\bP_k\bgamma_k \).

(b) Differentiating @eq-qnt-gamlss-score once more, the contribution of observation
\( i \) to the second derivative within block \( k \) is
\( \partial^2\ell_i/\partial\eta_{ik}^2 \) times the outer product of the
\( i \)th row of \( \Z_k \) with itself; taking expectations and summing over the
independent responses gives \( -\Z_k\T\W_k\Z_k-\lambda_k\bP_k \). A Fisher
scoring step therefore reads
\[
\bgamma_k^{\text{new}}
=\bgamma_k+\bigl(\Z_k\T\W_k\Z_k+\lambda_k\bP_k\bigr)^{-1}
   \bigl(\Z_k\T\bu_k-\lambda_k\bP_k\bgamma_k\bigr).
\]
Adding and subtracting \( \Z_k\T\W_k\Z_k\bgamma_k \) inside the bracket turns the right
side into \( (\Z_k\T\W_k\Z_k+\lambda_k\bP_k)^{-1}\Z_k\T\W_k
(\Z_k\bgamma_k+\W_k^{-1}\bu_k) \), which is @eq-qnt-gamlss-step.

(c) With \( \ell_i=-\log\sigma_i-(y_i-\mu_i)^2/(2\sigma_i^2)-\tfrac12\log(2\pi) \),
\( \mu_i=\eta_{i1} \) and \( \sigma_i=e^{\eta_{i2}} \), direct differentiation gives
\( \partial\ell_i/\partial\eta_{i1}=(y_i-\mu_i)/\sigma_i^2 \) and
\[
\frac{\partial\ell_i}{\partial\eta_{i2}}
=\sigma_i\frac{\partial\ell_i}{\partial\sigma_i}
=\sigma_i\Bigl\{-\frac1{\sigma_i}+\frac{(y_i-\mu_i)^2}{\sigma_i^{3}}\Bigr\}
=\frac{(y_i-\mu_i)^2}{\sigma_i^{2}}-1 .
\]
Differentiating again, \( \partial^2\ell_i/\partial\eta_{i1}^2=-1/\sigma_i^2 \),
which is already non-random, so \( w_{i1}=1/\sigma_i^2 \); and
\( \partial^2\ell_i/\partial\eta_{i2}^2=-2(y_i-\mu_i)^2/\sigma_i^2 \), whose
expectation is \( -2 \), so \( w_{i2}=2 \). Finally
\( \partial^2\ell_i/\partial\eta_{i1}\partial\eta_{i2}=-2(y_i-\mu_i)/\sigma_i^{2} \),
with expectation zero because \( \E(Y_i-\mu_i)=0 \). The working responses follow from
\( \bz_k=\boldsymbol{\upeta}_k+\W_k^{-1}\bu_k \): for \( k=1 \),
\( z_{i1}=\mu_i+\sigma_i^2(y_i-\mu_i)/\sigma_i^2=y_i \).
:::

::: {.idea}
A GAMLSS fit is a stack of weighted smoothers. Each parameter's block is updated by
the same penalized weighted least squares step that fitted a generalized additive
model, with weights and working responses supplied by the log-likelihood of the
chosen family. The rest is bookkeeping, which is why the framework accommodates any
family whose derivatives can be written down.
:::

Orthogonality in part (c) is special to this parameterization of the normal family;
for most families the cross-information is not zero, and cycling over blocks then
converges more slowly or must give way to a joint Newton step. The effective degrees
of freedom of a block is the trace
\( \tr\{(\Z_k\T\W_k\Z_k+\lambda_k\bP_k)^{-1}\Z_k\T\W_k\Z_k\} \)
that @thm-shr-ridge attached to ridge regression and @def-reg-penalized to penalized
least squares; the total of these traces goes into AIC (@prp-sel-aic-bic).

## Choosing the family

The distribution is a modelling choice on a par with the link, made from what the
response is. Positive with a right tail: gamma, inverse Gaussian, lognormal. Positive
with adjustable skewness: the Box–Cox family, whose third parameter is a power
estimated rather than assumed, generalizing @def-tr-box-cox. Heavy tailed: a
\( t \) family with the degrees of freedom as a third parameter. Counts with too many
zeros: a zero-inflated or hurdle family (@def-cnt-zero-models). Proportions: a beta
family reparameterized into a mean and a precision.

The ancestor is the LMS method of Cole and Green (1992) for growth charts: a Box–Cox
power, a median and a coefficient of variation, each a smooth function of age — a
three-parameter GAMLSS written a decade before the name, and still how paediatric
reference charts are built.

## Warnings

::: {.warning}
[Identifiability]

Each predictor in @eq-qnt-gamlss-predictor contains an intercept and several smooth
terms, and a constant moves freely between them. The remedy is that of
[Chapter 44](../ch44-additive-models/index.html): constrain every smooth term to sum
to zero over the data (@def-add-model), separately in every predictor. Without the
constraints the penalized information is singular in as many directions as there are
smooth terms. When the *same* covariate enters two predictors, the concurvity
of @prp-add-concurvity appears between parameters, and no ceteris paribus reading of a
single parameter's curve survives.
:::

::: {.warning}
[The likelihood is unbounded]

Let the mean predictor be flexible enough to interpolate one observation, say
\( \mu(x_1)=y_1 \), and the scale predictor flexible enough to be small near
\( x_1 \) and unchanged elsewhere. The contribution of observation \( 1 \)
to @eq-qnt-penalized-loglik is then \( -\log\sigma(x_1) \), which grows without bound
as \( \sigma(x_1)\to0 \) while the others stay finite: the unpenalized likelihood of
a flexible location-scale model has no maximum. What stops the fit running away is the
penalty, so the smoothing parameters are what make the problem well posed — and they
cannot be chosen by maximizing the same likelihood.
:::

::: {#exm-qnt-unbounded}
[Watching it happen]

Fitting @eq-qnt-locscale-gamlss to Engel's budgets with the scale penalty set to
\( 0.01 \) instead of the chosen \( 10.0 \) raises the in-sample log-likelihood from
\( 277.4 \) to \( 282.9 \), a gain a naive search would take. The ten-fold
cross-validated log score — minus the mean predictive log density, so that smaller is
better (@def-qnt-score and the remark in
[Section 45.6](06-comparing.html)) — tells the truth: it moves from \( -1.1102 \)
to \( 59.9 \). The wilder scale curve dips towards zero where the mean spline
interpolates, and a held-out household near such a dip is assigned a density
indistinguishable from zero.
:::

## A location-scale fit to Engel's budgets

::: {#exm-qnt-engel-gamlss}
[Mean and scale together]

Take @eq-qnt-locscale-gamlss with a cubic B-spline basis of \( 12 \) columns for each
predictor — knots at income quantiles, second-difference penalties (@def-smo-pspline)
— and measure income and food expenditure in thousands of francs. The two smoothing
parameters are chosen by ten-fold cross-validation of the log score over a grid,
giving \( 31.6 \) for the mean and \( 10.0 \) for the scale, with effective degrees
of freedom \( 8.75 \) and \( 5.59 \), \( 14.33 \) in all.

The fitted scale rises from \( 40.5 \) francs at an income of \( 500 \) to
\( 184.5 \) at an income of \( 2000 \), a ratio of \( 4.55 \); the single constant
the homoscedastic model fits is \( 107.4 \) francs, wrong in both directions at once.
At an income of \( 1000 \) francs the fitted mean budget is \( 645.3 \) francs and
the fitted ninetieth percentile \( 746.5 \). With the response in thousands of francs
the AIC is \( -526.1 \) against \( -368.0 \) for the same mean spline with a constant
scale, and the cross-validated log score \( -1.1102 \) against \( -0.6938 \): most of
what a homoscedastic model gets wrong here is the scale. The implied quantile curves
\( \hat\mu(x)+\hat\sigma(x)z_{\tau} \) cannot cross, one clear advantage over the
separate fits of [Section 45.3](03-inference.html).
:::

::: {when-format="html"}
![**Figure 45.5.1.** A normal location-scale GAMLSS for Engel's budgets, both
predictors penalized splines in income. (a) The implied \( 0.1 \), \( 0.5 \) and
\( 0.9 \) quantile curves. (b) The fitted scale against the homoscedastic
constant.](gamlss_engel.svg){#fig-qnt-gamlss width=100%}
:::

::: {when-format="pdf"}
![A normal location-scale GAMLSS for Engel's budgets, both predictors penalized
splines in income. (a) The implied \( 0.1 \), \( 0.5 \) and \( 0.9 \) quantile
curves. (b) The fitted scale against the homoscedastic
constant.](gamlss_engel.pdf){width=100%}
:::

```{.python .run #cell-gamlss-engel-fit}
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.interpolate import BSpline
data = sm.datasets.engel.load_pandas().data
DEGREE, N_INNER = 3, 8
CHOSEN = (10.0 ** 1.5, 10.0 ** 1.0)                # the pair the grid search picks out
WILD = 0.01                                        # a scale penalty small enough to run away

income = data["income"].to_numpy() / 1000.0        # thousands of francs
food = data["foodexp"].to_numpy() / 1000.0
n = len(food)

lo, hi = income.min(), income.max()
inner = np.quantile(income, np.linspace(0, 1, N_INNER + 2)[1:-1])
KNOTS = np.r_[[lo] * (DEGREE + 1), inner, [hi] * (DEGREE + 1)]

def basis(z):
    """The cubic B-spline basis of Section 43.3, with knots at income quantiles."""
    z = np.clip(np.atleast_1d(z), lo, hi)
    return np.asarray(BSpline.design_matrix(z, KNOTS, DEGREE).todense())

def difference_penalty(q, order=2):
    """K = D'D for the order-th difference matrix D: the penalty of Section 43.4."""
    D = np.diff(np.eye(q), order, axis=0)
    return D.T @ D

def gamlss_normal(Bm, Bs, y, lam_mu, lam_sig, iters=100, tol=1e-9):
    """Penalized Fisher scoring for y ~ N(mu, sigma^2), mu = Bm b_mu, log sigma = Bs b_sig.

    One block at a time: the mean step is a penalized weighted least squares fit with
    weights 1 / sigma^2, and the scale step is a penalized Fisher scoring step whose
    expected information is the constant 2.
    """
    Km, Ks = difference_penalty(Bm.shape[1]), difference_penalty(Bs.shape[1])
    b_mu = np.linalg.solve(Bm.T @ Bm + lam_mu * Km, Bm.T @ y)
    b_sig = np.linalg.solve(Bs.T @ Bs, Bs.T @ np.full(len(y), np.log(y.std())))
    for _ in range(iters):
        sigma = np.exp(np.clip(Bs @ b_sig, -8, 4))
        w = 1.0 / sigma ** 2                                  # working weights for the mean
        new_mu = np.linalg.solve(Bm.T @ (w[:, None] * Bm) + lam_mu * Km, Bm.T @ (w * y))
        u = (y - Bm @ new_mu) ** 2 / sigma ** 2 - 1.0         # score in the log-sigma direction
        new_sig = b_sig + np.linalg.solve(2 * Bs.T @ Bs + lam_sig * Ks,
                                          Bs.T @ u - lam_sig * Ks @ b_sig)
        step = max(np.max(np.abs(new_mu - b_mu)), np.max(np.abs(new_sig - b_sig)))
        b_mu, b_sig = new_mu, new_sig
        if step < tol:
            break
    return b_mu, b_sig

def neg_log_score(mu, sigma, y):
    """Minus the mean log density: the logarithmic score of Section 45.6."""
    return -np.mean(stats.norm.logpdf(y, mu, sigma))

folds = np.random.default_rng(7).permutation(n) % 10

def cross_validate(design, lam_mu, lam_sig):
    """Ten-fold cross-validated log score of one specification."""
    total = 0.0
    for k in range(10):
        train = folds != k
        Bm, Bs = design(income[train])
        b_mu, b_sig = gamlss_normal(Bm, Bs, food[train], lam_mu, lam_sig)
        Bm, Bs = design(income[~train])
        total += neg_log_score(Bm @ b_mu, np.exp(np.clip(Bs @ b_sig, -8, 4)), food[~train])
    return total / 10

spline_scale = lambda z: (basis(z), basis(z))
Bm, Bs = spline_scale(income)
b_mu, b_sig = gamlss_normal(Bm, Bs, food, *CHOSEN)
mu, sigma = Bm @ b_mu, np.exp(Bs @ b_sig)
cv_chosen = cross_validate(spline_scale, *CHOSEN)
print(f"smoothing parameters: {CHOSEN[0]:.3g} for the mean, {CHOSEN[1]:.3g} for the scale")
print(f"cross-validated log score {cv_chosen:.4f}")
```

```{.python .run #cell-gamlss-engel-unbounded}
# barely penalize the scale: the in-sample likelihood climbs, the honest score collapses
b_mu_w, b_sig_w = gamlss_normal(Bm, Bs, food, CHOSEN[0], WILD)
wild_loglik = np.sum(stats.norm.logpdf(food, Bm @ b_mu_w, np.exp(np.clip(Bs @ b_sig_w, -8, 4))))
loglik = np.sum(stats.norm.logpdf(food, mu, sigma))
wild_cv = cross_validate(spline_scale, CHOSEN[0], WILD)
print(f"log-likelihood: {loglik:.1f} at the chosen penalty, {wild_loglik:.1f} at"
      f" lambda_sigma = {WILD:g}")
print(f"cross-validated log score: {cv_chosen:.4f} against {wild_cv:.4f}")
```

## Variations

**Bayesian GAMLSS.** Every quadratic penalty is a normal prior on the
coefficients (@prp-prc-priors), so @eq-qnt-penalized-loglik is a log posterior up to
the smoothing-parameter priors, and
[Section 39.4](../ch39-glms-in-practice-bayes/04-bayesian-glms-mcmc.html) transfers
whole, with those parameters sampled rather than tuned (Klein, Kneib, Lang and Sohn
2015). **Boosting.** Feeding a GAMLSS log-likelihood to the componentwise boosting of
[Section 30.4](../ch30-regularization-boosting/04-componentwise-boosting.html) gives
an algorithm that also selects which terms enter which predictor (Mayr, Fenske,
Hofner, Kneib and Schmid 2012). **Multivariate responses.** With a multivariate
family the dependence parameters become regression problems too;
[Section 45.7](07-summary-and-notes.html) points at the literature.

## Exercises

### A. Check your understanding

::: {#exr-qnt-gamlss-edf}
[A1]

In @exm-qnt-engel-gamlss the mean block has \( 12 \) basis columns and effective
degrees of freedom \( 8.75 \). What are the two extreme values the effective degrees
of freedom of that block could take, and which smoothing parameters produce them?
:::

::: {.solution}
As \( \lambda_{\mu}\to0 \) the trace tends to \( 12 \), the number of columns; as
\( \lambda_{\mu}\to\infty \) it tends to the dimension of the null space of the
second-difference penalty, namely \( 2 \) (constants and straight lines). The fitted
value \( 8.75 \) is much nearer the flexible end.
:::

### B. Practice

::: {#exr-qnt-gamlss-gamma}
[B1]

For a gamma response with mean \( \mu \) and shape \( \nu \), density
\( d(y)=\{\nu/\mu\}^{\nu}y^{\nu-1}e^{-\nu y/\mu}/\Gamma(\nu) \), with log links on
both parameters, compute \( u_{i1} \) and \( w_{i1} \) of @thm-qnt-gamlss(b) for the
mean block. Is the model orthogonal in the sense of part (c)?
:::

::: {.solution}
With \( \eta_1=\log\mu \), \( \partial\ell/\partial\mu=\nu(y-\mu)/\mu^2 \), so
\( u_{1}=\mu\,\partial\ell/\partial\mu=\nu(y-\mu)/\mu \). Then
\( \partial^2\ell/\partial\eta_1^2=-\nu y/\mu \), with expectation \( -\nu \), so
\( w_{1}=\nu \). The cross derivative
\( \partial^2\ell/\partial\eta_1\partial\eta_2 \) involves
\( \nu(y-\mu)/\mu \) times \( \nu \), whose expectation is zero, so the gamma model
with these links is also orthogonal — a fact the algorithm exploits, and one that
does not extend to three-parameter families.
:::

::: {#exr-qnt-gamlss-working}
[B2]

Verify from @eq-qnt-normal-scores that the working response for the scale block of
the normal location-scale model is
\( z_{i2}=\log\sigma_i+\tfrac12\{(y_i-\mu_i)^2/\sigma_i^2-1\} \), and show that the
scale step is therefore an ordinary penalized least squares smooth of the squared
standardized residuals, shifted and halved.
:::

### C. Going deeper

::: {#exr-qnt-gamlss-orthogonality}
[C1]

Show that the normal family parameterized by \( (\mu,\sigma) \) is orthogonal for
*any* increasing link on \( \sigma \), not just the logarithm, and find a
parameterization of the normal family in which it is not. What does orthogonality buy
in @thm-qnt-gamlss(b)?
:::

