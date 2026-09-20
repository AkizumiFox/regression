# Covariance models for longitudinal data

A longitudinal study follows subjects over time. Its data are clustered by
subject, but with a structure the previous sections have not used: observations
inside a cluster are *ordered*, and the distance between two of them matters.
Compound symmetry ignores the ordering; an unstructured covariance uses it but
pays \( m(m+1)/2 \) parameters. Between the extremes lies a catalogue of models
with two to four parameters, each a different story about where the correlation
comes from.

## Three sources of correlation

A longitudinal measurement is usefully read as the sum of three parts, a
decomposition due to Diggle (1988):
\[
y_{ij}=\x_{(ij)}\T\bbeta+\bz_{(ij)}\T\bu_i+W_i(t_{ij})+e_{ij}.
\]{#eq-cls-three-sources}

The vector \( \bu_i \) holds **random effects**: a subject's own level or slope,
anything persisting for the whole study. The process \( W_i(\cdot) \) is **serial
correlation**, making nearby times alike and distant times nearly independent. And
\( e_{ij} \) is **measurement error**. The three leave different fingerprints —
correlation that does not decay with the gap, correlation that does, and a
discontinuity at gap zero — and most models in use are @eq-cls-three-sources with
one or two parts dropped.

::: {#def-cls-covariance}
[Covariance models for a subject's observations]

Let subject \( i \) be observed at times \( t_{i1}<\dots<t_{im_i} \), and write
\( \bSigma_i \) for the covariance matrix of its response vector, with entries
\( \sigma_{jk} \) and time gap \( d_{jk}=\lvert t_{ij}-t_{ik}\rvert \). The standard
models are

| Model | \( \sigma_{jk} \) | Parameters |
|:---|:---|---:|
| independence | \( \sigma^2\delta_{jk} \) | \( 1 \) |
| compound symmetry | \( \sigma_u^2+\sigma^2\delta_{jk} \) | \( 2 \) |
| first-order autoregressive | \( \sigma^2\rho^{\lvert j-k\rvert} \) | \( 2 \) |
| exponential (spatial power) | \( \sigma^2\rho^{d_{jk}} \) | \( 2 \) |
| autoregressive plus random intercept | \( \sigma_u^2+\sigma^2\rho^{\lvert j-k\rvert} \) | \( 3 \) |
| Toeplitz of order \( q \) | \( \sigma^2c_{\lvert j-k\rvert} \), with \( c_0=1 \) and \( c_u=0 \) for \( u\ge q \) | \( q \) |
| random intercept and slope | \( (1,t_{ij})\bD(1,t_{ik})\T+\sigma^2\delta_{jk} \) | \( 4 \) |
| unstructured | \( \sigma_{jk} \), unrestricted | \( m(m+1)/2 \) |

Each model must also be nonnegative definite; the constraint is
\( \sigma_u^2,\sigma^2\ge0 \) and \( \lvert\rho\rvert<1 \) for the first five,
\( \bD \) nonnegative definite for the random-coefficient model, and a genuine
restriction on \( (c_1,\dots,c_{q-1}) \) for the Toeplitz
model (@exr-cls-toeplitz-pd).
:::

Compound symmetry is @def-cls-cluster, what a random intercept alone produces: it
makes the correlation between the first and last visit equal to that between
consecutive visits, which over five years it rarely is. The autoregressive model
is @def-het-ar1 within a subject and needs equally spaced occasions; the
exponential model is its continuous-time version. A random intercept added to
either puts a floor under the decaying correlation, and that two-source version of
@eq-cls-three-sources fits much real longitudinal data. The Toeplitz model lets
the decay take an arbitrary shape, truncated after \( q-1 \) lags. Unstructured
assumes nothing, and is right when the occasions are few, common and complete and
the subjects comfortably outnumber \( m(m+1)/2 \): the covariance implicit
in @prp-cls-multivariate.

::: {.warning}
The covariance model is a choice with consequences. It changes the standard errors
of the fixed effects, and — unlike in the balanced cases of the previous
sections — the estimates themselves whenever the design is unbalanced or the
regressors vary within subjects.
:::

## Choosing, and what goes wrong if you choose badly

::: {#prp-cls-selection}
[Working covariances, likelihood comparisons and empirical standard errors]

Let the \( N \) subjects be independent, with \( \E(\Y_i)=\X_i\bbeta \) and true
covariance \( \bSigma_i^0 \). Let \( \bSigma_i(\boldsymbol{\uptheta}) \) be a working
covariance model whose eigenvalues are bounded above and away from zero
uniformly in \( i \) and, for fixed \( \boldsymbol{\uptheta} \), let
\[
\tilde{\bbeta}(\boldsymbol{\uptheta})=\A(\boldsymbol{\uptheta})^{-1}
 \sum_i\X_i\T\bSigma_i(\boldsymbol{\uptheta})^{-1}\Y_i,\qquad
\A(\boldsymbol{\uptheta})=\sum_i\X_i\T\bSigma_i(\boldsymbol{\uptheta})^{-1}\X_i,
\]
assumed nonsingular.

::: {.enumerate options="label=(\alph*)"}
1. \( \tilde{\bbeta}(\boldsymbol{\uptheta}) \) is unbiased for every
   \( \boldsymbol{\uptheta} \), whether or not the working model is correct, with exact
   covariance \( \A^{-1}\B\A^{-1} \) where
   \( \B=\sum_i\X_i\T\bSigma_i^{-1}\bSigma_i^0\bSigma_i^{-1}\X_i \). If the working
   model is correct, \( \B=\A \) and the covariance is \( \A^{-1} \).

2. The **empirical** covariance estimate \( \hat{\bU}_{\text{E}} \) replaces
   \( \bSigma_i^0 \) by \( \br_i\br_i\T \), with
   \( \br_i=\Y_i-\X_i\tilde{\bbeta} \). Under the conditions (K1)–(K3) of
   [Section 33.2](02-clustered-data.html), applied with
   \( \bSigma_i(\boldsymbol{\uptheta})^{-1}\X_i \) in place of \( \X_i \), it
   satisfies \( n(\hat{\bU}_{\text{E}}-\A^{-1}\B\A^{-1})\to\mathbf{0} \) in
   probability as \( N\to\infty \), whatever \( \bSigma_i^0 \) is.

3. Let \( \mathbf{K} \) be any \( n\times(n-p) \) matrix whose columns form an orthonormal
   basis of \( \C(\X)\perpc \). The restricted log likelihood of a covariance model is
   the log likelihood of \( \mathbf{K}\T\Y \) (@def-mix-reml) and does not depend on the
   choice of \( \mathbf{K} \) (@thm-mix-reml). Hence two covariance models fitted with
   the *same* \( \X \) have restricted likelihoods that are densities of one and the
   same random vector, and their difference is a genuine log likelihood ratio. If
   the models differ in \( \X \), the two are densities of different random vectors
   and cannot be compared.
:::

:::

::: {.proof}
(a) \( \tilde{\bbeta}=\A^{-1}\sum_i\X_i\T\bSigma_i^{-1}\Y_i \) is a fixed linear
function of \( \Y \), so \( \E\tilde{\bbeta}=\A^{-1}\sum_i\X_i\T\bSigma_i^{-1}\X_i\bbeta
 =\bbeta \) and, by independence of the subjects,
\( \Cov(\tilde{\bbeta})
 =\A^{-1}\{\sum_i\X_i\T\bSigma_i^{-1}\bSigma_i^0\bSigma_i^{-1}\X_i\}\A^{-1} \)
by @thm-rv-linear. If \( \bSigma_i^0=\bSigma_i \), the middle sum telescopes to
\( \A \).

(b) Put \( \tilde{\X}_i=\bSigma_i^{-1}\X_i \), so that
\( \tilde{\bbeta}-\bbeta=(\sum_i\tilde{\X}_i\T\X_i)^{-1}\sum_i\tilde{\X}_i\T\be_i \)
and the empirical estimate is the same expression with \( \br_i\br_i\T \) in the
middle. The cluster-sandwich argument of
[Section 33.2](02-clustered-data.html) — that is, the proof of @thm-het-sandwich
with the cluster in place of the observation — goes through with
\( \tilde{\X}_g \) in the meat and \( \X_g \) in the bread: the assumed bounds on
the eigenvalues of \( \bSigma_i \) make \( \norm{\tilde{\X}_i} \) satisfy (K2)
whenever \( \norm{\X_i} \) does, and (K1) and (K3) are unchanged.

(c) The first two sentences are @def-mix-reml and @thm-mix-reml, and \( \mathbf{K} \)
is determined up to an orthogonal rotation by \( \C(\X) \). Two models with the same
\( \X \) therefore attach densities to the same \( (n-p) \)-dimensional random vector,
so their ratio is a likelihood ratio; two models with different \( \X \) attach them
to different random vectors, and the ratio is not even dimensionless when the two
\( \X \) have different ranks.
:::

Part (c) has a practical corollary: **fix the mean model first and compare
covariance models by restricted likelihood; compare mean models by ordinary
likelihood, or by a Wald or \( F \) test, holding the covariance fixed.** The
criteria of [Chapter 29](../ch29-model-selection/index.html) then apply to the
restricted likelihood as to any likelihood (@def-sel-aic-bic):
\[
\text{AIC}_R=-2\ell_R+2q,\qquad \text{BIC}_R=-2\ell_R+q\log(n-p),
\]{#eq-cls-reml-aic}

with \( q \) the number of covariance parameters and \( n-p \) the number of error
contrasts, the dimension of the data whose likelihood \( \ell_R \) is.

Part (a) is why all this is less dangerous than it sounds: a wrong covariance
model costs efficiency, not validity, since \( \bbeta \) stays unbiased and part
(b) supplies a standard error that is right anyway — the longitudinal version of
the sandwich argument of
[Section 21.4](../ch21-nonnormality-heteroscedasticity-serial/04-sandwich-estimators.html),
and the basis of the generalized estimating equations of Chapter 40.

::: {.warning}
The empirical standard error protects against the covariance model and nothing
else. If \( \E(\Y_i)\ne\X_i\bbeta \), part (a) fails at the first line and there
is nothing left to be robust about; it does not protect against missing data that
depend on the unobserved responses; and it needs \( N \) large in the sense of
(K1)–(K3), being worse than the model-based estimate with a dozen subjects.
:::

## Checking the fit of a covariance model

Two displays do most of the work, both built from residuals
\( \hat{r}_{ij}=y_{ij}-\x_{(ij)}\T\hbeta \) after a fit of the mean model. The
**empirical variogram** plots, against the time gap \( u \), the average of
\( \tfrac12(\hat{r}_{ij}-\hat{r}_{ik})^2 \) over same-subject pairs at that gap. Its
population version \( \gamma(u)=\tfrac12\Var(r_{ij}-r_{ik}) \) is, for
@eq-cls-three-sources, \( \sigma_e^2+\sigma_s^2\{1-\rho(u)\} \), rising from
\( \sigma_e^2 \) to \( \sigma_e^2+\sigma_s^2 \); the gap between that plateau and
the total residual variance \( \sigma_u^2+\sigma_s^2+\sigma_e^2 \) measures
\( \sigma_u^2 \), the part of the correlation that never decays. The **residual
correlogram** plots the average correlation of residuals at gap \( u \): compound
symmetry predicts a horizontal line, the autoregressive model a geometric decay to
zero, the combination a decay to a positive floor. Comparing empirical with fitted
curves shows quickly which model does violence to the data.

::: {#exm-cls-grunfeld-covariance}
[Investment in eleven firms]

Grunfeld's panel records gross investment, market value and capital stock for
\( 11 \) US firms in each of the \( 20 \) years from 1935 to 1954, a complete
rectangle with \( n=220 \). Take log investment as the response and the mean model
\( \E(y)=\beta_0+\beta_1t+\beta_2\log(\text{value}) \), with \( t \) in decades from
the middle of the record, so \( p=3 \). An unstructured covariance would need
\( 210 \) parameters from \( 11 \) subjects, which is hopeless. Fitting the
alternatives by restricted maximum likelihood gives

| Covariance model | \( q \) | \( -2\ell_R \) | \( \text{AIC}_R \) | \( \text{BIC}_R \) |
|:---|---:|---:|---:|---:|
| independence | 1 | 476.95 | 478.95 | 482.33 |
| compound symmetry | 2 | 114.83 | 118.83 | 125.59 |
| autoregressive | 2 | 49.89 | 53.89 | 60.65 |
| autoregressive + random intercept | 3 | 29.87 | 35.87 | 46.01 |
| random intercept and slope | 4 | 98.66 | 106.66 | 120.18 |
| Toeplitz of order 4 | 4 | 147.97 | 155.97 | 169.49 |

Independence is hopeless: the firms differ enormously in size, and log investment
is nearly constant within a firm relative to the spread between firms. Compound
symmetry captures that and improves \( -2\ell_R \) by \( 362 \), the
autoregressive model on the same two parameters does better still, and the
combination is best on both criteria, by \( 82.96 \) of
\( \text{AIC}_R \) over compound symmetry. The random-coefficient model, which
gives each firm its own trend, describes these data worse than a persistent level
plus decaying serial correlation.

[Figure 33.4.1](#fig-cls-covariance) shows why. The residual correlations are
\( 0.953 \) at one year, \( 0.918 \) at two,
\( 0.884 \) at five and \( 0.839 \) at ten: high
everywhere, because the firm's level dominates, but visibly decaying, which
compound symmetry cannot reproduce. The variogram rises to a plateau of about
\( 0.099 \), far below the total residual variance of
\( 0.4839 \); the gap is the between-firm variance. The selected model
is not innocent either: its fitted correlogram flattens onto its own floor after
three or four years while the data go on decaying. That is what the diagnostic is
for, and a reason to read the empirical standard error of @prp-cls-selection(b)
alongside the model-based one.

With the selected model the coefficient of \( \log(\text{value}) \) is
\( 0.6295 \), with a model-based standard error of
\( 0.0712 \) and an empirical one of \( 0.0608 \), and the
coefficient of \( t \) is \( 0.4261 \) with standard errors
\( 0.0600 \) and \( 0.0601 \). Ordinary least squares,
treating the \( 220 \) observations as independent, reports
\( 0.9833 \) and \( 0.3768 \) with standard errors
\( 0.0351 \) and \( 0.0821 \). The two disagree substantially
about the elasticity with respect to market value, because least squares reads it
mostly from differences between firms while generalized least squares, having
absorbed the firm levels into the covariance, reads it from movements within
firms — the between-versus-within distinction of
[Section 6.6](../ch06-projections/06-fwl.html) again, and a question about the
subject matter.
:::

```{.python .run #cell-covariance-models-data}
import numpy as np
import statsmodels.api as sm

df = sm.datasets.grunfeld.load_pandas().data
firms = sorted(df["firm"].unique())
years = np.sort(df["year"].unique())
N, m = len(firms), len(years)                        # 11 firms, 20 years
t = (years - years.mean()) / 10                      # centred and scaled time

Ylist, Xlist = [], []
for f in firms:
    d = df[df["firm"] == f].sort_values("year")
    Ylist.append(np.log(d["invest"].to_numpy()))
    Xlist.append(np.column_stack([np.ones(m), t, np.log(d["value"].to_numpy())]))
p = Xlist[0].shape[1]
n = N * m
```

With a common \( \bSigma \), \( -2\ell_R \) is
\( N\log\lvert\bSigma\rvert+\sum_i\br_i\T\bSigma^{-1}\br_i
 +\log\lvert\sum_i\X_i\T\bSigma^{-1}\X_i\rvert \) plus a constant, evaluated at
the generalized least squares \( \tilde{\bbeta} \), the third term being what
distinguishes it from the ordinary likelihood (@thm-mix-reml).

```{.python .run #cell-covariance-models-reml}
from scipy.optimize import minimize

def reml(Sigma, Ylist, Xlist):
    """-2 x the restricted log likelihood, and the GLS estimate, for a common Sigma."""
    L = np.linalg.cholesky(Sigma)
    logdet = 2 * np.sum(np.log(np.diag(L)))
    solve = lambda B: np.linalg.solve(Sigma, B)
    XtWX = sum(X.T @ solve(X) for X in Xlist)
    XtWy = sum(X.T @ solve(y) for X, y in zip(Xlist, Ylist))
    beta = np.linalg.solve(XtWX, XtWy)
    quad = sum((y - X @ beta) @ solve(y - X @ beta) for X, y in zip(Xlist, Ylist))
    s, ld_info = np.linalg.slogdet(XtWX)
    n, p = len(Ylist) * Sigma.shape[0], len(beta)
    return len(Ylist) * logdet + quad + ld_info + (n - p) * np.log(2 * np.pi), beta

lag = np.abs(np.subtract.outer(np.arange(m), np.arange(m)))
Zrc = np.column_stack([np.ones(m), t])               # random intercept and slope basis

def build(name, th):
    """Sigma for each covariance model, from an unconstrained parameter vector."""
    if name == "independence":
        return np.exp(th[0]) * np.eye(m)
    if name == "compound symmetry":
        return np.exp(th[0]) * np.eye(m) + np.exp(th[1]) * np.ones((m, m))
    if name == "AR(1)":
        return np.exp(th[0]) * np.tanh(th[1])**lag
    if name == "AR(1) + random intercept":
        return (np.exp(th[0]) * np.tanh(th[1])**lag + np.exp(th[2]) * np.ones((m, m)))
    if name == "random intercept and slope":
        D = np.array([[np.exp(th[1]), th[3]], [th[3], np.exp(th[2])]])
        return np.exp(th[0]) * np.eye(m) + Zrc @ D @ Zrc.T
    if name == "Toeplitz(4)":
        a = np.concatenate([[1.0], th[1:4]])         # a moving average of order 3 ...
        g = np.array([a[:4 - u] @ a[u:] for u in range(4)])
        c = np.concatenate([g / g[0], np.zeros(m - 4)])   # ... so the band is nonneg. definite
        return np.exp(th[0]) * c[lag]
    raise ValueError(name)

MODELS = {"independence": 1, "compound symmetry": 2, "AR(1)": 2,
          "AR(1) + random intercept": 3, "random intercept and slope": 4,
          "Toeplitz(4)": 4}

def fit(name, start):
    def obj(th):
        S = build(name, th)
        if np.min(np.linalg.eigvalsh(S)) <= 1e-9:
            return 1e6
        return reml(S, Ylist, Xlist)[0]
    best = min((minimize(obj, s, method="Nelder-Mead",
                         options={"maxiter": 20000, "xatol": 1e-9, "fatol": 1e-9})
                for s in start), key=lambda r: r.fun)
    S = build(name, best.x)
    d2l, beta = reml(S, Ylist, Xlist)
    q = MODELS[name]
    return {"m2ll": d2l, "aic": d2l + 2 * q, "bic": d2l + q * np.log(n - p),
            "q": q, "Sigma": S, "beta": beta}

starts = {"independence": [[-1.0]],
          "compound symmetry": [[-1.0, -1.0]],
          "AR(1)": [[-1.0, 1.0]],
          "AR(1) + random intercept": [[-1.5, 1.0, -1.5]],
          "random intercept and slope": [[-2.0, -1.0, -2.0, 0.0], [-3.0, -2.0, -3.0, 0.0]],
          "Toeplitz(4)": [[-1.0, 1.5, 1.0, 0.5], [-1.0, 2.5, 2.0, 1.0]]}
res = {name: fit(name, starts[name]) for name in MODELS}
for name in MODELS:
    r = res[name]
    print(f"{name:28s} q={r['q']}  -2logL_R={r['m2ll']:9.3f}  "
          f"AIC={r['aic']:9.3f}  BIC={r['bic']:9.3f}")
```

::: {when-format="html"}
![**Figure 33.4.1.** Checking the covariance model for the investment panel.
(a) The empirical variogram of the residuals, with those implied by two fitted
models; the gap between the grey total-variance line and the plateau is the
between-firm variance. (b) The average residual correlation at each gap, with
three fitted models.](covariance_models.svg){#fig-cls-covariance width=100%}
:::

::: {when-format="pdf"}
![Checking the covariance model for the investment panel.
(a) The empirical variogram of the residuals, with those implied by two fitted
models; the gap between the grey total-variance line and the plateau is the
between-firm variance. (b) The average residual correlation at each gap, with
three fitted models.](covariance_models.pdf){width=100%}
:::

```{.python .run #cell-covariance-models-standard-errors}
best_name = min(MODELS, key=lambda k: res[k]["aic"])
Xall = np.vstack(Xlist)
yall = np.concatenate(Ylist)
Sbest = res[best_name]["Sigma"]
A = np.linalg.inv(sum(X.T @ np.linalg.solve(Sbest, X) for X in Xlist))
beta = res[best_name]["beta"]
resid = [y - X @ beta for X, y in zip(Xlist, Ylist)]
meat = sum(np.outer(X.T @ np.linalg.solve(Sbest, r), X.T @ np.linalg.solve(Sbest, r))
           for X, r in zip(Xlist, resid))
cov_robust = A @ meat @ A * N / (N - 1)              # empirical ("sandwich") covariance
XtXinv = np.linalg.inv(Xall.T @ Xall)
beta_ols = XtXinv @ Xall.T @ yall
s2_ols = np.sum((yall - Xall @ beta_ols)**2) / (n - p)
for j, lab in enumerate(["intercept", "time/10", "log value"]):
    print(f"{lab:10s} GLS {beta[j]: .4f}  model se {np.sqrt(A[j, j]):.4f}  "
          f"empirical se {np.sqrt(cov_robust[j, j]):.4f}  "
          f"| OLS {beta_ols[j]: .4f} se {np.sqrt(s2_ols * XtXinv[j, j]):.4f}")
```

## Missing data

Longitudinal data are rarely complete. The mixed-model machinery handles ragged
data without complaint — @prp-cls-selection never assumed a common \( m_i \) — an
advantage over the classical analysis of
[Section 33.3](03-repeated-measures.html), which needs a complete rectangle. But
"handles" means "computes". If the probability of a missing value depends only on
the covariates and the responses that *were* observed, the observed-data
likelihood is the right one: the data are missing at random. If it depends on the
value that would have been observed — dropout because a patient is not improving —
the likelihood is wrong and no covariance model repairs it.

::: {.warning}
Software that fits a mixed model to incomplete longitudinal data makes a
missing-at-random assumption on the analyst's behalf, silently, and the empirical
standard error of @prp-cls-selection(b) does *not* protect against its failure.
Chapter 41 treats missing data properly.
:::

## Exercises

### A. Check your understanding

::: {#exr-cls-model-count}
[A1]

A trial measures \( 60 \) patients at \( 8 \) visits. How many parameters does each
model of @def-cls-covariance need? Which could not be fitted with only \( 10 \)
patients, and why?
:::

::: {#exr-cls-variogram-shapes}
[A2]

Sketch the population variogram \( \gamma(u) \) for compound symmetry, for the
autoregressive model and for their combination, marking the total variance on
each. Which feature identifies \( \sigma_u^2 \), and which identifies measurement
error?
:::

::: {.solution}
Compound symmetry gives \( \gamma(u)=\sigma^2 \) for every \( u>0 \), flat and well
below the total variance \( \sigma^2+\sigma_u^2 \). The autoregressive model gives
\( \gamma(u)=\sigma^2(1-\rho^u) \), rising from \( 0 \) to its total variance
\( \sigma^2 \). The combination rises from \( 0 \) to \( \sigma^2 \) and stops,
leaving the gap \( \sigma_u^2 \) to the total. The gap between plateau and total
variance is \( \sigma_u^2 \); the intercept of the variogram as \( u\to0 \) is the
measurement-error variance.
:::

### B. Practice

::: {#exr-cls-reml-comparison}
[B1]

Two analysts fit the investment data of @exm-cls-grunfeld-covariance. The first
uses the mean model of the example with an autoregressive covariance; the second
adds a quadratic term in \( t \) and uses compound symmetry. They compare
\( \text{AIC}_R \). Explain what is wrong, and describe a correct sequence of
comparisons.
:::

::: {.solution}
By @prp-cls-selection(c) the two restricted likelihoods are densities of different
random vectors, of dimensions \( n-3 \) and \( n-4 \), so their difference is not a
log likelihood ratio and \( \text{AIC}_R \) cannot compare them. Correct sequence:
fix a rich mean model and compare covariance models by \( \text{AIC}_R \); then,
with the chosen covariance fixed, compare mean models by a Wald or \( F \) test, or
by AIC from the *full* likelihood; then refit the covariance.
:::

::: {#exr-cls-toeplitz-pd}
[B2]

Show that the Toeplitz model of @def-cls-covariance with \( q=3 \) is nonnegative
definite for every \( m \) iff \( (1,c_1,c_2,0,\dots) \) is the autocorrelation
sequence of a moving-average process of order two, and explain why the
parameterization \( c_u\propto\sum_la_la_{l+u} \) used in the listing above satisfies
this automatically.
:::

### C. Going deeper

::: {#exr-cls-reml-boundary}
[C1]

The autoregressive model is nested in the autoregressive-plus-random-intercept
model at \( \sigma_u^2=0 \), a boundary of the parameter space. Explain why the
likelihood ratio statistic comparing the two does not have a \( \chi^2(1) \) null
distribution, state what it does have (@prp-mix-inference), and say what this
means for reading \( \text{AIC}_R \) differences as informal evidence.
:::
