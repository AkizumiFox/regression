# Inference when the covariance is estimated

Every formula in the last three sections assumed \( \V \) known. In practice
\( \V(\hat{\boldsymbol{\uptheta}}) \) is substituted, and the substitution is not
free: it perturbs the estimator, it perturbs the estimator's covariance, and it
destroys the exactness of the \( F \) and \( t \) distributions Part III relied on.
Testing a variance component raises a separate difficulty, since the null value sits
on the boundary of the parameter space and the usual chi-square calibration is simply
wrong; @thm-mix-boundary computes the right one exactly in the balanced layout.

## Feasible generalized least squares is still unbiased

::: {#prp-mix-inference}
[Fixed-effect inference with an estimated covariance]

Consider the normal linear mixed model with \( \rank(\X)=p \), and write
\[
\hbeta(\boldsymbol{\uptheta})=(\X\T\V^{-1}\X)^{-1}\X\T\V^{-1}\Y
\]
for the generalized least squares estimator.

::: {.enumerate options="label=(\alph*)"}
1. *(Unbiasedness.)* Let \( \hat{\boldsymbol{\uptheta}}
   =\hat{\boldsymbol{\uptheta}}(\Y) \) be **translation invariant**, meaning
   \( \hat{\boldsymbol{\uptheta}}(\y+\X\bb)=\hat{\boldsymbol{\uptheta}}(\y) \) for all
   \( \bb \), and **even**, meaning
   \( \hat{\boldsymbol{\uptheta}}(-\y)=\hat{\boldsymbol{\uptheta}}(\y) \). If
   \( \E\norm{\hbeta(\hat{\boldsymbol{\uptheta}})}<\infty \), then
   \( \E\hbeta(\hat{\boldsymbol{\uptheta}})=\bbeta \). Both REML and maximum
   likelihood estimators of \( \boldsymbol{\uptheta} \) have these two properties.

2. *(A design where everything is exact.)* Consider a balanced cluster-randomized
   trial: \( g \) clusters, \( g/2 \) assigned to each arm, \( m \) units per cluster,
   with a random cluster effect of variance \( \sigma_a^2 \) and independent errors of
   variance \( \sigma^2 \). Let \( \hat\beta \) be the difference of the two arms'
   means of cluster means, and let \( \text{MS}_C \) be the mean square between
   clusters within arms, on \( g-2 \) degrees of freedom, and \( \text{MSE} \) the
   mean square within clusters. Then, with \( \lambda_1=\sigma^2+m\sigma_a^2 \),
   \( \Var(\hat\beta)=4\lambda_1/(gm) \) and
   \[
   \frac{\hat\beta-\beta}{\sqrt{4\,\text{MS}_C/(gm)}}\sim t(g-2)
   \]{#eq-mix-cluster-t}

   exactly. If \( \text{MS}_C\ge\text{MSE} \) then \( \text{MS}_C \) is also the REML
   estimate of \( \lambda_1 \); otherwise the REML estimate lies on the boundary
   \( \hat\sigma_a^2=0 \), as in @prp-mix-oneway-reml(c).

3. *(Satterthwaite's approximation.)* Let \( U=\sum_ic_i\text{MS}_i>0 \), where the
   \( \text{MS}_i \) are independent and
   \( \nu_i\text{MS}_i/\E\text{MS}_i\sim\chi^2(\nu_i) \). The two-moment match of
   [Section 4.3](../ch04-quadratic-forms/03-chisq.html), fitting
   \( \nu U/\E U\sim\chi^2(\nu) \) to the mean and variance of \( U \), has
   \[
   \nu=\frac{(\E U)^2}{\sum_i c_i^2(\E\text{MS}_i)^2/\nu_i},
   \]{#eq-mix-satterthwaite}

   estimated in practice by replacing \( \E\text{MS}_i \) by \( \text{MS}_i \).
:::

The companion statement for the variance components themselves, where the null value
sits on the boundary of the parameter space, is @thm-mix-boundary below.
:::

::: {.proof}
(a) Write \( \be=\Y-\X\bbeta \) and
\( \A(\boldsymbol{\uptheta})=(\X\T\V^{-1}\X)^{-1}\X\T\V^{-1} \), so that
\( \A(\boldsymbol{\uptheta})\X=\I \) and
\( \hbeta(\hat{\boldsymbol{\uptheta}})-\bbeta
=\A(\hat{\boldsymbol{\uptheta}}(\Y))\be \). Translation invariance gives
\( \hat{\boldsymbol{\uptheta}}(\Y)=\hat{\boldsymbol{\uptheta}}(\X\bbeta+\be)
=\hat{\boldsymbol{\uptheta}}(\be) \), so
\( \hbeta(\hat{\boldsymbol{\uptheta}})-\bbeta=\mathbf{h}(\be) \) with
\( \mathbf{h}(\mathbf{e})=\A(\hat{\boldsymbol{\uptheta}}(\mathbf{e}))\mathbf{e} \).
Evenness gives \( \mathbf{h}(-\mathbf{e})=-\mathbf{h}(\mathbf{e}) \), and \( \be \)
has the same distribution as \( -\be \), so \( \E\mathbf{h}(\be)=-\E\mathbf{h}(\be) \)
and the expectation is \( \bzero \). This is the argument of @thm-het-wls(c) with the
weights replaced by a covariance. REML estimators are functions of
\( \mathbf{K}\T\y \) alone, hence translation invariant, and @eq-mix-reml-loglik is
unchanged by \( \y\mapsto-\y \), hence even; the same holds for maximum likelihood,
whose profile @eq-mix-profile depends on \( \y \) only through \( \y\T\bP\y \) and is
translation invariant because \( \bP\X=\bzero \).

(b) Both columns of \( \X=[\bone,\mathbf{T}] \) are constant within clusters, so by
@exm-mix-ols-blue ordinary and generalized least squares agree and
\( \hat\beta \) is the difference of arm means of cluster means. The \( g \) cluster
means \( \bar Y_k \) are independent \( \Normal(\beta_0+\beta T_k,\lambda_1/m) \),
because \( \Var(\bar Y_k)=\sigma_a^2+\sigma^2/m=\lambda_1/m \). Applying the exact normal theory of
[Chapter 11](../ch11-general-linear-hypothesis/index.html) to these \( g \)
independent normal observations gives
\( \hat\beta\sim\Normal(\beta,4\lambda_1/(gm)) \), since the two arm means are
independent with variance \( \lambda_1/(m\cdot g/2) \) each. The quadratic form
\( m\sum_k(\bar Y_k-\bar Y_{\text{arm}(k)})^2/\lambda_1=(g-2)\text{MS}_C/\lambda_1 \)
is the residual sum of squares of a two-group model fitted to the \( g \) cluster
means, so it is \( \chi^2(g-2) \) and is independent of the fitted arm means
by @thm-qf-indep-linear. Dividing gives @eq-mix-cluster-t.

For the REML claim, the error contrast space \( \C(\X)\perpc \) decomposes
orthogonally into \( \C(\Z)\cap\C(\X)\perpc \), of dimension \( g-2 \), on which
\( \V \) acts as \( \lambda_1\I \), and \( \C(\Z)\perpc \), of dimension \( n-g \),
on which \( \V \) acts as \( \sigma^2\I \). Choosing \( \mathbf{K} \) adapted to this
decomposition, @eq-mix-reml-loglik becomes
\( -2\ell_R=(g-2)\log\lambda_1+(n-g)\log\sigma^2+\text{SS}_C/\lambda_1
+\text{SSE}/\sigma^2 \) up to constants, to be minimized over
\( \lambda_1\ge\sigma^2>0 \), exactly as in @prp-mix-oneway-reml. Its minimizers are
\( \hat\lambda_1=\text{MS}_C \) and \( \hat\sigma^2=\text{MSE} \) when
\( \text{MS}_C\ge\text{MSE} \); otherwise the constraint binds, by the argument
of @prp-mix-oneway-reml(c).

(c) The two-moment match of [Section 4.3](../ch04-quadratic-forms/03-chisq.html)
gives \( \nu=2(\E U)^2/\Var(U) \), and here
\( \Var(U)=\sum_ic_i^2\Var(\text{MS}_i)=2\sum_ic_i^2(\E\text{MS}_i)^2/\nu_i \),
which is @eq-mix-satterthwaite — the calculation already made for a linear combination
of independent variance estimates in
[Section 15.5](../ch15-anova-subspaces/05-balance.html).
:::

Part (b) isolates the problem. The standard error a mixed-model program reports for
this design is *exactly right*: it is \( \sqrt{4\,\text{MS}_C/(gm)} \) whenever the
REML estimate of \( \lambda_1 \) is interior, that is, unless
\( \text{MS}_C<\text{MSE} \) — probability \( 0.157 \) in the
design simulated below, and \( 0.558 \) when \( \sigma_a^2 \) is
in fact zero. What is wrong is the reference distribution. Treating
\( \hat\beta/\text{se} \) as standard normal,
which is what a Wald test does, ignores that \( \text{MS}_C \) came from \( g-2 \)
degrees of freedom. With \( g=8 \) clusters that is six, and the
\( t \) critical value \( 2.447 \) exceeds
\( 1.96 \) by \( 25 \) per cent.

::: {#exm-mix-cluster-coverage}
[Eight clusters]

Simulate the trial with \( g=8 \) clusters of \( m=5 \), an
intraclass correlation of \( 0.2 \), and no treatment effect. Over
\( 200\,000 \) trials, the nominal 95% normal Wald interval covers
\( 0.901 \) of the time, while the
\( t(g-2) \) interval covers \( 0.949 \).
[Figure 32.5.1](05-inference.html#fig-mix-wald) traces the coverage against the number
of clusters: \( 0.812 \) with four clusters, still
\( 0.943 \) with forty. The normal interval's coverage is exactly
\( 2F_{t(g-2)}(1.96)-1 \), which is why it approaches but never reaches
\( 0.95 \).

The coverage study is the null one; take now a single data set simulated from the
same design but with a true effect of \( 0.8 \). On it the estimated effect is
\( 1.154 \) with
standard error \( 0.479 \). The normal interval is
\( 1.154\pm0.940 \), which excludes zero; the exact
interval is \( 1.154\pm1.173 \), which does not. The
whole difference between "significant" and "not" is the reference distribution.
:::

::: {when-format="html"}
![**Figure 32.5.1.** Coverage of a nominal 95% interval for the treatment effect in a
balanced cluster-randomized trial with \( m=5 \) units per cluster and intraclass
correlation \( 0.2 \), against the number of clusters. The normal Wald interval treats
the estimated covariance as known; the \( t(g-2) \) interval is exact
by @eq-mix-cluster-t.](wald_coverage.svg){#fig-mix-wald width=62%}
:::

::: {when-format="pdf"}
![Coverage of a nominal 95% interval for the treatment effect in a balanced
cluster-randomized trial with \( m=5 \) units per cluster and intraclass correlation
\( 0.2 \), against the number of clusters. The normal Wald interval treats the
estimated covariance as known; the \( t(g-2) \) interval is exact
by @eq-mix-cluster-t.](wald_coverage.pdf){width=62%}
:::

```{.python .run #cell-fixed-effects-trial}
import numpy as np
from scipy import stats

def simulate(g, m, sigma_a, sigma, effect, B, seed):
    """B cluster-randomized trials; returns the estimate and its estimated variance."""
    rng = np.random.default_rng(seed)
    treat = np.tile([0.0, 1.0], g // 2)                 # cluster k is treated if treat[k]=1
    Y = (effect * treat[:, None] + sigma_a * rng.standard_normal((B, g, 1))
         + sigma * rng.standard_normal((B, g, m)))
    cluster_mean = Y.mean(axis=2)
    diff = (cluster_mean[:, treat == 1].mean(axis=1)
            - cluster_mean[:, treat == 0].mean(axis=1))
    # mean square between clusters within arms, on g - 2 degrees of freedom
    centred = cluster_mean - np.where(treat == 1,
                                      cluster_mean[:, treat == 1].mean(axis=1, keepdims=True),
                                      cluster_mean[:, treat == 0].mean(axis=1, keepdims=True))
    lam1 = m * np.sum(centred**2, axis=1) / (g - 2)     # estimates sigma^2 + m sigma_a^2
    return diff, 4 * lam1 / (g * m)

g, m = 8, 5
# a shorter run than the 200 000 trials quoted in the text, so the cell finishes quickly
diff, var_hat = simulate(g, m, 0.5, 1.0, 0.0, 20000, seed=5150)
t = diff / np.sqrt(var_hat)
print(f"normal interval covers {np.mean(np.abs(t) <= stats.norm.ppf(0.975)):.3f}")
print(f"t(g-2) interval covers {np.mean(np.abs(t) <= stats.t.ppf(0.975, g - 2)):.3f}")
```

## Degrees of freedom in general

Outside designs as tidy as @prp-mix-inference(b) there is no exact answer, and two
approximations are in common use. Neither is proved here; both rest on simulation
evidence rather than exact distribution theory.

**Satterthwaite.** The estimated variance
\( \blambda\T(\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\X)^{-1}\blambda \) is a
smooth function of \( \hat{\boldsymbol{\uptheta}} \). Approximating its distribution
by a scaled chi-square and matching two moments, with the variance from the delta
method and the asymptotic covariance of \( \hat{\boldsymbol{\uptheta}} \), gives a
denominator degrees of freedom that is in general not an integer. Giesbrecht and
Burns (1985) developed this for mixed models; in the design of
@prp-mix-inference(b) it returns the exact \( g-2 \) (@exr-mix-satterthwaite-exact).

**Kenward and Roger.** Satterthwaite corrects the reference distribution but not the
variance estimate, and the variance estimate is biased downwards: the plug-in
covariance \( (\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\X)^{-1} \) is what
\( \hbeta \) would have if \( \hat{\boldsymbol{\uptheta}} \) were the truth, and the
estimator actually used varies for a second reason, because
\( \hat{\boldsymbol{\uptheta}} \) itself varies, adding a nonnegative definite term.
Kackar and Harville (1984) approximated it; Kenward and Roger (1997) add an estimate
of it, then scale the \( F \) statistic and match two moments for its degrees of
freedom. The correction is the default in several packages and performs well in small
samples; @exm-mix-gibbs exhibits the inflation.

::: {.warning}
Unlike @prp-mix-inference(a), these are approximations with no finite-sample
guarantee. When the design is small, unbalanced and the covariance structure rich,
the safest calibration is a parametric bootstrap: simulate from the fitted model,
refit, and use the simulated distribution of the statistic (@prp-bs-parametric and
[Section 23.3](../ch23-resampling-inference/03-confidence-intervals.html)). It costs
a few hundred fits and needs no asymptotics.
:::

## Comparing models

Two rules govern likelihood comparisons in mixed models, and both follow from
[Section 32.4](04-likelihood.html).

**Different fixed effects: use maximum likelihood.** Since \( \mathbf{K} \) is
determined by \( \C(\X) \), two fixed-effects structures give restricted likelihoods
of different random vectors, and their difference means nothing. Likelihood ratio
tests and information criteria across fixed-effects structures must use the ordinary
likelihood @eq-mix-loglik.

**Same fixed effects, different covariance: REML is preferred.** Here \( \mathbf{K} \)
is the same, the restricted likelihoods are likelihoods of the same data, and the
additive constants of @thm-mix-reml(c) cancel. Restricted likelihood ratio tests, and
the AIC and BIC of @def-sel-aic-bic computed from \( \ell_R \) with \( d_k \) the
number of covariance parameters, are legitimate. For BIC the sample size to use is
\( n-p \), the dimension of the data actually modelled, not \( n \), following
Verbeke and Molenberghs (2000).

::: {#exm-mix-grunfeld-slopes}
[Do the firms' slopes differ?]

Continue @exm-mix-grunfeld. Adding a firm-specific slope on market value, with
variance \( \sigma_s^2 \) independent of the firm intercept, raises the restricted
log-likelihood from \( -1176.53 \) to
\( -1159.01 \). The two models have the same fixed
effects, so the comparison is legitimate: the statistic is
\( 35.05 \) on one boundary parameter, far beyond any reasonable
critical value, and AIC falls from \( 2357.1 \) to
\( 2324.0 \), BIC from \( 2363.8 \)
to \( 2334.2 \). The firms' responses to market value genuinely
differ: the estimated standard deviation of the firm-specific slopes is
\( 0.1492 \).

The consequence for inference is dramatic. In the random-intercept model the
coefficient of market value is \( 0.1094 \) with standard error
\( 0.0100 \); in the random-slope model it is
\( 0.0127 \) with standard error
\( 0.0562 \), five and a half times larger. Nothing about the
data changed. Two things did. The standard error grew because the effective sample
size for a *population-average* slope is \( 11 \) firms, not
\( 220 \) firm-years. The estimate moved because the two models
estimate different functionals. Firm-specific intercepts and slopes with a common
capital coefficient give eleven slopes of mean
\( -0.0009 \) and standard deviation
\( 0.198 \): their spread dwarfs their average. The random-slope
model estimates that average; the random-intercept model, forced to give every firm
the same slope, reports a precision-weighted blend instead. Panel (a) of
[Figure 32.5.2](05-inference.html#fig-mix-grunfeld) shows the four analyses side by
side; panel (b) shows that with \( 20 \) years per firm the predicted
firm effects are indistinguishable from the firm-indicator estimates.
:::

::: {when-format="html"}
![**Figure 32.5.2.** Grunfeld's panel. (a) The coefficient of market value with 95%
intervals under four analyses: one common intercept, firm indicators, a random firm
intercept, and random firm intercepts and slopes. (b) Firm effects from the
firm-indicator fit against the predicted firm effects of the random-intercept model;
with twenty years per firm there is nothing left to
shrink.](grunfeld_mixed.svg){#fig-mix-grunfeld width=100%}
:::

::: {when-format="pdf"}
![Grunfeld's panel. (a) The coefficient of market value with 95% intervals under four
analyses: one common intercept, firm indicators, a random firm intercept, and random
firm intercepts and slopes. (b) Firm effects from the firm-indicator fit against the
predicted firm effects of the random-intercept model; with twenty years per firm there
is nothing left to shrink.](grunfeld_mixed.pdf){width=100%}
:::

```{.python .run #cell-grunfeld-mixed-reml}
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import optimize

df = sm.datasets.grunfeld.load_pandas().data.sort_values(["firm", "year"])
firms = df["firm"].unique()
y = df["invest"].to_numpy()
X = np.column_stack([np.ones(len(y)), df["value"], df["capital"]])
blocks = [np.where(df["firm"].to_numpy() == f)[0] for f in firms]
n, p, g = len(y), X.shape[1], len(firms)

def restricted_loglik(theta, Zcols):
    """Restricted log-likelihood; theta holds log sigma^2 and the log variances of Zcols."""
    s2, dvar = np.exp(theta[0]), np.exp(theta[1:])
    logdet, W, XtVy, yVy = 0.0, np.zeros((p, p)), np.zeros(p), 0.0
    for ix in blocks:
        Zi = Zcols[ix]
        Vi = (Zi * dvar) @ Zi.T + s2 * np.eye(len(ix))
        Vinv = np.linalg.inv(Vi)
        logdet += np.linalg.slogdet(Vi)[1]
        W += X[ix].T @ Vinv @ X[ix]
        XtVy += X[ix].T @ Vinv @ y[ix]
        yVy += y[ix] @ Vinv @ y[ix]
    beta = np.linalg.solve(W, XtVy)
    return -0.5 * (logdet + yVy - XtVy @ beta + np.linalg.slogdet(W)[1]
                   - np.linalg.slogdet(X.T @ X)[1] + (n - p) * np.log(2 * np.pi)), beta, W

Z_int = np.ones((n, 1))                                   # random intercept only
fit = optimize.minimize(lambda t: -restricted_loglik(t, Z_int)[0],
                        np.log([2000.0, 5000.0]), method="Nelder-Mead",
                        options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 5000})
ll_int, beta_int, W_int = restricted_loglik(fit.x, Z_int)
s2_hat, s2_firm = np.exp(fit.x)
se_int = np.sqrt(np.diag(np.linalg.inv(W_int)))

print(f"sigma^2 {s2_hat:.1f}   firm variance {s2_firm:.1f}")
print(f"value   {beta_int[1]:.4f} ({se_int[1]:.4f})")
print(f"capital {beta_int[2]:.4f} ({se_int[2]:.4f})")
```

## Testing a variance component

The hypothesis \( \sigma_a^2=0 \) sits on the edge of the parameter space, so the
usual chi-square calibration of a likelihood ratio statistic does not apply: the
standard argument expands the log-likelihood about an *interior* maximum, and here the
unconstrained maximum falls outside the parameter space about half the time. In the
balanced one-way model the correct answer can be computed exactly.

::: {#thm-mix-boundary}
[Exact null distribution of the restricted likelihood ratio]

In the balanced normal one-way random effects model with \( n=gm \), let
\[
L=2\bigl\{\sup_{\sigma_a^2\ge0}\ell_R-\sup_{\sigma_a^2=0}\ell_R\bigr\}
\]
be the restricted likelihood ratio statistic for \( H_0:\sigma_a^2=0 \), and write
\( F=\text{MSB}/\text{MSE} \), \( a=(g-1)/(n-1) \), \( b=(n-g)/(n-1) \). Then
\[
L=\begin{cases}
0,& F\le1,\\[2pt]
(n-1)\log(aF+b)-(g-1)\log F,& F>1 .
\end{cases}
\]{#eq-mix-ratio-exact}

Under \( H_0 \), \( F\sim F(g-1,n-g) \), so the null distribution of \( L \) is free
of \( \sigma^2 \) and exactly computable. In particular \( L \) is a nondecreasing
function of \( F \) and
\[
\Pr(L=0)=\Pr\{F(g-1,n-g)\le1\},
\]
a probability that tends to \( \tfrac12 \) as \( g \) and \( n-g \) grow.
:::

::: {.proof}
By @eq-mix-oneway-reml, \( -2\ell_R=(g-1)\log\lambda_1+(n-g)\log\lambda_2
+\text{SSB}/\lambda_1+\text{SSE}/\lambda_2 \) up to a constant, over the region
\( \lambda_1\ge\lambda_2>0 \). Under \( H_0 \) the region collapses to
\( \lambda_1=\lambda_2=\lambda \), where the objective is
\( (n-1)\log\lambda+(\text{SSB}+\text{SSE})/\lambda \), minimized at
\( \hat\lambda=(\text{SSB}+\text{SSE})/(n-1)=\text{MST} \), with minimum
\( (n-1)\log\text{MST}+(n-1) \).

If \( \text{MSB}\ge\text{MSE} \) the unrestricted maximizers of
@prp-mix-oneway-reml(a) are feasible, and the minimum of \( -2\ell_R \) is
\( (g-1)\log\text{MSB}+(n-g)\log\text{MSE}+(n-1) \). Subtracting,
\[
L=(n-1)\log\text{MST}-(g-1)\log\text{MSB}-(n-g)\log\text{MSE}.
\]
Dividing inside every logarithm by \( \text{MSE} \), which is legitimate because the
coefficients \( (n-1) \), \( -(g-1) \) and \( -(n-g) \) sum to zero, and using
\( \text{MST}/\text{MSE}=aF+b \), gives @eq-mix-ratio-exact. If
\( \text{MSB}<\text{MSE} \) then by @prp-mix-oneway-reml(c) the constrained maximizer
is on the boundary, which is the null model, so \( L=0 \).

Under \( H_0 \) the distribution of \( F \) is \( F(g-1,n-g) \) by
@prp-mix-oneway-reml and @prp-mix-oneway(d) with \( \gamma=0 \), and \( L \) is a
function of \( F \) alone. It is nondecreasing in \( F \) because, using
\( (n-1)a=g-1 \) and then \( 1-a=b \),
\[
\begin{aligned}
\frac{dL}{dF}&=\frac{(n-1)a}{aF+b}-\frac{g-1}{F}
=\frac{(g-1)\{F-(aF+b)\}}{F(aF+b)}\\
&=\frac{(g-1)b(F-1)}{F(aF+b)}\ \ge\ 0
\end{aligned}
\]
for \( F\ge1 \). Hence \( \{L=0\}=\{F\le1\} \). Finally, as \( r=g-1 \) and
\( s=n-g \) grow, \( \log F \) is asymptotically normal with variance
\( 2/r+2/s\to0 \) and mean of order \( 1/r+1/s \), which is small compared with the
standard deviation; so \( \Pr(\log F\le0)\to1/2 \).
:::

::: {.remark}
[Sketch of the chi-bar-squared limit]

The exact result explains the shape of the general answer. Put \( \delta=F-1 \) and
expand @eq-mix-ratio-exact for small \( \delta>0 \):
\( \log(1+a\delta)=a\delta-a^2\delta^2/2+O(\delta^3) \) and
\( \log(1+\delta)=\delta-\delta^2/2+O(\delta^3) \). Using \( (n-1)a=g-1 \) and
\( (n-1)a^2=(g-1)a \), the first-order terms cancel and
\[
L=\tfrac12(g-1)(1-a)\,\delta^2+O\bigl((n-1)\delta^3\bigr).
\]
Under \( H_0 \), \( \delta=F-1 \) has, to first order, the variance of
\( \log F \), namely \( \tau^2=2/(g-1)+2/(n-g)=2(n-1)/\{(g-1)(n-g)\} \), and
\( \delta/\tau \) is asymptotically standard normal. Since
\( (g-1)(1-a)=(g-1)(n-g)/(n-1) \), the product \( (g-1)(1-a)\tau^2 \) equals
\( 2 \) exactly, so \( L\approx Z^2 \) on \( \{\delta>0\} \) for a standard normal
\( Z \), while \( L=0 \) on \( \{\delta\le0\} \), an event of limiting probability
\( 1/2 \). So
\( L \) converges in distribution to
\( \tfrac12\chi^2(0)+\tfrac12\chi^2(1) \), the **chi-bar-squared** mixture: a point
mass at zero of weight one half, plus one half of a chi-square on one degree of
freedom. This is a sketch, not a proof; the expansions above are heuristic. The
general theorem, covering several variance components on the boundary simultaneously
with weights given by the geometry of the constraint cone, is due to Self and Liang
(1987), with the mixed-model case worked out by Stram and Lee (1994).
:::

::: {#exm-mix-boundary-size}
[What the calibration costs]

Take \( g=6 \) groups of \( m=4 \). Exactly
\( 0.554 \) of the null distribution of \( L \) sits at zero, and a
simulation of \( 200\,000 \) data sets puts
\( 0.554 \) there. Referring \( L \) to the chi-bar-squared
mixture, whose \( 5 \) per cent point is the \( 90 \) per cent point
\( 2.706 \) of \( \chi^2(1) \), gives an actual size of
\( 0.0411 \): slightly conservative. Referring it to
\( \chi^2(1) \), with critical value \( 3.841 \), gives
\( 0.0204 \), badly conservative — the test throws away about
half its power. With \( g=40 \) groups the mixture is nearly exact,
at \( 0.0462 \).
[Figure 32.5.3](05-inference.html#fig-mix-boundary) compares the three tail functions.
:::

::: {when-format="html"}
![**Figure 32.5.3.** Null tail probabilities of the restricted likelihood ratio
statistic for \( \sigma_a^2=0 \) in a balanced one-way layout, on a logarithmic scale:
exact (@eq-mix-ratio-exact), the chi-bar-squared mixture, and \( \chi^2(1) \). The
mixture is conservative for small \( g \) and accurate for large
\( g \).](boundary_lrt.svg){#fig-mix-boundary width=100%}
:::

::: {when-format="pdf"}
![Null tail probabilities of the restricted likelihood ratio statistic for
\( \sigma_a^2=0 \) in a balanced one-way layout, on a logarithmic scale:
exact (@eq-mix-ratio-exact), the chi-bar-squared mixture, and \( \chi^2(1) \). The mixture is
conservative for small \( g \) and accurate for large
\( g \).](boundary_lrt.pdf){width=100%}
:::

```{.python .run #cell-boundary-statistic}
import numpy as np
from scipy import stats

def lrt_from_f(f, g, m):
    """Restricted likelihood ratio statistic for sigma_a^2 = 0, from the F ratio."""
    n = g * m
    a, b = (g - 1) / (n - 1), (n - g) / (n - 1)
    return np.where(f <= 1, 0.0,
                    (n - 1) * np.log(a * np.maximum(f, 1) + b) - (g - 1) * np.log(np.maximum(f, 1)))

g, m = 6, 4
n = g * m
p_zero = stats.f.cdf(1.0, g - 1, n - g)             # P(L = 0) = P(F <= 1)
crit_mix = stats.chi2.ppf(0.90, 1)                  # 5% point of the 50:50 mixture
print(f"g = {g}, m = {m}:  P(L = 0) = {p_zero:.3f},  5% point of the mixture {crit_mix:.3f}")
```

Crainiceanu and Ruppert (2004) showed that for models with a single variance
component the exact null distribution can be computed in general, not only in the
balanced case, and that the chi-bar-squared approximation can be poor when the design
is unbalanced. Where their results do not apply, the parametric bootstrap is again the
practical answer.

## Intervals for variance components

Three routes, in decreasing order of reliability. In balanced designs,
@cor-mix-icc-interval gives exact intervals for ratios of variance components and for
the intraclass correlation. In general, the profile restricted likelihood of
@exm-mix-reml-lab gives an interval by cutting \( 1.92 \) units below the maximum; it
respects the boundary and is asymmetric, as the likelihood is. Last and least, the
Wald interval from the inverse information matrix is symmetric, frequently includes
negative values, and is badly calibrated whenever \( \hat\sigma_a^2 \) is near zero.
Burdick and Graybill (1992) is a book-length treatment of what can and cannot be done
exactly.

## Exercises

### A. Check your understanding

::: {#exr-mix-why-t}
[A1]

In @prp-mix-inference(b), the reported standard error is exactly right and the test is
still wrong. Explain in one or two sentences what "exactly right" means here and what
is wrong, and say why the problem gets worse as the intraclass correlation grows with
\( n \) fixed.
:::

::: {#exr-mix-reml-lrt-fixed}
[A2]

A colleague tests whether a regressor belongs in the fixed part of a mixed model by
comparing restricted log-likelihoods with and without it, referring twice the
difference to \( \chi^2(1) \). Give two separate reasons why this is wrong, and say
what to do instead.
:::

### B. Practice

::: {#exr-mix-satterthwaite-exact}
[B1]

Apply @eq-mix-satterthwaite to the estimated variance
\( 4\text{MS}_C/(gm) \) of @prp-mix-inference(b), a single multiple of one mean
square, and confirm that it returns \( \nu=g-2 \) exactly. Then apply it to
\( \hat\sigma_a^2=(\text{MSB}-\text{MSE})/m \) in the balanced one-way model and write
down the resulting degrees of freedom.
:::

::: {.solution}
With one term, \( U=c\,\text{MS}_1 \), @eq-mix-satterthwaite gives
\( \nu=c^2(\E\text{MS}_1)^2/\{c^2(\E\text{MS}_1)^2/\nu_1\}=\nu_1=g-2 \). For the
variance component, \( c_1=1/m \) with \( \nu_1=g-1 \) and \( c_2=-1/m \) with
\( \nu_2=n-g \), so
\[
\hat\nu=\frac{(\hat\sigma_a^2)^2}
{(\text{MSB}/m)^2/(g-1)+(\text{MSE}/m)^2/(n-g)} .
\]
Note that this can be far below \( g-1 \) when \( \text{MSB} \) and \( \text{MSE} \)
are close, which is the usual situation when \( \sigma_a^2 \) is small.
:::

::: {#exr-mix-boundary-monotone}
[B2]

Verify directly that the function \( F\mapsto(n-1)\log(aF+b)-(g-1)\log F \) in
@eq-mix-ratio-exact vanishes at \( F=1 \), is increasing for \( F>1 \) and tends to
infinity as \( F\to\infty \). What does the last fact say about the power of the test
when \( \sigma_a^2 \) is large?
:::

::: {#exr-mix-wald-variance-component}
[B3]

In the balanced one-way model, compute the Wald interval
\( \hat\sigma_a^2\pm1.96\,\text{se} \) for the proficiency study using
@exr-mix-variance-of-estimate for the standard error with the estimates plugged in,
and compare it with the exact interval implied by @cor-mix-icc-interval and with the
profile interval \( (0.0086,0.1275) \) of @exm-mix-reml-lab. Comment.
:::

### C. Going deeper

::: {#exr-mix-boundary-two}
[C1]

Consider testing \( \sigma_a^2=0 \) when a second variance component
\( \sigma_b^2>0 \) is present and is *not* on the boundary. Argue informally that the
limiting null distribution is still
\( \tfrac12\chi^2(0)+\tfrac12\chi^2(1) \), and explain why testing *both* components
simultaneously gives a mixture over \( \chi^2(0),\chi^2(1),\chi^2(2) \) with weights
determined by the probability that the unconstrained maximizer lands in each orthant.
Which weights arise if the two estimates are asymptotically independent?
:::

::: {.solution}
Locally the log-likelihood is approximately quadratic in the parameter, so the
likelihood ratio behaves like the squared distance from the unconstrained maximizer to
the constrained set, measured in the information metric. With one parameter on the
boundary, the constraint acts only in that coordinate: half the time the unconstrained
maximizer already satisfies it and the ratio is zero, half the time the constrained
maximum is on the face and the ratio is a squared standard normal. With two
parameters constrained, the ratio is zero when the maximizer lies in the nonnegative
orthant, a \( \chi^2(1) \) when it lies in one of the two adjacent orthants, and a
\( \chi^2(2) \) when it lies in the fully negative orthant. If the two coordinates are
asymptotically independent, those four orthants have probabilities
\( 1/4,1/4,1/4,1/4 \), giving weights \( 1/4 \) for \( \chi^2(0) \), \( 1/2 \) for
\( \chi^2(1) \) and \( 1/4 \) for \( \chi^2(2) \).
:::

::: {#exr-mix-kackar-counterexample}
[C2]

@prp-mix-inference(a) needs the error distribution to be symmetric about zero. Give an
example of a translation-invariant, even \( \hat{\boldsymbol{\uptheta}} \) and a
skewed error distribution for which \( \hbeta(\hat{\boldsymbol{\uptheta}}) \) is
biased, and explain why symmetry cannot simply be dropped. (Compare the discussion
after @thm-het-wls.)
:::
