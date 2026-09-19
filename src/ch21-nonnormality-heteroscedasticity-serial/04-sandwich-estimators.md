# Heteroscedasticity-consistent covariance estimators

Weighting needs a model for the variances. Often there is none worth trusting, and the analyst would rather
keep ordinary least squares, which is unbiased whatever the variances are, and fix its standard errors. The
covariance of \( \hbeta \) under \( \Cov(\be)=\bSigma=\diag(\sigma_1^2,\dots,\sigma_n^2) \) is
\[
\Cov(\hbeta)=(\X\T\X)^{-1}\Bigl(\sum_{i=1}^n\sigma_i^2\x_{(i)}\x_{(i)}\T\Bigr)(\X\T\X)^{-1}
\]{#eq-het-sandwich-true}

This is the covariance given by @thm-dep-covariance. The \( n \) variances cannot be estimated, one observation each. The idea of Eicker (1963),
Huber (1967) and White (1980) is that they do not need to be. Only the \( p\times p \) matrix in the middle is
needed, and a sum of \( n \) terms can be estimated consistently even when its terms cannot. The matrix is
sandwiched between two copies of \( (\X\T\X)^{-1} \), which gives the estimators their name.

## The estimators

::: {#def-het-hc}
[Heteroscedasticity-consistent covariance estimators]

Let \( \hat{\varepsilon}_i \) and \( h_{ii} \) be the least squares residuals and leverages, with \( h_{ii}<1 \). For
weights \( \omega_i \), put
\[
\hat{\V}_\omega=(\X\T\X)^{-1}\Bigl(\sum_{i=1}^n\omega_i\x_{(i)}\x_{(i)}\T\Bigr)(\X\T\X)^{-1}.
\]
The estimators **HC0** to **HC3** take
\[
\omega_i=\hat{\varepsilon}_i^2,\qquad \frac{n}{n-p}\,\hat{\varepsilon}_i^2,\qquad \frac{\hat{\varepsilon}_i^2}{1-h_{ii}},\qquad
\frac{\hat{\varepsilon}_i^2}{(1-h_{ii})^2}.
\]
:::

HC0 is White's estimator. HC1 applies the degrees-of-freedom correction of \( s^2 \). HC2 and HC3, proposed by
MacKinnon and White (1985), correct each squared residual for its own leverage. Their finite-sample
properties follow from @eq-het-residual-second-moment and from the deletion formulas of
[Chapter 10](../ch10-computation/index.html).

::: {#prp-het-hc-facts}
[Finite-sample properties]

Let \( \blambda\in\Real^p \) and \( \mathbf{c}=\X(\X\T\X)^{-1}\blambda \), so that \( \blambda\T\hbeta-\blambda\T\bbeta=\mathbf{c}\T\be \) and
\( \blambda\T\hat{\V}_\omega\blambda=\sum_ic_i^2\omega_i \).

::: {.enumerate options="label=(\alph*)"}
1. Under constant variance \( \sigma^2 \), \( \E(\blambda\T\hat{\V}_{\text{HC2}}\blambda)=\Var(\blambda\T\hbeta) \), while
   \( \E(\blambda\T\hat{\V}_{\text{HC0}}\blambda)=\Var(\blambda\T\hbeta)-\sigma^2\sum_ic_i^2h_{ii}\le\Var(\blambda\T\hbeta) \) and
   \( \E(\blambda\T\hat{\V}_{\text{HC3}}\blambda)\ge\Var(\blambda\T\hbeta) \).

2. In general, \( \E(\blambda\T\hat{\V}_{\text{HC0}}\blambda)=\sum_ic_i^2\sum_j(\I-\M)_{ij}^2\sigma_j^2 \).

3. \( \hat{\V}_{\text{HC0}}\preceq\hat{\V}_{\text{HC2}}\preceq\hat{\V}_{\text{HC3}} \) in the order of nonnegative definiteness.

4. \( \hat{\V}_{\text{HC3}}=\sum_{i=1}^n(\hbeta_{(i)}-\hbeta)(\hbeta_{(i)}-\hbeta)\T \), where \( \hbeta_{(i)} \) is the estimate with
   observation \( i \) deleted.
:::

:::

::: {.proof}
(a) and (b). \( \E\hat{\varepsilon}_i^2=\sum_j(\I-\M)_{ij}^2\sigma_j^2 \) by @eq-het-residual-second-moment, which proves (b).
Under constant variance it is \( \sigma^2(1-h_{ii}) \), so the HC2 weights have expectation \( \sigma^2 \) and
\( \E\sum_ic_i^2\omega_i=\sigma^2\norm{\mathbf{c}}^2=\Var(\mathbf{c}\T\be) \). The HC0 weights have expectation
\( \sigma^2(1-h_{ii}) \) and the HC3 weights \( \sigma^2/(1-h_{ii})\ge\sigma^2 \).

(c) The HC weights satisfy \( \hat{\varepsilon}_i^2\le\hat{\varepsilon}_i^2/(1-h_{ii})\le\hat{\varepsilon}_i^2/(1-h_{ii})^2 \) term by
term, and each difference \( \hat{\V}_{\omega'}-\hat{\V}_\omega \) with \( \omega'_i\ge\omega_i \) is
\( (\X\T\X)^{-1}\{\sum_i(\omega'_i-\omega_i)\x_{(i)}\x_{(i)}\T\}(\X\T\X)^{-1} \), a nonnegative definite matrix.

(d) By @prp-cmp-loo(a), \( \hbeta_{(i)}-\hbeta=-(\X\T\X)^{-1}\x_{(i)}\hat{\varepsilon}_i/(1-h_{ii}) \). The sum of the outer
products of these vectors is \( \hat{\V}_{\text{HC3}} \).
:::

Part (d) says that HC3 measures how much the estimate moves when single observations are left out. This is
the jackknife idea: the variability of an estimate is judged from its sensitivity to the data, not from a
model for the errors. Part (a) makes HC0 biased downward when the model is right, and the bias is
concentrated where \( c_i^2h_{ii} \) is large, at high-leverage points that dominate the coefficient. Part (b)
shows the same effect under heteroscedasticity: a high-leverage observation's squared residual understates its
own variance, precisely when that variance matters most.

```{.python .run #cell-sandwich-engel}
import numpy as np
import statsmodels.api as sm
from scipy import stats

df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])
Q, R = np.linalg.qr(X)
beta = np.linalg.solve(R, Q.T @ food)
e = food - X @ beta
h = np.sum(Q**2, axis=1)                              # leverages
Rinv = np.linalg.inv(R)                               # (X'X)^{-1} = Rinv Rinv'


def sandwich(omega):
    """(X'X)^{-1} X' diag(omega) X (X'X)^{-1}, computed through Q."""
    meat = (Q * omega[:, None]).T @ Q
    return Rinv @ meat @ Rinv.T


p = X.shape[1]
covs = {
    "classical": (e @ e / (n - p)) * (Rinv @ Rinv.T),
    "HC0": sandwich(e**2),
    "HC1": sandwich(e**2) * n / (n - p),
    "HC2": sandwich(e**2 / (1 - h)),
    "HC3": sandwich(e**2 / (1 - h) ** 2),
}
for name, V in covs.items():
    print(f"{name:9s} slope se {np.sqrt(V[1, 1]):.4f}")
```

```{.python .run #cell-sandwich-jackknife}
loo = np.array([np.linalg.lstsq(np.delete(X, i, 0), np.delete(food, i),
                                rcond=None)[0]
                for i in range(n)])
# sum over i of (b_(i) - b)(b_(i) - b)'
jack = (loo - beta).T @ (loo - beta)
print("HC3 equals the sum of squared deletion changes:",
      np.allclose(jack, covs["HC3"]))
```

::: {#exm-het-engel-hc}
[Standard errors for Engel's slope]

For the least squares slope of food expenditure on income, the classical standard error
\( s\{(\X\T\X)^{-1}\}_{22}^{1/2} \) is \( 0.0144 \). The heteroscedasticity-consistent standard errors are \( 0.0518 \)
(HC0), \( 0.0520 \) (HC1), \( 0.0585 \) (HC2) and \( 0.0664 \) (HC3). The classical value is too small by a factor of
\( 4.6 \) relative to HC3. Evidently the high-income households, which pin down the slope, are the noisy ones, so
the classical formula, which averages the noise over all households, understates the uncertainty. The second
listing above checks the jackknife identity of @prp-het-hc-facts(d).
:::

## Consistency

The estimators are justified by a large-sample argument. We give it under simple conditions that make every
step elementary: bounded regressors, errors with bounded fourth moments, and variances bounded away from
zero. The theorem covers every \( \hat{\V}_\omega \) of @def-het-hc.

::: {#thm-het-sandwich}
[Consistency of the sandwich estimators]

Let \( \Y=\X\bbeta+\be \) with \( p \) fixed. Assume that

- (H1) the errors are independent with mean zero and variances \( \sigma_i^2 \), with \( \sigma_i^2\ge\sigma_{\min}^2>0 \) and
  \( \E\varepsilon_i^4\le K \) for all \( i \);
- (H2) the rows satisfy \( \norm{\x_{(i)}}\le C \) for all \( i \);
- (H3) \( n^{-1}\X\T\X\to\A \), a positive definite matrix.

Let \( \hat{\V} \) be any of HC0–HC3, and \( \V=\Cov(\hbeta) \) as in @eq-het-sandwich-true. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( n(\hat{\V}-\V)\to\mathbf{0} \) in probability;

2. for every fixed \( \blambda\ne\bzero \),
   \( (\blambda\T\hbeta-\blambda\T\bbeta)/(\blambda\T\hat{\V}\blambda)^{1/2}\to\Normal(0,1) \) in distribution;

3. for every fixed \( p\times q \) matrix \( \bLambda \) of rank \( q \),
   \[
W=(\bLambda\T\hbeta-\bLambda\T\bbeta)\T\bigl(\bLambda\T\hat{\V}\bLambda\bigr)^{-1}(\bLambda\T\hbeta-\bLambda\T\bbeta)
\]
   tends to \( \chi^2(q) \) in distribution.
:::

:::

::: {.proof}
Write \( \mathbf{G}_n=n^{-1}\X\T\X \), and let \( \mathbf{B}_n=n^{-1}\sum_i\sigma_i^2\x_{(i)}\x_{(i)}\T \) and
\( \hat{\mathbf{B}}_n=n^{-1}\sum_i\omega_i\x_{(i)}\x_{(i)}\T \) be the true and estimated middle matrices, so that
\( n\V=\mathbf{G}_n^{-1}\mathbf{B}_n\mathbf{G}_n^{-1} \) and \( n\hat{\V}=\mathbf{G}_n^{-1}\hat{\mathbf{B}}_n\mathbf{G}_n^{-1} \). Note
\( \sigma_i^2\le(\E\varepsilon_i^4)^{1/2}\le K^{1/2} \). By (H3), \( \mathbf{G}_n^{-1}\to\A^{-1} \), and
\( \lambda_{\min}(\X\T\X)\ge cn \) for some \( c>0 \) and all large \( n \).

*Step 1: \( \hbeta \) is consistent.* \( \E\norm{\hbeta-\bbeta}^2=\tr\V\le K^{1/2}\tr(\X\T\X)^{-1}\le K^{1/2}p/(cn)\to0 \).

*Step 2: the errors estimate the middle.* Each entry of \( n^{-1}\sum_i(\varepsilon_i^2-\sigma_i^2)\x_{(i)}\x_{(i)}\T \) is a
sum of independent mean-zero terms with variance at most \( n^{-2}\sum_iKC^4=KC^4/n \), so it tends to zero in
probability by Chebyshev's inequality.

*Step 3: the residuals estimate the middle.* Let \( d_i=\x_{(i)}\T(\hbeta-\bbeta) \), so
\( \hat{\varepsilon}_i=\varepsilon_i-d_i \) and \( \lvert d_i\rvert\le C\norm{\hbeta-\bbeta} \). Then, in the spectral norm,
\[
\begin{aligned}
\Bigl\lVert n^{-1}\sum_i(\hat{\varepsilon}_i^2-\varepsilon_i^2)\x_{(i)}\x_{(i)}\T\Bigr\rVert
&\le\frac{C^2}{n}\sum_i\bigl(d_i^2+2\lvert d_i\rvert\lvert\varepsilon_i\rvert\bigr)\\
&\le C^4\norm{\hbeta-\bbeta}^2+2C^3\norm{\hbeta-\bbeta}\,\frac1n\sum_i\lvert\varepsilon_i\rvert ,
\end{aligned}
\]
which tends to zero because \( \hbeta-\bbeta=o_p(1) \) and \( \E n^{-1}\sum_i\lvert\varepsilon_i\rvert\le K^{1/4} \). With
Step 2, \( \hat{\mathbf{B}}_n-\mathbf{B}_n\to\mathbf{0} \) for HC0.

*Step 4: the leverage corrections vanish.* \( h_{ii}\le C^2/\lambda_{\min}(\X\T\X)\le C^2/(cn) \), so
\( \delta_n=(1-\max_ih_{ii})^{-2}-1\to0 \). For HC1–HC3, \( \hat{\varepsilon}_i^2\le\omega_i\le(1+\delta_n)\hat{\varepsilon}_i^2 \)
for large \( n \) (with \( \delta_n \) replaced by \( p/(n-p) \) for HC1), so
\( \mathbf{0}\preceq\hat{\mathbf{B}}_n-\hat{\mathbf{B}}_n^{\text{HC0}}\preceq\delta_n\hat{\mathbf{B}}_n^{\text{HC0}} \), which tends to
zero because \( \delta_n\to0 \) and \( \hat{\mathbf{B}}_n^{\text{HC0}}=\mathbf{B}_n+o_p(1) \) is bounded in probability,
since \( \mathbf{0}\preceq\mathbf{B}_n\preceq K^{1/2}\mathbf{G}_n \) and \( \mathbf{G}_n \) converges. This proves (a), since
\( n(\hat{\V}-\V)=\mathbf{G}_n^{-1}(\hat{\mathbf{B}}_n-\mathbf{B}_n)\mathbf{G}_n^{-1} \).

(c) Let \( \V_\Lambda=\bLambda\T\V\bLambda \) and \( \mathbf{T}_n=\V_\Lambda^{-1/2}\bLambda\T(\hbeta-\bbeta) \), so that
\( \Cov(\mathbf{T}_n)=\I_q \). For fixed \( \mathbf{a}\in\Real^q \), \( \mathbf{a}\T\mathbf{T}_n=\mathbf{c}\T\be \) with
\( \mathbf{c}=\X(\X\T\X)^{-1}\blambda \) and \( \blambda=\bLambda\V_\Lambda^{-1/2}\mathbf{a} \). By the Cauchy–Schwarz inequality
in the metric \( (\X\T\X)^{-1} \), \( c_i^2\le h_{ii}\norm{\mathbf{c}}^2 \), so
\( \max_ic_i^2/\norm{\mathbf{c}}^2\le\max_ih_{ii}\to0 \). @lem-het-clt(b), with \( \delta=2 \), shows that
\( \mathbf{a}\T\mathbf{T}_n/\norm{\mathbf{a}}\to\Normal(0,1) \) for \( \mathbf{a}\ne\bzero \). By the Cramér–Wold device,
\( \mathbf{T}_n\to\Normal_q(\bzero,\I_q) \). Next, \( n\V_\Lambda=\bLambda\T\mathbf{G}_n^{-1}\mathbf{B}_n\mathbf{G}_n^{-1}\bLambda \succeq
\sigma_{\min}^2\bLambda\T\mathbf{G}_n^{-1}\bLambda \), whose limit \( \sigma_{\min}^2\bLambda\T\A^{-1}\bLambda \) is positive definite; and
\( n\V_\Lambda \) is bounded because \( \mathbf{B}_n\preceq K^{1/2}\mathbf{G}_n \). So \( (n\V_\Lambda)^{-1/2} \) is bounded, and by (a)
\[
\begin{aligned}
\V_\Lambda^{-1/2}\bigl(\bLambda\T\hat{\V}\bLambda\bigr)\V_\Lambda^{-1/2}
&=\I_q+(n\V_\Lambda)^{-1/2}\,\bLambda\T\\
&\quad\times n(\hat{\V}-\V)\bLambda\,(n\V_\Lambda)^{-1/2}\to\I_q .
\end{aligned}
\]
Since \( W=\mathbf{T}_n\T\bigl\{\V_\Lambda^{-1/2}(\bLambda\T\hat{\V}\bLambda)\V_\Lambda^{-1/2}\bigr\}^{-1}\mathbf{T}_n \), continuous mapping
and Slutsky's lemma give \( W\to\mathbf{T}\T\mathbf{T} \) with \( \mathbf{T}\sim\Normal_q(\bzero,\I_q) \), which is \( \chi^2(q) \). Part (b)
is the case \( q=1 \), with the sign kept.
:::

The theorem does not assume that the errors are normal, that their variances follow any model, or that they
are identically distributed. It does assume independence, which the next two sections relax. The bounded
regressors of (H2) are stronger than needed. Eicker (1967) and White (1980) replace them by conditions that
amount to \( \max_ih_{ii}\to0 \) together with moment bounds, and White's version allows random regressors
drawn independently from a common distribution, which covers most observational data. We do not prove these
extensions.

In practice the statistic in (b) is compared with the \( t(n-p) \) distribution rather than the normal. That
choice has no exact justification, but it costs nothing asymptotically and helps in small samples. A
heteroscedasticity-robust test of \( \bLambda\T\bbeta=\mathbf{d} \) compares \( W/q \) with \( F(q,n-p) \).

## Finite-sample behaviour

Consistency says nothing about how large \( n \) must be. The answer depends on the leverages, as
@prp-het-hc-facts suggests, and it can be sobering.

::: {#exm-het-hc-coverage}
[Coverage in a skewed design]

Let the regressor values be the \( n \) quantiles \( x_i=\exp\{\Phi^{-1}((i-1/2)/n)\} \) of a standard lognormal
distribution, a design with a long right tail like incomes. The largest leverage is \( 0.58 \) for \( n=25 \) and
still \( 0.15 \) for \( n=1000 \). [Figure 21.4.1](#fig-het-hc-coverage) shows the coverage of nominal 95% intervals
for the slope, with normal errors, over \( 20{,}000 \) simulated data sets for each \( n \).

With constant variance (panel a) the classical interval is exact, and the robust intervals pay a price for
their robustness: at \( n=25 \) they cover \( 0.852 \) (HC0), \( 0.866 \) (HC1), \( 0.895 \) (HC2) and \( 0.932 \) (HC3).
With standard deviation proportional to \( x \) (panel b) the classical interval fails badly. Its coverage is
\( 0.456 \) at \( n=25 \) and falls to \( 0.287 \) at \( n=1000 \), because the growing right tail of the design gives ever
more influence to the noisiest observations. The robust intervals improve with \( n \), but slowly. At \( n=25 \)
they cover \( 0.631 \) (HC0), \( 0.650 \) (HC1), \( 0.760 \) (HC2) and \( 0.853 \) (HC3); at \( n=1000 \), \( 0.912 \),
\( 0.912 \), \( 0.924 \) and \( 0.935 \). HC3 is the best of the four throughout, as Long and Ervin (2000) found in
a wide range of designs.
:::

::: {when-format="html"}
![**Figure 21.4.1.** Coverage of nominal 95% intervals for the slope in a design with lognormal quantiles as
regressor values. (a) Constant error variance. (b) Error standard deviation proportional to \( x \). HC1 is
omitted; it lies just above HC0.](sandwich_coverage.svg){#fig-het-hc-coverage width=100%}
:::

::: {when-format="pdf"}
![Coverage of nominal 95% intervals for the slope in a design with lognormal quantiles as
regressor values. (a) Constant error variance. (b) Error standard deviation proportional to \( x \). HC1 is
omitted; it lies just above HC0.](sandwich_coverage.pdf){width=100%}
:::

```{.python .run #cell-sandwich-coverage}
rng = np.random.default_rng(2107)


def coverage(n, hetero, reps):
    """Coverage of nominal 95% slope intervals, classical and HC0-HC3."""
    # lognormal quantiles
    x = np.exp(stats.norm.ppf((np.arange(1, n + 1) - 0.5) / n))
    Xs = np.column_stack([np.ones(n), x])
    Qs, Rs = np.linalg.qr(Xs)
    hs = np.sum(Qs**2, axis=1)
    # slope estimate - slope = c @ errors
    c = Qs @ np.linalg.inv(Rs).T[:, 1]
    E = rng.normal(size=(reps, n)) * (x if hetero else 1.0)
    err = E @ c
    res = E - (E @ Qs) @ Qs.T
    var = {"classical": (res**2).sum(axis=1) / (n - 2) * (c @ c),
           "HC0": (res**2) @ c**2,
           "HC1": (res**2) @ c**2 * n / (n - 2),
           "HC2": (res**2 / (1 - hs)) @ c**2,
           "HC3": (res**2 / (1 - hs) ** 2) @ c**2}
    q = stats.t.ppf(0.975, n - 2)
    cover = {k: np.mean(np.abs(err) <= q * np.sqrt(v)) for k, v in var.items()}
    return cover, hs.max()


for n_sim in (25, 100, 400):
    cov, hmax = coverage(n_sim, hetero=True, reps=2000)
    print(f"n = {n_sim:3d} (max leverage {hmax:.2f}):",
          {k: round(v, 3) for k, v in cov.items()})
```

## Leverage and the sandwich

The poor small-sample behaviour has a simple source. The variance estimate \( \sum_ic_i^2\omega_i \) is a weighted
sum of squared residuals, and when a few weights \( c_i^2 \) dominate, it behaves like a variance estimate with
very few degrees of freedom, however large \( n \) is. A \( t(n-p) \) reference distribution then has tails that
are far too light.

::: {#exm-het-engel-leverage}
[One household carries the standard error]

In Engel's data the household with the largest income has leverage \( 0.255 \). Its term
\( c_i^2\hat{\varepsilon}_i^2/(1-h_{ii})^2 \) makes up \( 0.85 \) of the HC3 variance of the slope. The robust standard
error of @exm-het-engel-hc is, in effect, an estimate of one household's error variance. That is not a
reason to go back to the classical standard error, which is far worse. It is a reason to look at that
household (the influence measures of [Chapter 20](../ch20-residuals-leverage-influence/index.html),
such as @def-res-cooks, are the tools), to report how the conclusions change without it, and to prefer
weighting, which here gives that household much less influence (@exm-het-engel-fwls).
:::

Several refinements address the degrees-of-freedom problem. Bell and McCaffrey's (2002) adjustment, as recommended
by Imbens and Kolesár (2016), combines HC2 with a Satterthwaite approximation to the degrees of freedom of
\( \sum_ic_i^2\omega_i \). The wild bootstrap of [Chapter 23](../ch23-resampling-inference/index.html) (@def-bs-wild)
resamples the residuals with random signs, keeping each one attached to its own row, and often gives better
calibrated intervals than any of the HC formulas.

## Exercises

### A. Check your understanding

::: {#exr-het-welch}
[A1]

Take the two-group model \( y_{gj}=\mu_g+\varepsilon_{gj} \), \( g=1,2 \), \( j=1,\dots,n_g \). Show that HC2 for
\( \hat{\mu}_1-\hat{\mu}_2 \) equals \( s_1^2/n_1+s_2^2/n_2 \), the variance estimate of Welch's test, where \( s_g^2 \) is the
sample variance of group \( g \). What does HC0 give?
:::

::: {.solution}
Here \( h_{ii}=1/n_g \), \( \hat{\varepsilon}_{gj}=y_{gj}-\bar y_g \), and \( \mathbf{c} \) has entries \( 1/n_1 \) in group 1 and
\( -1/n_2 \) in group 2. HC2 is
\( \sum_g\sum_j\hat{\varepsilon}_{gj}^2/\{(1-1/n_g)n_g^2\}=\sum_g\sum_j(y_{gj}-\bar y_g)^2/\{(n_g-1)n_g\}=\sum_gs_g^2/n_g \). HC0 gives
\( \sum_g(n_g-1)s_g^2/n_g^2 \), which is smaller.
:::

::: {#exr-het-leverage-one}
[A2]

What happens to HC0, HC2 and HC3 if some observation has \( h_{ii}=1 \)? Explain why no estimator built from
residuals can estimate that observation's contribution to \( \Var(\blambda\T\hbeta) \).
:::

### B. Practice

::: {#exr-het-hc-homoscedastic-bias}
[B1]

In simple regression with an intercept and constant variance, show that the relative bias of HC0 for the slope
variance is \( -\sum_i(x_i-\bar x)^2h_{ii}/S_{xx} \). Evaluate it for \( x=(1,2,3,4,10) \) of @exr-het-residual-variance.
:::

::: {.solution}
For the slope, \( c_i=(x_i-\bar x)/S_{xx} \) and \( \norm{\mathbf{c}}^2=1/S_{xx} \), so @prp-het-hc-facts(a) gives relative bias
\( -\sum_ic_i^2h_{ii}/\norm{\mathbf{c}}^2=-\sum_i(x_i-\bar x)^2h_{ii}/S_{xx} \). With \( S_{xx}=50 \) and the leverages
\( (0.38,0.28,0.22,0.20,0.92) \), the sum is \( 9(0.38)+4(0.28)+1(0.22)+0+36(0.92)=37.88 \), and the relative bias is
\( -37.88/50=-0.76 \). HC0 estimates the slope variance at a quarter of its true value, on average, because the
point \( x=10 \) dominates the slope and has a tiny residual.
:::

::: {#exr-het-hc-wls}
[B2]

Using @thm-het-wls(b), write down a sandwich estimator of the covariance of a weighted least squares estimate
with fixed weights \( \mathbf{W} \), and show that it reduces to HC0 applied to the rescaled data
\( (\mathbf{W}^{1/2}\X,\mathbf{W}^{1/2}\y) \).
:::

::: {#exr-het-white-test}
[B3]

White's test compares \( s^2\X\T\X \) with \( \sum_i\hat{\varepsilon}_i^2\x_{(i)}\x_{(i)}\T \). Show that the difference of these
matrices, divided by \( n \), has entries \( n^{-1}\sum_i(\hat{\varepsilon}_i^2-s^2)x_{ij}x_{ik} \), and explain why
Koenker's statistic with \( \bz_i \) the distinct products \( x_{ij}x_{ik} \) tests whether these are jointly zero.
:::

### C. Going deeper

::: {#exr-het-hc-random}
[C1]

Suppose the rows \( (\x_{(i)},\varepsilon_i) \) are independent copies of a random pair with \( \E(\varepsilon\mid\x)=0 \),
\( \E\norm{\x}^4<\infty \), \( \E\varepsilon^4<\infty \) and \( \E\x\x\T \) positive definite. Adapt the proof of
@thm-het-sandwich to show that HC0 is consistent in the sense of part (a), replacing Chebyshev's inequality by the
weak law of large numbers.
:::

::: {#exr-het-hc-jackknife-centred}
[C2]

The usual jackknife variance is \( \frac{n-1}{n}\sum_i(\hbeta_{(i)}-\mathbf{m})(\hbeta_{(i)}-\mathbf{m})\T \),
where \( \mathbf{m} \) is the mean of the \( \hbeta_{(i)} \). Express it in terms of \( \hat{\V}_{\text{HC3}} \) and
the vector \( \sum_i\x_{(i)}\hat{\varepsilon}_i/(1-h_{ii}) \), and show that the two agree to first order when
\( \max_ih_{ii}\to0 \).
:::

