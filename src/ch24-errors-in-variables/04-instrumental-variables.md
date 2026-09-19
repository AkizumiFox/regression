# Instrumental variables

The corrections so far used a *number* from outside the regression. Instrumental variables use a *variable*,
related to the mismeasured regressor and unrelated to everything that makes least squares fail. The idea goes
back to Philip Wright (1928), who used it to estimate supply and demand curves.

## Endogeneity

By @eq-eiv-composite, \( Y=\alpha+\beta W+e \) with \( e=\varepsilon-\beta U \) correlated with \( W \). An omitted common cause (@thm-dep-omitted)
or a regressor partly determined by the response has the same effect. Such regressors are
called **endogenous**. In every case the equation
\[
y_i=\x_{(i)}\T\bbeta+e_i,\qquad \E(\x_{(i)}e_i)\ne\bzero ,
\]{#eq-eiv-structural-equation}

describes a relation of scientific interest, not the best linear predictor that least squares estimates. Here
\( \x_{(i)} \) holds all the regressors, including the constant. In this section \( \x \) collects all regressors and
\( \bz \) all instruments; the exactly measured covariates of earlier sections appear in both.

::: {#def-eiv-instrument}
[Instruments]

A random \( q \)-vector \( \bz \) is a vector of **instruments** for @eq-eiv-structural-equation, with random
regressor \( p \)-vector \( \x \), if

::: {.enumerate options="label=(\roman*)"}
1. *exogeneity:* \( \E(\bz e)=\bzero \);

2. *relevance:* \( \bSigma_{zx}=\E(\bz\x\T) \) has rank \( p \), which requires \( q\ge p \);
:::

and \( \bSigma_{zz}=\E(\bz\bz\T) \) is positive definite.
:::

Regressors uncorrelated with \( e \), such as the constant, are their own instruments; the rank condition asks the
*excluded* instruments to carry extra information about the endogenous regressors. With \( q=p \) the equation is **just identified**, and
with \( q>p \) **overidentified**.

In the blood-pressure study the second reading instruments the first: \( \x=(1,W_1,\text{age}) \),
\( \bz=(1,W_2,\text{age}) \) and \( e=\varepsilon-\beta U_1 \), which is uncorrelated with \( W_2=X+U_2 \). Nothing about \( \sigma_u^2 \)
needs to be known, but the two errors must be *independent*, which rules out two readings that share a
transient cause (@exr-eiv-correlated-replicates).

## The estimators

With the \( n\times p \) matrix \( \X \) of regressors and the \( n\times q \) matrix \( \Z \) of instruments, the sample
version of \( \E\{\bz(y-\x\T\bbeta)\}=\bzero \) in the just-identified case is \( \Z\T(\y-\X\bbeta)=\bzero \), with solution
\[
\hbeta_{\text{IV}}=(\Z\T\X)^{-1}\Z\T\y .
\]
When \( q>p \) these are more equations than unknowns. **Two-stage least squares** (2SLS) combines them through the
projection \( \M_{\Z}=\Z(\Z\T\Z)^{-1}\Z\T \) onto \( \C(\Z) \):
\[
\hbeta_{\text{2SLS}}=(\X\T\M_{\Z}\X)^{-1}\X\T\M_{\Z}\y .
\]{#eq-eiv-2sls}

The name comes from a two-step computation.

::: {#lem-eiv-2sls}
[Two-stage least squares as projection]

Let \( \Z \) have full column rank and \( \X\T\M_{\Z}\X \) be nonsingular, and let \( \hat{\X}=\M_{\Z}\X \) be the matrix of
first-stage fitted values, obtained by regressing each column of \( \X \) on \( \Z \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta_{\text{2SLS}} \) is the least squares coefficient of \( \y \) on \( \hat{\X} \), and also the IV estimator with
   instruments \( \hat{\X} \): \( \hbeta_{\text{2SLS}}=(\hat{\X}\T\hat{\X})^{-1}\hat{\X}\T\y=(\hat{\X}\T\X)^{-1}\hat{\X}\T\y \).

2. Columns of \( \X \) that are also columns of \( \Z \) are unchanged by the first stage.

3. If \( q=p \) and \( \Z\T\X \) is nonsingular, \( \hbeta_{\text{2SLS}}=\hbeta_{\text{IV}} \).
:::

:::

::: {.proof}
*(a)* \( \M_{\Z} \) is the orthogonal projection onto \( \C(\Z) \) (@thm-proj-M-formula), so it is symmetric and idempotent (@thm-proj-sym-idem).
Hence \( \hat{\X}\T\hat{\X}=\X\T\M_{\Z}\M_{\Z}\X=\X\T\M_{\Z}\X \), and
\( \hat{\X}\T\X=\X\T\M_{\Z}\X \) and \( \hat{\X}\T\y=\X\T\M_{\Z}\y \) as well. Substitute these in @eq-eiv-2sls.
*(b)* A column in \( \C(\Z) \) is fixed by \( \M_{\Z} \).
*(c)* With \( \Z\T\X \) square and invertible,
\( (\X\T\Z(\Z\T\Z)^{-1}\Z\T\X)^{-1}=(\Z\T\X)^{-1}(\Z\T\Z)(\X\T\Z)^{-1} \), and multiplying by
\( \X\T\Z(\Z\T\Z)^{-1}\Z\T\y \) leaves \( (\Z\T\X)^{-1}\Z\T\y \).
:::

The first stage replaces each endogenous regressor by the part of it that the instruments explain, which is
uncorrelated with \( e \). As with the residualized regressions of [Section 6.6](../ch06-projections/06-fwl.html), the
second-stage *coefficients* are right but its *residuals* \( \y-\hat{\X}\hbeta \) are not; the equation's residuals are
\( \y-\X\hbeta \).

::: {#thm-eiv-iv}
[Consistency and asymptotic normality of IV]

Let \( (\x_{(i)},\bz_{(i)},e_i) \) be independent copies of \( (\x,\bz,e) \) with finite second moments, satisfying
@eq-eiv-structural-equation and @def-eiv-instrument. Let \( \mathbf{H}=\bSigma_{xz}\bSigma_{zz}^{-1}\bSigma_{zx} \), where
\( \bSigma_{xz}=\bSigma_{zx}\T \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta_{\text{2SLS}}\to\bbeta \) with probability one.

2. If \( \E(e^2\norm{\bz}^2)<\infty \), then \( \sqrt n(\hbeta_{\text{2SLS}}-\bbeta)\to\Normal_p(\bzero,\mathbf{V}) \) in distribution, where
   \[
   \mathbf{V}=\mathbf{H}^{-1}\bSigma_{xz}\bSigma_{zz}^{-1}\boldsymbol{\Omega}\bSigma_{zz}^{-1}\bSigma_{zx}\mathbf{H}^{-1},\qquad
   \boldsymbol{\Omega}=\E(e^2\bz\bz\T).
   \]
   If \( \E(e^2\mid\bz)=\sigma_e^2 \) is constant, then \( \mathbf{V}=\sigma_e^2\mathbf{H}^{-1} \).

3. \( \hat{\sigma}_e^2=n^{-1}\norm{\y-\X\hbeta_{\text{2SLS}}}^2\to\sigma_e^2=\E(e^2) \) with probability one.
:::

:::

::: {.proof}
*(a)* Write \( \hat{\mathbf{C}}=n^{-1}\Z\T\X \) and \( \hat{\bS}=n^{-1}\Z\T\Z \). Substituting \( \y=\X\bbeta+\be \) in
@eq-eiv-2sls,
\[
\hbeta_{\text{2SLS}}-\bbeta=\bigl(\hat{\mathbf{C}}\T\hat{\bS}^{-1}\hat{\mathbf{C}}\bigr)^{-1}\hat{\mathbf{C}}\T\hat{\bS}^{-1}\,\tfrac1n\Z\T\be .
\]
By the strong law, \( \hat{\mathbf{C}}\to\bSigma_{zx} \), \( \hat{\bS}\to\bSigma_{zz} \) and \( n^{-1}\Z\T\be\to\E(\bz e)=\bzero \). The
matrix \( \mathbf{H} \) is positive definite: for \( \mathbf{a}\ne\bzero \), \( \bSigma_{zx}\mathbf{a}\ne\bzero \) by the rank condition, so
\( \mathbf{a}\T\mathbf{H}\mathbf{a}=(\bSigma_{zx}\mathbf{a})\T\bSigma_{zz}^{-1}(\bSigma_{zx}\mathbf{a})>0 \). By continuity the right side tends to
\( \mathbf{H}^{-1}\bSigma_{xz}\bSigma_{zz}^{-1}\bzero=\bzero \).

*(b)* Multiply the display by \( \sqrt n \). The vectors \( \bz_{(i)}e_i \) are independent and identically distributed with mean
\( \bzero \) and covariance \( \boldsymbol{\Omega} \), finite by assumption, so \( n^{-1/2}\Z\T\be\to\Normal_q(\bzero,\boldsymbol{\Omega}) \) by the
central limit theorem. The matrix in front converges to \( \mathbf{L}=\mathbf{H}^{-1}\bSigma_{xz}\bSigma_{zz}^{-1} \), and the
splitting used for @thm-eiv-correction(b) gives the limit \( \mathbf{L}\,\Normal_q(\bzero,\boldsymbol{\Omega})=\Normal_p(\bzero,\mathbf{V}) \). If \( \E(e^2\mid\bz)=\sigma_e^2 \), then \( \boldsymbol{\Omega}=\E\{\E(e^2\mid\bz)\bz\bz\T\}=\sigma_e^2\bSigma_{zz} \), and
\( \mathbf{V} \) collapses to \( \sigma_e^2\mathbf{H}^{-1} \).

*(c)* With \( \boldsymbol{\updelta}=\hbeta_{\text{2SLS}}-\bbeta \), \( \y-\X\hbeta_{\text{2SLS}}=\be-\X\boldsymbol{\updelta} \), so
\[
\hat{\sigma}_e^2=\tfrac1n\be\T\be-2\boldsymbol{\updelta}\T\tfrac1n\X\T\be+\boldsymbol{\updelta}\T\bigl(\tfrac1n\X\T\X\bigr)\boldsymbol{\updelta} .
\]
The averages converge to \( \E(e^2) \), \( \E(\x e) \) and \( \E(\x\x\T) \), all finite, and \( \boldsymbol{\updelta}\to\bzero \) by (a).
:::

In practice \( \mathbf{V} \) is estimated with sample moments and
\( \hat{\boldsymbol{\Omega}}=n^{-1}\sum_i\hat{e}_i^2\bz_{(i)}\bz_{(i)}\T \), as in @thm-het-sandwich, or under homoscedasticity by
\( \hat{\sigma}_e^2(\X\T\M_{\Z}\X)^{-1} \). Among estimators using linear combinations of the instruments, 2SLS has the
smallest asymptotic covariance under homoscedasticity (@exr-eiv-2sls-optimal).

::: {#exm-eiv-second-reading}
[The second reading as an instrument]

Instrument the first reading with the second, and let age and the constant instrument themselves. The IV estimates
are \( 0.0708 \) for blood pressure, with standard error \( 0.0085 \)
(\( 0.0081 \) from the sandwich), and \( 0.0054 \) for age (standard error
\( 0.0081 \)). The spurious age effect has gone. The first-stage \( F \) statistic for the second reading is
\( 160.1 \), a strong instrument. The correct residuals give \( \hat{\sigma}_e^2=1.438 \), estimating
\( \sigma^2+\beta^2\sigma_u^2=1.5184 \); the second-stage residuals give \( 1.376 \). The moment correction of
@exm-eiv-replicates gave \( 0.0757 \) with standard error \( 0.0076 \). It uses both readings
symmetrically, while IV uses one only as an instrument and pays in precision. But IV needs no estimate of the error
variance, and the instrument need only be correlated with \( X \), not a second unbiased reading of it.
:::

```{.python .run #cell-iv-replicate-iv}
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x

X = np.column_stack([np.ones(n), w[:, 0], age])    # regressors: first reading, age
Z = np.column_stack([np.ones(n), w[:, 1], age])    # instruments: second reading, age
b_iv = np.linalg.solve(Z.T @ X, Z.T @ y)
e = y - X @ b_iv                                    # residuals use X, not the first-stage fit
ZXinv = np.linalg.inv(Z.T @ X)
cov_h = e @ e / n * ZXinv @ (Z.T @ Z) @ ZXinv.T               # homoscedastic
cov_s = ZXinv @ ((Z * e[:, None] ** 2).T @ Z) @ ZXinv.T       # heteroscedasticity-robust
print("IV (1, x, age):", b_iv.round(4))
print("se, homoscedastic:", np.sqrt(np.diag(cov_h)).round(4))
print("se, sandwich:     ", np.sqrt(np.diag(cov_s)).round(4))

# first stage: regress the first reading on the instruments; F for the excluded instrument
g, rss1 = np.linalg.lstsq(Z, w[:, 0], rcond=None)[:2]
rss0 = np.sum((w[:, 0] - np.column_stack([np.ones(n), age]) @
               np.linalg.lstsq(np.column_stack([np.ones(n), age]), w[:, 0], rcond=None)[0]) ** 2)
F1 = (rss0 - rss1[0]) / (rss1[0] / (n - 3))
print(f"first-stage F = {F1:.1f}")
```

Wald's (1940) grouping estimator ([Section 7.1](../ch07-optimality/01-linear-unbiased.html)) is an IV estimator
whose instrument is a group indicator (@exr-eiv-wald-grouping).

## Weak instruments

If the instruments explain little of the endogenous regressor, the limit theory of @thm-eiv-iv is a poor
guide at realistic \( n \). The simplest case can be solved exactly.

::: {#prp-eiv-weak}
[The exact distribution with one instrument]

Let \( y_i=\beta x_i+e_i \) and \( x_i=\pi z_i+v_i \), \( i=1,\dots,n \), with fixed \( z_i \), \( S_{zz}=\sum_iz_i^2>0 \), and \( (e_i,v_i) \)
independent bivariate normal with means zero, standard deviations \( \sigma_e,\sigma_v \) and correlation \( \rho \),
\( \lvert\rho\rvert<1 \). Let \( \hat{\beta}=\sum_iz_iy_i/\sum_iz_ix_i \) and \( \mu=\pi\sqrt{S_{zz}}/\sigma_v \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\beta}-\beta \) has the same distribution as
   \[
   \frac{\sigma_e}{\sigma_v}\cdot\frac{\rho\eta+\sqrt{1-\rho^2}\,\xi}{\mu+\eta},\qquad \eta,\xi\iid\Normal(0,1).
   \]{#eq-eiv-weak-representation}

2. \( \E\lvert\hat{\beta}\rvert=\infty \) for every \( \mu \).

3. As \( \mu\to\infty \), \( (\mu\sigma_v/\sigma_e)(\hat{\beta}-\beta)\to\Normal(0,1) \) in distribution. If \( \mu=0 \),
   \( \hat{\beta}-\beta \) has median \( \rho\sigma_e/\sigma_v \), which is the large-sample bias of least squares when \( \pi=0 \).

4. \( (\sum_iz_ix_i)^2/(\sigma_v^2S_{zz})\sim\chi^2(1,\mu^2) \), with mean \( 1+\mu^2 \).
:::

:::

::: {.proof}
*(a)* Put \( \eta=\sum_iz_iv_i/(\sigma_v\sqrt{S_{zz}}) \) and \( \zeta=\sum_iz_ie_i/(\sigma_e\sqrt{S_{zz}}) \). They are jointly normal (@thm-mvn-linear)
with unit variances and covariance \( \sum_iz_i^2\rho\sigma_e\sigma_v/(\sigma_e\sigma_vS_{zz})=\rho \). By
@thm-mvn-conditional, \( \xi=(\zeta-\rho\eta)/\sqrt{1-\rho^2} \) is standard normal and independent of \( \eta \). Now
\( \sum_iz_ix_i=\pi S_{zz}+\sigma_v\sqrt{S_{zz}}\,\eta=\sigma_v\sqrt{S_{zz}}(\mu+\eta) \), and
\( \hat{\beta}-\beta=\sum_iz_ie_i/\sum_iz_ix_i=\sigma_e\zeta/\{\sigma_v(\mu+\eta)\} \).

*(b)* Given \( \eta \), the numerator \( \rho\eta+\sqrt{1-\rho^2}\,\xi \) is normal with standard deviation \( \sqrt{1-\rho^2}>0 \), so its
conditional absolute mean is at least \( c=\sqrt{2(1-\rho^2)/\pi} \). Hence
\( \E\lvert\hat{\beta}-\beta\rvert\ge(\sigma_e/\sigma_v)\,c\,\E\lvert\mu+\eta\rvert^{-1}=\infty \), because the density of \( \mu+\eta \) is continuous
and positive at zero, as in the proof of @prp-ci-ratio-no-mean. Finally
\( \lvert\hat{\beta}\rvert\ge\lvert\hat{\beta}-\beta\rvert-\lvert\beta\rvert \).

*(c)* \( (\mu\sigma_v/\sigma_e)(\hat{\beta}-\beta) \) has the law of \( (\rho\eta+\sqrt{1-\rho^2}\,\xi)/(1+\eta/\mu) \). The numerator is
standard normal and \( \eta/\mu\to0 \) in probability, so Slutsky's lemma applies. For \( \mu=0 \) the representation is
\( (\sigma_e/\sigma_v)(\rho+\sqrt{1-\rho^2}\,\xi/\eta) \), and \( \xi/\eta \) is symmetric about zero, so the median is
\( \rho\sigma_e/\sigma_v \). When \( \pi=0 \), \( x=v \) and least squares converges to \( \beta+\Cov(v,e)/\sigma_v^2=\beta+\rho\sigma_e/\sigma_v \).

*(d)* \( \sum_iz_ix_i/(\sigma_v\sqrt{S_{zz}})=\mu+\eta\sim\Normal(\mu,1) \), and its square is \( \chi^2(1,\mu^2) \), whose mean is
\( 1+\mu^2 \) (the convention of [Chapter 4](../ch04-quadratic-forms/index.html), without a factor \( \tfrac12 \)).
:::

The design and \( \pi \) enter only through the **concentration parameter** \( \mu^2 \), which acts as a sample size,
and by (d) the first-stage statistic estimates \( 1+\mu^2 \). Hence the rule of thumb of Staiger and Stock (1997):
an instrument is weak if the first-stage \( F \) is below about \( 10 \). A useless instrument (\( \mu=0 \)) gives an
estimator centred on the least squares bias, with Cauchy tails.

::: {#exm-eiv-weak}
[How weak is too weak?]

With \( n=200 \), \( \beta=1 \), \( \sigma_e=\sigma_v=1 \) and \( \rho=0.9 \), least squares tends to about \( 1.9 \) for a weak
instrument. In \( 20000 \) samples for each concentration:

| \( \mu^2 \) | median | quartiles | coverage of nominal 95% Wald interval |
|---|---|---|---|
| 1 | 1.298 | 0.711 to 1.646 | 0.843 |
| 10 | 1.008 | 0.740 to 1.181 | 0.921 |
| 100 | 1.002 | 0.929 to 1.064 | 0.947 |

With \( \mu^2=1 \) the median moves a third of the way to the least squares limit, a proportion \( 0.200 \)
of estimates fall outside \( [-1,3] \), and the Wald interval undercovers, as it still does at \( \mu^2=10 \)
([Figure 24.4.1](#fig-eiv-iv-weak)).
:::

::: {when-format="html"}
![**Figure 24.4.1.** Sampling distributions of the IV estimate of \( \beta=1 \) for \( \mu^2=1 \), \( 10 \), \( 100 \) on
\( [-1,3] \); dotted, the least squares limit.](iv_weak.svg){#fig-eiv-iv-weak width=78%}
:::

::: {when-format="pdf"}
![Sampling distributions of the IV estimate of \( \beta=1 \) for \( \mu^2=1 \), \( 10 \), \( 100 \) on \( [-1,3] \); dotted,
the least squares limit.](iv_weak.pdf){width=78%}
:::

```{.python .run #cell-iv-weak}
def iv_draws(mu2, reps, seed, n=200, rho=0.9):
    """IV estimates of beta = 1 and Wald-interval coverage when the concentration is mu2."""
    r = np.random.default_rng(seed)
    z = r.normal(size=n)
    z = z * np.sqrt(n / (z @ z))                  # fixed instrument with sum z^2 = n
    pi = np.sqrt(mu2 / n)
    ev = r.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=(reps, n))
    xx = pi * z + ev[:, :, 1]
    yy = 1.0 * xx + ev[:, :, 0]
    b = (yy @ z) / (xx @ z)
    res = yy - b[:, None] * xx
    se = np.sqrt((res ** 2).mean(axis=1) * n / (xx @ z) ** 2)
    return b, np.abs(b - 1) <= 1.96 * se


for mu2 in [1, 10, 100]:
    b, hit = iv_draws(mu2, 2000, 2440)
    print(f"mu^2 = {mu2:3d}: median {np.median(b):.3f}, "
          f"quartiles {np.percentile(b, 25):.3f} to {np.percentile(b, 75):.3f}, coverage {hit.mean():.3f}")
```

With many weak instruments, 2SLS is biased towards least squares, and Bound, Jaeger and Baker (1995) showed that
it can nearly reproduce least squares; we do not prove this. For inference, invert a test whose null distribution
does not depend on instrument strength, such as the Anderson–Rubin test (@exr-eiv-anderson-rubin). Its confidence
sets are of Fieller type (@prp-ci-fieller) and can be unbounded, as they must sometimes be (@prp-ci-calibration-unbounded).

::: {.warning}
Exogeneity cannot be checked in a just-identified model: all the equations \( \E(\bz e)=\bzero \) are used to
estimate \( \bbeta \). With more instruments, tests of overidentifying restrictions (Sargan, 1958) check only that
the instruments *agree*. An instrument is justified by knowledge of how the data arose; for replicates, by the
independence of their errors. Instruments for causal effects are taken up in
[Chapter 25](../ch25-causal-interpretation/index.html).
:::

## Exercises

### A. Check your understanding

::: {#exr-eiv-iv-simple}
[A1]

For \( \x=(1,W) \) and \( \bz=(1,Z) \), show that the IV slope is \( s_{zy}/s_{zw} \), the ratio of two sample covariances.
Show that the population ratio equals \( \beta \) when \( Z \) is correlated with \( X \) and
uncorrelated with \( U \) and \( \varepsilon \).
:::

::: {#exr-eiv-correlated-replicates}
[A2]

If the errors of two readings are correlated, \( \Cov(U_1,U_2)=\kappa \), find the limit of the IV slope in simple
regression that instruments \( W_1 \) with \( W_2 \). In which direction is it biased when \( \kappa>0 \)?
:::

::: {.solution}
\( \Cov(W_2,Y)=\beta\sigma_x^2 \) and \( \Cov(W_2,W_1)=\sigma_x^2+\kappa \), so the limit is \( \beta\sigma_x^2/(\sigma_x^2+\kappa) \). With
\( \kappa>0 \) it is attenuated, though less than least squares unless \( \kappa=\sigma_u^2 \). The instrument is not exogenous:
\( \Cov(W_2,\varepsilon-\beta U_1)=-\beta\kappa \).
:::

### B. Practice

::: {#exr-eiv-wald-grouping}
[B1]

Let \( g_i\in\{0,1\} \) be a group indicator, and let \( \bz_{(i)}=(1,g_i) \) instrument \( \x_{(i)}=(1,w_i) \). Show that the IV slope
is \( (\bar{y}_1-\bar{y}_0)/(\bar{w}_1-\bar{w}_0) \), the difference of group means of \( y \) divided by that of \( w \). This is Wald's
(1940) estimator. Show that it is consistent if \( g \) is independent of \( (U,\varepsilon) \) and \( \E(X\mid g) \) differs between the
groups. Explain why forming the groups by splitting at the median of \( w \) itself, as Wald proposed, does not satisfy
these conditions under classical error.
:::

::: {.solution}
By @exr-eiv-iv-simple the slope is \( s_{gy}/s_{gw} \). For a binary \( g \) with \( n_1 \) ones, \( \sum_i(g_i-\bar{g})y_i=n_1(\bar{y}_1-\bar{y}) \)
and \( \bar{y}_1-\bar{y}=(n_0/n)(\bar{y}_1-\bar{y}_0) \). The same holds for \( w \), and the factors cancel. Exogeneity needs
\( \Cov(g,\varepsilon-\beta U)=0 \), which holds when \( g \) is independent of \( (U,\varepsilon) \). Relevance needs \( \Cov(g,W)=\Cov(g,X)\ne0 \),
that is, different group means of \( X \). If the groups are formed from \( w \), membership depends on \( U \): readings with
large positive errors tend to land in the upper group, so \( \Cov(g,U)>0 \) and exogeneity fails. The grouping must come
from information independent of the measurement error, such as a second reading or a design variable.
:::

### C. Going deeper

::: {#exr-eiv-anderson-rubin}
[C1]

In @eq-eiv-structural-equation with \( \x=(1,W) \), let \( \bz=(1,\bz_1) \) with \( \bz_1 \) a \( k \)-vector of excluded instruments,
and suppose \( e \) is independent of \( \bz \) and \( \Normal(0,\sigma_e^2) \). To test \( H_0:\beta=\beta_0 \), regress
\( y-\beta_0w \) on \( \bz \) and compute the \( F \) statistic for \( \bz_1 \). Show that under \( H_0 \) it has an exact \( F(k,n-k-1) \)
distribution, whatever the strength of the instruments. This is the test of Anderson and Rubin (1949). Inverting
it gives a confidence set for \( \beta \). Show that for \( k=1 \) the set is of the form in @prp-ci-fieller.
:::

::: {.solution}
Under \( H_0 \), \( y-\beta_0w=\alpha+e \) with \( e \) normal and independent of \( \bz \): given \( \Z \), a normal linear model in which the
coefficients of \( \bz_1 \) are zero. The \( F \) statistic of @thm-glh-f-test is exactly \( F(k,n-k-1) \) given \( \Z \), hence
unconditionally, and \( w \) plays no role. For \( k=1 \), let \( \hat{a} \) and
\( \hat{b} \) be the slopes of \( y \) and of \( w \) on \( (1,z) \), let \( v=1/\sum_i(z_i-\bar{z})^2 \), and let \( s_{11} \), \( s_{12} \), \( s_{22} \)
be the residual mean squares and mean cross-product of \( y \) and \( w \) from these two regressions. The regression of
\( y-\beta_0w \) on \( (1,z) \) has slope \( \hat{a}-\beta_0\hat{b} \) and residual mean square
\( s_{11}-2\beta_0s_{12}+\beta_0^2s_{22} \), so \( F=t^2 \) and the acceptance region is
\[
(\hat{a}-\beta_0\hat{b})^2\le t_{n-2,\alpha/2}^2\,v\,\bigl(s_{11}-2\beta_0s_{12}+\beta_0^2s_{22}\bigr).
\]
This is @eq-ci-fieller with \( s^2v_{ij} \) replaced by \( v\,s_{ij} \): a quadratic inequality in \( \beta_0 \). Like Fieller's
set, it can be an interval, the complement of an interval, or the whole line.
:::

::: {#exr-eiv-2sls-optimal}
[C2]

Let \( \mathbf{F} \) be a fixed \( q\times p \) matrix with \( \bSigma_{zx}\T\mathbf{F} \) nonsingular, and consider the IV estimator with instruments
\( \Z\mathbf{F} \). Under homoscedasticity its asymptotic covariance is
\( \sigma_e^2(\mathbf{F}\T\bSigma_{zx})^{-1}\mathbf{F}\T\bSigma_{zz}\mathbf{F}(\bSigma_{zx}\T\mathbf{F})^{-1} \). Show that this is never smaller,
in the nonnegative definite ordering, than \( \sigma_e^2\mathbf{H}^{-1} \), with equality when \( \mathbf{F}=\bSigma_{zz}^{-1}\bSigma_{zx} \),
which is the population version of 2SLS.
:::
