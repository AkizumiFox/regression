# Sampling distributions of the estimates

A point estimate is of little use without a measure of its uncertainty. Under the second-moment assumptions we
know the mean and variance of \( \blambda\T\hbeta \) and the mean of \( s^2 \), but not their distributions. Under
normality we know everything. The estimates of the estimable functions are jointly normal, the residual
sum of squares is a scaled chi-squared variable, and the two are independent. These three facts are what
the whole of Part III rests on. Tests, confidence intervals, confidence ellipsoids and prediction intervals
are all built from them. The machinery was assembled in [Chapter 3](../ch03-multivariate-normal/index.html)
and [Chapter 4](../ch04-quadratic-forms/index.html). Here we apply it.

## The main theorem

::: {#thm-opt-sampling}
[Sampling distributions in the normal linear model]

Assume @eq-opt-normal-model with \( r=\rank(\X)<n \). Let \( \bLambda \) be \( p\times q \) with
\( \C(\bLambda)\subseteq\C(\X\T) \), let \( (\X\T\X)\ginv \) be any generalized inverse, and let \( s^2=\text{SSE}/(n-r) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \bLambda\T\hbeta\sim\Normal_q\bigl(\bLambda\T\bbeta,\ \sigma^2\bLambda\T(\X\T\X)\ginv\bLambda\bigr) \). In particular
   \( \hY=\M\Y\sim\Normal_n(\X\bbeta,\sigma^2\M) \), and if \( \rank(\X)=p \),
   \( \hbeta\sim\Normal_p\bigl(\bbeta,\sigma^2(\X\T\X)^{-1}\bigr) \).

2. \( \text{SSE}/\sigma^2=(n-r)s^2/\sigma^2\sim\chi^2(n-r) \).

3. \( \M\Y \) and \( \text{SSE} \) are independent. Hence \( \bLambda\T\hbeta \), and every function of the fitted values, is
   independent of \( s^2 \).

4. \( \E(s^2)=\sigma^2 \) and \( \Var(s^2)=2\sigma^4/(n-r) \).
:::
:::

::: {.proof}
Write \( \bLambda=\X\T\mathbf{P} \) for an \( n\times q \) matrix \( \mathbf{P} \). As in the proof of @thm-opt-gauss-markov(d),
\( \bLambda\T\hbeta=\mathbf{P}\T\M\Y \) for every least squares estimate \( \hbeta \).

(a) By @thm-mvn-linear, \( \mathbf{P}\T\M\Y \) is normal with mean \( \mathbf{P}\T\M\X\bbeta=\mathbf{P}\T\X\bbeta=\bLambda\T\bbeta \) and
covariance \( \sigma^2\mathbf{P}\T\M\M\T\mathbf{P}=\sigma^2\mathbf{P}\T\M\mathbf{P}=\sigma^2\bLambda\T(\X\T\X)\ginv\bLambda \), by @thm-proj-M-formula.
The special cases are \( \bLambda=\X\T \), with \( \mathbf{P}=\I \), and \( \bLambda=\I_p \) in full rank.

(b) \( \text{SSE}=\Y\T(\I-\M)\Y \), and \( \I-\M \) is symmetric and idempotent of rank \( n-r \)
(@thm-proj-sym-idem and @prp-proj-trace-rank). By @thm-qf-chisq(a),
\( \text{SSE}/\sigma^2\sim\chi^2\bigl(n-r,\norm{(\I-\M)\X\bbeta}^2/\sigma^2\bigr) \), and the noncentrality is zero
because \( (\I-\M)\X=\mathbf{0} \).

(c) Apply @thm-qf-indep-linear with \( \B=\M \), \( \A=\I-\M \) and \( \bSigma=\sigma^2\I \). Then
\( \B\bSigma\A=\sigma^2\M(\I-\M)=\mathbf{0} \), so \( \M\Y \) and \( \Y\T(\I-\M)\Y \) are independent. Functions of independent
random vectors are independent.

(d) A \( \chi^2(k) \) variable has mean \( k \) and variance \( 2k \) (@thm-qf-ncchisq(b) with \( \gamma=0 \)). So
\( \E\,\text{SSE}=\sigma^2(n-r) \) and \( \Var(\text{SSE})=2\sigma^4(n-r) \). Divide by \( n-r \) and \( (n-r)^2 \).
:::

The proof is short because the work was done earlier. It is worth seeing where each ingredient comes from
geometrically. The fitted vector \( \M\Y \) and the residual vector \( (\I-\M)\Y \) are projections of \( \Y \) onto
orthogonal subspaces. For a spherical normal vector, orthogonal projections are independent
(@thm-qf-orthogonal-projections). Everything about \( \bbeta \) lives in the first projection. The second is pure
noise, and its squared length, measured in units of \( \sigma^2 \), counts \( n-r \) independent squared standard
normals.

The independence in (c) has no analogue under the second-moment assumptions. If the errors are
independent with zero third moment, \( \blambda\T\hbeta \) and \( s^2 \) are still uncorrelated, but in general they
are not independent. The simulation at the end of the section shows the difference.

## Standardizing: the \( t \) statistic

The distribution in (a) involves the unknown \( \sigma \). Since \( \E(s^2)=\sigma^2 \), the matrix
\( s^2\bLambda\T(\X\T\X)\ginv\bLambda \) is an unbiased estimate of \( \Cov(\bLambda\T\hbeta) \), and in full rank
\( s^2(\X\T\X)^{-1} \) is unbiased for \( \Cov(\hbeta) \); this needs only the second-moment assumptions. The square
roots of its diagonal entries are the standard errors reported by regression software. Replacing \( \sigma \) by
\( s \) changes the normal distribution into a \( t \) distribution, and independence is what makes the change exact.

::: {#cor-opt-t}
[The \( t \) statistic for an estimable function]

Assume @eq-opt-normal-model with \( r<n \), let \( \blambda\in\C(\X\T) \) with \( \blambda\ne\bzero \), and put
\( \text{se}(\blambda\T\hbeta)=s\sqrt{\blambda\T(\X\T\X)\ginv\blambda} \). For every constant \( d \),
\[
T_d=\frac{\blambda\T\hbeta-d}{\text{se}(\blambda\T\hbeta)}\sim t(n-r,\delta),\qquad
\delta=\frac{\blambda\T\bbeta-d}{\sigma\sqrt{\blambda\T(\X\T\X)\ginv\blambda}} .
\]
In particular \( T=(\blambda\T\hbeta-\blambda\T\bbeta)/\text{se}(\blambda\T\hbeta)\sim t(n-r) \), and
\[
\Pr\bigl\{\lvert\blambda\T\hbeta-\blambda\T\bbeta\rvert\le t_{n-r,\alpha/2}\,\text{se}(\blambda\T\hbeta)\bigr\}=1-\alpha,
\]{#eq-opt-t-interval}

where \( t_{n-r,\alpha/2} \) is the upper \( \alpha/2 \) point of \( t(n-r) \).
:::

::: {.proof}
Put \( c^2=\blambda\T(\X\T\X)\ginv\blambda \). Then \( c>0 \), since \( c^2=\norm{\M\boldsymbol{\uprho}}^2 \) with
\( \X\T\M\boldsymbol{\uprho}=\blambda\ne\bzero \). By @thm-opt-sampling(a),
\( Z=(\blambda\T\hbeta-d)/(\sigma c)\sim\Normal(\delta,1) \). By (b), \( V=(n-r)s^2/\sigma^2\sim\chi^2(n-r) \), and by
(c) \( Z \) and \( V \) are independent. Now \( T_d=Z/\sqrt{V/(n-r)} \), which is \( t(n-r,\delta) \) by @def-qf-noncentral-t.
With \( d=\blambda\T\bbeta \), \( \delta=0 \), and @eq-opt-t-interval restates
\( \Pr\{\lvert T\rvert\le t_{n-r,\alpha/2}\}=1-\alpha \).
:::

@eq-opt-t-interval is the familiar confidence interval \( \blambda\T\hbeta\pm t_{n-r,\alpha/2}\,\text{se} \), and
\( T_d \) is the statistic for testing \( \blambda\T\bbeta=d \). Its noncentrality \( \delta \) governs the power. [Chapter 11](../ch11-general-linear-hypothesis/index.html)
and [Chapter 12](../ch12-intervals-and-bands/index.html) develop both. For the air-flow coefficient in the stack loss regression this interval was
\( (0.431,\ 1.000) \) (@exm-opt-stackloss-profile). Two features of the corollary are easy to
overlook. The degrees of freedom are \( n-r \), the dimension of the residual space. They are not \( n \), and when
\( \X \) is rank deficient they are not \( n-p \) either. And the interval is exact for every sample size, which no
large-sample argument delivers.

## Several functions at once

For a vector of estimable functions, the natural summary of @thm-opt-sampling(a) is a quadratic form. The
covariance matrix \( \bLambda\T(\X\T\X)\ginv\bLambda \) may be singular when the columns of \( \bLambda \) are linearly
dependent, so a generalized inverse is needed.

::: {#cor-opt-quadratic}
[A chi-squared pivot for several functions]

In the setting of @thm-opt-sampling, let \( \mathbf{W}=\bLambda\T(\X\T\X)\ginv\bLambda \), let \( \mathbf{W}\ginv \) be any generalized
inverse of \( \mathbf{W} \), and put \( \mathbf{D}=\bLambda\T\hbeta-\bLambda\T\bbeta \). Then
\[
Q=\frac{\mathbf{D}\T\mathbf{W}\ginv\mathbf{D}}{\sigma^2}\sim\chi^2\bigl(\rank(\bLambda)\bigr),
\]
\( Q \) does not depend on the choice of \( \mathbf{W}\ginv \) (with probability one), and \( Q \) is independent of SSE. Consequently
\[
\frac{\mathbf{D}\T\mathbf{W}\ginv\mathbf{D}/\rank(\bLambda)}{s^2}\sim F\bigl(\rank(\bLambda),\ n-r\bigr).
\]
:::

::: {.proof}
By @thm-opt-sampling(a), \( \mathbf{D}\sim\Normal_q(\bzero,\sigma^2\mathbf{W}) \). The quadratic form is unchanged if
\( \mathbf{W}\ginv \) is replaced by its symmetric part \( \frac12(\mathbf{W}\ginv+\mathbf{W}^{-\top}) \), which is again a generalized
inverse of the symmetric matrix \( \mathbf{W} \) (@prp-mat-ginverse-props(d)). So we may assume \( \mathbf{W}\ginv \) symmetric.
Apply @thm-qf-chisq(c) with
\( \bSigma=\sigma^2\mathbf{W} \), \( \bmu=\bzero \) and \( \A=\mathbf{W}\ginv/\sigma^2 \). Conditions (ii) and (iii) hold trivially because
\( \bmu=\bzero \). Condition (i) is
\( \bSigma\A\bSigma\A\bSigma=\sigma^2\mathbf{W}\mathbf{W}\ginv\mathbf{W}\mathbf{W}\ginv\mathbf{W}=\sigma^2\mathbf{W}=\bSigma\A\bSigma \), by two uses of
\( \mathbf{W}\mathbf{W}\ginv\mathbf{W}=\mathbf{W} \). So \( Q\sim\chi^2(\tr(\A\bSigma)) \), where \( \tr(\A\bSigma)=\tr(\mathbf{W}\ginv\mathbf{W})=\rank(\mathbf{W}) \),
because \( \mathbf{W}\ginv\mathbf{W} \) is idempotent with the same rank as \( \mathbf{W} \) (@prp-mat-ginverse-props(a)), and the trace
of an idempotent matrix is its rank (@prp-mat-idempotent-basic).

It remains to show \( \rank(\mathbf{W})=\rank(\bLambda) \). With \( \bLambda=\X\T\mathbf{P} \),
\( \mathbf{W}=\mathbf{P}\T\M\mathbf{P}=(\M\mathbf{P})\T(\M\mathbf{P}) \), so \( \rank\mathbf{W}=\rank(\M\mathbf{P}) \). Also \( \bLambda=\X\T\mathbf{P}=\X\T\M\mathbf{P} \). If
\( \X\T\M\mathbf{P}\bu=\bzero \), then \( \M\mathbf{P}\bu\in\C(\X)\cap\C(\X)\perpc=\{\bzero\} \), so \( \X\T \) is one-to-one on \( \C(\M\mathbf{P}) \)
and \( \rank\bLambda=\rank(\M\mathbf{P}) \).

For invariance: \( \mathbf{D}\in\C(\mathbf{W}) \) with probability one by @thm-rv-cov-nnd(d). If \( \mathbf{D}=\mathbf{W}\bu \), then
\( \mathbf{D}\T\mathbf{W}\ginv\mathbf{D}=\bu\T\mathbf{W}\mathbf{W}\ginv\mathbf{W}\bu=\bu\T\mathbf{W}\bu \), whatever \( \mathbf{W}\ginv \) is. Independence follows from
@thm-opt-sampling(c), since \( \mathbf{D} \) is a function of \( \M\Y \). The \( F \) ratio is then a ratio of independent chi-squared
variables divided by their degrees of freedom (@def-qf-noncentral-f).
:::

The set of \( \bLambda\T\bbeta \) values for which this \( F \) ratio stays below an \( F \) quantile is a confidence
ellipsoid. [Chapter 12](../ch12-intervals-and-bands/index.html) studies these regions, and [Chapter 11](../ch11-general-linear-hypothesis/index.html) turns the same statistic into the test of the
general linear hypothesis. Here we note only that the degrees of freedom in the numerator are the number
of *linearly independent* functions being estimated, not the number of rows of \( \bLambda\T \).

::: {#exm-opt-sampling-simulation}
[A rank-deficient model, simulated]

Consider a one-way layout with three groups of sizes \( 4 \), \( 5 \), \( 6 \), an intercept, one indicator column per group,
and a covariate. So \( n=15 \) and \( p=5 \), but the rank is \( r=4 \), because the indicators add up to the
intercept. The contrast \( \blambda\T\bbeta=\alpha_1-\alpha_2 \), with \( \blambda=(0,1,-1,0,0)\T \), is estimable. With
\( \sigma=2 \), @thm-opt-sampling gives
\( \Var(\blambda\T\hbeta)=\sigma^2\blambda\T(\X\T\X)\ginv\blambda=1.9007 \), computed with the Moore–Penrose
inverse. In \( 100{,}000 \) simulated data sets the variance of the estimated contrast is
\( 1.9017 \), and the average of \( \text{SSE}/\sigma^2 \) is \( 11.001 \), against
\( n-r=11 \). The sample correlation between the contrast estimate and \( s^2 \) is
\( 0.0036 \). The standardized contrast \( T \) falls within \( \pm2.201 \), the
\( 97.5\% \) point of \( t(11) \), in a proportion \( 0.9498 \) of the samples.
[Figure 7.5.1](#fig-opt-sampling) compares the simulated distributions with the theory. Near the centre the
\( t(11) \) and \( \Normal(0,1) \) densities are hard to tell apart, so panel (b) uses a logarithmic density scale, on
which the tails separate: the simulated points follow \( t(11) \), and the normal density falls away too fast. The
tail mass says the same in numbers. A proportion \( 0.0758 \) of the samples have \( |T|>1.96 \), against \( 0.05 \)
under \( \Normal(0,1) \) and \( 0.0758 \) under \( t(11) \). A nominal \( 95\% \) interval built with the normal
point \( 1.96 \) would therefore undercover.
:::

```{.python .run #cell-sampling-distributions-design}
import numpy as np
from scipy import stats

groups = np.repeat([0, 1, 2], [4, 5, 6])
z = np.array([2.1, 3.4, 1.8, 4.0, 2.9, 3.3, 5.1, 2.2, 3.8, 1.5, 4.4, 2.7, 3.9, 5.0, 3.1])
n = len(z)
X = np.column_stack([np.ones(n), np.eye(3)[groups], z])   # p = 5 columns, rank 4
r = np.linalg.matrix_rank(X)
G = np.linalg.pinv(X.T @ X)                                # one generalized inverse
M = X @ G @ X.T
lam = np.array([0.0, 1.0, -1.0, 0.0, 0.0])                 # group 1 minus group 2: estimable
beta, sigma = np.array([1.0, 0.5, -0.5, 0.0, 0.8]), 2.0
var_lam = sigma ** 2 * lam @ G @ lam                       # Var(lambda^T beta_hat)
print(f"n = {n}, p = {X.shape[1]}, rank = {r}, Var(lambda^T beta_hat) = {var_lam:.4f}")
```

```{.python .run #cell-sampling-distributions-simulate}
rng = np.random.default_rng(75)
reps = 100_000
Y = X @ beta + sigma * rng.normal(size=(reps, n))          # one data set per row
est = Y @ (X @ G @ lam)                                    # lambda^T beta_hat = lambda^T G X^T y
sse = np.sum((Y - Y @ M) ** 2, axis=1)                     # y^T (I - M) y
s2 = sse / (n - r)
t = (est - lam @ beta) / np.sqrt(s2 * (lam @ G @ lam))
print(f"mean {est.mean():.4f}, variance {est.var():.4f}")
print(f"E(SSE)/sigma^2 = {np.mean(sse) / sigma**2:.3f}   (n - r = {n - r})")
print(f"corr(estimate, s^2) = {np.corrcoef(est, s2)[0, 1]:.4f}")
print("KS p-value, SSE/sigma^2 vs chi^2(n-r):", round(stats.kstest(sse / sigma**2, "chi2", args=(n - r,)).pvalue, 3))
print("KS p-value, t vs t(n-r):              ", round(stats.kstest(t, "t", args=(n - r,)).pvalue, 3))
```

::: {when-format="html"}
![**Figure 7.5.1.** Sampling distributions in a rank-deficient normal linear model (\( n=15 \), \( r=4 \); 100,000
simulated data sets). (a) \( \text{SSE}/\sigma^2 \) with the \( \chi^2(11) \) density. (b) The standardized contrast (binned density,
points) with the \( t(11) \) and standard normal densities, on a logarithmic scale. (c) The contrast estimate against \( s^2 \): no relationship, as independence
requires.](sampling_distributions.svg){#fig-opt-sampling width=100%}
:::

::: {when-format="pdf"}
![Sampling distributions in a rank-deficient normal linear model (\( n=15 \), \( r=4 \); 100,000
simulated data sets). (a) \( \text{SSE}/\sigma^2 \) with the \( \chi^2(11) \) density. (b) The standardized contrast (binned density,
points) with the \( t(11) \) and standard normal densities, on a logarithmic scale. (c) The contrast estimate against \( s^2 \): no relationship, as independence
requires.](sampling_distributions.pdf){width=100%}
:::

## Without normality

Which parts of @thm-opt-sampling survive when the errors are merely uncorrelated with common variance?
The first two moments of \( \blambda\T\hbeta \) do, by @thm-opt-gauss-markov, and so does \( \E(s^2)=\sigma^2 \). The
normality of \( \blambda\T\hbeta \) is lost, but it often returns approximately. \( \blambda\T\hbeta \) is a weighted sum of
independent errors, and a central limit theorem applies when no single observation dominates the weights.
[Chapter 19](../ch19-theory-of-departures/index.html) makes this precise (@thm-dep-nonnormal). The distribution of SSE is less forgiving, because its variance depends on the
fourth moment of the errors.

::: {#prp-opt-var-sse}
[Variance of SSE without normality]

Let the errors be independent with mean zero, variance \( \sigma^2 \), third moment \( \mu_3 \) and fourth moment
\( \mu_4 \), and let \( h_{ii} \) be the leverages. Then
\[
\Var(\text{SSE})=2\sigma^4(n-r)+(\mu_4-3\sigma^4)\sum_{i=1}^n(1-h_{ii})^2 .
\]
:::

::: {.proof}
This is @prp-lm-var-s2 with \( r \) in place of \( p \). That proof applies
@thm-rv-quadform-variance to \( \text{SSE}=\Y\T(\I-\M)\Y \) and uses only that
\( \A=\I-\M \) is idempotent with \( \tr(\A^2)=\tr(\A)=n-r \) and diagonal entries
\( 1-h_{ii} \), and that \( \A\X\bbeta=\bzero \), which kills the terms carrying
\( \boldsymbol{\uptheta} \) and \( \mu_3 \). None of it needs full column rank.
:::

The factor \( \mu_4-3\sigma^4 \) is \( \sigma^4 \) times the excess kurtosis. For normal errors it vanishes and we recover
\( 2\sigma^4(n-r) \). For heavy-tailed errors it is positive, and SSE is more variable than the chi-squared theory says.
Intervals for \( \sigma^2 \) built from @thm-opt-sampling(b) are then too short, and this does not improve with
sample size: the ratio of the true variance to the nominal one tends to
\( 1+\tfrac12(\mu_4/\sigma^4-3) \) when the leverages are small. The \( t \) interval for \( \blambda\T\bbeta \) is much
less affected, because \( s^2 \) enters only through \( s \), and \( s \) converges to \( \sigma \) whatever the kurtosis.

::: {#exm-opt-laplace}
[Laplace errors]

Rerun the simulation of @exm-opt-sampling-simulation with Laplace (double-exponential) errors of the same variance
\( \sigma^2=4 \). For these \( \mu_4=6\sigma^4 \). The contrast still has variance \( 1.9007 \). For SSE,
\( \sum_i(1-h_{ii})^2=8.138 \), and @prp-opt-var-sse gives
\( \Var(\text{SSE})=742.6 \), more than twice the normal-theory value
\( 352.0 \). The simulated variance is \( 741.5 \). Independence fails too.
The contrast and SSE are still uncorrelated, because the errors are symmetric, but the *squared* error of the
contrast has correlation \( 0.211 \) with SSE. Samples with a few large errors inflate both.
:::

```{.python .run #cell-sampling-distributions-laplace}
E_lap = rng.laplace(scale=sigma / np.sqrt(2), size=(reps, n))   # variance sigma^2, mu_4 = 6 sigma^4
Y_lap = X @ beta + E_lap
est_lap = Y_lap @ (X @ G @ lam)
sse_lap = np.sum((Y_lap - Y_lap @ M) ** 2, axis=1)
h = np.diag(np.eye(n) - M)
var_sse_lap = (6 - 3) * sigma**4 * h @ h + 2 * sigma**4 * (n - r)
print(f"Var(SSE): normal theory {2 * sigma**4 * (n - r):.1f}, Laplace formula {var_sse_lap:.1f}, "
      f"simulated {sse_lap.var():.1f}")
print(f"corr(squared error of estimate, SSE) under Laplace errors: "
      f"{np.corrcoef((est_lap - lam @ beta) ** 2, sse_lap)[0, 1]:.3f}")
```

In the listing, `h` holds the diagonal of \( \I-\M \), that is, \( 1-h_{ii} \).

## Exercises

### A. Check your understanding

::: {#exr-opt-fit-residual-independent}
[A1]

Under @eq-opt-normal-model, find the joint distribution of \( (\hY,\he) \). Show that \( \hY \) and \( \he \) are independent,
and that \( \he\sim\Normal_n(\bzero,\sigma^2(\I-\M)) \) whatever \( \bbeta \) is. Is \( \he \) independent of SSE?
:::

::: {#exr-opt-slope-t}
[A2]

In simple linear regression, write out \( T \) of @cor-opt-t for the slope, in terms of \( S_{xx} \) and \( s \), and state its
distribution. What is \( r \)?
:::

### B. Practice

::: {#exr-opt-independent-contrasts}
[B1]

Let \( \blambda_1,\blambda_2\in\C(\X\T) \). Show that
\( \Cov(\blambda_1\T\hbeta,\blambda_2\T\hbeta)=\sigma^2\blambda_1\T(\X\T\X)\ginv\blambda_2 \), and that under normality the two
estimates are independent iff \( \blambda_1\T(\X\T\X)\ginv\blambda_2=0 \). In a one-way layout with \( a \) groups of \( m \)
observations each and cell-means parameterization \( \E(Y_{ij})=\mu_i \), show that the estimates of
\( \sum_ic_i\mu_i \) and \( \sum_id_i\mu_i \) are independent iff \( \sum_ic_id_i=0 \).
:::

::: {.solution}
The covariance is the off-diagonal entry of @thm-opt-sampling(a) with \( \bLambda=[\blambda_1,\blambda_2] \). Under
normality the pair is bivariate normal, so zero covariance is equivalent to independence
(@thm-mvn-independence). In the cell-means model \( \X\T\X=m\I_a \), so
\( \blambda_1\T(\X\T\X)^{-1}\blambda_2=\sum_ic_id_i/m \), which vanishes iff \( \sum_ic_id_i=0 \).
:::

::: {#exr-opt-prediction-t}
[B2]

A new observation \( Y_0=\x_0\T\bbeta+\varepsilon_0 \) will be taken, with \( \varepsilon_0\sim\Normal(0,\sigma^2) \)
independent of \( \Y \) and \( \x_0\in\C(\X\T) \). Show that
\[
\frac{Y_0-\x_0\T\hbeta}{s\sqrt{1+\x_0\T(\X\T\X)\ginv\x_0}}\sim t(n-r),
\]
and deduce a prediction interval for \( Y_0 \). Why is the extra \( 1 \) in the denominator never negligible?
:::

::: {.solution}
\( Y_0 \) is independent of \( \Y \), so \( Y_0-\x_0\T\hbeta \) is normal with mean \( \x_0\T\bbeta-\x_0\T\bbeta=0 \) and variance
\( \sigma^2+\sigma^2\x_0\T(\X\T\X)\ginv\x_0 \). It is a function of \( (Y_0,\M\Y) \), which is independent of SSE. Dividing
by \( \sigma\sqrt{1+\x_0\T(\X\T\X)\ginv\x_0} \) gives a standard normal, independent of \( (n-r)s^2/\sigma^2\sim\chi^2(n-r) \),
and the ratio is \( t(n-r) \) by @def-qf-noncentral-t. The interval is
\( \x_0\T\hbeta\pm t_{n-r,\alpha/2}\,s\sqrt{1+\x_0\T(\X\T\X)\ginv\x_0} \). The \( 1 \) is the variance of the new error
itself, which no amount of data about \( \bbeta \) can reduce. The second term tends to zero as the design grows. The
first does not.
:::

::: {#exr-opt-sigma-interval}
[B3]

Use @thm-opt-sampling(b) to construct an exact \( 1-\alpha \) confidence interval for \( \sigma^2 \). Then, using
@prp-opt-var-sse with small leverages, explain why its coverage does not tend to \( 1-\alpha \) as \( n\to\infty \) when the
errors have nonzero excess kurtosis.
:::

### C. Going deeper

::: {#exr-opt-cochran-route}
[C1]

Give a second proof of @thm-opt-sampling(b) and (c) using @thm-qf-orthogonal-projections with the two projections
\( \M \) and \( \I-\M \). Extend it: if \( \C(\X_0)\subseteq\C(\X) \), show that \( \norm{(\M-\Mo)\Y}^2 \) and SSE are independent, and
find the distribution of \( \norm{(\M-\Mo)\Y}^2/\sigma^2 \) when \( \E(\Y)=\X\bbeta \) is arbitrary.
:::

::: {#exr-opt-t-power}
[C2]

For the test that rejects \( \blambda\T\bbeta=d \) when \( \lvert T_d\rvert>t_{n-r,\alpha/2} \), show that the power depends on the
parameters only through \( \lvert\delta\rvert \) of @cor-opt-t. For a fixed design, how does the power change when the sample is
replicated \( m \) times (\( \X \) replaced by \( m \) stacked copies)?
:::
