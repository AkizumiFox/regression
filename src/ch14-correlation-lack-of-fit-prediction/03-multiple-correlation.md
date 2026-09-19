# Multiple correlation

The population squared multiple correlation
\[
\rho^2_{Y\cdot X}=\frac{\boldsymbol{\sigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\sigma}_{XY}}{\sigma_Y^2}
\]
is the largest squared correlation between \( Y \) and a linear combination of \( \mathbf{X} \) (@prp-rv-multiple-correlation). [Section 9.5](../ch09-sums-of-squares/05-r-squared.html) studied its sample
counterpart \( R^2 \) for a fixed design (@thm-ss-r2-max-correlation, @thm-ss-r2-null). With random rows,
\( R^2 \) estimates \( \rho^2_{Y\cdot X} \). This section finds its exact distribution under multivariate normality,
and uses it to test \( \rho^2_{Y\cdot X}=0 \), to build an interval and to measure the bias.

Throughout, \( (\mathbf{X}_i\T,Y_i) \), \( i=1,\dots,n \), are independent draws from a \( (k+1) \)-variate normal
distribution with \( \bSigma_{XX} \) positive definite and \( \sigma_Y^2>0 \), and \( n\ge k+2 \). By
@eq-mvn-random-x-model, \( Y_i=\alpha+\bbeta\T\mathbf{X}_i+e_i \) with independent \( e_i\sim\Normal(0,\sigma^2) \) independent of the
regressors, where
\[
\sigma^2=\sigma_Y^2\bigl(1-\rho^2_{Y\cdot X}\bigr),\quad
\bbeta\T\bSigma_{XX}\bbeta=\sigma_Y^2\rho^2_{Y\cdot X},\quad\text{so}\quad
\frac{\bbeta\T\bSigma_{XX}\bbeta}{\sigma^2}=\frac{\rho^2_{Y\cdot X}}{1-\rho^2_{Y\cdot X}} .
\]{#eq-cor-signal-noise}

The ratio on the right is the population *signal-to-noise ratio*. We write \( \rho^2 \) for
\( \rho^2_{Y\cdot X} \) when no confusion can arise.

## The sample multiple correlation as an estimator

Let \( \bS_{XX} \), \( \mathbf{s}_{XY} \) and \( s_Y^2 \) be the sample covariances, with any common divisor. The fit with an
intercept has \( \text{SSR}=\hbeta\T\W\hbeta \), where \( \W \) is the centred cross-product matrix of the
regressors and \( \hbeta=\W^{-1}\sum_i(\mathbf{X}_i-\bar{\mathbf{X}})(Y_i-\bar Y) \) (@thm-proj-fwl with \( \X_2=\bone \)). Dividing by
\( \text{SST} \) gives
\[
R^2=\frac{\mathbf{s}_{XY}\T\bS_{XX}^{-1}\mathbf{s}_{XY}}{s_Y^2},
\]{#eq-cor-r2-plugin}

the population formula with sample moments in place of population moments. Under normality it is
more than a plug-in.

::: {#prp-cor-mvn-mle}
[Maximum likelihood with normal rows]

Suppose the sample covariance matrix of the rows is nonsingular. The maximum likelihood estimates
of \( (\alpha,\bbeta,\sigma^2) \) are the least squares estimates and \( \text{SSE}/n \), and the maximum likelihood
estimate of \( \rho^2_{Y\cdot X} \) is \( R^2 \).
:::

::: {.proof}
The joint density of a row factorizes as \( f(\mathbf{x};\bmu_X,\bSigma_{XX})\,f(y\mid\mathbf{x};\alpha,\bbeta,\sigma^2) \), a normal
density for the regressors times the conditional normal density of @thm-mvn-conditional. The map
from \( (\bmu,\bSigma) \) to \( (\bmu_X,\bSigma_{XX},\alpha,\bbeta,\sigma^2) \) is one-to-one onto the set where \( \bSigma_{XX} \) is positive
definite and \( \sigma^2>0 \), with inverse \( \boldsymbol{\sigma}_{XY}=\bSigma_{XX}\bbeta \), \( \sigma_Y^2=\sigma^2+\bbeta\T\bSigma_{XX}\bbeta \),
\( \mu_Y=\alpha+\bbeta\T\bmu_X \). So the log-likelihood is a sum of two terms with separate parameters, and
each can be maximized on its own. The conditional term is the likelihood of a normal linear model
with model matrix \( [\bone,\X_1] \), maximized by least squares and \( \text{SSE}/n \) (@thm-opt-mle).

In the marginal term the quadratic form is minimized over \( \bmu_X \) at \( \bar{\mathbf{x}} \), leaving
\( -\tfrac n2\bigl[\log\det\bSigma_{XX}+\tr(\bSigma_{XX}^{-1}\bS)\bigr] \) with \( \bS=\W/n \) positive definite. If
\( \lambda_1,\dots,\lambda_k>0 \) are the eigenvalues of \( \bS^{1/2}\bSigma_{XX}^{-1}\bS^{1/2} \), then
\[
\log\det\bSigma_{XX}+\tr(\bSigma_{XX}^{-1}\bS)-\log\det\bS-k=\sum_{j=1}^k(\lambda_j-\log\lambda_j-1)\ge0,
\]
with equality iff every \( \lambda_j=1 \), that is, iff \( \bSigma_{XX}=\bS \). By invariance of maximum likelihood
under reparameterization, the estimate of
\( \rho^2=\bbeta\T\bSigma_{XX}\bbeta/(\sigma^2+\bbeta\T\bSigma_{XX}\bbeta) \) is
\( \hbeta\T\bS\hbeta/(\text{SSE}/n+\hbeta\T\bS\hbeta)=\text{SSR}/(\text{SSE}+\text{SSR})=R^2 \).
:::

The rest of the section concerns the sampling distribution of this estimator.

## The distribution of R-squared

::: {#thm-cor-r2-distribution}
[Distribution of \( R^2 \) with normal rows]

Under the assumptions above, let \( F=\dfrac{n-k-1}{k}\,\dfrac{R^2}{1-R^2} \) and
\( \lambda=\rho^2/(1-\rho^2) \).

::: {.enumerate options="label=(\alph*)"}
1. Given the regressors, \( F\sim F\bigl(k,n-k-1,\gamma\bigr) \) with \( \gamma=\bbeta\T\W\bbeta/\sigma^2 \), where \( \W \) is the
   centred cross-product matrix of the regressors.

2. If \( \rho^2=0 \), then \( F\sim F(k,n-k-1) \) and \( R^2\sim\mathrm{Beta}\bigl(\tfrac k2,\tfrac{n-k-1}2\bigr) \), independently of
   the regressors. The overall \( F \) test of the regression is an exact size-\( \alpha \) test of
   \( \rho^2=0 \), which under joint normality is the hypothesis that \( Y \) is independent of \( \mathbf{X} \).

3. If \( \rho^2>0 \), then \( \gamma=\lambda V \) with \( V\sim\chi^2(n-1) \), and for \( 0<c<1 \)
   \[
\Pr(R^2\le c)=\E\Bigl[\Pr\bigl\{F(k,n-k-1,\lambda V)\le\tfrac{(n-k-1)c}{k(1-c)}\bigr\}\Bigr].
\]{#eq-cor-r2-mixture}

4. For each \( c\in(0,1) \), \( \Pr(R^2>c) \) is a strictly increasing function of \( \rho^2\in[0,1) \).
:::

:::

::: {.proof}
(a) The centred regressor matrix has rank \( k \) with probability one, because \( \W\sim W_k(n-1,\bSigma_{XX}) \) (@lem-cor-centred-wishart) is nonsingular when \( n-1\ge k \) (the nonsingularity step in the proof of
@lem-cor-inverse-wishart needs only \( m\ge k \)). Given the regressors, the model is a normal linear model
with \( r=k+1 \) and mean \( \alpha\bone+\X_1\bbeta \). By @thm-ss-r2-null(a), \( F\sim F(k,n-k-1,\gamma) \) with
\( \gamma=\norm{(\M-\bP_1)(\alpha\bone+\X_1\bbeta)}^2/\sigma^2 \). Since \( (\M-\bP_1)\bone=\bzero \) and
\( (\M-\bP_1)\X_1\bbeta=(\I-\bP_1)\X_1\bbeta \), this is \( \bbeta\T\X_1\T(\I-\bP_1)\X_1\bbeta/\sigma^2=\bbeta\T\W\bbeta/\sigma^2 \).

(b) \( \rho^2=0 \) means \( \boldsymbol{\sigma}_{XY}=\bzero \), so \( \bbeta=\bzero \) and \( \gamma=0 \) for every design. By
@lem-cor-conditioning, \( F\sim F(k,n-k-1) \) independently of the regressors, and so is
\( R^2=kF/(kF+n-k-1) \), which is \( \mathrm{Beta}(k/2,(n-k-1)/2) \) by @thm-ss-r2-null(b). For jointly normal
variables, zero covariance is independence (@thm-mvn-independence).

(c) By the proof of @lem-cor-centred-wishart, \( \W=(\mathbf{K}\X_1)\T(\mathbf{K}\X_1) \), where the \( n-1 \) rows of
\( \mathbf{K}\X_1 \) are independent \( \Normal_k(\bzero,\bSigma_{XX}) \). So \( \mathbf{K}\X_1\bbeta \) has independent
\( \Normal(0,\bbeta\T\bSigma_{XX}\bbeta) \) entries, and \( \bbeta\T\W\bbeta=\norm{\mathbf{K}\X_1\bbeta}^2=(\bbeta\T\bSigma_{XX}\bbeta)\,V \) with
\( V\sim\chi^2(n-1) \). By @eq-cor-signal-noise, \( \gamma=\lambda V \). Since \( R^2\le c \) iff
\( F\le(n-k-1)c/(k(1-c)) \), taking expectations of the conditional probability from (a) gives @eq-cor-r2-mixture.

(d) Fix \( v>0 \) and \( f>0 \). By @thm-qf-f-power, \( \Pr\{F(k,n-k-1,\lambda v)>f\} \) is strictly increasing in
\( \lambda \), and \( \lambda=\rho^2/(1-\rho^2) \) is strictly increasing in \( \rho^2 \). If \( \rho_1^2<\rho_2^2 \), the conditional
tail probabilities therefore differ by a positive amount for every value of \( V>0 \), and so do their
expectations.
:::

Part (c) is Fisher's (1928) distribution of the multiple correlation coefficient, reached here
without any integration. Given the design, \( R^2 \) is a transformed noncentral \( F \) variable, and the
design enters only through the single number \( \bbeta\T\W\bbeta \), which is a scaled chi-squared variable
because the regressors are normal. [Figure 14.3.1](#fig-cor-r2-distribution) shows the resulting
densities for \( n=25 \) and \( k=4 \), computed from @eq-cor-r2-mixture by averaging over a grid of
quantiles of \( V \).

::: {#exm-cor-r2-bias}
[How biased is R-squared?]

With \( n=25 \), \( k=4 \) and \( 100000 \) simulated samples, the averages of \( R^2 \) are
\( 0.1669 \), \( 0.4026 \) and \( 0.6500 \) for \( \rho^2=0 \), \( 0.3 \) and \( 0.6 \), against exact values
\( 0.1667 \), \( 0.4026 \) and \( 0.6498 \). The exact values come from the formula
\( \E(R^2)=1-\frac{n-k-1}{n-1}(1-\rho^2)\,{}_2F_1\bigl(1,1;\tfrac{n+1}2;\rho^2\bigr) \), which follows from Fisher's (1928) distribution of \( R^2 \) (see Olkin and Pratt 1958), where
\( {}_2F_1 \) is the hypergeometric function; it can also be computed from @eq-cor-r2-mixture. The bias
shrinks as \( \rho^2 \) grows. The adjusted \( R^2 \) of
@eq-ss-adjusted-r2 averages \( 0.0003 \), \( 0.2831 \) and \( 0.5800 \): unbiased under the null (@prp-ss-adjusted-r2(e)), but biased *downwards* when \( \rho^2>0 \). Neither is unbiased everywhere.
The standard deviations of \( R^2 \) are \( 0.104 \), \( 0.146 \) and \( 0.118 \): with
\( 25 \) cases, a single \( R^2 \) says little about \( \rho^2 \).
:::

::: {when-format="html"}
![**Figure 14.3.1.** The law of \( R^2 \) with normal rows, \( n=25 \), \( k=4 \), for \( \rho^2=0 \), \( 0.3 \) and
\( 0.6 \) (dashed). Step outlines: simulation. Curves: @eq-cor-r2-mixture.](r2_distribution.svg){#fig-cor-r2-distribution width=72%}
:::

::: {when-format="pdf"}
![The law of \( R^2 \) with normal rows, \( n=25 \), \( k=4 \), for \( \rho^2=0 \), \( 0.3 \) and
\( 0.6 \) (dashed). Step outlines: simulation. Curves: @eq-cor-r2-mixture.](r2_distribution.pdf){width=72%}
:::

```{.python .run #cell-multiple-correlation-mixture}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

def r2_cdf(r, rho2, n, k, m=4000):
    """P(R^2 <= r) for iid normal rows: average the conditional noncentral F over V."""
    df2 = n - k - 1
    f = df2 * r / (k * (1 - r))
    if rho2 == 0:
        return stats.f.cdf(f, k, df2)
    V = stats.chi2.ppf((np.arange(m) + 0.5) / m, n - 1)      # equal-probability nodes
    return np.mean(stats.ncf.cdf(f, k, df2, rho2 / (1 - rho2) * V))
```

Part (d) says that the overall \( F \) test is unbiased and has power increasing in \( \rho^2 \), as a test of
\( \rho^2=0 \) against random-regressor alternatives. For \( n=25 \) and \( k=4 \) it rejects when
\( R^2>0.364 \), and its power is \( 0.182 \), \( 0.377 \), \( 0.599 \) and
\( 0.922 \) at \( \rho^2=0.1 \), \( 0.2 \), \( 0.3 \) and \( 0.5 \). A population in which the regressors
account for a fifth of the variance of the response is detected less than half the time.

## An interval for the squared multiple correlation

Part (d) also makes an exact confidence interval possible. For observed \( R^2=r \), let
\( \rho^2_L \) be the value at which \( \Pr(R^2\le r)=0.975 \) and \( \rho^2_U \) the value at which it equals
\( 0.025 \) (with \( \rho^2_L=0 \) if \( \Pr(R^2\le r)<0.975 \) already at \( \rho^2=0 \)). Because the probability
decreases in \( \rho^2 \), \( [\rho^2_L,\rho^2_U] \) covers the true \( \rho^2 \) exactly when \( R^2 \) lies between the
\( 2.5\% \) and \( 97.5\% \) points of its law at that \( \rho^2 \), which has probability \( 0.95 \). This is the usual inversion of a family of tests (@prp-ci-duality).

::: {#exm-cor-states-r2}
[Violent crime in the states]

For the \( 50 \) states in 2009, regress the violent crime rate on the high-school graduation rate,
the poverty rate, the percentage of single-parent households, the percentage of white residents
and the urbanization rate. Then \( R^2=0.4607 \), the adjusted value is \( 0.3994 \), and
\( F=7.52 \) on \( 5 \) and \( 44 \) degrees of freedom, with \( p \)-value \( 3.6\times 10^{-5} \). Read with random
regressors, the \( 95\% \) interval for \( \rho^2 \) is \( (0.175, 0.610) \). The interval is wide: five regressors and fifty cases leave much room for chance. It also treats
the states as a random sample, which they are not; it describes hypothetical replications of a
process generating rows like these, and here, as in @exm-cor-states-correlation, no more can be claimed.
:::

```{.python .run #cell-multiple-correlation-states}
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
regs = ["hs_grad", "poverty", "single", "white", "urban"]
y = data["violent"].to_numpy()
X = np.column_stack([np.ones(len(y))] + [data[c].to_numpy() for c in regs])
fit = sm.OLS(y, X).fit()
n_s, k_s = len(y), len(regs)
R2_obs = fit.rsquared
print(f"R^2 = {R2_obs:.4f}, adjusted {fit.rsquared_adj:.4f}, F = {fit.fvalue:.2f},"
      f" p = {fit.f_pvalue:.2g}")

# 95% interval for rho^2: the values at which the observed R^2 is not extreme
def limit(prob):
    g = lambda rho2: r2_cdf(R2_obs, rho2, n_s, k_s) - prob
    return 0.0 if g(1e-12) < 0 else optimize.brentq(g, 1e-12, 0.999)

lower, upper = limit(0.975), limit(0.025)
print(f"95% interval for rho^2: ({lower:.3f}, {upper:.3f})")
```

::: {.remark}
[Fixed or random?]

Under a fixed design the expectation of \( R^2 \) depends on the chosen regressor values (@exr-ss-design-r2).
With random rows \( R^2 \) estimates \( \rho^2_{Y\cdot X} \), which still depends on the regressor distribution:
a narrower range gives a smaller \( \rho^2 \). The null law of the \( F \) test is the same in both readings.
:::

## Exercises

### A. Check your understanding

::: {#exr-cor-r2-to-f}
[A1]

A regression on \( k=3 \) random regressors from \( n=40 \) cases has \( R^2=0.25 \). Compute the overall \( F \)
statistic and its degrees of freedom. Under the null hypothesis \( \rho^2=0 \), what are the mean and the
distribution of \( R^2 \)?
:::

### B. Practice

::: {#exr-cor-r2-inverse-correlation}
[B1]

Let \( \R \) be the sample correlation matrix of \( (Y,X_1,\dots,X_k) \), assumed nonsingular, and let \( r^{YY} \) be
the \( (1,1) \) entry of \( \R^{-1} \). Show that \( R^2=1-1/r^{YY} \), and that
\( 1-R^2=\det\R/\det\R_{XX} \), where \( \R_{XX} \) is the correlation matrix of the regressors.
:::

::: {.solution}
@eq-cor-r2-plugin is unchanged when every variable is standardized, so
\( R^2=\mathbf{r}_{XY}\T\R_{XX}^{-1}\mathbf{r}_{XY} \). By @thm-mat-partitioned-inverse, the leading entry of \( \R^{-1} \) is the inverse of
the Schur complement \( 1-\mathbf{r}_{XY}\T\R_{XX}^{-1}\mathbf{r}_{XY}=1-R^2 \), so \( r^{YY}=1/(1-R^2) \). By
@thm-mat-block-determinant(b), \( \det\R=\det\R_{XX}\,(1-R^2) \).
:::

::: {#exr-cor-r2-invariance}
[B2]

Show that \( R^2 \) is unchanged if \( Y \) is replaced by \( aY+b \) with \( a\ne0 \), or the regressors by
\( \B\T\mathbf{X}+\mathbf{c} \) with \( \B \) nonsingular. Show that \( \rho^2_{Y\cdot X} \) has the same invariance.
:::

::: {#exr-cor-r2-consistent}
[B3]

Show that \( R^2\to\rho^2_{Y\cdot X} \) with probability one as \( n\to\infty \), without assuming normality, provided the
rows are independent and identically distributed with finite second moments and \( \bSigma_{XX} \) positive
definite.
:::

### C. Going deeper

::: {#exr-cor-r2-lrt}
[C1]

With normal rows, show that the likelihood ratio statistic for \( \rho^2_{Y\cdot X}=0 \) against an unrestricted
alternative is \( (1-R^2)^{n/2} \), so the likelihood ratio test is the overall \( F \) test. *Hint:* under the
null the covariance matrix is block diagonal; use the proof of @prp-cor-mvn-mle and @exr-cor-r2-inverse-correlation.
:::

::: {.solution}
With \( \bS \) the covariance matrix of the rows with divisor \( n \), the maximized log-likelihood of a
\( (k+1) \)-variate normal sample is \( -\tfrac n2[\log\det\bS+(k+1)(1+\log2\pi)] \), by the argument in the proof of
@prp-cor-mvn-mle applied to all \( k+1 \) variables. Under the null, \( Y \) and \( \mathbf{X} \) are independent, the
likelihood factorizes, and the maximum is the same expression with \( \log(s_Y^2\det\bS_{XX}) \) in place of
\( \log\det\bS \). The likelihood ratio is therefore
\( \bigl(\det\bS/(s_Y^2\det\bS_{XX})\bigr)^{n/2}=(1-R^2)^{n/2} \), by @thm-mat-block-determinant(b) as in @exr-cor-r2-inverse-correlation. It is a decreasing function of \( R^2 \), hence of \( F \).
:::

