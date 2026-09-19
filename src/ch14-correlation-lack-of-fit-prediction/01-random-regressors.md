# Random regressors and conditional inference

The theory of Part III was built for a model matrix \( \X \) of known constants. In observational
data the regressors are drawn with the response, and a second sample would bring a different \( \X \). [Section 5.1](../ch05-model-and-least-squares/01-what-a-model-claims.html)
already argued that the assumptions of the linear model should then be read *conditionally on*
\( \X \), and @cor-lm-random-regressors carried unbiasedness over to random regressors. This section
does the same for the exact distribution theory. Tests and intervals stay exact; the estimates
acquire heavier-tailed distributions, and power becomes a random quantity that depends on the design
the data deliver.

## The model with random regressors

::: {#def-cor-random-regressors}
[Normal linear model with random regressors]

Let \( \X \) be a random \( n\times p \) matrix and \( \Y \) a random \( n \)-vector. They follow the
**normal linear model with random regressors** if \( \rank(\X)=p \) with probability one and,
conditionally on \( \X \),
\[
\Y\mid\X\ \sim\ \Normal_n(\X\bbeta,\sigma^2\I_n)
\]
for constants \( \bbeta\in\Real^p \) and \( \sigma^2>0 \). Equivalently, \( \Y=\X\bbeta+\be \) with
\( \be\sim\Normal_n(\bzero,\sigma^2\I_n) \) independent of \( \X \).
:::

The leading case has independent and identically distributed rows \( (\x_{(i)}\T,Y_i) \) with
\( Y_i=\x_{(i)}\T\bbeta+\varepsilon_i \), where \( \varepsilon_i\sim\Normal(0,\sigma^2) \) is independent of
\( \x_{(i)} \). Independence of the rows then makes the errors independent given the whole of
\( \X \). Jointly normal rows satisfy the definition (@eq-mvn-random-x-model), but the regressors may
equally be skewed, discrete or dependent on one another; only the errors must be normal,
homoscedastic and independent of the regressors.

## Conditioning on the design

::: {#lem-cor-conditioning}
[Pivots given the design]

Let \( S \) be a random vector and \( U \) any random element. Suppose there is a distribution \( Q \) such
that, for almost every value \( u \), the conditional distribution of \( S \) given \( U=u \) is \( Q \). Then
\( S\sim Q \), and \( S \) is independent of \( U \).
:::

::: {.proof}
For a Borel set \( B \) and an event \( \{U\in A\} \), the tower rule gives
\[
\begin{aligned}
\Pr(S\in B,\ U\in A)&=\E\bigl[\mathbf 1\{U\in A\}\Pr(S\in B\mid U)\bigr]\\
&=\E\bigl[\mathbf 1\{U\in A\}\,Q(B)\bigr]=Q(B)\Pr(U\in A).
\end{aligned}
\]
Taking \( A \) to be the whole space gives \( \Pr(S\in B)=Q(B) \), and then the display says that the joint
probability factorizes, which is independence.
:::

A statistic whose conditional law given \( \X \) is the same for every design is a *pivot given
the design*. The \( t \) and \( F \) statistics of Part III are of this kind.

::: {#thm-cor-conditional}
[Conditional inference is valid unconditionally]

Assume the model of @def-cor-random-regressors, with \( n>p \).
Let \( \hbeta=(\X\T\X)^{-1}\X\T\Y \) and \( s^2=\norm{\Y-\X\hbeta}^2/(n-p) \).

::: {.enumerate options="label=(\alph*)"}
1. For every fixed \( \blambda\ne\bzero \), the statistic
   \( T=(\blambda\T\hbeta-\blambda\T\bbeta)/\bigl(s\sqrt{\blambda\T(\X\T\X)^{-1}\blambda}\bigr) \) has the
   \( t(n-p) \) distribution and is independent of \( \X \).

2. Let \( \bLambda \) be a fixed \( p\times q \) matrix of rank \( q \), let \( \mathbf{d}\in\Real^q \), and put
   \[
F=\frac{(\bLambda\T\hbeta-\mathbf{d})\T\bigl[\bLambda\T(\X\T\X)^{-1}\bLambda\bigr]^{-1}(\bLambda\T\hbeta-\mathbf{d})}{q\,s^2}.
\]
   Given \( \X \), \( F\sim F\bigl(q,n-p,\gamma(\X)\bigr) \) with
   \[
\gamma(\X)=\frac{(\bLambda\T\bbeta-\mathbf{d})\T\bigl[\bLambda\T(\X\T\X)^{-1}\bLambda\bigr]^{-1}(\bLambda\T\bbeta-\mathbf{d})}{\sigma^2}.
\]{#eq-cor-random-noncentrality}

   If \( \bLambda\T\bbeta=\mathbf{d} \), then \( F\sim F(q,n-p) \) unconditionally, independently of \( \X \).

3. A confidence set whose coverage probability given \( \X \) equals \( 1-\alpha \) for almost every \( \X \)
   has unconditional coverage \( 1-\alpha \). A test whose rejection probability under the null
   hypothesis equals \( \alpha \) given almost every \( \X \) has unconditional size \( \alpha \).

4. The unconditional power of the level-\( \alpha \) test that rejects when \( F>F_\alpha(q,n-p) \) is
   \( \E\bigl[\pi(\gamma(\X))\bigr] \), where \( \pi(\gamma)=\Pr\{F(q,n-p,\gamma)>F_\alpha(q,n-p)\} \).
:::

:::

::: {.proof}
Fix a value of \( \X \) of rank \( p \). Given this value, \( \Y \) follows the normal linear model with a
fixed model matrix, so the results of [Chapter 7](../ch07-optimality/index.html) apply to the
conditional law given \( \X \).

(a) By @cor-opt-t with \( d=\blambda\T\bbeta \), the conditional law of \( T \) is \( t(n-p) \). It is the same
for every design, so @lem-cor-conditioning gives the unconditional law and the independence.

(b) Write \( \mathbf{W}=\bLambda\T(\X\T\X)^{-1}\bLambda \), which is positive definite because \( \bLambda \) has rank
\( q \). By @thm-opt-sampling(a), given \( \X \), \( \mathbf{D}=\bLambda\T\hbeta-\mathbf{d}\sim\Normal_q(\bLambda\T\bbeta-\mathbf{d},\sigma^2\mathbf{W}) \),
so by @thm-qf-chisq(b) with \( \A=\mathbf{W}^{-1}/\sigma^2 \) (for which \( \A\cdot\sigma^2\mathbf{W}=\I \) is idempotent) the
quadratic form \( \mathbf{D}\T\mathbf{W}^{-1}\mathbf{D}/\sigma^2 \) is \( \chi^2(q,\gamma(\X)) \). It is independent of
\( (n-p)s^2/\sigma^2\sim\chi^2(n-p) \) by @thm-opt-sampling(b) and (c). Hence the conditional law of \( F \)
is \( F(q,n-p,\gamma(\X)) \) (@def-qf-noncentral-f). Under \( \bLambda\T\bbeta=\mathbf{d} \) the noncentrality is
zero for every design, and @lem-cor-conditioning applies.

(c) Let \( A \) be the event that the set covers the true value, or that the test rejects. Then
\( \Pr(A)=\E\bigl[\Pr(A\mid\X)\bigr]=\E(1-\alpha)=1-\alpha \), and likewise for the size.

(d) Given \( \X \), the rejection probability is \( \pi(\gamma(\X)) \) by (b). Take expectations.
:::

The statistic in (b) is the \( F \) statistic of the general linear hypothesis (@thm-glh-general-f) for a full-rank design. So every test of
[Chapter 11](../ch11-general-linear-hypothesis/index.html) and every interval of
[Chapter 12](../ch12-intervals-and-bands/index.html) that is exact for fixed \( \X \) is exact for
random \( \X \) too, whatever the distribution of the regressors, even without finite moments.

::: {#exm-cor-skewed-design}
[Skewed regressors]

Take \( n=10 \) cases with an intercept and two regressors drawn afresh for every sample, one
exponential and one lognormal, and \( \E(Y)=1+0\cdot x_1+0.5\,x_2 \). The \( t \) test of
\( \beta_1=0 \) at level \( 0.05 \) uses the \( t(7) \) quantile. Over \( 100000 \) samples with standard
normal errors independent of the regressors, the test rejects in a proportion \( 0.0508 \) of
them, as @thm-cor-conditional requires; panel (b) of [Figure 14.1.1](#fig-cor-random-regressors)
shows the simulated statistics against the \( t(7) \) density. If instead the error of each case is
multiplied by its own value of \( x_1 \), so that the spread depends on the regressor, the same test
rejects in a proportion \( 0.2798 \). The conditional model no longer holds, and nothing in the
theorem protects the test.
:::

```{.python .run #cell-random-regressors-skewed}
import numpy as np
from scipy import stats
reps = 100_000

rng = np.random.default_rng(1402)
n2, p2 = 10, 3
beta = np.array([1.0, 0.0, 0.5])                     # H0: beta_1 = 0 is true
tq = stats.t.ppf(0.975, n2 - p2)

def t_stats(errors_from):
    """t statistics for beta_1 = 0 over many samples with skewed random regressors."""
    out = np.empty(reps)
    for s in range(reps):
        X = np.column_stack([np.ones(n2), rng.exponential(size=n2), rng.lognormal(size=n2)])
        y = X @ beta + errors_from(X)
        XtX_inv = np.linalg.inv(X.T @ X)
        b = XtX_inv @ X.T @ y
        s2 = np.sum((y - X @ b) ** 2) / (n2 - p2)
        out[s] = b[1] / np.sqrt(s2 * XtX_inv[1, 1])
    return out

t_indep = t_stats(lambda X: rng.normal(size=n2))           # errors independent of X
t_hetero = t_stats(lambda X: X[:, 1] * rng.normal(size=n2))  # spread grows with x_1
for label, t in [("independent errors", t_indep), ("spread depends on X", t_hetero)]:
    print(f"{label:20s} rejection rate {np.mean(np.abs(t) > tq):.4f}")
```

::: {.warning}
The conditional model can fail in ways no amount of data repairs. If the error variance depends
on the regressors, the \( t \) and \( F \) statistics lose their exact laws ([Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html)). If the regressors are
correlated with the errors, as with measurement error (@thm-eiv-attenuation), an omitted common cause
([Chapter 25](../ch25-causal-interpretation/index.html)) or a lagged response, least squares estimates something other than \( \bbeta \). And if the mean
is not linear, the population projection of [Section 6.11](../ch06-projections/11-population.html) is
still estimated consistently, but its errors depend on the regressors and the usual standard errors are wrong.
:::

## What changes: the law of the estimates

Pivots keep their laws, but estimates do not. Given the design, \( \hbeta\sim\Normal_p(\bbeta,\sigma^2(\X\T\X)^{-1}) \).
Unconditionally \( \hbeta \) is a mixture of normal distributions over the random covariance
\( \sigma^2(\X\T\X)^{-1} \). Its mean is still \( \bbeta \), and its covariance is
\( \sigma^2\E\bigl[(\X\T\X)^{-1}\bigr] \) by @cor-lm-random-regressors, but it is not normal. For a single
normal regressor the mixture can be written down exactly. We need one moment of the chi-squared
distribution.

::: {#lem-cor-inverse-chisq}
[Mean of an inverse chi-squared variable]

If \( V\sim\chi^2(\nu) \) with \( \nu>2 \), then \( \E(1/V)=1/(\nu-2) \), the fact used without proof for @eq-ss-f-mean. For \( \nu\le2 \) the mean is infinite.
:::

::: {.proof}
With the density @eq-qf-chisq-density,
\[
\E\Bigl(\frac1V\Bigr)=\int_0^\infty\frac{v^{\nu/2-2}e^{-v/2}}{2^{\nu/2}\Gamma(\nu/2)}\,dv
=\frac{2^{\nu/2-1}\Gamma(\nu/2-1)}{2^{\nu/2}\Gamma(\nu/2)}=\frac{1}{2(\nu/2-1)}=\frac1{\nu-2},
\]
using \( \Gamma(a+1)=a\Gamma(a) \). The integral converges at zero iff \( \nu/2-2>-1 \), that is, \( \nu>2 \).
:::

::: {#prp-cor-slope-unconditional}
[The slope under a random normal regressor]

Let \( (X_i,Y_i) \), \( i=1,\dots,n \) with \( n\ge3 \), be independent draws from a bivariate normal
distribution with \( \Var(X_i)=\sigma_X^2 \), regression slope \( \beta_1 \) and conditional variance
\( \sigma^2=\Var(Y_i\mid X_i) \). Let \( \hat{\beta}_1 \) be the least squares slope. Then
\[
\frac{\sigma_X\sqrt{n-1}}{\sigma}\,(\hat{\beta}_1-\beta_1)\ \sim\ t(n-1),
\]
so \( \hat{\beta}_1 \) has heavier tails than any normal distribution, and for \( n>3 \)
\[
\Var(\hat{\beta}_1)=\frac{\sigma^2}{\sigma_X^2\,(n-3)} .
\]{#eq-cor-slope-variance}

:::

::: {.proof}
Given the \( X \) values, \( \hat{\beta}_1\sim\Normal(\beta_1,\sigma^2/S_{xx}) \) with
\( S_{xx}=\sum_i(X_i-\bar X)^2 \) (@thm-opt-sampling(a) and @cor-lm-simple-moments). So
\( Z=\sqrt{S_{xx}}(\hat{\beta}_1-\beta_1)/\sigma \) is standard normal given every design, and by
@lem-cor-conditioning it is standard normal and independent of the \( X \) values, hence of
\( V=S_{xx}/\sigma_X^2 \). The \( X_i \) are independent \( \Normal(\mu_X,\sigma_X^2) \), so
\( V\sim\chi^2(n-1) \) by @cor-qf-sample-variance(b). Therefore
\[
\frac{\sigma_X\sqrt{n-1}}{\sigma}(\hat{\beta}_1-\beta_1)=\frac{Z}{\sqrt{V/(n-1)}}\sim t(n-1)
\]
by @def-qf-noncentral-t. For the variance, \( \E(Z^2/V)=\E(Z^2)\E(1/V)=1/(n-3) \) by independence and
@lem-cor-inverse-chisq, so \( \Var(\hat{\beta}_1)=\sigma^2\sigma_X^{-2}/(n-3) \).
:::

The naive guess, the fixed-design variance \( \sigma^2/S_{xx} \) with \( S_{xx} \) replaced by its mean
\( (n-1)\sigma_X^2 \), is too small by the factor \( (n-3)/(n-1) \). This is Jensen's inequality at work, as in
@exm-lm-fixed-random: designs that happen to bunch the regressor values together cost more
precision than spread-out designs gain. [Section 14.2](02-best-linear-prediction.html) extends
@eq-cor-slope-variance to several regressors.

::: {#exm-cor-slope-simulation}
[A heavy-tailed slope]

With \( n=12 \), \( \sigma_X=2 \), \( \sigma=1 \) and slope \( 0.5 \), the variance of the slope over
\( 100000 \) simulated samples, each with its own design, is \( 0.02751 \), against
\( 0.02778 \) from @eq-cor-slope-variance and \( 0.02273 \) from the naive
guess. The excess kurtosis of the simulated slopes is \( 0.85 \), close to the value
\( 0.86 \) of the \( t(11) \) distribution. Panel (a) of
[Figure 14.1.1](#fig-cor-random-regressors) shows the tails on a logarithmic scale. The scaled
\( t(11) \) density follows the simulated slopes, and the normal density with the same variance
falls away too fast.
:::

```{.python .run #cell-random-regressors-slope}
rng = np.random.default_rng(1401)
n, beta1, sigma, sigma_x = 12, 0.5, 1.0, 2.0
x = sigma_x * rng.normal(size=(reps, n))            # a new design for every sample
y = 1.0 + beta1 * x + sigma * rng.normal(size=(reps, n))
xc = x - x.mean(axis=1, keepdims=True)
Sxx = np.sum(xc ** 2, axis=1)
b1 = np.sum(xc * y, axis=1) / Sxx                    # least squares slope, one per sample

scale = sigma / (sigma_x * np.sqrt(n - 1))           # b1 - beta1 = scale * t(n-1)
print(f"variance of the slope   {b1.var():.5f}")
print(f"sigma^2/(sigma_x^2 (n-3)) {sigma**2 / (sigma_x**2 * (n - 3)):.5f}")
print(f"sigma^2/E(Sxx)           {sigma**2 / (sigma_x**2 * (n - 1)):.5f}")
print("KS p-value against the scaled t(n-1):",
      round(stats.kstest((b1 - beta1) / scale, "t", args=(n - 1,)).pvalue, 3))
```

::: {when-format="html"}
![**Figure 14.1.1.** Random regressors. (a) The standardized slope with a random normal regressor
(\( n=12 \)) follows the scaled \( t(11) \) law of @prp-cor-slope-unconditional. (b) With skewed regressors
the \( t \) statistic is exactly \( t(7) \). (c) Conditional power varies with the design; the line marks the
unconditional power.](random_regressors.svg){#fig-cor-random-regressors width=100%}
:::

::: {when-format="pdf"}
![Random regressors. (a) The standardized slope with a random normal regressor
(\( n=12 \)) follows the scaled \( t(11) \) law of @prp-cor-slope-unconditional. (b) With skewed regressors
the \( t \) statistic is exactly \( t(7) \). (c) Conditional power varies with the design; the line marks the
unconditional power.](random_regressors.pdf){width=100%}
:::

Which of the two distributions should a report describe? For the data in hand, the conditional
one. The design is *ancillary*: its distribution does not involve \( (\bbeta,\sigma^2) \), and it tells us how
precise this particular sample is. A slope estimated from widely spread regressor values deserves a
shorter interval than one estimated from bunched values, which the conditional standard error
\( s/\sqrt{S_{xx}} \) delivers, and by @thm-cor-conditional the interval is also exact on average over designs.
The unconditional law answers questions asked before the data arrive, such as the power of a planned study.

## Power depends on the design

Part (d) of @thm-cor-conditional says that power, unlike size, is not free of the design. For a
single normal regressor the noncentrality of the test of a zero slope is
\( \gamma(\X)=\beta_1^2S_{xx}/\sigma^2 \), which is proportional to a \( \chi^2(n-1) \) variable.

::: {#exm-cor-conditional-power}
[How much power has this design?]

Take \( n=12 \), \( \sigma=\sigma_X=1 \) and slope \( 0.6 \). The test of a zero slope at level \( 0.05 \) has
conditional power between \( 0.321 \) and \( 0.518 \) for the middle half of the designs, with
median \( 0.414 \) (panel (c) of [Figure 14.1.1](#fig-cor-random-regressors)). The
unconditional power, the average over designs, is \( 0.425 \); a direct simulation of the
whole procedure gives \( 0.426 \). Plugging in the average design, \( S_{xx}=n-1 \), gives
\( 0.436 \), which overstates the power a study should expect.
:::

```{.python .run #cell-random-regressors-power}
rng = np.random.default_rng(1403)
n3, b3, df3 = 12, 0.6, 10                            # slope 0.6, sigma = sigma_x = 1
Sxx3 = stats.chi2.rvs(n3 - 1, size=reps, random_state=rng)   # Sxx / sigma_x^2 ~ chi^2(n-1)
gamma = b3 ** 2 * Sxx3                               # noncentrality given the design
fq = stats.f.ppf(0.95, 1, df3)
cond_power = stats.ncf.sf(fq, 1, df3, gamma)         # power of the F = t^2 test given X
print(f"conditional power: quartiles {np.percentile(cond_power, [25, 50, 75]).round(3)}")
print(f"unconditional power {cond_power.mean():.3f}")
print(f"power at the average design {stats.ncf.sf(fq, 1, df3, b3**2 * (n3 - 1)):.3f}")
```

## Exercises

### A. Check your understanding

::: {#exr-cor-slope-numbers}
[A1]

In the setting of @prp-cor-slope-unconditional take \( n=20 \), \( \sigma=2 \) and \( \sigma_X=1 \). Compute the
unconditional variance of \( \hat{\beta}_1 \) and compare it with \( \sigma^2/\E(S_{xx}) \). Which quantile would you
use for a \( 95\% \) interval for \( \beta_1 \) computed from one sample?
:::

::: {.solution}
\( \Var(\hat{\beta}_1)=4/17=0.235 \), against \( 4/19=0.211 \) for the naive value. For an interval from
one sample, use the conditional standard error \( s/\sqrt{S_{xx}} \) and the \( t(18) \) quantile, as for a
fixed design; by @thm-cor-conditional its coverage is exactly \( 95\% \) over random designs too.
The unconditional \( t(19) \) law of the proposition involves the unknown \( \sigma/\sigma_X \) and is not
used for intervals.
:::

::: {#exr-cor-which-exact}
[A2]

Under @def-cor-random-regressors, which of the following hold unconditionally: (i)
\( \E(s^2)=\sigma^2 \); (ii) \( \Cov(\hbeta)=\sigma^2(\X\T\X)^{-1} \); (iii) \( (n-p)s^2/\sigma^2\sim\chi^2(n-p) \),
independent of \( \X \); (iv) the \( 95\% \) confidence ellipsoid for \( \bbeta \) has coverage \( 0.95 \)?
:::

### B. Practice

::: {#exr-cor-fitted-random}
[B1]

In the setting of @prp-cor-slope-unconditional, let \( x_0 \) be a fixed value and
\( \hat\mu_0=\hat{\beta}_0+\hat{\beta}_1x_0 \) the fitted mean there. Show that
\[
\Var(\hat\mu_0)=\sigma^2\Bigl[\frac1n+\frac{(x_0-\mu_X)^2+\sigma_X^2/n}{\sigma_X^2(n-3)}\Bigr].
\]
*Hint:* \( \bar X \) and \( S_{xx} \) are independent.
:::

::: {.solution}
Given the design, \( \hat\mu_0 \) is unbiased with variance \( \sigma^2[1/n+(x_0-\bar X)^2/S_{xx}] \) (@eq-lm-var-mean-response).
By @prp-rv-total-covariance, the unconditional variance is the expectation of this, since the conditional mean
\( \beta_0+\beta_1x_0 \) is constant. By @cor-qf-sample-variance, \( \bar X\sim\Normal(\mu_X,\sigma_X^2/n) \) is
independent of \( S_{xx}=\sigma_X^2V \) with \( V\sim\chi^2(n-1) \). Hence
\( \E[(x_0-\bar X)^2/S_{xx}]=\E(x_0-\bar X)^2\,\E(1/V)/\sigma_X^2=[(x_0-\mu_X)^2+\sigma_X^2/n]/[\sigma_X^2(n-3)] \) by @lem-cor-inverse-chisq.
:::

::: {#exr-cor-adaptive-design}
[B2]

Suppose a pilot sample is used to choose the regressor values of the main study, so that the
distribution of \( \X \) in the main study depends on \( \bbeta \). Assume that, given \( \X \), the main-study
responses follow the normal linear model with errors independent of the pilot. Show that the \( t \) test
of @thm-cor-conditional(a), computed from the main study alone, still has exact size. Why is the
design no longer ancillary, and what information about \( \bbeta \) does the \( t \) test ignore?
:::

### C. Going deeper

::: {#exr-cor-average-variance-interval}
[C1]

In the setting of @prp-cor-slope-unconditional, with \( \sigma \) known, consider the interval
\( \hat{\beta}_1\pm1.96\,\sigma/\sqrt{(n-1)\sigma_X^2} \), which uses the average design instead of the observed
one. Show that its coverage given the design is \( 2\Phi(1.96\sqrt{V/(n-1)})-1 \) with \( V\sim\chi^2(n-1) \), and
that its unconditional coverage is below \( 0.95 \). Which designs make it badly wrong, and in which
direction?
:::
