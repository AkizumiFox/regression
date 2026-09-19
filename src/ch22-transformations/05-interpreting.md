# Interpreting a transformed model

A model fitted on a transformed scale answers questions on that scale, but readers ask about the original
one: how much will a household with income \( 2000 \) spend on food, within what range, and how much more if its
income rises? The dividing line is simple. Quantiles pass through a monotone transformation unchanged;
means do not.

Throughout this section \( g \) is strictly increasing with inverse \( h \), and the model on the transformed scale is
\[
g(Y_i)=\x_{(i)}\T\bbeta+\varepsilon_i,\qquad i=1,\dots,n,
\]{#eq-tr-transformed-model}

where \( \varepsilon_1,\dots,\varepsilon_n \) are independent with distribution \( F \), mean \( 0 \) and variance \( \sigma^2 \).

## Medians and means

::: {#prp-tr-retransformation}
[Retransformation bias and the smearing estimate]

Assume @eq-tr-transformed-model.

::: {.enumerate options="label=(\alph*)"}
1. **Quantiles.** If \( \xi_q \) is a \( q \)-quantile of \( F \), then \( h(\x\T\bbeta+\xi_q) \) is a \( q \)-quantile of \( Y \) at \( \x \).
   In particular, if \( F \) has median \( 0 \) (for instance, if it is symmetric), the median of \( Y \) at \( \x \) is \( h(\x\T\bbeta) \).

2. **Means.** \( \E(Y\mid\x)=\int h(\x\T\bbeta+e)\,dF(e) \) whenever this is finite. If \( h \) is strictly convex and
   \( \sigma^2>0 \), then \( \E(Y\mid\x)>h(\x\T\bbeta) \). For \( h=\exp \), \( \E(Y\mid\x)=\exp(\x\T\bbeta)\,\E e^{\varepsilon} \), and for normal
   errors \( \E e^{\varepsilon}=\exp(\sigma^2/2) \).

3. **Smearing.** Let \( h=\exp \), let \( \X \) have full column rank \( p \) for every \( n \), with leverages \( h_{ii} \), and let the
   \( \varepsilon_i \) be the first \( n \) terms of an independent and identically distributed sequence with
   \( \E e^{\varepsilon}<\infty \). If \( \max_ih_{ii}\to0 \) as \( n\to\infty \), then the **smearing factor** satisfies
   \[
\frac1n\sum_{i=1}^n\exp(\hat{\varepsilon}_i)\to\E e^{\varepsilon}\quad\text{in probability}.
\]{#eq-tr-smearing}

   If moreover \( \x_0\T\hbeta\to\x_0\T\bbeta \) in probability, the **smearing estimate**
   \( \exp(\x_0\T\hbeta)\,n^{-1}\sum_i\exp(\hat{\varepsilon}_i) \) converges in probability to \( \E(Y\mid\x_0) \).
:::

:::

::: {.proof}
(a) \( Y\le h(\x\T\bbeta+\xi_q) \) iff \( g(Y)\le\x\T\bbeta+\xi_q \) iff \( \varepsilon\le\xi_q \), because \( g \) is strictly increasing, and
likewise with strict inequalities. So the two events have the probabilities that define a \( q \)-quantile.
(b) \( Y=h(\x\T\bbeta+\varepsilon) \). If \( h \) is strictly convex and \( \varepsilon \) is not degenerate, Jensen's inequality is strict:
\( \E h(\x\T\bbeta+\varepsilon)>h(\x\T\bbeta+\E\varepsilon)=h(\x\T\bbeta) \). For \( h=\exp \) the factor \( \exp(\x\T\bbeta) \) comes out of the
expectation, and the normal moment generating function (@thm-mvn-mgf) at \( 1 \) is \( \exp(\sigma^2/2) \).
(c) Put \( \mathbf{d}=\X(\hbeta-\bbeta)=\M\be \), so that \( \hat{\varepsilon}_i=\varepsilon_i-d_i \). With \( \vect{e}_i \) the \( i \)th standard basis
vector, \( d_i=\vect{e}_i\T\M\be=(\M\vect{e}_i)\T(\M\be) \), so by the Cauchy–Schwarz inequality and @prp-proj-leverage(a),
\[
\lvert d_i\rvert\le\norm{\M\vect{e}_i}\,\norm{\M\be}=\sqrt{h_{ii}}\,\norm{\M\be}.
\]
Since \( \E\norm{\M\be}^2=\sigma^2\tr(\M)=p\sigma^2 \), Markov's inequality gives
\( \Pr(\norm{\M\be}>K)\le p\sigma^2/K^2 \) for every \( K \), and therefore
\( \delta_n=\max_i\lvert d_i\rvert\le(\max_ih_{ii})^{1/2}\norm{\M\be}\to0 \) in probability. Let
\( A_n=n^{-1}\sum_i\exp(\varepsilon_i) \), which tends to \( \E e^{\varepsilon} \) in probability by the weak law of large
numbers. Because \( \hat{\varepsilon}_i=\varepsilon_i-d_i \), we have
\[
e^{-\delta_n}A_n\le\frac1n\sum_{i=1}^n\exp(\hat{\varepsilon}_i)\le e^{\delta_n}A_n ,
\]
and both bounds tend to \( \E e^{\varepsilon} \) in probability. This proves @eq-tr-smearing. The last statement follows
because \( \exp(\x_0\T\hbeta)\to\exp(\x_0\T\bbeta) \) in probability, and products of convergent sequences converge.
:::

By part (a), a transformed model is naturally read through medians: \( h(\x\T\hbeta) \) estimates the conditional
median when the errors are symmetric. Part (b) is the **retransformation bias**: for the exponential,
\( h(\x\T\hbeta) \) underestimates the mean by the factor \( \E e^{\varepsilon} \). Under normal errors this suggests
\( \exp(\x_0\T\hbeta+s^2/2) \), but the factor depends on the whole error distribution, and with non-normal errors
that estimate is in general inconsistent: it converges to \( \exp(\sigma^2/2) \), not to \( \E e^{\varepsilon} \). Part (c) is Duan's (1983) remedy, the empirical distribution of the residuals
standing in for \( F \), under the familiar condition that no observation keep a fixed share of the fit
([Section 6.8](../ch06-projections/08-leverage.html)). For general \( h \) the smearing estimate is
\( n^{-1}\sum_ih(\x_0\T\hbeta+\hat{\varepsilon}_i) \) (@exr-tr-smearing-general).

::: {#exm-tr-skewed-errors}
[Skewed errors]

The script simulates the log-linear model \( \log Y=1+0.5x+\varepsilon \) with \( n=100 \) equally spaced values of \( x \)
in \( [0,2] \) and right-skewed errors \( \varepsilon=0.6(E-1) \), \( E \) exponential with mean \( 1 \). The errors have mean \( 0 \)
and variance \( 0.36 \). The true factor is \( \E e^{\varepsilon}=e^{-0.6}/0.4 \), about \( 1.372 \), while the normal-theory
factor \( \exp(0.36/2) \) is about \( 1.197 \). At \( x_0=1 \), averaged over \( 4000 \) samples, the ratio of the estimate to the
true mean is \( 0.730 \) for the naive \( \exp(\x_0\T\hbeta) \), \( 0.877 \) for the normal-theory correction and \( 0.986 \) for
the smearing estimate. Its remaining bias of about one percent reflects that residuals vary slightly less
than errors.
:::

```{.python .run #cell-retransformation-skewed}
import numpy as np


def simulate(reps, n_sim, seed):
    """Log-linear model with right-skewed errors e = 0.6 (E - 1), E ~ Exp(1): estimate / true mean."""
    rng = np.random.default_rng(seed)
    xs = np.linspace(0, 2, n_sim)
    Xs = np.column_stack([np.ones(n_sim), xs])
    beta, sig = np.array([1.0, 0.5]), 0.6
    x0 = np.array([1.0, 1.0])
    true_mean = np.exp(x0 @ beta) * np.exp(-sig) / (1 - sig)       # E exp(e) = e^(-sig)/(1 - sig)
    out = np.empty((reps, 3))
    for r in range(reps):
        e = sig * (rng.exponential(size=n_sim) - 1)
        z = Xs @ beta + e
        coef, *_ = np.linalg.lstsq(Xs, z, rcond=None)
        res = z - Xs @ coef
        s2_sim = res @ res / (n_sim - 2)
        naive = np.exp(x0 @ coef)
        out[r] = [naive, naive * np.exp(s2_sim / 2), naive * np.mean(np.exp(res))]
    return out.mean(axis=0) / true_mean


ratios = simulate(reps=300, n_sim=100, seed=2205)
print("mean of estimate / true mean:  naive %.3f   lognormal %.3f   smearing %.3f" % tuple(ratios))
```

The browser cell uses \( 300 \) samples, so its ratios differ slightly from those quoted.

::: {.warning}
Smearing needs identically distributed errors on the transformed scale. If their variance depends on the
regressors, so does \( \E e^{\varepsilon} \), and a single average corrects too much where the variance is small and
too little where it is large (@exr-tr-hetero-smearing); Manning (1998) showed how large the errors can be.
Check the residual spread on the transformed scale
([Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html)) before retransforming a mean.
:::

## Intervals on the original scale

::: {#prp-tr-back-intervals}
[Intervals after a transformation]

Assume @eq-tr-transformed-model, and let \( Y_0 \) be a new response at \( \x_0 \), independent of \( \Y \) and following the
same model.

::: {.enumerate options="label=(\alph*)"}
1. If \( L\le U \) are statistics with values in the range of \( g \) such that \( \Pr\{L\le g(Y_0)\le U\}=1-\alpha \), then
   \( \Pr\{h(L)\le Y_0\le h(U)\}=1-\alpha \). Under normal errors the prediction interval of @thm-ci-prediction-interval
   therefore transforms into an exact prediction interval for \( Y_0 \).

2. If \( [L,U] \) is a confidence interval for \( \x_0\T\bbeta \) with coverage \( 1-\alpha \) and \( F \) has median \( 0 \), then \( [h(L),h(U)] \)
   is a confidence interval for the median of \( Y_0 \) with coverage \( 1-\alpha \).

3. Let \( g=\log \), let the errors be normal, and let \( \X \) have rank \( r<n \) with \( \x_0\in\C(\X\T) \) and
   \( h_0=\x_0\T(\X\T\X)\ginv\x_0 \). The logarithm of the mean, \( \theta=\x_0\T\bbeta+\sigma^2/2 \), is estimated by
   \( \hat{\theta}=\x_0\T\hbeta+s^2/2 \), which is unbiased with
   \[
\Var(\hat{\theta})=\sigma^2h_0+\frac{\sigma^4}{2(n-r)} .
\]{#eq-tr-mean-variance}

:::

:::

::: {.proof}
(a) and (b): since \( g \) is strictly increasing with inverse \( h \), the events \( \{L\le g(Y_0)\le U\} \) and
\( \{h(L)\le Y_0\le h(U)\} \) are equal, and so are \( \{L\le\x_0\T\bbeta\le U\} \) and \( \{h(L)\le h(\x_0\T\bbeta)\le h(U)\} \); by
@prp-tr-retransformation(a), \( h(\x_0\T\bbeta) \) is the median of \( Y_0 \). (c) By @thm-opt-sampling, \( \x_0\T\hbeta \) is normal
with mean \( \x_0\T\bbeta \) and variance \( \sigma^2h_0 \), \( (n-r)s^2/\sigma^2\sim\chi^2(n-r) \), and the two are independent. So
\( \E\hat{\theta}=\theta \), \( \Var(s^2)=2\sigma^4/(n-r) \), and the variances add.
:::

Part (c) gives the approximate interval \( \exp\{\hat{\theta}\pm z_{\alpha/2}(s^2h_0+s^4/\{2(n-r)\})^{1/2}\} \) for the mean,
approximate because the variance is estimated and \( \hat{\theta} \) is only approximately normal, and dependent on the
lognormal assumption. The intervals of parts (a) and (b) are exact under the normal model on the log scale.

::: {#exm-tr-engel-back}
[Engel's data on the original scale]

With the log–log fit of @exm-tr-engel-scales, [Figure 22.5.1](#fig-tr-back) shows the data on the original
scale with the fitted median curve, the smearing estimate of the mean and a \( 95\% \) prediction band. For a
household with income \( 500 \) the estimated median food expenditure is \( 352.2 \), with \( 95\% \) confidence interval
\( [342.2,\ 362.5] \), and a \( 95\% \) prediction interval for a single household (@thm-ci-prediction-interval and @prp-tr-back-intervals(a)) is \( [269,\ 462] \). At income \( 2000 \) the
median is \( 1153.7 \), with interval \( [1111.7,\ 1197.3] \), and the prediction interval is \( [879,\ 1514] \). A fraction \( 0.455 \) of the households lie below the fitted median curve.

The log-scale residual variance \( s^2=0.01871 \) is small, and so is the retransformation bias: the smearing factor
is \( 1.0091 \) and \( \exp(s^2/2) \) is \( 1.0094 \). The estimated mean at income \( 2000 \) is \( 1164.2 \) by smearing and
\( 1164.5 \) under normality, about one percent above the median. The approximate \( 95\% \) interval for the mean from @eq-tr-mean-variance is \( [1122.4,\ 1208.3] \) at
income \( 2000 \) and \( [345.4,\ 365.9] \) at income \( 500 \).
:::

::: {when-format="html"}
![**Figure 22.5.1.** Engel's data on the original scale with the log–log fit transformed back: the median
curve \( \exp(\hat{\beta}_0+\hat{\beta}_1\log x) \), the smearing estimate of the mean (dashed), nearly indistinguishable here
because \( s^2 \) is small, and the \( 95\% \) prediction band, which widens and becomes asymmetric as income
grows.](retransformation.svg){#fig-tr-back width=70%}
:::

::: {when-format="pdf"}
![Engel's data on the original scale with the log–log fit transformed back: the median
curve \( \exp(\hat{\beta}_0+\hat{\beta}_1\log x) \), the smearing estimate of the mean (dashed), nearly indistinguishable here
because \( s^2 \) is small, and the \( 95\% \) prediction band, which widens and becomes asymmetric as income
grows.](retransformation.pdf){width=70%}
:::

```{.python .run #cell-retransformation-back}
import numpy as np
import statsmodels.api as sm
from scipy import stats

df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), np.log(income)])
fit = sm.OLS(np.log(food), X).fit()
b, s2, resid = fit.params, fit.scale, fit.resid
XtX_inv = np.linalg.inv(X.T @ X)
t = stats.t.ppf(0.975, n - 2)

smear = np.mean(np.exp(resid))                    # Duan's smearing factor
lognormal = np.exp(s2 / 2)                        # the normal-theory factor
print(f"smearing factor {smear:.4f}, lognormal factor {lognormal:.4f}")

for x0 in (500, 2000):
    v0 = np.array([1, np.log(x0)])
    eta, h0 = v0 @ b, v0 @ XtX_inv @ v0
    median = np.exp(eta)
    ci_med = np.exp(eta + np.array([-1, 1]) * t * np.sqrt(s2 * h0))          # for the median, exact
    pi = np.exp(eta + np.array([-1, 1]) * t * np.sqrt(s2 * (1 + h0)))        # for a new household, exact
    print(f"income {x0}: median {median:.1f}, CI [{ci_med[0]:.1f}, {ci_med[1]:.1f}];"
          f" mean {median * smear:.1f} (smearing), {median * lognormal:.1f} (lognormal);"
          f" 95% prediction interval [{pi[0]:.0f}, {pi[1]:.0f}]")
```

## Elasticities and marginal effects

In a transformed model a coefficient is not a rate of change of \( y \) per unit of \( x \); rates and elasticities
follow by differentiating through both transformations.

::: {#prp-tr-elasticity}
[Elasticities in a transformed model]

Let \( x_1>0 \), and suppose that the median \( m(\x) \) of \( Y \) satisfies
\( m(\x)^{(\lambda)}=\beta_0+\beta_1x_1^{(\alpha)}+c(\x_{-1}) \), where \( c \) does not involve \( x_1 \) and the Box–Cox
transformation @eq-tr-box-cox is applied to both \( m \) and \( x_1 \). Then
\[
\frac{\partial m}{\partial x_1}=\beta_1\,x_1^{\alpha-1}m^{1-\lambda},\qquad
\frac{x_1}{m}\,\frac{\partial m}{\partial x_1}=\beta_1\,x_1^{\alpha}m^{-\lambda}.
\]{#eq-tr-elasticity}

In particular the elasticity is the constant \( \beta_1 \) in the log–log model (\( \lambda=\alpha=0 \)), \( \beta_1x_1 \) in the
log–level model (\( \lambda=0 \), \( \alpha=1 \)), \( \beta_1/m \) in the level–log model (\( \lambda=1 \), \( \alpha=0 \)) and
\( \beta_1x_1/m \) in the untransformed model. Under @eq-tr-multiplicative the mean is a constant multiple of the median,
so for \( \lambda=0 \) the elasticity of the mean is the same.
:::

::: {.proof}
Differentiate the identity with respect to \( x_1 \). By @prp-tr-family(b) the derivative of \( m^{(\lambda)} \) with respect
to \( m \) is \( m^{\lambda-1} \) and that of \( x_1^{(\alpha)} \) is \( x_1^{\alpha-1} \), so \( m^{\lambda-1}\,\partial m/\partial x_1=\beta_1x_1^{\alpha-1} \).
Multiply by \( m^{1-\lambda} \), and then by \( x_1/m \). If the mean is \( k\,m(\x) \) for a constant \( k \), the factor \( k \) cancels
in the elasticity.
:::

The model fixes the *shape* of the elasticity profile before the data have a say. For Engel's data the log–log
slope \( 0.8559 \) (interval \( [0.816,\ 0.896] \)) is one elasticity for all households. The untransformed linear fit
implies an elasticity \( \hat{\beta}_1x/\hat{m} \) rising from \( 0.62 \) at income \( 500 \) to \( 0.87 \) at income \( 2000 \), only
because a line with a positive intercept forces that shape. The Box–Tidwell estimate \( \hat{\alpha}=-0.182 \) of
@exm-tr-engel-box-tidwell, with \( \lambda=0 \), implies an elasticity proportional to \( x^{-0.182} \), which *falls* between
the same incomes by the factor \( 4^{-0.182}=\exp(-0.182\log4) \), about \( 0.78 \). Which shape the data support is a
question about the mean function ([Section 22.3](03-transforming-regressors.html)).

## Which scale to report

A transformation chosen for statistical reasons need not be the reporting scale. The logarithm is often worth
keeping, since its coefficients are ratios and elasticities (@prp-tr-multiplicative). For other powers, report
on the original scale: medians and prediction intervals, which transform exactly, and means, which need a
correction factor and a constant variance. Compare models fitted on different scales on a common scale,
through the Box–Cox likelihood ([Section 22.2](02-box-cox.html)) or the accuracy of predictions of \( y \)
(Chapter 29), not through \( R^2 \) (@exr-tr-r2-scales).

## Exercises

### A. Check your understanding

::: {#exr-tr-smear-vs-normal}
[A1]

For Engel's log–log fit \( s^2=0.01871 \). Check that \( \exp(s^2/2) \) is \( 1.0094 \). The smearing factor is \( 1.0091 \).
Why are the two so close here, and in what kind of data would they differ substantially?
:::

::: {#exr-tr-asymmetric-pi}
[A2]

Explain why the prediction interval of @exm-tr-engel-back, obtained from @thm-ci-prediction-interval, is longer above the median than below it,
and why its coverage is nevertheless exactly \( 95\% \) under the normal model on the log scale.
:::

### B. Practice

::: {#exr-tr-median-bias}
[B1]

Under normal errors on the log scale, show that \( \E\exp(\x_0\T\hbeta)=\exp(\x_0\T\bbeta+\sigma^2h_0/2) \), so the naive median
estimate is biased upwards by a factor that tends to one as \( h_0\to0 \). Suggest an approximately mean-unbiased
alternative, and explain why \( \exp(\x_0\T\hbeta) \) is exactly median-unbiased.
:::

::: {.solution}
\( \x_0\T\hbeta\sim\Normal(\x_0\T\bbeta,\sigma^2h_0) \) (@thm-opt-sampling), and the normal moment generating function (@thm-mvn-mgf) at \( 1 \) gives
the mean. The median of \( \x_0\T\hbeta \) is \( \x_0\T\bbeta \), and \( \exp \) preserves medians, so \( \exp(\x_0\T\hbeta) \) is exactly
median-unbiased for \( \exp(\x_0\T\bbeta) \). For mean-unbiasedness one can use \( \exp(\x_0\T\hbeta-s^2h_0/2) \), which is
approximately unbiased. The factor \( \exp(\sigma^2h_0/2) \) is tiny when \( h_0 \) is small, so the choice rarely matters.
:::

::: {#exr-tr-sqrt-smearing}
[B2]

Let \( \sqrt{Y_i}=\x_{(i)}\T\bbeta+\varepsilon_i \) with \( \bone\in\C(\X) \). Show that \( \E(Y\mid\x_0)=(\x_0\T\bbeta)^2+\sigma^2 \) and that the
smearing estimate \( n^{-1}\sum_i(\x_0\T\hbeta+\hat{\varepsilon}_i)^2 \) equals \( (\x_0\T\hbeta)^2+\text{SSE}/n \). How does this compare with
the plug-in estimate \( (\x_0\T\hbeta)^2+s^2 \)?
:::

::: {#exr-tr-delta-check}
[B3]

For Engel's data at income \( 500 \), \( s^2=0.01871 \), \( h_0=0.01147 \) and \( n-r=233 \). Evaluate the two terms of
@eq-tr-mean-variance with \( \sigma^2 \) replaced by \( s^2 \), and check the standard error \( 0.0147 \) used for the
interval \( [345.4,\ 365.9] \). Which term dominates, and why?
:::

### C. Going deeper

::: {#exr-tr-smearing-general}
[C1]

Extend @prp-tr-retransformation(c) to a general strictly increasing \( h \) that is continuously differentiable with
\( \lvert h'(t)\rvert\le K(1+e^{\kappa\lvert t\rvert}) \) for some constants, assuming \( \E e^{\kappa'\lvert\varepsilon\rvert}<\infty \) for some \( \kappa'>\kappa \).
Show that \( n^{-1}\sum_ih(\x_0\T\hbeta+\hat{\varepsilon}_i) \) converges in probability to \( \E(Y\mid\x_0) \). (Use the mean value theorem
and the bound on \( \max_i\lvert\hat{\varepsilon}_i-\varepsilon_i\rvert \) from the proof.)
:::

::: {#exr-tr-hetero-smearing}
[C2]

Suppose that on the log scale \( \varepsilon_i\sim\Normal(0,\sigma^2(\x_{(i)})) \), independent, with the variance depending on the
regressors. Show that \( \E(Y\mid\x_0)=\exp\{\x_0\T\bbeta+\sigma^2(\x_0)/2\} \), and that, when the leverages tend to zero and the
average \( n^{-1}\sum_i\exp\{\sigma^2(\x_{(i)})/2\} \) converges to a limit \( \bar{c} \), the smearing factor tends to \( \bar{c} \). For which
\( \x_0 \) does the smearing estimate overestimate the mean, and for which does it underestimate it?
:::

::: {.solution}
Given \( \x_0 \), \( \log Y \) is normal with mean \( \x_0\T\bbeta \) and variance \( \sigma^2(\x_0) \), and the moment generating
function (@thm-mvn-mgf) gives the mean. The argument of @prp-tr-retransformation(c) goes through with \( A_n=n^{-1}\sum_ie^{\varepsilon_i} \),
whose mean is \( n^{-1}\sum_i\exp\{\sigma^2(\x_{(i)})/2\} \); if the variances are bounded, its variance tends to zero, so
\( A_n\to\bar{c} \). The smearing estimate therefore targets \( \exp(\x_0\T\bbeta)\,\bar{c} \). It overestimates the mean where
\( \exp\{\sigma^2(\x_0)/2\}<\bar{c} \), that is, where the log-scale variance is below its average in this sense, and
underestimates it where the variance is above.
:::
