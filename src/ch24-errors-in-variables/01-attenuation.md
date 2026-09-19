# Attenuation bias

When a regressor is recorded with error, conditioning on the recorded value is still legitimate, but it
describes how the response varies with the *reading*, not with the underlying *quantity*. This section measures
the difference.

## A running example

One simulated study is used throughout the chapter. For \( n=300 \) patients, let \( x_i \) be long-run systolic blood pressure, the average that many visits would give, with
mean \( 130 \) mmHg and standard deviation \( 12 \). Age has mean \( 55 \) years, standard deviation \( 10 \), and
correlation \( 0.5 \) with \( x_i \). Each patient has two clinic readings,
\( w_{ij}=x_i+u_{ij} \), with independent errors \( u_{ij}\sim\Normal(0,9^2) \). The outcome is
\[
y_i=-4+0.08\,x_i+0\cdot\text{age}_i+\varepsilon_i,\qquad \varepsilon_i\sim\Normal(0,1),
\]
so age has no effect given long-run pressure. In a real study \( x_i \) is unseen; here it lets every method be
checked.

## Classical error

Throughout the chapter, a regression model is stated in terms of *true* regressors \( \x \) (the
ones measured with error) and \( \bz \) (the ones measured exactly):
\[
Y=\alpha+\bbeta_x\T\x+\bbeta_z\T\bz+\varepsilon ,
\]{#eq-eiv-true-model}

where \( \varepsilon \) has mean zero and variance \( \sigma^2 \) and is uncorrelated with \( (\x,\bz) \). We write a
generic observation without an index. The data are \( n \) independent copies of it unless stated otherwise.

::: {#def-eiv-classical}
[Classical measurement error]

Suppose \( \x \) is observed only through \( \bw=\x+\bu \). The error \( \bu \) is **classical** if
\( \E(\bu)=\bzero \), \( \Cov(\bu)=\bSigma_{uu} \), and \( \bu \) is independent of \( (\x,\bz,\varepsilon) \).
:::

The error is *additive and unrelated to the true value*, and *nondifferential*: independent of \( \varepsilon \), so
the reading says nothing about the response beyond what \( \x \) says. (People with a disease may recall past
exposures differently, a differential error.) This section uses only that \( \bu \) has mean zero and is *uncorrelated* with \( \x \),
\( \bz \) and \( \varepsilon \).

::: {#def-eiv-reliability}
[Reliability ratio]

The **reliability ratio** of \( W=X+U \) is
\[
\lambda=\frac{\sigma_x^2}{\sigma_x^2+\sigma_u^2}=\frac{\Var(X)}{\Var(W)} .
\]
If \( \bz \) is measured exactly, let \( \sigma^2_{x\cdot z}=\sigma_x^2-\boldsymbol{\upsigma}_{zx}\T\bSigma_{zz}^{-1}\boldsymbol{\upsigma}_{zx} \)
be the residual variance of \( X \) after its best linear prediction from \( \bz \) (@thm-proj-blp),
where \( \boldsymbol{\upsigma}_{zx}=\Cov(\bz,X) \). The **reliability given
\( \bz \)** is
\[
\lambda_{x\mid z}=\frac{\sigma^2_{x\cdot z}}{\sigma^2_{x\cdot z}+\sigma_u^2} .
\]
:::

In the running example a single reading has \( \lambda=144/(144+81)=0.64 \). Age explains a quarter of the
variance of long-run pressure, leaving \( \sigma^2_{x\cdot z}=108 \), so \( \lambda_{x\mid z}=108/189=0.571 \). Covariates
remove part of the signal in \( W \) and none of the noise, so \( \lambda_{x\mid z}\le\lambda \).

## Least squares on the reading

Substituting \( X=W-U \) in \( Y=\alpha+\beta X+\varepsilon \) gives
\[
Y=\alpha+\beta W+(\varepsilon-\beta U).
\]{#eq-eiv-composite}

Its error is correlated with the regressor, \( \Cov(W,\varepsilon-\beta U)=-\beta\sigma_u^2 \), so least squares does not
estimate \( \beta \). It estimates the population projection of \( Y \) on the recorded regressors (@prp-proj-consistency),
which the next theorem computes.

::: {#thm-eiv-attenuation}
[Attenuation]

Let \( (\x_{(i)},\bz_{(i)},\bu_{(i)},\varepsilon_i) \), \( i=1,2,\dots \), be independent copies of
\( (\x,\bz,\bu,\varepsilon) \) with finite second moments, satisfying @eq-eiv-true-model and @def-eiv-classical.
Write \( \mathbf{t}=(\x\T,\bz\T)\T \) for the true regressors,
\( \mathbf{s}=(\bw\T,\bz\T)\T \) for the observed ones, \( \boldsymbol{\uptheta}=(\bbeta_x\T,\bbeta_z\T)\T \),
\( \bSigma_{tt}=\Cov(\mathbf{t}) \), assumed positive definite, and
\( \tilde{\bSigma}=\diag(\bSigma_{uu},\bzero) \) for the covariance matrix of the error in \( \mathbf{s} \).
Let \( \hat{\alpha}_n \) and \( \hat{\boldsymbol{\uptheta}}_n \) be the least squares intercept and slopes from regressing
\( Y \) on \( \mathbf{s} \) in the first \( n \) observations. Then, with probability one:

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\boldsymbol{\uptheta}}_n\to\boldsymbol{\uptheta}_w=(\bSigma_{tt}+\tilde{\bSigma})^{-1}\bSigma_{tt}\boldsymbol{\uptheta}
   =\boldsymbol{\uptheta}-(\bSigma_{tt}+\tilde{\bSigma})^{-1}\tilde{\bSigma}\boldsymbol{\uptheta} \), and
   \( \hat{\alpha}_n\to\alpha+\bmu_t\T(\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_w) \), where \( \bmu_t=\E(\mathbf{t}) \);

2. with a single regressor and no \( \bz \), the slope converges to \( \lambda\beta \);

3. with a single mismeasured scalar \( X \) and exactly measured \( \bz \), the coefficient of \( W \) converges
   to \( \lambda_{x\mid z}\beta_x \) and the coefficients of \( \bz \) converge to
   \( \bbeta_z+(1-\lambda_{x\mid z})\beta_x\boldsymbol{\uppi} \), where \( \boldsymbol{\uppi}=\bSigma_{zz}^{-1}\boldsymbol{\upsigma}_{zx} \)
   is the vector of coefficients from the best linear prediction of \( X \) from \( \bz \) (@thm-rv-blp);

4. the residual sum of squares divided by \( n \) (or by \( n \) minus the number of coefficients) converges to
   \( \sigma^2+\boldsymbol{\uptheta}\T\bSigma_{tt}(\bSigma_{tt}+\tilde{\bSigma})^{-1}\tilde{\bSigma}\boldsymbol{\uptheta} \),
   which in case (b) equals \( \sigma^2+\lambda\beta^2\sigma_u^2 \).
:::

:::

::: {.proof}
Because \( \bu \) is uncorrelated with \( \mathbf{t} \), \( \Cov(\mathbf{s})=\bSigma_{tt}+\tilde{\bSigma} \), which is
positive definite as the sum of a positive definite and a nonnegative definite matrix. By
@prp-proj-consistency, the least squares coefficients converge with probability one to those of the best
linear predictor of \( Y \) from \( \mathbf{s} \), which by @thm-proj-blp has slopes
\( \Cov(\mathbf{s})^{-1}\Cov(\mathbf{s},Y) \) and intercept \( \E Y-\boldsymbol{\uptheta}_w\T\E\mathbf{s} \). Now
\[
\Cov(\mathbf{s},Y)=\Cov\bigl(\mathbf{t}+(\bu\T,\bzero\T)\T,\ \mathbf{t}\T\boldsymbol{\uptheta}+\varepsilon\bigr)=\bSigma_{tt}\boldsymbol{\uptheta},
\]
since \( \bu \) is uncorrelated with \( \mathbf{t} \) and \( \varepsilon \), and \( \varepsilon \) with \( \mathbf{t} \). This is the first
form in (a). The second follows from \( \bSigma_{tt}=(\bSigma_{tt}+\tilde{\bSigma})-\tilde{\bSigma} \). For the
intercept, \( \E\mathbf{s}=\bmu_t \) and \( \E Y=\alpha+\boldsymbol{\uptheta}\T\bmu_t \).

*(b)* With scalars, \( \theta_w=\sigma_x^2\beta/(\sigma_x^2+\sigma_u^2)=\lambda\beta \).

*(c)* Here \( \tilde{\bSigma}\boldsymbol{\uptheta}=\sigma_u^2\beta_x\mathbf{e}_1 \), so by the second form in (a) the
bias \( \boldsymbol{\uptheta}-\boldsymbol{\uptheta}_w \) is \( \sigma_u^2\beta_x \) times the first column of
\( \Cov(\mathbf{s})^{-1} \). Order the observed regressors as \( (\bz,W) \) for the moment, so that
\[
\Cov(\mathbf{s})=\begin{pmatrix}\bSigma_{zz}&\boldsymbol{\upsigma}_{zx}\\ \boldsymbol{\upsigma}_{zx}\T&\sigma_x^2+\sigma_u^2\end{pmatrix}.
\]
The Schur complement of \( \bSigma_{zz} \) is \( \sigma_x^2+\sigma_u^2-\boldsymbol{\upsigma}_{zx}\T\bSigma_{zz}^{-1}\boldsymbol{\upsigma}_{zx}
=\sigma^2_{x\cdot z}+\sigma_u^2 \). By @thm-mat-partitioned-inverse the last column of the inverse is
\( (-\boldsymbol{\uppi}\T,1)\T/(\sigma^2_{x\cdot z}+\sigma_u^2) \). So the bias in the coefficient of \( W \) is
\( \sigma_u^2\beta_x/(\sigma^2_{x\cdot z}+\sigma_u^2)=(1-\lambda_{x\mid z})\beta_x \), and the bias in the
coefficients of \( \bz \) is \( -(1-\lambda_{x\mid z})\beta_x\boldsymbol{\uppi} \).

*(d)* The residual sum of squares over \( n \) is
\( n^{-1}\sum_i(Y_i-\hat{\alpha}_n-\hat{\boldsymbol{\uptheta}}_n\T\mathbf{s}_i)^2 \). Expanded, it is a continuous
function of the coefficients and of second sample moments, which converge by (a) and the strong law, so the limit is
\( \E(Y-\alpha_w-\boldsymbol{\uptheta}_w\T\mathbf{s})^2 \), the mean squared error of the best linear predictor. By
@eq-proj-population-pythagoras it equals \( \Var(Y)-\boldsymbol{\uptheta}_w\T\Cov(\mathbf{s})\boldsymbol{\uptheta}_w \), which is
\[
\sigma^2+\boldsymbol{\uptheta}\T\bSigma_{tt}\boldsymbol{\uptheta}-\boldsymbol{\uptheta}\T\bSigma_{tt}(\bSigma_{tt}+\tilde{\bSigma})^{-1}\bSigma_{tt}\boldsymbol{\uptheta}
=\sigma^2+\boldsymbol{\uptheta}\T\bSigma_{tt}(\bSigma_{tt}+\tilde{\bSigma})^{-1}\tilde{\bSigma}\boldsymbol{\uptheta} .
\]
In case (b) the last term is \( \beta^2\sigma_x^2\sigma_u^2/(\sigma_x^2+\sigma_u^2)=\lambda\beta^2\sigma_u^2 \).
:::

Part (b) is **attenuation**: the slope is shrunk towards zero by the reliability ratio, whatever the
sample size. More data give a more precise estimate of the wrong number \( \lambda\beta \).

::: {#exm-eiv-attenuation}
[One reading instead of long-run pressure]

In the running example, regressing \( y \) on the true \( x \) gives slope \( 0.0755 \)
(standard error \( 0.0050 \)). Regressing on the first clinic reading gives
\( 0.0499 \) (standard error \( 0.0046 \)), a ratio of \( 0.661 \),
close to \( \lambda=0.64 \); the limit is \( \lambda\beta=0.0512 \). The usual 95% interval,
\( [0.0409,\,0.0589] \), misses \( 0.08 \), and in \( 4000 \) repeated samples the observed
coverage was \( 0.000 \). The standard error measures variability around \( \lambda\beta \), not the
distance from \( \lambda\beta \) to \( \beta \). The residual variance tends to \( 1+0.64\cdot0.08^2\cdot81=1.3318 \)
by @thm-eiv-attenuation(d).
:::

## Several regressors: bias spreads

Part (c) is more troubling. Correctly measured covariates lower the relevant reliability to
\( \lambda_{x\mid z} \), and they inherit part of the effect of \( X \). The part of \( X \) that the reading misses acts
like an omitted regressor (@thm-dep-omitted), and covariates correlated with \( X \) pick up its effect in
proportion to \( \boldsymbol{\uppi} \). This is **residual confounding**: adjusting for a noisy measurement of a
confounder only partly adjusts for it.

::: {#exm-eiv-leak}
[Age looks important]

Add age to the regression on the first reading. Theory predicts a coefficient of
\( \lambda_{x\mid z}\beta=0.0457 \) for the reading and
\( (1-0.571)\cdot0.08\cdot0.6=0.0206 \) per year for age, which truly has no effect. The
sample gives \( 0.0455 \) for the reading and \( 0.0182 \) for age, with standard
error \( 0.0070 \) and \( p \)-value \( 0.0098 \). In \( 4000 \) repeated samples the \( t \) test
of "no age effect" rejected at the 5% level in a proportion \( 0.812 \) of them.
[Figure 24.1.1](#fig-eiv-attenuation)(b) shows both limits as functions of the reliability.
:::

::: {when-format="html"}
![**Figure 24.1.1.** (a) The outcomes against true pressure (grey) and one reading (blue), with their least
squares lines. (b) Limits of the coefficients of the reading and of age, in units of \( \beta_x \), against the
reliability; dashed, \( \lambda \) itself.](attenuation.svg){#fig-eiv-attenuation width=100%}
:::

::: {when-format="pdf"}
![(a) The outcomes against true pressure (grey) and one reading (blue), with their least squares lines.
(b) Limits of the coefficients of the reading and of age, in units of \( \beta_x \), against the reliability;
dashed, \( \lambda \) itself.](attenuation.pdf){width=100%}
:::

```{.python .run #cell-attenuation-setup}
import numpy as np
import statsmodels.api as sm

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x

oracle = sm.OLS(y, sm.add_constant(x)).fit()                   # uses the true x
naive = sm.OLS(y, sm.add_constant(w[:, 0])).fit()              # uses one reading
lam = 144 / (144 + 81)                                         # reliability of one reading
print(f"slope on x:          {oracle.params[1]:.4f}")
print(f"slope on reading:    {naive.params[1]:.4f}   limit lambda*beta = {lam * 0.08:.4f}")
print("naive 95% interval:", naive.conf_int()[1].round(4))

both = sm.OLS(y, sm.add_constant(np.column_stack([w[:, 0], age]))).fit()
print(f"with age: reading {both.params[1]:.4f}, age {both.params[2]:.4f} (p = {both.pvalues[2]:.4f})")
```

The repeated-sample figures come from this loop (with fewer samples in the cell).

```{.python .run #cell-attenuation-repeat}
def repeat(reps, seed):
    """Coverage of the naive 95% interval for the slope, and the rate at which age is 'significant'."""
    r = np.random.default_rng(seed)
    cover = reject = 0
    for _ in range(reps):
        a = r.normal(55, 10, n)
        xx = 130 + 0.6 * (a - 55) + r.normal(0, np.sqrt(108), n)
        ww = xx + r.normal(0, 9, n)
        yy = -4 + 0.08 * xx + r.normal(0, 1.0, n)
        lo, hi = sm.OLS(yy, sm.add_constant(ww)).fit().conf_int()[1]
        cover += lo <= 0.08 <= hi
        reject += sm.OLS(yy, sm.add_constant(np.column_stack([ww, a]))).fit().pvalues[2] < 0.05
    return cover / reps, reject / reps

print("coverage of naive interval, rate of rejecting 'no age effect':", repeat(200, 2404))
```

With several mismeasured regressors, \( (\bSigma_{tt}+\tilde{\bSigma})^{-1}\bSigma_{tt} \) need not shrink each
coefficient. A coefficient can grow, change sign, or appear from nothing (@exr-eiv-sign-reversal).

## What survives

Two things survive classical error intact.

**Error in the response.** Suppose instead that the response is measured with error, \( Y^*=Y+V \), where \( V \) has mean zero and variance \( \sigma_v^2 \) and is independent of the regressors and of \( \varepsilon \).

::: {#prp-eiv-response}
[Error in the response]

Under these assumptions \( Y^*=\alpha+\bbeta\T\x+(\varepsilon+V) \), a linear model with the same coefficients and error
variance \( \sigma^2+\sigma_v^2 \). If \( \varepsilon \) and \( V \) are normal, so is their sum, and all the exact theory of
Parts II and III applies to \( Y^* \) with \( \sigma^2 \) replaced by \( \sigma^2+\sigma_v^2 \).
:::

The proof is the display. Error in the response costs precision, not bias, unless it depends on the
regressors (heavier people under-reporting weight) or is correlated with a regressor's error (both from the
same self-report).

**Tests of no effect.** Since \( \lambda\beta=0 \) iff \( \beta=0 \), attenuation can hide an effect of the mismeasured
regressor but cannot create one. The protection does not extend to the covariates.

::: {#prp-eiv-tests}
[Which tests survive classical error]

Assume @eq-eiv-true-model and @def-eiv-classical with a single mismeasured scalar \( X \),
\( \varepsilon\sim\Normal(0,\sigma^2) \) independent of \( (X,\bz) \), and the covariance matrix of \( (X,\bz\T)\T \) positive definite. Consider the least squares regression of \( Y \) on
\( (1,W,\bz) \).

::: {.enumerate options="label=(\alph*)"}
1. If \( \beta_x=0 \), the usual \( t \) test of "the coefficient of \( W \) is zero" has exactly its nominal level.

2. If \( \beta_x\ne0 \) and \( \boldsymbol{\uppi}\ne\bzero \), the usual \( F \) test of "the coefficients of \( \bz \) equal
   \( \bbeta_z \)" rejects with probability tending to one as \( n\to\infty \).
:::

:::

::: {.proof}
*(a)* If \( \beta_x=0 \), then \( Y=\alpha+\bbeta_z\T\bz+\varepsilon \), and \( \varepsilon \) is independent of \( (W,\bz) \) because it
is independent of \( (X,\bz,U) \): it is independent of \( (X,\bz) \) by hypothesis, and \( U \) is independent of
\( (X,\bz,\varepsilon) \) by @def-eiv-classical. So, conditionally on the observed regressors, \( Y \) follows the normal linear model
of @def-cor-random-regressors with regressors \( (1,W,\bz) \), in which the coefficient of \( W \) is zero.
@thm-cor-conditional gives the exact \( t \) distribution.

*(b)* By @thm-eiv-attenuation(c), \( \hat{\bbeta}_z-\bbeta_z\to\mathbf{c}=(1-\lambda_{x\mid z})\beta_x\boldsymbol{\uppi}\ne\bzero \). The
\( F \) statistic is \( n(\hat{\bbeta}_z-\bbeta_z)\T\mathbf{K}_n(\hat{\bbeta}_z-\bbeta_z) \), where
\( \mathbf{K}_n=[n(\mathbf{S}\T\mathbf{S})^{-1}]_{zz}^{-1}/(q\,s^2) \), \( \mathbf{S} \) is the observed model matrix and \( q \) the
dimension of \( \bz \). By the strong law and @thm-eiv-attenuation(d), \( \mathbf{K}_n \) converges to a positive definite
limit, so the statistic grows like \( n \), while the critical value tends to a finite limit.
:::

So "is there any association with the underlying quantity?" can be tested with the reading, with reduced
power. "Does this covariate matter once the underlying quantity is accounted for?" cannot.

## Exercises

### A. Check your understanding

::: {#exr-eiv-reliability-correlation}
[A1]

Show that \( \lambda=\operatorname{Corr}(W,X)^2 \). If \( W_1=X+U_1 \) and \( W_2=X+U_2 \) are two readings with
independent classical errors of the same variance, show that \( \operatorname{Corr}(W_1,W_2)=\lambda \). This
is why test–retest correlations are used as reliability estimates.
:::

::: {.solution}
\( \Cov(W,X)=\sigma_x^2 \), so \( \operatorname{Corr}(W,X)^2=\sigma_x^4/\{\sigma_x^2(\sigma_x^2+\sigma_u^2)\}=\lambda \). For two readings,
\( \Cov(W_1,W_2)=\Var(X)=\sigma_x^2 \) because the errors are independent of each other and of \( X \), and
\( \Var(W_1)=\Var(W_2)=\sigma_x^2+\sigma_u^2 \). The correlation is \( \sigma_x^2/(\sigma_x^2+\sigma_u^2)=\lambda \).
:::

### B. Practice

::: {#exr-eiv-spearman-brown}
[B1]

Let \( \bar{W} \) be the mean of \( m \) readings with independent classical errors of variance \( \sigma_u^2 \), and let \( \lambda \) be the
reliability of one reading. Show that the reliability of \( \bar{W} \) is
\[
\lambda_m=\frac{m\lambda}{1+(m-1)\lambda}.
\]
How many readings are needed in the running example for the attenuation factor to reach \( 0.9 \)?
:::

::: {.solution}
\( \lambda_m=\sigma_x^2/(\sigma_x^2+\sigma_u^2/m) \). Dividing through by \( \sigma_x^2+\sigma_u^2 \) and using
\( \sigma_u^2/(\sigma_x^2+\sigma_u^2)=1-\lambda \) gives \( \lambda/\{\lambda+(1-\lambda)/m\} \), the formula. With \( \lambda=0.64 \),
\( \lambda_m\ge0.9 \) iff \( 0.64m\ge0.9+0.576(m-1) \), so \( m\ge5.0625 \): six readings (without age in the model).
:::

::: {#exr-eiv-r2}
[B2]

In the setting of @thm-eiv-attenuation(b), show that the population \( R^2 \) of \( Y \) on \( W \) equals \( \lambda \) times the
population \( R^2 \) of \( Y \) on \( X \). Measurement error lowers the apparent strength of a relation twice: once
through the slope and once through the residual variance.
:::

::: {#exr-eiv-sign-reversal}
[B3]

Let \( \mathbf{t}=(X_1,X_2) \) have unit variances and correlation \( \rho \), let \( X_1 \) be measured exactly and \( X_2 \)
with classical error of variance \( \sigma_u^2 \), and let \( \boldsymbol{\uptheta}=(\beta_1,\beta_2) \). Show that the limit
of the coefficient of \( X_1 \) is \( \beta_1+(1-\lambda_{2\mid1})\rho\beta_2 \), where
\( \lambda_{2\mid1}=(1-\rho^2)/(1-\rho^2+\sigma_u^2) \). With \( \rho=0.5 \), \( \sigma_u^2=1 \), \( \beta_2=1 \) and
\( \beta_1=-0.2 \), show that the sign of the first coefficient is reversed.
:::

::: {.solution}
This is @thm-eiv-attenuation(c) with \( \bz=X_1 \), \( X=X_2 \) and \( \boldsymbol{\uppi}=\rho \), since the coefficient of the best
linear prediction of \( X_2 \) from \( X_1 \) is \( \rho \) and \( \sigma^2_{x\cdot z}=1-\rho^2 \). With the numbers given,
\( \lambda_{2\mid1}=0.75/1.75=3/7 \), so the limit is \( -0.2+(4/7)(0.5)=-0.2+2/7\approx0.086>0 \). A
regressor with a negative effect acquires a positive coefficient.
:::

### C. Going deeper

::: {#exr-eiv-rounding}
[C1]

Regressors are often recorded after rounding. Model this as \( \mathbf{S}=\mathbf{T}+\boldsymbol{\Delta} \) with
\( \mathbf{T} \) the true model matrix and \( \boldsymbol{\Delta} \) a *fixed* matrix of rounding errors, and suppose
\( \E(\Y)=\mathbf{T}\boldsymbol{\uptheta} \), \( \Cov(\Y)=\sigma^2\I \). Use @prp-lm-misspecified to show that
least squares on \( \mathbf{S} \) has bias \( -(\mathbf{S}\T\mathbf{S})^{-1}\mathbf{S}\T\boldsymbol{\Delta}\boldsymbol{\uptheta} \), that its
covariance matrix is exactly \( \sigma^2(\mathbf{S}\T\mathbf{S})^{-1} \), and that the expected residual mean square
is at least \( \sigma^2 \).
:::

::: {.solution}
Here \( \E(\Y)=\mathbf{S}\boldsymbol{\uptheta}+\boldsymbol{\updelta} \) with \( \boldsymbol{\updelta}=-\boldsymbol{\Delta}\boldsymbol{\uptheta} \). By
@prp-lm-misspecified(a) with model matrix \( \mathbf{S} \), the expectation is
\( \boldsymbol{\uptheta}+(\mathbf{S}\T\mathbf{S})^{-1}\mathbf{S}\T\boldsymbol{\updelta} \), which gives the bias. The covariance is
\( \sigma^2(\mathbf{S}\T\mathbf{S})^{-1} \) by @prp-lm-misspecified(b) with \( \bSigma=\sigma^2\I \): nothing random
has been added. For the residual mean square, @thm-rv-quadform-mean gives
\( \E\norm{(\I-\M_S)\Y}^2=\sigma^2(n-p)+\boldsymbol{\updelta}\T(\I-\M_S)\boldsymbol{\updelta} \), where \( \M_S \) projects onto
\( \C(\mathbf{S}) \) and \( (\I-\M_S)\mathbf{S}\boldsymbol{\uptheta}=\bzero \). The last term is nonnegative. Unlike random error,
fixed rounding error does not inflate the variance of the estimator. It only biases it, and the bias
is small when rounding is fine compared with the spread of the regressors.
:::
