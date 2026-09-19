# Regression calibration and SIMEX

The moment correction and IV are tied to linear moment equations. Two more portable ideas were designed for
logistic, Cox and other nonlinear regressions, where they are approximations. **Regression calibration** replaces
the unobserved regressor by its best prediction from what was observed. **SIMEX** watches the naive estimate
deteriorate as error is added, and extrapolates back to no error. In the linear model both can be exact.

## Regression calibration

Error is **nondifferential** if the distribution of \( Y \) given \( (\x,\bz,\bw) \) does not depend on \( \bw \), as for
classical error.

::: {#prp-eiv-calibration}
[Regression calibration]

Suppose \( \E(Y\mid\x,\bz)=\alpha+\bbeta_x\T\x+\bbeta_z\T\bz \) and the error is nondifferential, and let
\( \mathbf{m}(\bw,\bz)=\E(\x\mid\bw,\bz) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \E(Y\mid\bw,\bz)=\alpha+\bbeta_x\T\mathbf{m}(\bw,\bz)+\bbeta_z\T\bz \).

2. If moreover \( \Var(Y\mid\x,\bz)=\sigma^2 \) and \( \Cov(\x\mid\bw,\bz)=\mathbf{C} \) does not depend on \( (\bw,\bz) \), as when
   \( (\x,\bz,\bu) \) is jointly normal, then \( \Var(Y\mid\bw,\bz)=\sigma^2+\bbeta_x\T\mathbf{C}\bbeta_x \). The regression of \( Y \) on
   \( (1,\mathbf{m}(\bw,\bz),\bz) \) is then a homoscedastic linear model with the original coefficients.

3. Without distributional assumptions, let the error be classical with finite second moments. With
   \( \mathbf{t}=(1,\x\T,\bz\T)\T \), \( \mathbf{s}=(1,\bw\T,\bz\T)\T \), \( \A=\E(\mathbf{t}\mathbf{t}\T) \) positive definite and \( \bS=\E(\mathbf{s}\mathbf{s}\T) \),
   the best linear predictor of \( \mathbf{t} \) from \( \mathbf{s} \) is \( \hat{\mathbf{t}}=\A\bS^{-1}\mathbf{s} \), and the best linear predictor of \( Y \)
   from \( \hat{\mathbf{t}} \) is \( \boldsymbol{\uptheta}\T\hat{\mathbf{t}} \), where \( \boldsymbol{\uptheta}=(\alpha,\bbeta_x\T,\bbeta_z\T)\T \).

4. The sample version of (c) reproduces the method-of-moments estimator exactly. If \( \hat{\bS}=n^{-1}\sum_i\mathbf{s}_i\mathbf{s}_i\T \),
   \( \hat{\A}=\hat{\bS}-\bD \) for a symmetric matrix \( \bD \) with \( \hat{\A} \) nonsingular, and \( \hat{\mathbf{t}}_i=\hat{\A}\hat{\bS}^{-1}\mathbf{s}_i \), then the least
   squares coefficients of \( y \) on the \( \hat{\mathbf{t}}_i \) are \( \hat{\A}^{-1}n^{-1}\sum_i\mathbf{s}_iy_i \), which is @eq-eiv-corrected with \( \bD_i=\bD \).
:::

:::

::: {.proof}
*(a)* By the tower property and nondifferential error,
\( \E(Y\mid\bw,\bz)=\E\{\E(Y\mid\x,\bz,\bw)\mid\bw,\bz\}=\E\{\alpha+\bbeta_x\T\x+\bbeta_z\T\bz\mid\bw,\bz\} \), which is the claim.

*(b)* The law of total variance, conditioning on \( \x \) inside \( (\bw,\bz) \), gives
\( \Var(Y\mid\bw,\bz)=\E\{\Var(Y\mid\x,\bz,\bw)\mid\bw,\bz\}+\Var\{\E(Y\mid\x,\bz,\bw)\mid\bw,\bz\}=\sigma^2+\bbeta_x\T\mathbf{C}\bbeta_x \).
Under joint normality the conditional covariance of \( \x \) given \( (\bw,\bz) \) does not depend on the conditioning
values (@thm-mvn-conditional).

*(c)* The best linear predictor (@def-proj-blp, @thm-rv-blp), with the constant among the predictors, is the
projection in \( L^2 \) onto the span of the entries of \( \mathbf{s} \), so it is \( \E(\mathbf{t}\mathbf{s}\T)\bS^{-1}\mathbf{s} \). Under classical error
\( \E(\mathbf{t}\mathbf{s}\T)=\E\{\mathbf{t}(\mathbf{t}+\tilde{\bu})\T\}=\A \), and \( \A\bS^{-1}\mathbf{s} \) follows. The entries of \( \hat{\mathbf{t}} \) span the same space as those
of \( \mathbf{s} \), since \( \A\bS^{-1} \) is nonsingular. Their second-moment matrix is \( \A\bS^{-1}\bS\bS^{-1}\A=\A\bS^{-1}\A \), and
\( \E(\hat{\mathbf{t}}Y)=\A\bS^{-1}\E(\mathbf{s}Y)=\A\bS^{-1}\A\boldsymbol{\uptheta} \), because \( \E(\mathbf{s}Y)=\A\boldsymbol{\uptheta} \) as in
[Section 24.3](03-moment-correction.html). The projection coefficients are
\( (\A\bS^{-1}\A)^{-1}\A\bS^{-1}\A\boldsymbol{\uptheta}=\boldsymbol{\uptheta} \).

*(d)* With the calibrated regressors \( \hat{\mathbf{t}}_i=\hat{\A}\hat{\bS}^{-1}\mathbf{s}_i \), the normal equations (@thm-proj-normal-equations) have matrix
\( n^{-1}\sum_i\hat{\mathbf{t}}_i\hat{\mathbf{t}}_i\T=\hat{\A}\hat{\bS}^{-1}\hat{\A} \) and right side
\( \hat{\A}\hat{\bS}^{-1}n^{-1}\sum_i\mathbf{s}_iy_i \). Solving,
\( (\hat{\A}\hat{\bS}^{-1}\hat{\A})^{-1}\hat{\A}\hat{\bS}^{-1}n^{-1}\sum_i\mathbf{s}_iy_i=\hat{\A}^{-1}n^{-1}\sum_i\mathbf{s}_iy_i \).
:::

The coordinates of \( \hat{\mathbf{t}} \) for the constant and \( \bz \) reproduce them, so only the mismeasured regressors are
replaced. By (c) and (d), linear regression calibration needs no distributional assumption and *is* the moment
correction, with the standard errors of @thm-eiv-correction. Second-stage standard errors, which treat the
calibrated values as data, are too small. A **validation subsample** in which \( \x \) is measured identifies the
calibration directly, by regressing \( \x \) on \( (\bw,\bz) \) there, and needs only (a): nondifferential error and a correct
model for \( \E(\x\mid\bw,\bz) \). In a logistic regression, (a) gives \( \E\{G(\alpha+\beta X)\mid w\} \), not
\( G(\alpha+\beta\,m(w)) \), so calibration is an approximation, good when \( \Var(X\mid W) \) is small (Carroll, Ruppert,
Stefanski and Crainiceanu, 2006).

::: {.remark}
[Two meanings of "calibration"]

In [Section 12.5](../ch12-intervals-and-bands/05-calibration.html) *calibration* meant inverse prediction of
\( x_0 \) from a new response. Regression calibration instead predicts the regressor from its measurement.
:::

::: {#exm-eiv-calibration}
[Calibrating the mean reading]

The replicates estimate the error variance of the mean reading as \( 37.04 \). The calibrated pressure
\( \hat{x}_i \) has standard deviation \( 10.36 \), against \( 13.11 \) for the mean reading: each
reading is shrunk towards what age predicts. Regressing \( y \) on \( (1,\hat{x},\text{age}) \) gives \( 0.0757 \)
and \( 0.0018 \), exactly the moment correction of @exm-eiv-replicates, by @prp-eiv-calibration(d).
:::

```{.python .run #cell-simex-calibration}
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x

wbar = w.mean(axis=1)
s2_ubar = np.mean((w[:, 0] - w[:, 1]) ** 2) / 4               # error variance of the mean reading
W = np.column_stack([np.ones(n), wbar, age])
S = W.T @ W / n                                               # raw second moments of (1, wbar, age)
Suu = np.diag([0.0, s2_ubar, 0.0])

# calibration: best linear predictor of (1, x, age) from (1, wbar, age), using E(w x^T) = S - Suu
Gamma = np.linalg.solve(S, S - Suu)
x_hat = (W @ Gamma)[:, 1]                                     # predicted long-run pressure
b_rc = np.linalg.lstsq(np.column_stack([np.ones(n), x_hat, age]), y, rcond=None)[0]
b_mom = np.linalg.solve(S - Suu, W.T @ y / n)                 # method of moments, Section 24.3
print("regression calibration:", b_rc.round(4))
print("method of moments:     ", b_mom.round(4))
```

## SIMEX

Cook and Stefanski (1994) treat the naive estimator as a black box. Data cannot be made more accurate, but they
can be made *less* accurate by adding simulated error; watching the naive estimate move as the error variance
grows, we extrapolate back to zero. Let \( \sigma_u^2 \) be the known or estimated error variance of \( W \).

1. **Simulation.** For each \( \zeta \) in a grid such as \( \{0.5,1,1.5,2\} \) and \( b=1,\dots,B \), form remeasured data
   \( W_{b,i}(\zeta)=W_i+\sqrt{\zeta}\,\sigma_uV_{b,i} \), with \( V_{b,i}\iid\Normal(0,1) \) independent of the data, and compute the
   naive estimate \( \hat{\boldsymbol{\uptheta}}_b(\zeta) \) from \( (W_b(\zeta),\bz,Y) \). Average over \( b \) to get
   \( \hat{\boldsymbol{\uptheta}}(\zeta) \), and put \( \hat{\boldsymbol{\uptheta}}(0) \) equal to the naive estimate.

2. **Extrapolation.** Fit a function \( g(\zeta;\boldsymbol{\upgamma}) \) of \( \zeta \) to the points \( (\zeta,\hat{\boldsymbol{\uptheta}}(\zeta)) \), one coefficient at a
   time, and report \( g(-1;\hat{\boldsymbol{\upgamma}}) \).

The remeasured error \( U+\sqrt{\zeta}\sigma_uV \) is still classical, with variance \( (1+\zeta)\sigma_u^2 \), so \( \zeta=-1 \)
means no error. Averaging over \( b \) removes simulation noise only.

::: {#prp-eiv-simex}
[SIMEX in the linear model]

Assume the conditions of @thm-eiv-attenuation(c), with one mismeasured scalar regressor of error variance \( \sigma_u^2 \), and fix
\( \zeta\ge0 \). For each \( b \), the naive least squares coefficients from the remeasured data converge with probability
one to \( \mathbf{g}(\zeta)=(g_x(\zeta),\mathbf{g}_z(\zeta)) \), where
\[
g_x(\zeta)=\frac{\beta_x\sigma^2_{x\cdot z}}{\sigma^2_{x\cdot z}+(1+\zeta)\sigma_u^2},\qquad
\mathbf{g}_z(\zeta)=\bbeta_z+\frac{(1+\zeta)\sigma_u^2}{\sigma^2_{x\cdot z}+(1+\zeta)\sigma_u^2}\,\beta_x\boldsymbol{\uppi} ,
\]{#eq-eiv-simex-limit}

and so does their average over \( b \). Moreover:

::: {.enumerate options="label=(\alph*)"}
1. \( g_x(-1)=\beta_x \) and \( \mathbf{g}_z(-1)=\bbeta_z \);

2. every coefficient has the form \( a+b/(c+\zeta) \) with \( c=1+\sigma^2_{x\cdot z}/\sigma_u^2 \); in particular
   \( 1/g_x(\zeta) \) is linear in \( \zeta \) when \( \beta_x\ne0 \);

3. if a polynomial of fixed degree is fitted by least squares to the points \( (\zeta_k,\hat{\boldsymbol{\uptheta}}(\zeta_k)) \) on a fixed grid, the
   SIMEX estimate converges to the value at \( -1 \) of the polynomial fitted in the same way to \( (\zeta_k,\mathbf{g}(\zeta_k)) \), which in
   general differs from \( \mathbf{g}(-1) \).
:::

:::

::: {.proof}
The remeasured data satisfy @def-eiv-classical with error \( U+\sqrt\zeta\sigma_uV \), whose variance is
\( (1+\zeta)\sigma_u^2 \), and the observations \( (X_i,\bz_i,U_i,V_{b,i},\varepsilon_i) \) are independent over \( i \). By
@thm-eiv-attenuation(c), with \( \sigma_u^2 \) replaced by \( (1+\zeta)\sigma_u^2 \), the coefficients converge to
\( \lambda_\zeta\beta_x \) and \( \bbeta_z+(1-\lambda_\zeta)\beta_x\boldsymbol{\uppi} \), where
\( \lambda_\zeta=\sigma^2_{x\cdot z}/\{\sigma^2_{x\cdot z}+(1+\zeta)\sigma_u^2\} \). This is @eq-eiv-simex-limit. A finite average of
sequences converging with probability one converges to the average of the limits.
*(a)* At \( \zeta=-1 \), \( \lambda_\zeta=1 \).
*(b)* The denominator is \( \sigma_u^2(c+\zeta) \). So \( g_x(\zeta)=(\beta_x\sigma^2_{x\cdot z}/\sigma_u^2)/(c+\zeta) \), and
\( (1+\zeta)/(c+\zeta)=1-(c-1)/(c+\zeta) \) gives \( \mathbf{g}_z(\zeta)=\bbeta_z+\beta_x\boldsymbol{\uppi}-(c-1)\beta_x\boldsymbol{\uppi}/(c+\zeta) \).
*(c)* On a fixed grid the fitted polynomial's coefficients are a fixed linear function of the fitted points (the least squares
solution with a fixed model matrix), so its value at \( -1 \) is a continuous function of the points, and the points converge.
:::

So the rational form matches the limit exactly, and the rational SIMEX estimate is consistent provided the nonlinear
fit is unique and continuous in the points, which we do not verify. The quadratic, a common default because it is linear in its
parameters, recovers only part of the attenuation (@exr-eiv-simex-quadratic). With several mismeasured regressors,
or in nonlinear models, no simple extrapolant is exact and SIMEX is only approximately consistent. Its appeal is that
it needs nothing but the naive estimator.

::: {#exm-eiv-simex}
[SIMEX for the blood-pressure study]

Apply SIMEX to the regression on the mean reading and age, with error variance \( 37.04 \), grid
\( \{0,0.5,1,1.5,2\} \) and \( B=200 \) ([Figure 24.5.1](#fig-eiv-simex)). Extrapolating to \( \zeta=-1 \):

| Coefficient | naive | SIMEX, quadratic | SIMEX, \( a+b/(c+\zeta) \) | moment correction | truth |
|---|---|---|---|---|---|
| blood pressure | 0.0562 | 0.0709 | 0.0744 | 0.0757 | 0.08 |
| age | 0.0120 | 0.0042 | 0.0024 | 0.0018 | 0 |

The rational extrapolant nearly reproduces the moment correction; the quadratic stops short, and by
@prp-eiv-simex(c) its large-sample limit is \( 0.0749 \), only \( 0.77 \) of the way from the naive limit
\( 0.0582 \) to \( 0.08 \), and \( 0.0030 \) for age. A pairs bootstrap (@def-bs-bootstrap) over
patients, re-estimating the error variance and rerunning SIMEX (\( B=50 \), \( 200 \) resamples), gives standard
error \( 0.0076 \) for the rational SIMEX slope, as for the moment correction.
:::

::: {when-format="html"}
![**Figure 24.5.1.** SIMEX: averaged naive coefficients (dots), quadratic and rational extrapolants with their
values at \( \zeta=-1 \) (diamonds), the limits @eq-eiv-simex-limit (dotted) and the truth (dashed).](simex.svg){#fig-eiv-simex width=100%}
:::

::: {when-format="pdf"}
![SIMEX: averaged naive coefficients (dots), quadratic and rational extrapolants with their values at
\( \zeta=-1 \) (diamonds), the limits @eq-eiv-simex-limit (dotted) and the truth (dashed).](simex.pdf){width=100%}
:::

```{.python .run #cell-simex-simex}
from scipy.optimize import curve_fit

def naive_fits(wcol, age, y):
    """Least squares coefficients of y on (1, wcol, age) for a stack of remeasured columns."""
    m, k = wcol.shape                                         # m remeasured data sets of size k
    Xs = np.stack([np.ones((m, k)), wcol, np.broadcast_to(age, (m, k))], axis=2)
    XtX = np.einsum("mki,mkj->mij", Xs, Xs)
    Xty = np.einsum("mki,k->mi", Xs, y)
    return np.linalg.solve(XtX, Xty[:, :, None])[:, :, 0]


def simex(wbar, age, y, s2_ubar, zetas, B, rng):
    """Average naive coefficients after adding error of variance zeta * s2_ubar."""
    out = [naive_fits(wbar[None, :], age, y)[0]]
    for z in zetas[1:]:
        extra = rng.normal(0, np.sqrt(z * s2_ubar), (B, len(y)))
        out.append(naive_fits(wbar + extra, age, y).mean(axis=0))
    return np.array(out)                                      # one row per zeta


def rational(z, a, b, c):
    return a + b / (c + z)


zetas = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
G = simex(wbar, age, y, s2_ubar, zetas, 200, np.random.default_rng(2450))
est_quad, est_rat = [], []
for j in [1, 2]:                                              # coefficient of x, of age
    est_quad.append(np.polyval(np.polyfit(zetas, G[:, j], 2), -1.0))
    p, _ = curve_fit(rational, zetas, G[:, j], p0=[0.0, G[0, j], 1.5], maxfev=20000)
    est_rat.append(rational(-1.0, *p))
print("naive (zeta = 0):   ", G[0, 1:].round(4))
print("SIMEX, quadratic:   ", np.round(est_quad, 4))
print("SIMEX, rational:    ", np.round(est_rat, 4))
```

## Which correction?

In the linear model the moment correction and regression calibration coincide, rational SIMEX agrees up to
simulation noise, and IV trades precision for not needing an error variance. Beyond it, only regression calibration
and SIMEX carry over. All rest on information from outside the regression of \( Y \) on \( W \). Without it,
@prp-eiv-nonidentified is the honest summary: a range of slopes, with the naive estimate at one end.

## Exercises

### A. Check your understanding

::: {#exr-eiv-rc-simple}
[A1]

In simple regression without \( \bz \), show directly that regression calibration (@prp-eiv-calibration) with the linear calibration
\( \hat{x}=\bar{w}+\hat{\lambda}(w-\bar{w}) \), \( \hat{\lambda}=(s_{ww}-\sigma_u^2)/s_{ww} \), gives the slope \( \hat{\beta}/\hat{\lambda} \), where
\( \hat{\beta} \) is the naive slope.
:::

### B. Practice

::: {#exr-eiv-simex-quadratic}
[B1]

In simple regression with reliability \( \lambda \), write \( g(\zeta)=\lambda\beta/\{1+\zeta(1-\lambda)\} \). Show that the quadratic
through \( (\zeta,g(\zeta)) \) for \( \zeta=0,1,2 \) takes the value \( 3g(0)-3g(1)+g(2) \) at \( \zeta=-1 \). Evaluate it for \( \lambda=0.5 \)
and \( \lambda=0.64 \), as fractions of \( \beta \).
:::

::: {.solution}
The quadratic through three equally spaced points is determined by the forward differences. Written as
\( p(\zeta)=g_0+\zeta\Delta g_0+\tfrac12\zeta(\zeta-1)\Delta^2g_0 \), it gives
\( p(-1)=g_0-\Delta g_0+\Delta^2g_0=g_0-(g_1-g_0)+(g_2-2g_1+g_0)=3g_0-3g_1+g_2 \). With \( \lambda=0.5 \):
\( g_0=\beta/2 \), \( g_1=\beta/3 \), \( g_2=\beta/4 \), so \( p(-1)=\beta(3/2-1+1/4)=0.75\beta \). With \( \lambda=0.64 \):
\( g_0=0.64\beta \), \( g_1=0.64\beta/1.36 \), \( g_2=0.64\beta/1.72 \), so
\( p(-1)=\beta(1.92-1.4118+0.3721)\approx0.880\beta \). The quadratic removes most of the attenuation but not
all of it, and it does worse as the reliability falls.
:::

::: {#exr-eiv-validation}
[B2]

A validation subsample of \( n_v \) patients has the true \( x \) measured along with \( (w,z) \), and the main sample has only
\( (w,z,y) \). Describe regression calibration using the validation data. Which of the assumptions of
@prp-eiv-calibration does it need, and which does it avoid? Why does using the validation patients' \( x \) directly
in the outcome regression not help if their \( y \) is not recorded?
:::

### C. Going deeper

::: {#exr-eiv-simex-logistic}
[C1]

Simulate a logistic regression \( \Pr(Y=1\mid X)=1/(1+e^{-(\alpha+\beta X)}) \) with \( X\sim\Normal(0,1) \), \( \beta=1 \), \( n=2000 \), and a
reading \( W=X+U \) with \( \sigma_u^2=0.5 \). Compare the naive estimate, regression calibration (@prp-eiv-calibration: replace \( W \) by
\( \E(X\mid W) \)), and SIMEX with quadratic and rational extrapolants, over repeated samples. Which is closest to \( \beta \)?
:::
