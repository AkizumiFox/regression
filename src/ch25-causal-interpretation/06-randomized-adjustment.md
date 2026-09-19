# Regression adjustment in randomized experiments

In an experiment covariates are not needed to remove bias (@thm-cau-randomized); they are used to reduce variance, which
analysis of covariance does under a normal linear model with a common slope (@prp-dsn-precision). Freedman (2008a) analysed
adjustment in Neyman's framework of [Section 25.2](02-potential-outcomes.html), with no model, and found that the adjusted
estimator is biased in finite samples, that its usual standard error can be wrong, and that it can be *less* precise than
the difference in means. Lin (2013) showed that interacting the covariates with the treatment repairs the precision
problem. This section proves both halves.

## The setting

A completely randomized experiment assigns \( n_1 \) of \( n \) units to treatment and \( n_0=n-n_1 \) to control. Each unit has
potential outcomes \( y_i(1),y_i(0) \) and a vector \( \x_i\in\Real^q \) of covariates measured before randomization, hence
fixed numbers unaffected by treatment, which we centre: \( \sum_i\x_i=\bzero \). Write
\[
\begin{aligned}
\bS_{xx}&=\frac1{n-1}\sum_i\x_i\x_i\T,\qquad
\mathbf{s}_{t}=\frac1{n-1}\sum_i\x_i\bigl(y_i(t)-\bar{y}(t)\bigr),\\
\bbeta_t&=\bS_{xx}^{-1}\mathbf{s}_{t}\qquad(t=0,1),
\end{aligned}
\]
assuming \( \bS_{xx} \) positive definite. The vector \( \bbeta_t \) is the least squares slope of the potential outcomes
\( y_i(t) \) on \( \x_i \) over *all* \( n \) units, a finite-population projection coefficient in the sense of
[Section 6.11](../ch06-projections/11-population.html). No linear model is assumed.

With \( \bar{Y}_t \) and \( \bar{\x}_t \) the arm means, three estimators are in use: the **difference in means**
\( \hat{\tau}_{\mathrm{D}}=\bar{Y}_1-\bar{Y}_0 \); **Fisher's** \( \hat{\tau}_{\mathrm{F}} \), the coefficient of \( T \) in the least squares
regression of \( Y \) on \( 1,T,\x \) (the analysis of covariance of [Section 18.1](../ch18-covariance-and-design/01-ancova.html));
and **Lin's** \( \hat{\tau}_{\mathrm{L}} \), the coefficient of \( T \) in the regression on \( 1,T,\x,T\x \). All three belong to one
family. For fixed vectors \( \bb_1,\bb_0\in\Real^q \) put
\[
\hat{\tau}(\bb_1,\bb_0)=\bigl(\bar{Y}_1-\bb_1\T\bar{\x}_1\bigr)-\bigl(\bar{Y}_0-\bb_0\T\bar{\x}_0\bigr).
\]{#eq-cau-fixed-slope-estimator}

Each arm mean is corrected for the chance imbalance of its covariate mean. Then \( \hat{\tau}_{\mathrm{D}}=\hat{\tau}(\bzero,\bzero) \). By @thm-dsn-ancova(b,c), the analysis of covariance
estimates the slope by the pooled within-arm slope \( \hat{\bbeta}_{\mathrm{p}} \) and the treatment effect by the difference of
the arm means adjusted with it, so \( \hat{\tau}_{\mathrm{F}}=\hat{\tau}(\hat{\bbeta}_{\mathrm{p}},\hat{\bbeta}_{\mathrm{p}}) \). The regression with
the interaction fits a separate line in each arm, so by @exr-dsn-interacted (whose
argument is unchanged for a vector covariate) \( \hat{\tau}_{\mathrm{L}}=\hat{\tau}(\hat{\bbeta}_1,\hat{\bbeta}_0) \), with \( \hat{\bbeta}_t \) the least squares slope within arm \( t \).

## Fixed coefficients: an exact result

If the coefficients in @eq-cau-fixed-slope-estimator are fixed in advance, everything is exact.

::: {#lem-cau-fixed-coefficients}
[Adjustment with fixed coefficients]

For fixed \( \bb_1,\bb_0 \), put \( \mathbf{c}=(n_0\bb_1+n_1\bb_0)/n \) and \( \mathbf{c}^*=(n_0\bbeta_1+n_1\bbeta_0)/n \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \E\,\hat{\tau}(\bb_1,\bb_0)=\tau \);

2. the variance is
   \[
   \begin{aligned}
   \Var\hat{\tau}(\bb_1,\bb_0)&=\Var\hat{\tau}_{\mathrm{D}}\\
   &\quad+\frac{n}{n_0n_1}\Bigl[(\mathbf{c}-\mathbf{c}^*)\T\bS_{xx}(\mathbf{c}-\mathbf{c}^*)-\mathbf{c}^{*\top}\bS_{xx}\mathbf{c}^*\Bigr].
   \end{aligned}
   \]
:::

In particular the variance depends on \( (\bb_1,\bb_0) \) only through \( \mathbf{c} \), and it is smallest, over all fixed
coefficients, exactly when \( \mathbf{c}=\mathbf{c}^* \), for example when \( (\bb_1,\bb_0)=(\bbeta_1,\bbeta_0) \).
:::

::: {.proof}
Define adjusted potential outcomes \( \tilde{y}_i(t)=y_i(t)-\bb_t\T\x_i \). Because the treated units are exactly those whose
\( \tilde{y}_i(1) \) is observed, \( \hat{\tau}(\bb_1,\bb_0) \) is the difference in means for the \( \tilde{y} \)'s, and their average
effect is \( \tau-(\bb_1-\bb_0)\T\bar{\x}=\tau \) since the covariates are centred. Part (a) is @thm-cau-randomized(a). By
@eq-cau-u-representation and @lem-cau-srs, the variance is \( \{n_0/(nn_1)\}S^2_{\tilde{u}} \), where
\[
\tilde{u}_i=\tilde{y}_i(1)+\frac{n_1}{n_0}\tilde{y}_i(0)=u_i-\frac{n}{n_0}\mathbf{c}\T\x_i ,\qquad u_i=y_i(1)+\frac{n_1}{n_0}y_i(0).
\]
Expanding,
\[
S^2_{\tilde{u}}=S_u^2-2(n/n_0)\mathbf{c}\T\mathbf{s}_{xu}+(n/n_0)^2\mathbf{c}\T\bS_{xx}\mathbf{c},
\]
with
\[
\mathbf{s}_{xu}=\mathbf{s}_1+(n_1/n_0)\mathbf{s}_0=\bS_{xx}(\bbeta_1+(n_1/n_0)\bbeta_0)=(n/n_0)\bS_{xx}\mathbf{c}^*.
\]
Completing the square,
\[
S^2_{\tilde{u}}=S_u^2+\frac{n^2}{n_0^2}\Bigl[(\mathbf{c}-\mathbf{c}^*)\T\bS_{xx}(\mathbf{c}-\mathbf{c}^*)-\mathbf{c}^{*\top}\bS_{xx}\mathbf{c}^*\Bigr].
\]
Multiplying by \( n_0/(nn_1) \), and noting that \( \{n_0/(nn_1)\}S_u^2=\Var\hat{\tau}_{\mathrm{D}} \), gives (b). The bracket is
smallest iff \( \mathbf{c}=\mathbf{c}^* \), because \( \bS_{xx} \) is positive definite.
:::

The optimal \( \mathbf{c}^* \) weights the treated slope by the *control* fraction and the control slope by the *treated*
fraction: by @eq-cau-u-representation, \( \mathbf{c}^* \) is, up to the factor \( n_0/n \), the slope of
\( u_i=y_i(1)+(n_1/n_0)y_i(0) \) on \( \x_i \).

## Estimated coefficients: Freedman's critique and Lin's repair

Estimated coefficients cost exact unbiasedness but, in large experiments, nothing else. Consider a sequence of finite populations and experiments indexed by \( n \), with covariates
centred in each, and assume:

::: {.enumerate options="label=(A\arabic*)"}
1. \( n_1/n\to p\in(0,1) \);

2. the mean and the covariance matrix, over the \( n \) units, of \( \mathbf{w}_i=(\x_i\T,y_i(1),y_i(0))\T \) converge to finite
   limits, and the limit of \( \bS_{xx} \) is positive definite;

3. \( n^{-1}\sum_i\norm{\mathbf{w}_i-\bar{\mathbf{w}}}^4 \) is bounded.
:::

A sequence of random vectors \( \mathbf{g}_n \) is **bounded in probability** if for every \( \epsilon>0 \) there is an \( M \) with
\( \Pr(\norm{\mathbf{g}_n}>M)<\epsilon \) for all \( n \). If \( \mathbf{g}_n \) is bounded in probability and \( \mathbf{h}_n\to\bzero \) in probability,
then \( \mathbf{h}_n\T\mathbf{g}_n\to0 \) in probability.

::: {#thm-cau-adjustment}
[Regression adjustment in a completely randomized experiment]

Under (A1)–(A3), let \( \mathbf{b}_n=(n_1\bbeta_1+n_0\bbeta_0)/n \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\bbeta}_t-\bbeta_t\to\bzero \) for \( t=0,1 \), and \( \hat{\bbeta}_{\mathrm{p}}-\mathbf{b}_n\to\bzero \), in probability.

2. \( \sqrt n\,\bigl\{\hat{\tau}_{\mathrm{L}}-\hat{\tau}(\bbeta_1,\bbeta_0)\bigr\}\to0 \) and
   \( \sqrt n\,\bigl\{\hat{\tau}_{\mathrm{F}}-\hat{\tau}(\mathbf{b}_n,\mathbf{b}_n)\bigr\}\to0 \) in probability. So, to order \( n^{-1/2} \),
   Lin's estimator behaves like the optimal fixed-coefficient estimator and Fisher's like the fixed-coefficient estimator
   with the common coefficient \( \mathbf{b}_n \).

3. The variances of these approximating estimators satisfy
   \[
   \Var\hat{\tau}(\bbeta_1,\bbeta_0)\le\Var\hat{\tau}(\bb,\bb)\quad\text{for every }\bb,
   \]
in particular \( \Var\hat{\tau}(\bbeta_1,\bbeta_0)\le\Var\hat{\tau}_{\mathrm{D}} \); and
   \[
   \Var\hat{\tau}(\mathbf{b}_n,\mathbf{b}_n)-\Var\hat{\tau}(\bbeta_1,\bbeta_0)
   =\frac{(n_1-n_0)^2}{nn_0n_1}(\bbeta_1-\bbeta_0)\T\bS_{xx}(\bbeta_1-\bbeta_0),
   \]{#eq-cau-fisher-excess}

   which is zero iff \( n_1=n_0 \) or \( \bbeta_1=\bbeta_0 \). Fisher's approximating estimator has *larger* variance than the
   difference in means iff \( \mathbf{b}_n\T\bS_{xx}\mathbf{b}_n>2\,\mathbf{b}_n\T\bS_{xx}\mathbf{c}^* \).
:::

:::

::: {.proof}
(a) Fix components \( j,l \) of \( \mathbf{w} \), and let \( z_i \) be \( w_{ij}-\bar{w}_j \) or
\( (w_{ij}-\bar{w}_j)(w_{il}-\bar{w}_l) \). The treated units are a simple random sample of size \( n_1 \), so by @lem-cau-srs the
mean of the \( z_i \) over the treated arm has mean equal to their mean over all units and variance at most \( S_z^2/n_1 \), where
\( S_z^2\le(n-1)^{-1}\sum_iz_i^2\le(n-1)^{-1}\sum_i(1+\norm{\mathbf{w}_i-\bar{\mathbf{w}}}^4) \), which is bounded by (A3). By (A1) and
Chebyshev's inequality, every such arm mean differs from the corresponding mean over all units by a quantity tending to zero in
probability. The within-arm covariance matrix of \( \mathbf{w} \) is \( n_1/(n_1-1) \) times the arm mean of the products minus the
product of the arm means, so it differs from the covariance matrix over all units by a quantity tending to zero in probability.
In particular the arm's \( \x \)-covariance matrix and its covariance vector of \( \x \) with \( y(1) \) are close to \( \bS_{xx} \) and \( \mathbf{s}_1 \).
Since the limit of \( \bS_{xx} \) is positive definite and matrix inversion is continuous there, and \( \bbeta_1 \) converges by (A2),
\( \hat{\bbeta}_1-\bbeta_1\to\bzero \) in probability. The control arm is the same. The pooled slope is
\[
\begin{aligned}
\hat{\bbeta}_{\mathrm{p}}={}&\bigl[(n_1-1)\hat{\bS}_{xx,1}+(n_0-1)\hat{\bS}_{xx,0}\bigr]^{-1}\\
&\times\bigl[(n_1-1)\hat{\mathbf{s}}_1+(n_0-1)\hat{\mathbf{s}}_0\bigr]
\end{aligned}
\]
in terms of the arm moments, which by the same argument is close to
\( \bigl[(n-2)\bS_{xx}\bigr]^{-1}\bS_{xx}\bigl[(n_1-1)\bbeta_1+(n_0-1)\bbeta_0\bigr] \), and this differs from \( \mathbf{b}_n \) by \( O(1/n) \).

(b) By @eq-cau-fixed-slope-estimator,
\[
\hat{\tau}_{\mathrm{L}}-\hat{\tau}(\bbeta_1,\bbeta_0)=-(\hat{\bbeta}_1-\bbeta_1)\T\bar{\x}_1+(\hat{\bbeta}_0-\bbeta_0)\T\bar{\x}_0.
\]
By
@lem-cau-srs, \( \E\bar{\x}_1=\bzero \) and \( \E\norm{\bar{\x}_1}^2=(1/n_1-1/n)\tr\bS_{xx} \), so \( n\E\norm{\bar{\x}_1}^2 \) is bounded, and
\( \sqrt n\,\bar{\x}_1 \) is bounded in probability by Markov's inequality; likewise \( \sqrt n\,\bar{\x}_0 \). With (a), both products tend to
zero in probability. The argument for \( \hat{\tau}_{\mathrm{F}} \) is identical, with \( \hat{\bbeta}_{\mathrm{p}}-\mathbf{b}_n \) in both places.

(c) Apply @lem-cau-fixed-coefficients. The pair \( (\bbeta_1,\bbeta_0) \) has \( \mathbf{c}=\mathbf{c}^* \), the minimum, and
\( (\bb,\bb) \) has \( \mathbf{c}=\bb \). For \( \bb=\mathbf{b}_n \), \( \mathbf{b}_n-\mathbf{c}^*=\{(n_1-n_0)/n\}(\bbeta_1-\bbeta_0) \), which gives @eq-cau-fisher-excess. Comparing with \( \bb=\bzero \), the variance is larger iff
\[
(\mathbf{b}_n-\mathbf{c}^*)\T\bS_{xx}(\mathbf{b}_n-\mathbf{c}^*)>\mathbf{c}^{*\top}\bS_{xx}\mathbf{c}^*,
\]
which simplifies to the stated condition.
:::

Part (b) concerns the leading term only. Because the slopes and \( \bar{\x}_t \) come from the same assignment,
\( \hat{\tau}_{\mathrm{L}} \) and \( \hat{\tau}_{\mathrm{F}} \) have a bias of order \( 1/n \) under moment conditions, negligible beside a standard
error of order \( n^{-1/2} \) (Freedman 2008a; Lin 2013). With a finite-population central limit theorem (Li
and Ding 2017), part (b) makes both estimators asymptotically normal with the variances of part (c).

**Freedman's critique** is part (c) for Fisher's estimator: its common slope converges to \( \mathbf{b}_n \), which weights each
arm's slope by the arm's *own* size, while the optimal \( \mathbf{c}^* \) uses the *other* arm's size. With an unbalanced design and an effect that varies with the covariate, the analysis of covariance can lose to no
adjustment at all. **Lin's repair** is part (c) for the interacted regression, which converges to the optimal
coefficients and is asymptotically at least as precise as both alternatives *whether or not the linear model is true*. When
\( n_1=n_0 \) the two adjusted estimators agree to first order.

::: {#exm-cau-freedman-lin}
[When adjustment hurts]

A synthetic population of \( n=400 \) units has one standardized covariate, on which the response depends under treatment
(\( \beta_1=2.007 \)) but not under control (\( \beta_0=-0.009 \)); \( \tau=1.003 \). With
\( \beta_0\approx0 \) and \( S_{xx}=1 \), the condition of part (c) reads \( p\beta_1(p\beta_1-2(1-p)\beta_1)>0 \), that is \( p>2/3 \).
[Figure 25.6.1](#fig-cau-adjustment-variance) compares \( n \) times the variance of each estimator over \( 20000 \)
randomizations with the exact variances of the approximating fixed-coefficient estimators:

| fraction treated | difference in means | Fisher (ANCOVA) | Lin (interacted) |
|---|---|---|---|
| \( 0.2 \) | \( 20.88 \) (\( 20.81 \)) | \( 14.03 \) (\( 13.87 \)) | \( 4.80 \) (\( 4.73 \)) |
| \( 0.5 \) | \( 6.01 \) (\( 6.02 \)) | \( 2.04 \) (\( 2.03 \)) | \( 2.03 \) (\( 2.03 \)) |
| \( 0.8 \) | \( 5.20 \) (\( 5.20 \)) | \( 13.50 \) (\( 13.37 \)) | \( 4.26 \) (\( 4.23 \)) |

(simulated; values from @lem-cau-fixed-coefficients in brackets). At \( p=0.8 \) the analysis of covariance has more than twice
the variance of the unadjusted estimate while Lin's estimator still gains; at \( p=0.5 \) both cut it by two thirds.
:::

::: {when-format="html"}
![**Figure 25.6.1.** \( n \) times the variance of the three estimators in @exm-cau-freedman-lin as a function of the
fraction treated. Curves: exact variances of the fixed-coefficient estimators that the estimators approach (@thm-cau-adjustment). Points: simulation over \( 20000 \) randomizations. The vertical line is \( p=2/3 \), beyond which
the analysis of covariance loses to the difference in means.](adjustment_variance.svg){#fig-cau-adjustment-variance width=72%}
:::

::: {when-format="pdf"}
![\( n \) times the variance of the three estimators in @exm-cau-freedman-lin as a function of the
fraction treated. Curves: exact variances of the fixed-coefficient estimators that the estimators approach (@thm-cau-adjustment). Points: simulation over \( 20000 \) randomizations. The vertical line is \( p=2/3 \), beyond which
the analysis of covariance loses to the difference in means.](adjustment_variance.pdf){width=72%}
:::

```{.python .run #cell-adjustment-population}
import numpy as np

rng = np.random.default_rng(2506)
n = 400
x = rng.normal(size=n)
x = (x - x.mean()) / x.std(ddof=1)                    # centred, S_xx = 1
y0 = rng.normal(size=n)                               # control: no dependence on x
y1 = 1 + 2 * x + rng.normal(size=n)                   # treated: slope 2
tau = np.mean(y1 - y0)
beta1, beta0 = np.polyfit(x, y1, 1)[0], np.polyfit(x, y0, 1)[0]

def estimators(T):
    """Difference in means, ANCOVA and Lin's estimator for assignment rows T (reps x n)."""
    n1 = T.sum(axis=1, keepdims=True); n0 = n - n1
    yobs = np.where(T, y1, y0)
    xb1, xb0 = (T * x).sum(1, keepdims=True) / n1, (~T * x).sum(1, keepdims=True) / n0
    yb1, yb0 = (T * yobs).sum(1, keepdims=True) / n1, (~T * yobs).sum(1, keepdims=True) / n0
    dx = np.where(T, x - xb1, x - xb0)                # deviations from the arm means
    dy = np.where(T, yobs - yb1, yobs - yb0)
    b1 = (T * dx * dy).sum(1) / (T * dx ** 2).sum(1)  # within-arm slopes
    b0 = (~T * dx * dy).sum(1) / (~T * dx ** 2).sum(1)
    bp = (dx * dy).sum(1) / (dx ** 2).sum(1)          # pooled within-arm slope
    xb1, xb0, yb1, yb0 = xb1[:, 0], xb0[:, 0], yb1[:, 0], yb0[:, 0]
    unadj = yb1 - yb0
    ancova = unadj - bp * (xb1 - xb0)
    lin = (yb1 - b1 * xb1) - (yb0 - b0 * xb0)         # x is centred, so x-bar = 0
    return unadj, ancova, lin
```

```{.python .run #cell-adjustment-quick}
p, reps = 0.8, 2000
T = np.array([rng.permutation(n) < p * n for _ in range(reps)])
for name, est in zip(["difference in means", "ANCOVA", "Lin (interacted)"], estimators(T)):
    print(f"{name:20s} mean {est.mean():6.3f}  (tau = {tau:.3f})   n * variance {n * est.var():6.2f}")
```

## Standard errors

Freedman (2008a) showed that the classical standard error of the analysis of covariance can be inconsistent in either
direction. Lin (2013) showed that the sandwich standard errors of
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (@thm-het-sandwich) for the interacted regression
are consistent or asymptotically conservative, for the reason seen in @exr-cau-hc2: the residuals estimate the variances of
@exr-cau-lin-residuals arm by arm, and the one term that cannot be estimated enters with a negative sign.

In @exm-cau-freedman-lin at \( p=0.8 \), over \( 10000 \) randomizations, nominal \( 95\% \) intervals covered \( \tau \) with
frequency \( 0.997 \) for the difference in means with Neyman's standard error, \( 0.921 \) for the
analysis of covariance with its classical standard error, and \( 0.984 \) for Lin's estimator with the HC2 standard
error. The first and last are conservative because the effect varies across units.

::: {.idea}
In a randomized experiment, adjust only for covariates measured before randomization; include their interactions with the
treatment (with the covariates centred at their overall mean); and use sandwich standard errors. Then adjustment cannot hurt in
large samples, whatever the true relation between covariates and responses, and it helps as much as a linear adjustment can.
:::

Freedman (2008b) extends the critique to several treatments; the interacted regression, with one set of slopes per arm,
extends Lin's repair. With many covariates relative to \( n \) the \( O(1/n) \) terms are no
longer negligible.

## Exercises

### A. Check your understanding

::: {#exr-cau-common-b}
[A1]

With one covariate and a fixed common coefficient \( b \), show from @lem-cau-fixed-coefficients that \( \hat{\tau}(b,b) \) has smaller
variance than the difference in means iff \( b \) lies strictly between \( 0 \) and \( 2c^* \). What does this say about guessing a slope
from a previous study?
:::

::: {.solution}
With \( \mathbf{c}=b \) and \( S_{xx}>0 \), the variance difference is \( \{n/(n_0n_1)\}S_{xx}\bigl[(b-c^*)^2-c^{*2}\bigr]=\{n/(n_0n_1)\}S_{xx}\,b(b-2c^*) \),
which is negative iff \( b \) is strictly between \( 0 \) and \( 2c^* \). A slope guessed from earlier data, fixed before randomization,
keeps the estimator exactly unbiased, and it reduces the variance as long as it has the right sign and is less than twice the
optimal value: a crude guess is enough.
:::

### B. Practice

::: {#exr-cau-lin-residuals}
[B1]

Let \( e_i(t)=y_i(t)-\bar{y}(t)-\bbeta_t\T\x_i \) be the residuals of the finite-population regressions. Show that
\[
\Var\hat{\tau}(\bbeta_1,\bbeta_0)=\frac{S^2_{e(1)}}{n_1}+\frac{S^2_{e(0)}}{n_0}-\frac{S^2_{e(1)-e(0)}}{n},
\]
Neyman's formula with the potential outcomes replaced by their residuals. Deduce that the proportional gain over the difference in
means is large when the covariates explain much of both potential outcomes.
:::

::: {.solution}
As in the proof of @lem-cau-fixed-coefficients, \( \hat{\tau}(\bbeta_1,\bbeta_0) \) is the difference in means for the adjusted
outcomes \( y_i(t)-\bbeta_t\T\x_i=\bar{y}(t)+e_i(t) \). Adding constants changes nothing in @eq-cau-neyman-variance, so its variance is
Neyman's formula for the \( e_i(t) \), whose unit effects are \( e_i(1)-e_i(0) \) up to a constant. The first two terms are
\( (1-R_t^2) \) times those of the unadjusted formula, where \( R_t^2 \) is the finite-population \( R^2 \) of \( y(t) \) on \( \x \).
:::

::: {#exr-cau-harm-region}
[B2]

In the setting of @exm-cau-freedman-lin, let \( \beta_0=0 \), \( \beta_1=\beta\ne0 \) and \( S_{xx}=1 \), and write \( p=n_1/n \). Show that the
excess @eq-cau-fisher-excess, multiplied by \( n \), is \( (2p-1)^2\beta^2/\{p(1-p)\} \), and that the analysis of covariance loses to the
difference in means iff \( p>2/3 \). What happens if instead \( \beta_1=0 \) and \( \beta_0=\beta \)?
:::

### C. Going deeper

::: {#exr-cau-imputation}
[C1]

Show that Lin's estimator is the average over *all* \( n \) units of the difference between the two arms' fitted regression lines:
if \( \hat{\alpha}_t+\hat{\bbeta}_t\T\x \) is the least squares fit in arm \( t \), then
\[
\hat{\tau}_{\mathrm{L}}=n^{-1}\sum_i\bigl\{(\hat{\alpha}_1+\hat{\bbeta}_1\T\x_i)-(\hat{\alpha}_0+\hat{\bbeta}_0\T\x_i)\bigr\}.
\]
Relate each arm's term
to the regression estimator of a population mean in survey sampling.
:::
