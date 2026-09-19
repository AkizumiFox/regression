# Reliability and moment correction

This section supplies the missing information of [Section 24.2](02-measurement-error-model.html) directly, as a
reliability ratio, an error covariance or a variance ratio, and solves @thm-eiv-attenuation for \( \boldsymbol{\uptheta} \).

## Correcting with a known reliability

If \( \lambda \) is known, for instance from a test–retest correlation (@exr-eiv-reliability-correlation in
[Section 24.1](01-attenuation.html)), then
\( \hat{\beta}/\lambda \) is consistent in simple regression. For correlations this is Spearman's (1904) **correction
for attenuation**, \( \operatorname{Corr}(X_1,X_2)=\operatorname{Corr}(W_1,W_2)/\sqrt{\lambda_1\lambda_2} \) (@exr-eiv-spearman).
With covariates the right divisor is \( \lambda_{x\mid z} \) (@exr-eiv-conditional-reliability). Dividing by \( \lambda \)
undercorrects and leaves the covariates' bias untouched.

## A known error covariance: the corrected estimator

Work with raw second moments, so the intercept needs no special treatment. Collect the true regressors, *including* the constant, in \( \mathbf{t}=(1,\x\T,\bz\T)\T\in\Real^p \), and the observed
ones in \( \mathbf{s}=\mathbf{t}+\tilde{\bu} \), where \( \tilde{\bu}=(0,\bu\T,\bzero\T)\T \) has zeros for the constant and for the exactly
measured \( \bz \). Write the model as \( Y=\mathbf{t}\T\boldsymbol{\uptheta}+\varepsilon \). Under classical error,
\[
\E(\mathbf{s}\mathbf{s}\T)=\E(\mathbf{t}\mathbf{t}\T)+\tilde{\bSigma},\qquad
\E(\mathbf{s}Y)=\E(\mathbf{t}\mathbf{t}\T)\boldsymbol{\uptheta},
\qquad \tilde{\bSigma}=\E(\tilde{\bu}\tilde{\bu}\T),
\]
since the cross terms vanish. So \( \boldsymbol{\uptheta}=\{\E(\mathbf{s}\mathbf{s}\T)-\tilde{\bSigma}\}^{-1}\E(\mathbf{s}Y) \): subtract the
error's share of the moment matrix before inverting. To cover a known \( \tilde{\bSigma} \) and one estimated from
replicates at once, let observation \( i \) carry a matrix \( \bD_i \), computed from its own data, with expectation
\( \tilde{\bSigma} \). If \( \tilde{\bSigma} \) is known, \( \bD_i=\tilde{\bSigma} \).

::: {#thm-eiv-correction}
[Method-of-moments correction]

Let \( (\mathbf{t}_i,\tilde{\bu}_i,\varepsilon_i,\bD_i) \), \( i=1,2,\dots \), be independent copies of
\( (\mathbf{t},\tilde{\bu},\varepsilon,\bD) \), where \( \mathbf{s}_i=\mathbf{t}_i+\tilde{\bu}_i \) and \( y_i=\mathbf{t}_i\T\boldsymbol{\uptheta}+\varepsilon_i \)
are observed, and \( \bD_i \) is a symmetric \( p\times p \) matrix computed from the data of observation \( i \). Assume
that \( \E(\varepsilon)=0 \), \( \E(\mathbf{t}\varepsilon)=\bzero \), \( \E(\tilde{\bu}\varepsilon)=\bzero \), \( \E(\mathbf{t}\tilde{\bu}\T)=\bzero \),
\( \E(\bD)=\tilde{\bSigma}=\E(\tilde{\bu}\tilde{\bu}\T) \), and that \( \A=\E(\mathbf{t}\mathbf{t}\T) \) is positive definite. Let
\[
\hat{\A}_n=\frac1n\sum_{i=1}^n\bigl(\mathbf{s}_i\mathbf{s}_i\T-\bD_i\bigr),\qquad
\hat{\boldsymbol{\uptheta}}_{\text{mm}}=\hat{\A}_n^{-1}\,\frac1n\sum_{i=1}^n\mathbf{s}_iy_i .
\]{#eq-eiv-corrected}

::: {.enumerate options="label=(\alph*)"}
1. With probability one, \( \hat{\A}_n \) is nonsingular for all large \( n \) and \( \hat{\boldsymbol{\uptheta}}_{\text{mm}}\to\boldsymbol{\uptheta} \).

2. If \( \mathbf{t} \), \( \tilde{\bu} \) and \( \varepsilon \) have finite fourth moments and \( \bD \) has finite second moments, then
   \[
   \sqrt n(\hat{\boldsymbol{\uptheta}}_{\text{mm}}-\boldsymbol{\uptheta})\to\Normal_p\bigl(\bzero,\A^{-1}\B\A^{-1}\bigr)
   \quad\text{in distribution},
   \]
   where \( \B=\E(\boldsymbol{\uppsi}\boldsymbol{\uppsi}\T) \) and \( \boldsymbol{\uppsi}=\mathbf{s}(y-\mathbf{s}\T\boldsymbol{\uptheta})+\bD\boldsymbol{\uptheta} \).

3. Under the conditions of (b), \( \hat{\A}_n^{-1}\hat{\B}_n\hat{\A}_n^{-1}\to\A^{-1}\B\A^{-1} \) in probability, where
   \( \hat{\B}_n=n^{-1}\sum_i\hat{\boldsymbol{\uppsi}}_i\hat{\boldsymbol{\uppsi}}_i\T \) and \( \hat{\boldsymbol{\uppsi}}_i \) is \( \boldsymbol{\uppsi}_i \) with
   \( \hat{\boldsymbol{\uptheta}}_{\text{mm}} \) in place of \( \boldsymbol{\uptheta} \).
:::

:::

::: {.proof}
*(a)* By the strong law of large numbers, \( n^{-1}\sum_i\mathbf{s}_i\mathbf{s}_i\T\to\E(\mathbf{s}\mathbf{s}\T)=\A+\tilde{\bSigma} \)
(the cross terms vanish since \( \E(\mathbf{t}\tilde{\bu}\T)=\bzero \)) and \( n^{-1}\sum_i\bD_i\to\tilde{\bSigma} \). So
\( \hat{\A}_n\to\A \), which is nonsingular, and \( \hat{\A}_n \) is nonsingular for all large \( n \) because the determinant is
continuous. Also
\( n^{-1}\sum_i\mathbf{s}_iy_i\to\E\{(\mathbf{t}+\tilde{\bu})(\mathbf{t}\T\boldsymbol{\uptheta}+\varepsilon)\}=\A\boldsymbol{\uptheta} \). By continuity of the
inverse, \( \hat{\boldsymbol{\uptheta}}_{\text{mm}}\to\A^{-1}\A\boldsymbol{\uptheta}=\boldsymbol{\uptheta} \).

*(b)* Directly from the definitions,
\[
\frac1n\sum_i\mathbf{s}_iy_i-\hat{\A}_n\boldsymbol{\uptheta}=\frac1n\sum_i\bigl\{\mathbf{s}_iy_i-\mathbf{s}_i\mathbf{s}_i\T\boldsymbol{\uptheta}+\bD_i\boldsymbol{\uptheta}\bigr\}
=\frac1n\sum_i\boldsymbol{\uppsi}_i ,
\]
so \( \sqrt n(\hat{\boldsymbol{\uptheta}}_{\text{mm}}-\boldsymbol{\uptheta})=\hat{\A}_n^{-1}\,n^{-1/2}\sum_i\boldsymbol{\uppsi}_i \). The \( \boldsymbol{\uppsi}_i \) are independent and
identically distributed with mean
\[
\E(\boldsymbol{\uppsi})=\E(\mathbf{s}Y)-\E(\mathbf{s}\mathbf{s}\T)\boldsymbol{\uptheta}+\E(\bD)\boldsymbol{\uptheta}
=\A\boldsymbol{\uptheta}-(\A+\tilde{\bSigma})\boldsymbol{\uptheta}+\tilde{\bSigma}\boldsymbol{\uptheta}=\bzero .
\]
Each entry of \( \boldsymbol{\uppsi} \) is a sum of products of at most two entries of \( (\mathbf{t},\tilde{\bu},\varepsilon) \), plus a linear
function of \( \bD \), so the moment conditions give \( \B \) finite entries. By the central limit theorem
(the remark in [Section 14.4](../ch14-correlation-lack-of-fit-prediction/04-correlation-coefficient.html)),
\( \boldsymbol{\upxi}_n=n^{-1/2}\sum_i\boldsymbol{\uppsi}_i\to\Normal_p(\bzero,\B) \). Split
\[
\hat{\A}_n^{-1}\boldsymbol{\upxi}_n=\A^{-1}\boldsymbol{\upxi}_n+(\hat{\A}_n^{-1}-\A^{-1})\boldsymbol{\upxi}_n .
\]

By continuous mapping ([Section 19.3](../ch19-theory-of-departures/03-non-normal-errors.html)) the first term
converges in distribution to \( \A^{-1}\Normal_p(\bzero,\B)=\Normal_p(\bzero,\A^{-1}\B\A^{-1}) \). Each entry of the second is a
sum of products of an entry of \( \hat{\A}_n^{-1}-\A^{-1} \), tending to \( 0 \), and an entry of \( \boldsymbol{\upxi}_n \); by
the scalar Slutsky lemma it tends to \( 0 \) in probability. The sum form of Slutsky's lemma finishes the proof.

*(c)* Write \( \mathbf{G}_i=\mathbf{s}_i\mathbf{s}_i\T-\bD_i \) and \( \boldsymbol{\updelta}=\hat{\boldsymbol{\uptheta}}_{\text{mm}}-\boldsymbol{\uptheta} \). Then
\( \hat{\boldsymbol{\uppsi}}_i=\boldsymbol{\uppsi}_i-\mathbf{G}_i\boldsymbol{\updelta} \) and
\[
\hat{\B}_n=\frac1n\sum_i\boldsymbol{\uppsi}_i\boldsymbol{\uppsi}_i\T
-\frac1n\sum_i\bigl(\boldsymbol{\uppsi}_i\boldsymbol{\updelta}\T\mathbf{G}_i\T+\mathbf{G}_i\boldsymbol{\updelta}\boldsymbol{\uppsi}_i\T\bigr)
+\frac1n\sum_i\mathbf{G}_i\boldsymbol{\updelta}\boldsymbol{\updelta}\T\mathbf{G}_i\T .
\]
The first term tends to \( \B \) by the strong law. Every entry of the second is at most
\( 2\norm{\boldsymbol{\updelta}}\,n^{-1}\sum_i\norm{\boldsymbol{\uppsi}_i}\norm{\mathbf{G}_i} \) in absolute value (with the Frobenius norm for
\( \mathbf{G}_i \)). The average converges to \( \E(\norm{\boldsymbol{\uppsi}}\norm{\mathbf{G}}) \), which is finite by the Cauchy–Schwarz
inequality, and \( \norm{\boldsymbol{\updelta}}\to0 \) by (a). Similarly the third term is at most
\( \norm{\boldsymbol{\updelta}}^2\,n^{-1}\sum_i\norm{\mathbf{G}_i}^2\to0 \). Combine with \( \hat{\A}_n\to\A \).
:::

In simple regression with known error variance, @eq-eiv-corrected is \( s_{wy}/(s_{ww}-\sigma_u^2) \), with sample
moments taken with divisor \( n \), the naive slope
divided by the estimated reliability. The covariance \( \A^{-1}\B\A^{-1} \) has the sandwich form of
@thm-het-sandwich, to which it reduces when \( \tilde{\bu}=\bzero \) and \( \bD=\bzero \). No assumption about the
distribution of \( X \) or homoscedasticity is needed.

**Replicates.** Suppose each unit has \( m\ge2 \) readings \( \bw_{ij}=\x_i+\bu_{ij} \), \( j=1,\dots,m \), with errors independent
across \( j \), independent of \( (\x_i,\bz_i,\varepsilon_i) \), and with covariance \( \bSigma_{uu} \). Use the mean
\( \bar{\bw}_i \) as the observed regressor. Its error \( \bar{\bu}_i \) has covariance \( \bSigma_{uu}/m \), and an unbiased estimate of that covariance
from unit \( i \) alone is
\[
\frac1{m(m-1)}\sum_{j=1}^m(\bw_{ij}-\bar{\bw}_i)(\bw_{ij}-\bar{\bw}_i)\T ,
\]
padded with zeros to make \( \bD_i \) (@exr-eiv-replicate-D). For \( m=2 \) and a scalar regressor this is \( (w_{i1}-w_{i2})^2/4 \).
Because \( \bD_i \) enters \( \boldsymbol{\uppsi}_i \), the sandwich of (c) includes the uncertainty from estimating the error
variance. An external estimate treated as known would need that uncertainty added separately.

::: {#exm-eiv-replicates}
[Correcting with two readings]

In the running example, regress \( y \) on \( (1,\bar{w},\text{age}) \), where \( \bar{w} \) is the mean of the two
readings. The reliability of the mean is \( 144/(144+40.5)=0.780 \), and given age it is
\( 108/(108+40.5)=0.727 \). The naive coefficients are \( 0.0562 \) (standard error
\( 0.0054 \)) for blood pressure and \( 0.0120 \) (standard error
\( 0.0070 \)) for age. The replicates estimate \( \sigma_u^2 \) as \( 74.1 \) (true \( 81 \)), and
@eq-eiv-corrected gives
\[
\tilde{\beta}_x=0.0757\ (0.0076),\qquad \tilde{\beta}_{\text{age}}=0.0018\ (0.0076),
\]
with sandwich standard errors. The 95% interval for \( \beta_x \) is \( [0.0608,\,0.0905] \), and
the spurious age effect is gone. The standard error is \( 1.40 \) times the naive one: the bias is paid
for in variance. With \( \sigma_u^2=81 \) known the estimate would be \( 0.0782 \). In \( 10000 \)
repeated samples the corrected intervals covered \( 0.08 \) in a proportion \( 0.949 \), the naive ones in
\( 0.013 \).
:::

```{.python .run #cell-moment-correction-correct}
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x

def corrected(w, age, y):
    """Moment-corrected coefficients of (1, x, age) and their sandwich standard errors."""
    n = len(y)
    W = np.column_stack([np.ones(n), w.mean(axis=1), age])   # observed regressors
    D = np.zeros((n, 3, 3))
    D[:, 1, 1] = (w[:, 0] - w[:, 1]) ** 2 / 4                  # E D_i = Var(error of the mean)
    A = (W.T @ W - D.sum(axis=0)) / n
    theta = np.linalg.solve(A, W.T @ y / n)
    psi = W * (y - W @ theta)[:, None] + D @ theta               # estimating function, one row per patient
    B = psi.T @ psi / n
    Ainv = np.linalg.inv(A)
    cov = Ainv @ B @ Ainv / n
    return theta, np.sqrt(np.diag(cov))

theta, se = corrected(w, age, y)
W = np.column_stack([np.ones(n), w.mean(axis=1), age])
naive = np.linalg.lstsq(W, y, rcond=None)[0]
s2_u = np.mean((w[:, 0] - w[:, 1]) ** 2) / 2                   # estimate of sigma_u^2 (true 81)
print(f"sigma_u^2 estimate {s2_u:.1f}")
print("naive     (1, x, age):", naive.round(4))
print("corrected (1, x, age):", theta.round(4), " se:", se.round(4))
```

```{.python .run #cell-moment-correction-coverage}
def coverage(reps, seed):
    """Coverage of nominal 95% intervals for beta_x = 0.08: corrected (sandwich) and naive."""
    r = np.random.default_rng(seed)
    hits_c = hits_n = 0
    for _ in range(reps):
        a = r.normal(55, 10, n)
        xx = 130 + 0.6 * (a - 55) + r.normal(0, np.sqrt(108), n)
        ww = xx[:, None] + r.normal(0, 9, (n, 2))
        yy = -4 + 0.08 * xx + r.normal(0, 1.0, n)
        th, s = corrected(ww, a, yy)
        hits_c += abs(th[1] - 0.08) <= 1.96 * s[1]
        Wr = np.column_stack([np.ones(n), ww.mean(axis=1), a])
        b = np.linalg.lstsq(Wr, yy, rcond=None)[0]
        e = yy - Wr @ b
        sb = np.sqrt(e @ e / (n - 3) * np.linalg.inv(Wr.T @ Wr)[1, 1])
        hits_n += abs(b[1] - 0.08) <= 1.96 * sb
    return hits_c / reps, hits_n / reps

print("coverage (corrected, naive):", coverage(500, 2421))
```

::: {.warning}
The theorem is asymptotic, and the correction is fragile when the reliability is low. \( \hat{\A}_n \) can fail
to be positive definite, and under normality the corrected slope has no finite mean (@exr-eiv-no-mean), as
with @prp-ci-ratio-no-mean; Fuller (1987) gives modified estimators. Gleser and Hwang (1987) showed that, as in
@prp-ci-calibration-unbounded, any confidence set with guaranteed coverage has infinite expected length. Wald
intervals are reliable only when the estimated reliability is well away from zero.
:::

## Deming and orthogonal regression

Sometimes only a *ratio* of error variances is known. In a method-comparison study two instruments measure the
same specimens, \( W=X+U \) and \( Y=\alpha+\beta X+\varepsilon \), with \( \varepsilon \) the second instrument's error, and
repeatability experiments give
\[
\delta=\frac{\Var(\varepsilon)}{\Var(U)} .
\]
The question is whether \( \alpha=0 \) and \( \beta=1 \).

::: {#prp-eiv-deming}
[Deming regression]

Let \( (W,Y) \) have covariances \( \sigma_{ww},\sigma_{wy},\sigma_{yy} \) with \( \sigma_{wy}\ne0 \) and \( \sigma_{ww}\sigma_{yy}\ge\sigma_{wy}^2 \), and let
\( \delta>0 \) be known.

::: {.enumerate options="label=(\alph*)"}
1. The equations \( \sigma_{ww}=\sigma_x^2+\sigma_u^2 \), \( \sigma_{wy}=\beta\sigma_x^2 \), \( \sigma_{yy}=\beta^2\sigma_x^2+\delta\sigma_u^2 \) have exactly
   one solution with \( \sigma_x^2>0 \), namely
   \[
   \beta_\delta=\frac{\sigma_{yy}-\delta\sigma_{ww}+\sqrt{(\sigma_{yy}-\delta\sigma_{ww})^2+4\delta\sigma_{wy}^2}}{2\sigma_{wy}} ,
   \]{#eq-eiv-deming}

   and for it \( \sigma_u^2\ge0 \). The slope \( \beta_\delta \) lies in the interval of @prp-eiv-nonidentified. It tends to
   \( b_2 \) as \( \delta\to0 \) and to \( b_1 \) as \( \delta\to\infty \).

2. For \( \delta=1 \), \( \beta_1 \) is the slope of the line that minimizes the mean squared *perpendicular*
   distance from \( (W,Y) \), whose direction is an eigenvector of the covariance matrix of \( (W,Y) \) for its largest
   eigenvalue. For general \( \delta \), \( \beta_\delta \) is \( \sqrt\delta \) times the \( \delta=1 \) slope computed for
   \( (W,Y/\sqrt\delta) \).
:::

:::

::: {.proof}
*(a)* For \( \beta\ne0 \) the first two equations give \( \sigma_x^2=\sigma_{wy}/\beta \) and \( \sigma_u^2=\sigma_{ww}-\sigma_{wy}/\beta \).
Substituting in the third and multiplying by \( \beta \),
\[
q(\beta)=\sigma_{wy}\beta^2+(\delta\sigma_{ww}-\sigma_{yy})\beta-\delta\sigma_{wy}=0 .
\]
The product of the roots is \( -\delta<0 \), so there is one root of each sign, and \( \sigma_x^2>0 \) selects the root
with the sign of \( \sigma_{wy} \). That root is @eq-eiv-deming. Take \( \sigma_{wy}>0 \) (otherwise replace \( Y \) by \( -Y \)). Then
\[
q(b_1)=\frac{\sigma_{wy}}{\sigma_{ww}^2}\bigl(\sigma_{wy}^2-\sigma_{ww}\sigma_{yy}\bigr)\le0,\qquad
q(b_2)=\frac{\delta}{\sigma_{wy}}\bigl(\sigma_{ww}\sigma_{yy}-\sigma_{wy}^2\bigr)\ge0,
\]
and \( q \) is an upward parabola with \( q(0)<0 \), so its positive root lies in \( [b_1,b_2] \). Then
\( \sigma_u^2=\sigma_{ww}(1-b_1/\beta_\delta)\ge0 \). As \( \delta\to0 \) the formula tends to
\( (\sigma_{yy}+\sigma_{yy})/(2\sigma_{wy})=b_2 \). Multiplying the numerator and denominator of @eq-eiv-deming by
the conjugate of the numerator gives
\( \beta_\delta=2\delta\sigma_{wy}/\{\sqrt{(\sigma_{yy}-\delta\sigma_{ww})^2+4\delta\sigma_{wy}^2}+\delta\sigma_{ww}-\sigma_{yy}\} \), and the
denominator is \( 2\delta\sigma_{ww}+O(1) \), so \( \beta_\delta\to b_1 \) as \( \delta\to\infty \).

*(b)* Let \( \bSigma \) be the covariance matrix of \( \mathbf{v}=(W,Y)\T \). A line through the point \( \E\mathbf{v} \) with unit normal
vector \( \mathbf{n} \) has mean squared perpendicular distance \( \E\{\mathbf{n}\T(\mathbf{v}-\E\mathbf{v})\}^2=\mathbf{n}\T\bSigma\mathbf{n} \),
and a line through any other point \( \mathbf{c} \) adds \( \{\mathbf{n}\T(\E\mathbf{v}-\mathbf{c})\}^2\ge0 \). By
@thm-mat-extremal-rayleigh, \( \mathbf{n}\T\bSigma\mathbf{n} \) is minimized by an eigenvector for the smallest eigenvalue, so
the line's direction is an eigenvector \( (v_1,v_2) \) for the largest eigenvalue \( \mu \). From
\( (\sigma_{ww}-\mu)v_1+\sigma_{wy}v_2=0 \), the slope is \( v_2/v_1=(\mu-\sigma_{ww})/\sigma_{wy} \). Put \( b=v_2/v_1 \). Then
\( \mu-\sigma_{ww}=b\sigma_{wy} \), and the second eigenvector equation \( \sigma_{wy}v_1+(\sigma_{yy}-\mu)v_2=0 \) gives
\( \mu-\sigma_{yy}=\sigma_{wy}/b \). Subtracting, \( \sigma_{yy}-\sigma_{ww}=b\sigma_{wy}-\sigma_{wy}/b \), which is \( q(b)=0 \) with
\( \delta=1 \). Since \( \mu\ge\sigma_{ww} \), \( b \) has the sign of \( \sigma_{wy} \), so \( b=\beta_1 \). For the last claim, substitute
\( \sigma_{yy}/\delta \) and \( \sigma_{wy}/\sqrt\delta \) for \( \sigma_{yy} \) and \( \sigma_{wy} \) in @eq-eiv-deming with \( \delta=1 \) and multiply by
\( \sqrt\delta \).
:::

With sample moments, @eq-eiv-deming is **Deming regression**, consistent because it is continuous where
\( \sigma_{wy}\ne0 \). With \( \delta=1 \) it is **orthogonal regression** (total least squares), whose line is the first
principal axis of the scatter.

::: {.warning}
Orthogonal regression makes sense only when the errors have equal variances in the units used. Changing the
units of \( Y \) changes the fitted line by more than a change of scale. Deming regression, with \( \delta \) rescaled
along with the units, does not have this problem.
:::

::: {#exm-eiv-method-comparison}
[Comparing two assays]

Sixty specimens with true concentrations uniform on \( [2,20] \) are measured by method A with error standard deviation
\( 0.8 \) and by method B, which reads \( 0.3+1.05x \) with error standard deviation \( 1.0 \), so
\( \delta=1.5625 \). The four slopes are

| Method | Slope |
|---|---|
| least squares of B on A | 1.0042 |
| Deming, \( \delta=1.5625 \) | 1.0308 |
| orthogonal (\( \delta=1 \), unjustified) | 1.0382 |
| least squares of A on B, inverted | 1.0709 |

with Deming intercept \( 0.482 \). Over \( 20000 \) simulated studies the mean slopes were
\( 1.0262 \) for least squares, attenuated by \( 27/(27+0.64)=0.9768 \), and \( 1.0505 \)
for Deming. With B in units ten times smaller, orthogonal regression gives \( 1.0703 \) (converted back to the
original units) instead of
\( 1.0382 \).
:::

```{.python .run #cell-deming-deming}
import numpy as np


def deming(w, y, delta):
    """Slope and intercept of Deming regression; delta = Var(error in y) / Var(error in w)."""
    S = np.cov(w, y)
    s_ww, s_wy, s_yy = S[0, 0], S[0, 1], S[1, 1]
    d = s_yy - delta * s_ww
    b = (d + np.sqrt(d ** 2 + 4 * delta * s_wy ** 2)) / (2 * s_wy)
    return b, y.mean() - b * w.mean()


rng = np.random.default_rng(2430)
x = rng.uniform(2, 20, 60)                        # true concentrations
w = x + rng.normal(0, 0.8, 60)                    # method A
y = 0.3 + 1.05 * x + rng.normal(0, 1.0, 60)       # method B
delta = 1.0 ** 2 / 0.8 ** 2

b_ls = np.polyfit(w, y, 1)[0]                     # regress y on w
b_rev = 1 / np.polyfit(y, w, 1)[0]                # regress w on y, inverted
b_dem, a_dem = deming(w, y, delta)
b_orth, _ = deming(w, y, 1.0)                     # orthogonal regression
print(f"least squares {b_ls:.4f}, Deming {b_dem:.4f} (intercept {a_dem:.3f}), "
      f"orthogonal {b_orth:.4f}, reverse {b_rev:.4f}")
```

## Exercises

### A. Check your understanding

::: {#exr-eiv-spearman}
[A1]

Prove Spearman's correction: if \( W_k=X_k+U_k \) for \( k=1,2 \), with classical errors independent of each other and
of \( (X_1,X_2) \), then \( \operatorname{Corr}(W_1,W_2)=\sqrt{\lambda_1\lambda_2}\operatorname{Corr}(X_1,X_2) \).
:::

::: {#exr-eiv-replicate-D}
[A2]

Show that the matrix \( \bD_i \) defined from \( m \) replicates has expectation \( \bSigma_{uu}/m \), the covariance matrix of
\( \bar{\bu}_i \), padded with zeros.
:::

::: {.solution}
\( \bw_{ij}-\bar{\bw}_i=\bu_{ij}-\bar{\bu}_i \), so \( \x_i \) drops out. For independent vectors with common covariance
\( \bSigma_{uu} \), the sample covariance matrix \( (m-1)^{-1}\sum_j(\bu_{ij}-\bar{\bu}_i)(\bu_{ij}-\bar{\bu}_i)\T \) is unbiased for \( \bSigma_{uu} \) (@prp-rv-sample-moments).
Dividing by \( m \) gives expectation \( \bSigma_{uu}/m=\Cov(\bar{\bu}_i) \).
:::

### B. Practice

::: {#exr-eiv-conditional-reliability}
[B1]

For a scalar \( X \) and a scalar exactly measured \( Z \), show that
\[
\lambda_{x\mid z}=\frac{\lambda-\rho_{wz}^2}{1-\rho_{wz}^2},
\]
where \( \rho_{wz} \) is the correlation between \( W \) and \( Z \). Check the formula on the running example, and deduce that
\( \lambda_{x\mid z}<\lambda \) whenever \( \rho_{wz}\ne0 \) and \( \lambda<1 \).
:::

::: {.solution}
\( \sigma^2_{x\cdot z}=\sigma_x^2-\sigma_{xz}^2/\sigma_z^2 \), and since \( \sigma_{wz}=\sigma_{xz} \),
\( \sigma_{xz}^2/\sigma_z^2=\rho_{wz}^2\sigma_w^2 \). Divide the numerator and the denominator of
\( \lambda_{x\mid z}=\sigma^2_{x\cdot z}/(\sigma^2_{x\cdot z}+\sigma_u^2) \) by \( \sigma_w^2 \). The numerator becomes \( \lambda-\rho_{wz}^2 \) and the denominator
\( \lambda-\rho_{wz}^2+(1-\lambda)=1-\rho_{wz}^2 \). In the example \( \rho_{wz}^2=60^2/(225\cdot100)=0.16 \), so
\( \lambda_{x\mid z}=(0.64-0.16)/0.84=0.571 \), as in [Section 24.1](01-attenuation.html). Finally
\( (\lambda-\rho^2)/(1-\rho^2)<\lambda \) iff \( \lambda-\rho^2<\lambda-\lambda\rho^2 \) iff \( \rho^2(1-\lambda)>0 \).
:::

::: {#exr-eiv-normal-variance}
[B2]

In simple regression with centred variables and known \( \sigma_u^2 \), take \( \bD=\diag(0,\sigma_u^2) \) in @thm-eiv-correction
and suppose \( (X,U,\varepsilon) \) are jointly normal and independent. Show that the asymptotic variance of
\( \sqrt n(\tilde{\beta}-\beta) \) is
\[
\frac{\sigma_w^2\sigma_v^2+\beta^2\sigma_u^4}{\sigma_x^4},\qquad \sigma_v^2=\sigma^2+\beta^2\sigma_u^2 .
\]
Compare it with \( \sigma^2/\sigma_x^2 \), the asymptotic variance when \( X \) is observed. What does this give for the running
example with the mean of two readings?
:::

::: {.solution}
The slope's estimating function (the intercept decouples for centred variables) is
\( \psi=W(Y-W\beta)+\sigma_u^2\beta=Wv+\sigma_u^2\beta \) with \( v=\varepsilon-\beta U \), and the slope entry of \( \A \) is \( \sigma_x^2 \).
Now \( (W,v) \) is bivariate normal with \( \Cov(W,v)=-\beta\sigma_u^2 \), and for such a pair
\( \E(W^2v^2)=\sigma_w^2\sigma_v^2+2\Cov(W,v)^2 \). So
\( \Var(\psi)=\E(W^2v^2)-\{\E(Wv)\}^2=\sigma_w^2\sigma_v^2+\beta^2\sigma_u^4 \), and the sandwich is
\( \Var(\psi)/\sigma_x^4 \). It exceeds \( \sigma^2/\sigma_x^2 \) because \( \sigma_w^2\sigma_v^2>\sigma_x^2\sigma^2 \). For the example,
use \( \sigma_u^2=40.5 \) for the mean of two readings: \( \sigma_w^2=184.5 \), \( \sigma_v^2=1+0.0064\cdot40.5=1.2592 \),
and the variance is \( (184.5\cdot1.2592+0.0064\cdot1640.25)/144^2\approx0.01171 \), a standard error of about
\( \sqrt{0.01171/300}\approx0.0062 \), below the \( 0.0076 \) of @exm-eiv-replicates, where \( \sigma_u^2 \) is estimated and
age lowers the reliability.
:::

### C. Going deeper

::: {#exr-eiv-no-mean}
[C1]

In the normal structural model with known \( \sigma_u^2 \) and \( n\ge3 \), show that
\( \tilde{\beta}=s_{wy}/(s_{ww}-\sigma_u^2) \) has \( \E\lvert\tilde{\beta}\rvert=\infty \). (Condition on \( W_1,\dots,W_n \), and use the fact that
\( s_{ww} \) has a continuous density that is positive at \( \sigma_u^2 \).)
:::

::: {.solution}
Given the \( W_i \), \( (W_i,Y_i) \) is jointly normal and \( \E(Y\mid W)=a+\lambda\beta W \) with constant conditional variance
\( \tau^2>0 \) (@thm-mvn-conditional). So, given the \( W_i \), \( s_{wy}\sim\Normal\bigl(\lambda\beta s_{ww},\tau^2s_{ww}/(n-1)\bigr) \), and
\( \E(\lvert s_{wy}\rvert\mid W)\ge\tau\sqrt{2s_{ww}/\{\pi(n-1)\}} \), because a normal variable's absolute mean is at least that of its
centred version. Hence
\( \E\lvert\tilde{\beta}\rvert\ge c\,\E\bigl\{\sqrt{s_{ww}}/\lvert s_{ww}-\sigma_u^2\rvert\bigr\} \). Here \( (n-1)s_{ww}/\sigma_w^2\sim\chi^2(n-1) \), whose density
is continuous and positive on \( (0,\infty) \). So the density of \( s_{ww} \) is bounded below near \( \sigma_u^2 \), and the integral of
\( 1/\lvert s-\sigma_u^2\rvert \) diverges there.
:::
