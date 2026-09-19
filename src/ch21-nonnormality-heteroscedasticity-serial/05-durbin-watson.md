# Serial correlation and the Durbin–Watson test

When observations are taken in sequence, over time or along a transect, neighbouring errors are often alike.
A shock to the economy lasts more than one quarter, a drifting instrument affects consecutive readings, and a
variable missing from the model usually changes slowly, so its effect appears in the errors as a slow wave.
Errors that are correlated with their neighbours are **serially correlated** or **autocorrelated**. As with
heteroscedasticity, least squares stays unbiased and its usual standard errors become
wrong (@thm-dep-covariance). The error is typically in the dangerous direction: with positive autocorrelation and a
slowly varying regressor, the true variance of a slope can be many times its nominal value.
[Section 21.6](06-autoregressive-errors.html) measures this. This section is about detection.

## Looking at residuals in time order

Three displays are standard. The residuals plotted against time show runs of the same sign when the
autocorrelation is positive, and rapid alternation when it is negative. The **lag plot** of \( \hat{\varepsilon}_t \)
against \( \hat{\varepsilon}_{t-1} \) shows the lag-one dependence as a tilt. And the **residual autocorrelations**
\[
r_k=\frac{\sum_{t=k+1}^n\hat{\varepsilon}_t\hat{\varepsilon}_{t-k}}{\sum_{t=1}^n\hat{\varepsilon}_t^2},\qquad k=1,2,\dots,
\]
summarize dependence at each lag. Under the model the residuals are themselves slightly correlated, with
covariance \( \sigma^2(\I-\M) \), so the \( r_k \) are not exactly centred at zero, but for large \( n \) they are
approximately independent with standard deviation \( n^{-1/2} \).

::: {#exm-het-okun}
[A quarterly Okun's-law regression]

Okun's law, named after the economist Arthur Okun, relates changes in unemployment to output growth. The
statsmodels macroeconomic data give US quarterly series from 1959 to 2009. Regress the change in the
unemployment rate (percentage points per quarter) on real GDP growth (percent, at an annual rate) for the
\( n=202 \) quarters from the second quarter of 1959. The fitted line is
\[
\widehat{\Delta u}_t=0.227-0.0671\,g_t ,
\]
so unemployment is stable when output grows at about \( 0.227/0.0671\approx3.4 \) percent a year, and each
extra percentage point of growth lowers it by about \( 0.07 \) points per quarter.
[Figure 21.5.1](#fig-het-dw)(a,b) shows the residuals. They wander in long swings, and the lag plot tilts,
with \( r_1=0.300 \).
:::

```{.python .run #cell-durbin-watson-data}
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import brentq

macro = sm.datasets.macrodata.load_pandas().data
# change in unemployment rate
du = np.diff(macro["unemp"].to_numpy())
# GDP growth, % a year
growth = 400 * np.diff(np.log(macro["realgdp"].to_numpy()))
year = (macro["year"] + (macro["quarter"] - 1) / 4).to_numpy()[1:]
n = len(du)
X = np.column_stack([np.ones(n), growth])
beta, *_ = np.linalg.lstsq(X, du, rcond=None)
e = du - X @ beta
# Durbin-Watson statistic
d = np.sum(np.diff(e) ** 2) / np.sum(e**2)
# lag-one residual autocorrelation
r1 = np.sum(e[1:] * e[:-1]) / np.sum(e**2)
print(f"n = {n}, slope {beta[1]:.4f}, d = {d:.3f}, r1 = {r1:.3f}, "
      f"2(1 - r1) = {2 * (1 - r1):.3f}")
```

## The Durbin–Watson statistic

::: {#def-het-dw}
[Durbin–Watson statistic]

For residuals \( \hat{\varepsilon}_1,\dots,\hat{\varepsilon}_n \) in time order,
\[
d=\frac{\sum_{t=2}^n(\hat{\varepsilon}_t-\hat{\varepsilon}_{t-1})^2}{\sum_{t=1}^n\hat{\varepsilon}_t^2}
=\frac{\he\T\A\he}{\he\T\he},\qquad \A=\bD\T\bD ,
\]
where \( \bD \) is the \( (n-1)\times n \) first-difference matrix, with rows \( (\dots,0,-1,1,0,\dots) \).
:::

Expanding the numerator gives
\( d=2(1-r_1)-(\hat{\varepsilon}_1^2+\hat{\varepsilon}_n^2)/\sum_t\hat{\varepsilon}_t^2 \), so \( d\approx2(1-r_1) \). Values near \( 2 \) indicate
no lag-one correlation, values towards \( 0 \) positive correlation, and values towards \( 4 \) negative
correlation. For the Okun regression \( d=1.386 \), against \( 2(1-r_1)=1.401 \).

The matrix \( \A \) is tridiagonal with diagonal \( (1,2,\dots,2,1) \) and off-diagonal entries \( -1 \). Its spectrum is
known exactly, and that is what makes the Durbin–Watson theory work.

::: {#lem-het-dw-eigen}
[Eigenvalues of the difference matrix]

The eigenvalues of \( \A \) are
\[
\lambda_j=4\sin^2\Bigl(\frac{\pi(j-1)}{2n}\Bigr),\qquad j=1,\dots,n,
\]
with eigenvectors \( \bv_j=(\cos(\theta_j/2),\cos(3\theta_j/2),\dots,\cos((2n-1)\theta_j/2))\T \), where
\( \theta_j=\pi(j-1)/n \). In particular \( 0=\lambda_1<\lambda_2<\dots<\lambda_n<4 \), and \( \bv_1=\bone \).
:::

::: {.proof}
Let \( \theta=\theta_j \) and \( v_t=\cos((t-\tfrac12)\theta) \), and extend the sequence by \( v_0=\cos(-\theta/2)=v_1 \) and
\( v_{n+1}=\cos((n+\tfrac12)\theta) \). Since \( n\theta=\pi(j-1) \),
\( v_{n+1}=\cos(\pi(j-1)+\theta/2)=\cos(\pi(j-1)-\theta/2)=v_n \). With these boundary values every row of \( \A\bv \) has the
interior form \( (\A\bv)_t=2v_t-v_{t-1}-v_{t+1} \), because the first row \( v_1-v_2 \) equals \( 2v_1-v_0-v_2 \) and similarly
for the last. By the identity \( \cos(a-\theta)+\cos(a+\theta)=2\cos a\cos\theta \),
\( 2v_t-v_{t-1}-v_{t+1}=2(1-\cos\theta)v_t=4\sin^2(\theta/2)\,v_t \). The \( n \) values \( \theta_j\in[0,\pi) \) give \( n \)
distinct eigenvalues, so these are all of them.
:::

The distribution theory also needs a classical fact about eigenvalues of compressed matrices.

::: {#lem-het-interlace}
[Interlacing]

Let \( \B \) be \( N\times N \) symmetric with eigenvalues \( \beta_1\le\dots\le\beta_N \), and let \( \mathbf{U} \) be \( N\times m \) with
orthonormal columns. If \( \mu_1\le\dots\le\mu_m \) are the eigenvalues of \( \mathbf{U}\T\B\mathbf{U} \), then
\( \beta_k\le\mu_k\le\beta_{k+N-m} \) for \( k=1,\dots,m \).
:::

::: {.proof}
Let \( \mathbf{f}_1,\dots,\mathbf{f}_m \) be orthonormal eigenvectors of \( \mathbf{U}\T\B\mathbf{U} \) and \( \mathbf{q}_1,\dots,\mathbf{q}_N \) those of \( \B \).
The subspaces \( \mathcal S=\spn\{\mathbf{U}\mathbf{f}_1,\dots,\mathbf{U}\mathbf{f}_k\} \), of dimension \( k \), and
\( \mathcal T=\spn\{\mathbf{q}_k,\dots,\mathbf{q}_N\} \), of dimension \( N-k+1 \), have dimensions adding to \( N+1 \), so they share a
nonzero vector \( \mathbf{x}=\mathbf{U}\mathbf{f} \). Since \( \mathbf{U} \) preserves lengths,
\( \mathbf{x}\T\B\mathbf{x}/\mathbf{x}\T\mathbf{x}=\mathbf{f}\T\mathbf{U}\T\B\mathbf{U}\mathbf{f}/\mathbf{f}\T\mathbf{f}\le\mu_k \), because
\( \mathbf{f}\in\spn\{\mathbf{f}_1,\dots,\mathbf{f}_k\} \); and \( \mathbf{x}\T\B\mathbf{x}/\mathbf{x}\T\mathbf{x}\ge\beta_k \), because
\( \mathbf{x}\in\mathcal T \). Both bounds are @thm-mat-extremal-rayleigh applied to the restriction of the quadratic form to
the span of eigenvectors. So \( \beta_k\le\mu_k \). The upper bound follows in the same way from
\( \spn\{\mathbf{U}\mathbf{f}_k,\dots,\mathbf{U}\mathbf{f}_m\} \) and \( \spn\{\mathbf{q}_1,\dots,\mathbf{q}_{k+N-m}\} \), whose dimensions add to \( N+1 \).
:::

## The null distribution of d

::: {#thm-het-durbin-watson}
[Durbin–Watson]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \rank(\X)=p<n \), let the columns of the \( n\times(n-p) \) matrix
\( \Q_\perp \) be an orthonormal basis of \( \C(\X)\perpc \), and let \( \nu_1\le\dots\le\nu_{n-p} \) be the eigenvalues of
\( \Q_\perp\T\A\Q_\perp \). Put \( m=n-p \).

::: {.enumerate options="label=(\alph*)"}
1. \( d \) has the distribution of \( \sum_k\nu_k\xi_k^2/\sum_k\xi_k^2 \), where \( \xi_1,\dots,\xi_m \) are independent
   \( \Normal(0,1) \). It depends on \( \X \) but not on \( \bbeta \) or \( \sigma^2 \), and \( d \) is independent of
   \( \text{SSE}=\he\T\he \).

2. \( \E d=\bar\nu=m^{-1}\sum_k\nu_k=\tr\{(\I-\M)\A\}/m \) and
   \[
\Var d=\frac{2\sum_k(\nu_k-\bar\nu)^2}{m(m+2)} .
\]

3. \( \Pr(d\le c)=\Pr\bigl\{\sum_k(\nu_k-c)\xi_k^2\le0\bigr\} \) for every \( c \).

4. *(Bounds.)* If \( \bone\in\C(\X) \), then \( \lambda_{k+1}\le\nu_k\le\lambda_{k+p} \) for \( k=1,\dots,m \), with \( \lambda_j \) as
   in @lem-het-dw-eigen. Consequently \( d_L\le d\le d_U \) with probability one, where, for the same \( \xi_k \),
   \[
d_L=\frac{\sum_k\lambda_{k+1}\xi_k^2}{\sum_k\xi_k^2},\qquad d_U=\frac{\sum_k\lambda_{k+p}\xi_k^2}{\sum_k\xi_k^2},
\]
   whose distributions depend only on \( n \) and \( p \).
:::

:::

::: {.proof}
(a) \( \he=(\I-\M)\be=\Q_\perp\Q_\perp\T\be \), because \( \Q_\perp\Q_\perp\T \) is the projection onto \( \C(\X)\perpc \). Let
\( \boldsymbol{\upzeta}=\Q_\perp\T\be/\sigma\sim\Normal_m(\bzero,\I_m) \) (@thm-mvn-linear), and take a spectral decomposition
\( \Q_\perp\T\A\Q_\perp=\mathbf{P}\,\diag(\nu_k)\,\mathbf{P}\T \) with \( \mathbf{P} \) orthogonal. Then
\( \boldsymbol{\upxi}=\mathbf{P}\T\boldsymbol{\upzeta}\sim\Normal_m(\bzero,\I_m) \),
\( \he\T\A\he=\sigma^2\boldsymbol{\upzeta}\T\Q_\perp\T\A\Q_\perp\boldsymbol{\upzeta}=\sigma^2\sum_k\nu_k\xi_k^2 \) and
\( \he\T\he=\sigma^2\norm{\boldsymbol{\upzeta}}^2=\sigma^2\sum_k\xi_k^2 \). For independence, put \( G_k=\xi_k^2/2 \), which are
independent gamma variables with shape \( \tfrac12 \) and density proportional to \( g^{-1/2}e^{-g} \). Change variables
from \( (G_1,\dots,G_m) \) to \( S=\sum_kG_k \) and \( W_k=G_k/S \), \( k<m \); the Jacobian is \( S^{m-1} \), and the joint density
becomes proportional to \( S^{m/2-1}e^{-S}\prod_kw_k^{-1/2} \) (with \( w_m=1-\sum_{k<m}w_k \)), a product of a function of
\( S \) and a function of the \( W_k \). So \( (W_1,\dots,W_m) \) is independent of \( S \). Now \( d=\sum_k\nu_kW_k \) and
\( \text{SSE}=2\sigma^2S \).

(b) Let \( N=\sum_k\nu_k\xi_k^2 \) and \( S_2=\sum_k\xi_k^2 \), so \( N=dS_2 \) with \( d \) independent of \( S_2 \). Then
\( \E N=\E d\,\E S_2 \) and \( \E N^2=\E d^2\,\E S_2^2 \). Here \( \E N=\sum_k\nu_k \), \( \E S_2=m \),
\( \E N^2=\Var N+(\E N)^2=2\sum_k\nu_k^2+(\sum_k\nu_k)^2 \) and \( \E S_2^2=2m+m^2 \). Dividing gives
\( \E d=\bar\nu \) and \( \E d^2=\{2\sum_k\nu_k^2+(\sum_k\nu_k)^2\}/\{m(m+2)\} \), and \( \Var d=\E d^2-\bar\nu^2 \) simplifies to
the stated form. Finally \( \sum_k\nu_k=\tr(\Q_\perp\T\A\Q_\perp)=\tr(\A\Q_\perp\Q_\perp\T)=\tr\{\A(\I-\M)\} \).

(c) \( d\le c \) iff \( \sum_k\nu_k\xi_k^2\le c\sum_k\xi_k^2 \).

(d) Let \( \mathbf{P}_1 \) be the \( n\times(n-1) \) matrix whose columns are the normalized eigenvectors
\( \bv_2,\dots,\bv_n \) of \( \A \); they span \( \bone\perpc \). Since \( \bone\in\C(\X) \), \( \C(\Q_\perp)\subseteq\bone\perpc \), so
\( \Q_\perp=\mathbf{P}_1\mathbf{U} \) for an \( (n-1)\times m \) matrix \( \mathbf{U}=\mathbf{P}_1\T\Q_\perp \) with orthonormal columns. Then
\( \Q_\perp\T\A\Q_\perp=\mathbf{U}\T\diag(\lambda_2,\dots,\lambda_n)\mathbf{U} \), and @lem-het-interlace with \( N=n-1 \) gives
\( \lambda_{k+1}\le\nu_k\le\lambda_{k+1+(n-1-m)}=\lambda_{k+p} \). Multiplying by \( \xi_k^2 \), summing, and dividing by
\( \sum_k\xi_k^2 \) gives \( d_L\le d\le d_U \). The \( \lambda_j \) depend only on \( n \).
:::

Part (a) reduces the exact null distribution to that of a ratio of quadratic forms in normal variables, of the
kind studied in [Chapter 4](../ch04-quadratic-forms/index.html). Part (c) turns its distribution function into
that of a single indefinite quadratic form, whose distribution function Imhof (1961) showed how to compute by
numerical inversion of its characteristic function. So the exact p-value, given \( \X \), is a routine
computation, and the listing below does it. When Durbin and Watson (1950, 1951) proposed the test, that computation
was impractical, and part (d) was their way around it. The bounds \( d_L \) and \( d_U \) do not depend on \( \X \), so their
percentage points can be tabulated once for each \( n \) and \( p \). The **bounds test** of the hypothesis of no
autocorrelation against positive autocorrelation, at level \( \alpha \), rejects if \( d<d_{L,\alpha} \), accepts if
\( d>d_{U,\alpha} \), and is inconclusive in between. To test against negative autocorrelation, use \( 4-d \).

::: {#exm-het-okun-dw}
[The exact test for the Okun regression]

Under the hypothesis of independent normal errors, and given the growth series, \( d \) has mean \( 2.0031 \) and
standard deviation \( 0.1399 \), and its 5% point is \( 1.773 \). The bounds for \( n=202 \), \( p=2 \) are
\( d_{L,0.05}=1.760 \) and \( d_{U,0.05}=1.780 \), close together for a sample this large. The observed \( d=1.386 \) is far below
both, and the exact one-sided p-value is \( 3.3\times 10^{-6} \).

The bounds matter in small samples. For the \( 40 \) quarters from the second quarter of 1999 to the first quarter
of 2009, the same regression gives \( d=1.478 \), between the bounds \( 1.442 \) and \( 1.544 \): the bounds test is
inconclusive. The exact p-value, \( 0.037 \), settles it. [Figure 21.5.1](#fig-het-dw)(c) shows the three null
densities for this window.
:::

::: {when-format="html"}
![**Figure 21.5.1.** The Okun regression. (a) Residuals in time order. (b) Lag plot. (c) For the 40-quarter window
1999–2009: the exact null density of \( d \) given the regressors, between the densities of the bounds \( d_L \) and
\( d_U \); the vertical line marks the observed \( d \).](durbin_watson.svg){#fig-het-dw width=100%}
:::

::: {when-format="pdf"}
![The Okun regression. (a) Residuals in time order. (b) Lag plot. (c) For the 40-quarter window
1999–2009: the exact null density of \( d \) given the regressors, between the densities of the bounds \( d_L \) and
\( d_U \); the vertical line marks the observed \( d \).](durbin_watson.pdf){width=100%}
:::

```{.python .run #cell-durbin-watson-exact}
def dw_matrix(n):
    """A = D'D, where D takes first differences, so that d = e'Ae / e'e."""
    D = np.diff(np.eye(n), axis=0)
    return D.T @ D


def null_eigenvalues(X):
    """The n - p eigenvalues nu_k of the DW numerator on C(X)-perp."""
    n, p = X.shape
    Qfull, _ = np.linalg.qr(X, mode="complete")
    # orthonormal basis of C(X)-perp
    Qp = Qfull[:, p:]
    return np.linalg.eigvalsh(Qp.T @ dw_matrix(n) @ Qp)


def prob_negative(c):
    """Pr(sum_k c_k xi_k^2 < 0), xi_k iid standard normal (Imhof's formula)."""
    def integrand(u):
        theta = 0.5 * np.sum(np.arctan(c * u))
        rho = np.exp(0.25 * np.sum(np.log1p((c * u) ** 2)))
        return np.sin(theta) / (u * rho)
    return 0.5 - quad(integrand, 0, np.inf, limit=400)[0] / np.pi


nu = null_eigenvalues(X)
# Pr(d <= observed) under H0
p_value = prob_negative(nu - d)
m = len(nu)
mean_d = nu.mean()
var_d = 2 * (m * np.sum(nu**2) - np.sum(nu) ** 2) / (m**2 * (m + 2))
print(f"exact one-sided p-value {p_value:.2e}; "
      f"null mean {mean_d:.4f}, sd {np.sqrt(var_d):.4f}")
```

```{.python .run #cell-durbin-watson-window}
w = slice(160, 200)
Xw, yw = X[w], du[w]
ew = yw - Xw @ np.linalg.lstsq(Xw, yw, rcond=None)[0]
d_w = np.sum(np.diff(ew) ** 2) / np.sum(ew**2)
nw_, pw_ = Xw.shape
lam_w = 4 * np.sin(np.pi * np.arange(nw_) / (2 * nw_)) ** 2
dL_w = brentq(lambda t: prob_negative(lam_w[1:nw_ - pw_ + 1] - t) - 0.05,
              0.5, 2.0)
dU_w = brentq(lambda t: prob_negative(lam_w[pw_:] - t) - 0.05, 0.5, 2.0)
p_w = prob_negative(null_eigenvalues(Xw) - d_w)
print(f"1999Q2-2009Q1: d = {d_w:.3f}, bounds ({dL_w:.3f}, {dU_w:.3f}), "
      f"exact p = {p_w:.3f}")
```

Three limitations should be kept in mind. The exact theory needs normal errors and regressors fixed in the
sense of [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html): conditional on the regressors
the argument goes through as long as they are not functions of past responses. It fails when a lagged response
is among the regressors, and then \( d \) is biased towards \( 2 \), hiding autocorrelation; Durbin (1970) gave an
alternative for that case. The bounds need an intercept. And the statistic is aimed at lag-one dependence. It
is sensitive to other forms only insofar as they produce lag-one correlation, and it can miss, for example,
quarterly seasonal correlation at lag four.

## The Breusch–Godfrey test

A more flexible test comes from the score principle of [Section 21.2](02-detecting-heteroscedasticity.html).
Regress the residual \( \hat{\varepsilon}_t \) on the original regressors \( \x_{(t)} \) and on the lagged residuals
\( \hat{\varepsilon}_{t-1},\dots,\hat{\varepsilon}_{t-q} \), with lagged values before the first observation set to zero, and
let \( R^2 \) be the coefficient of determination. The statistic \( nR^2 \) is the score statistic for the hypothesis
of no autocorrelation against autoregressive, or moving-average, errors of order \( q \), and under the hypothesis
it is approximately \( \chi^2(q) \) in large samples (Breusch 1978; Godfrey 1978). Unlike \( d \), it remains valid
when lagged responses are among the regressors, and it looks at several lags at once. We state its null limit
without proof; the argument parallels that of @thm-het-breusch-pagan(c).

For the Okun regression with \( q=4 \), \( nR^2=22.43 \) on four degrees of freedom, with p-value
\( 1.6\times 10^{-4} \). Portmanteau statistics such as that of Ljung and Box (1978), based on
\( \sum_{k\le q}r_k^2 \), serve the same purpose in time-series modelling.

```{.python .run #cell-durbin-watson-bg}
lags = np.column_stack([np.concatenate([np.zeros(k), e[:-k]])
                        for k in range(1, 5)])
Z = np.column_stack([X, lags])
coef, *_ = np.linalg.lstsq(Z, e, rcond=None)
# n R^2
bg = n * (1 - np.sum((e - Z @ coef) ** 2) / np.sum((e - e.mean()) ** 2))
print(f"Breusch-Godfrey (4 lags) = {bg:.2f}, p = {stats.chi2.sf(bg, 4):.1e}")
```

## Exercises

### A. Check your understanding

::: {#exr-het-dw-range}
[A1]

Show that \( 0\le d<4 \) and that \( d=2(1-r_1)-(\hat{\varepsilon}_1^2+\hat{\varepsilon}_n^2)/\sum_t\hat{\varepsilon}_t^2 \). Which residual
patterns give values near \( 0 \) and near \( 4 \)?
:::

::: {#exr-het-dw-intercept}
[A2]

For the intercept-only model \( \X=\bone \), show that \( \E d=2 \) exactly under the hypothesis.
:::

::: {.solution}
By @thm-het-durbin-watson(b), \( \E d=\tr\{(\I-n^{-1}\bone\bone\T)\A\}/(n-1)=\{\tr\A-n^{-1}\bone\T\A\bone\}/(n-1) \). The diagonal of
\( \A \) sums to \( 2(n-1) \), and \( \A\bone=\bzero \). So \( \E d=2 \).
:::

### B. Practice

::: {#exr-het-dw-alternative}
[B1]

If \( \Cov(\be)=\sigma^2\V \), show that \( \E(\he\T\A\he)=\sigma^2\tr\{\A(\I-\M)\V(\I-\M)\} \). For \( n=3 \), \( \X=\bone \) and
\( v_{st}=\rho^{\lvert s-t\rvert} \), compute it and show that it decreases in \( \rho \).
:::

::: {#exr-het-dw-bounds-attained}
[B2]

Suppose \( \C(\X) \) is spanned by the eigenvectors \( \bv_1,\dots,\bv_p \) of \( \A \). Show that \( \nu_k=\lambda_{k+p} \) for all
\( k \), so that \( d=d_U \). Which columns of \( \X \) have this property, and what do they look like?
:::

::: {.solution}
Then \( \C(\X)\perpc \) is spanned by \( \bv_{p+1},\dots,\bv_n \), and \( \Q_\perp\T\A\Q_\perp \) has eigenvalues
\( \lambda_{p+1},\dots,\lambda_n \). The \( \bv_j \) are cosines of increasing frequency, starting with the constant \( \bv_1 \);
\( \bv_2 \) is a single slow half-wave, close to a centred linear trend. Regressors that are smooth functions of time
are close to this case, so \( d \) is close to \( d_U \) for trend regressions, and the upper bound is the relevant one.
:::

::: {#exr-het-dw-invariance}
[B3]

Show that \( d \) is unchanged if \( \y \) is replaced by \( c\y+\X\mathbf{a} \) for any \( c\ne0 \) and \( \mathbf{a} \). Use this to explain why its
null distribution cannot depend on \( \bbeta \) or \( \sigma \), without appeal to normality.
:::

### C. Going deeper

::: {#exr-het-dw-asymptotic}
[C1]

For the intercept-only model, show that under the hypothesis \( \Var d=4(n-2)/\{(n-1)(n+1)\} \) exactly, so that the
standard deviation of \( d \) is close to \( 2/\sqrt n \). (Use @thm-het-durbin-watson(b), \( \sum_j\lambda_j^2=\tr(\A^2) \), and
the entries of \( \A \).) Compare with the Okun regression of @exm-het-okun-dw.
:::

::: {.solution}
Here \( m=n-1 \) and the \( \nu_k \) are \( \lambda_2,\dots,\lambda_n \), so \( \sum_k\nu_k=\tr\A=2(n-1) \) and
\( \sum_k\nu_k^2=\tr(\A^2) \), the sum of the squared entries of \( \A \): \( 1+4(n-2)+1 \) from the diagonal and \( 2(n-1) \)
from the off-diagonal, in total \( 6n-8 \). Then \( m\sum_k\nu_k^2-(\sum_k\nu_k)^2=(n-1)(6n-8)-4(n-1)^2=2(n-1)(n-2) \),
and \( \Var d=4(n-1)(n-2)/\{(n-1)^2(n+1)\}=4(n-2)/\{(n-1)(n+1)\} \). For \( n=202 \), \( 2/\sqrt{202}=0.141 \), close to the
standard deviation \( 0.1399 \) computed for the Okun design, which has one more column.
:::

::: {#exr-het-dw-lagged}
[C2]

In the model \( y_t=\beta y_{t-1}+\varepsilon_t \) with independent errors, fitted by least squares on \( t=2,\dots,n \),
argue that the residuals are less autocorrelated than the errors would be under an autocorrelated alternative,
so that \( d \) is biased towards \( 2 \). Relate this to the requirement in @thm-het-durbin-watson that \( \X \) be fixed.
:::
