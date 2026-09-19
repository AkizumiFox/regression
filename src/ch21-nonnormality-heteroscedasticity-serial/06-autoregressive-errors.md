# Autoregressive errors

Serial correlation raises the same choice as heteroscedasticity: model it and use generalized least squares,
or keep least squares and estimate its covariance without a model. Either way dependence must be described by
few parameters, and the simplest description, the one behind the Durbin–Watson test, is the first-order
autoregression.

## The AR(1) model

::: {#def-het-ar1}
[First-order autoregressive errors]

The errors follow a **first-order autoregressive process**, AR(1), if
\[
\varepsilon_t=\rho\,\varepsilon_{t-1}+u_t,\qquad t=2,\dots,n,
\]
with \( \lvert\rho\rvert<1 \), where \( u_2,\dots,u_n \) are uncorrelated with mean zero and variance \( \sigma_u^2 \), and the
starting value \( \varepsilon_1 \) has mean zero, variance \( \sigma_u^2/(1-\rho^2) \) and is uncorrelated with the \( u_t \).
:::

The variance of \( \varepsilon_1 \) is the one that makes the process **stationary**: it is the variance an AR(1)
process reaches after running for a long time, and with it every \( \varepsilon_t \) has the same variance.

::: {#thm-het-ar1}
[Regression with AR(1) errors]

Let \( \Y=\X\bbeta+\be \) with \( \rank(\X)=p \) and AR(1) errors.

::: {.enumerate options="label=(\alph*)"}
1. \( \Cov(\varepsilon_s,\varepsilon_t)=\sigma_u^2\rho^{\lvert s-t\rvert}/(1-\rho^2) \). Write \( \Cov(\be)=\sigma_u^2\V_\rho \).

2. *(Prais–Winsten transformation.)* Let \( \bL_\rho \) be the \( n\times n \) lower bidiagonal matrix with
   \( (\bL_\rho)_{11}=\sqrt{1-\rho^2} \), \( (\bL_\rho)_{tt}=1 \) and \( (\bL_\rho)_{t,t-1}=-\rho \) for \( t\ge2 \). Then
   \( \Cov(\bL_\rho\be)=\sigma_u^2\I \), so \( \V_\rho^{-1}=\bL_\rho\T\bL_\rho \) and \( \det\V_\rho=1/(1-\rho^2) \). For known \( \rho \),
   the generalized least squares estimate is ordinary least squares on \( (\bL_\rho\y,\bL_\rho\X) \), and it is the best
   linear unbiased estimator.

3. *(Likelihood.)* If the errors are normal, the log-likelihood is
   \[
\begin{aligned}
\ell(\bbeta,\rho,\sigma_u^2)={}&-\frac n2\log(2\pi\sigma_u^2)+\frac12\log(1-\rho^2)\\
&-\frac{\norm{\bL_\rho(\y-\X\bbeta)}^2}{2\sigma_u^2}.
\end{aligned}
\]

4. *(Cochrane–Orcutt.)* Let \( S(\bbeta,\rho)=\sum_{t=2}^n\{(y_t-\rho y_{t-1})-(\x_{(t)}-\rho\x_{(t-1)})\T\bbeta\}^2 \), the sum of
   squares of the last \( n-1 \) rows of \( \bL_\rho(\y-\X\bbeta) \). Alternate the steps
   \[
\begin{gathered}
\rho_{k+1}=\frac{\sum_{t\ge2}e_te_{t-1}}{\sum_{t\ge2}e_{t-1}^2}\ \text{ with } \mathbf{e}=\y-\X\bbeta_k,\\
\bbeta_{k+1}=\argmin_{\bb}S(\bb,\rho_{k+1}).
\end{gathered}
\]
   Then \( S(\bbeta_{k+1},\rho_{k+1})\le S(\bbeta_k,\rho_{k+1})\le S(\bbeta_k,\rho_k) \), so \( S \) decreases along the iterations, and
   any limit point \( (\bbeta^*,\rho^*) \) with \( \lvert\rho^*\rvert<1 \) at which the steps are continuous is a stationary point of \( S \).
:::

:::

::: {.proof}
(a) By induction \( \Var\varepsilon_t=\sigma_u^2/(1-\rho^2) \) for all \( t \), because
\( \Var\varepsilon_t=\rho^2\Var\varepsilon_{t-1}+\sigma_u^2 \), \( \varepsilon_{t-1} \) being a linear combination of
\( \varepsilon_1,u_2,\dots,u_{t-1} \) and so uncorrelated with \( u_t \). For \( k\ge1 \), iterating the recursion gives
\( \varepsilon_t=\rho^k\varepsilon_{t-k}+\sum_{j=0}^{k-1}\rho^ju_{t-j} \), and the sum is uncorrelated with \( \varepsilon_{t-k} \), so
\( \Cov(\varepsilon_t,\varepsilon_{t-k})=\rho^k\Var\varepsilon_{t-k} \).

(b) \( \bL_\rho\be=(\sqrt{1-\rho^2}\,\varepsilon_1,u_2,\dots,u_n)\T \), whose entries are uncorrelated with variance \( \sigma_u^2 \).
So \( \sigma_u^2\bL_\rho\V_\rho\bL_\rho\T=\sigma_u^2\I \), and since \( \bL_\rho \) is nonsingular,
\( \V_\rho^{-1}=\bL_\rho\T\bL_\rho \). \( \bL_\rho \) is triangular with determinant \( \sqrt{1-\rho^2} \), so
\( \det\V_\rho=1/(1-\rho^2) \). The last statement is @thm-proj-A-projection and @cor-opt-aitken with the whitening
matrix \( \bL_\rho \), as in [Section 6.9](../ch06-projections/09-inner-products.html).

(c) The normal density of \( \be \) with covariance \( \sigma_u^2\V_\rho \) has
\( -\tfrac12\log\det(\sigma_u^2\V_\rho)=-\tfrac n2\log\sigma_u^2+\tfrac12\log(1-\rho^2) \) and quadratic form
\( (\y-\X\bbeta)\T\V_\rho^{-1}(\y-\X\bbeta)/\sigma_u^2=\norm{\bL_\rho(\y-\X\bbeta)}^2/\sigma_u^2 \).

(d) For fixed \( \bbeta \), \( S(\bbeta,\rho)=\sum_{t\ge2}(e_t-\rho e_{t-1})^2 \) is a quadratic in \( \rho \) minimized by
\( \rho_{k+1} \), which gives the second inequality. For fixed \( \rho \), \( S \) is a least squares criterion in \( \bb \),
minimized by \( \bbeta_{k+1} \), which gives the first. At a limit point both steps leave \( (\bbeta^*,\rho^*) \) unchanged,
so \( \partial S/\partial\rho=0 \) and \( \partial S/\partial\bbeta=\bzero \) there.
:::

The whitening in (b) is the transformation of Prais and Winsten (1954). Row \( t\ge2 \) of the transformed model is
the **quasi-difference** \( y_t-\rho y_{t-1}=(\x_{(t)}-\rho\x_{(t-1)})\T\bbeta+u_t \), and row one is the first observation
rescaled by \( \sqrt{1-\rho^2} \). The inverse covariance \( \V_\rho^{-1} \) is tridiagonal, with diagonal
\( (1,1+\rho^2,\dots,1+\rho^2,1) \) and off-diagonal entries \( -\rho \) (@exr-het-ar1-inverse), which is why AR(1) models are
cheap to fit however long the series. Cochrane and Orcutt (1949) dropped the first row, which leaves the
conditional sum of squares \( S \) of (d). Dropping it loses little when the regressors are stationary, but it can
lose a great deal when they trend, because the first observation is then the only one whose transformed regressor
is not nearly constant. Maximizing the likelihood of (c) keeps the first row and the term
\( \tfrac12\log(1-\rho^2) \), which keeps \( \hat\rho \) inside \( (-1,1) \).

## Feasible generalized least squares

In practice \( \rho \) is unknown. The two-step **feasible** estimator estimates it from the ordinary least squares
residuals, \( \hat\rho=\sum_{t\ge2}\hat{\varepsilon}_t\hat{\varepsilon}_{t-1}/\sum_{t\ge2}\hat{\varepsilon}_{t-1}^2 \), and then applies the
Prais–Winsten transformation with \( \hat\rho \). Iterating the two steps gives either the Cochrane–Orcutt estimate
or, keeping the first row, a close relative of the maximum likelihood estimate.

::: {.remark}
[What is known about feasible GLS with AR(1) errors]

Suppose the errors are a stationary AR(1) process with independent innovations of finite variance, and the
regressors are stationary, or have polynomial trends, with suitable limits for their sample moments. Then
\( \hat\rho\to\rho \) in probability, and \( \sqrt n(\hbeta_{\text{FGLS}}-\hbeta_{\text{GLS}})\to\bzero \) after the appropriate
normalization for trending regressors, so the feasible estimator has the limiting distribution of the
optimal one and the usual standard errors from the transformed regression are valid in large samples. We do not
prove this. A sketch: \( \hat\rho \) differs from the sample autocorrelation of the true errors by terms bounded by
multiples of \( \norm{\M\be}\,\norm{\be}/n \) and \( \norm{\M\be}^2/n \), which tend to zero because \( \norm{\M\be} \) is bounded in
probability (compare the proof of @thm-het-breusch-pagan(c)), and that autocorrelation is consistent by a law of
large numbers for stationary sequences; the comparison of \( \hbeta_{\text{FGLS}} \) with \( \hbeta_{\text{GLS}} \) then follows the pattern of
@thm-het-wls(d), because \( \bL_\rho \) depends smoothly on \( \rho \). Fuller (1996) gives complete statements and
proofs. The finite-sample behaviour must be checked by simulation.
:::

::: {#exm-het-okun-pw}
[Okun's law with AR(1) errors]

For the regression of @exm-het-okun, \( \hat\rho=0.303 \). The two-step Prais–Winsten estimate of the slope is
\( -0.0506 \) with standard error \( 0.0049 \), against the least squares \( -0.0671 \) with nominal standard error
\( 0.0050 \). Iterating Cochrane–Orcutt to convergence takes \( 17 \) iterations and ends at \( \hat\rho=0.548 \) and slope
\( -0.0380 \). The estimate of the effect of growth on unemployment thus changes by almost half depending on the
treatment of the errors, and the residuals of the transformed regression still have lag-one autocorrelation
\( 0.148 \). Both facts suggest that the AR(1) error model does not describe these data. We return to this
below.
:::

```{.python .run #cell-ar1-errors-pw}
import numpy as np
import statsmodels.api as sm
from scipy import stats

macro = sm.datasets.macrodata.load_pandas().data
du = np.diff(macro["unemp"].to_numpy())
growth = 400 * np.diff(np.log(macro["realgdp"].to_numpy()))
n = len(du)
X = np.column_stack([np.ones(n), growth])


def prais_winsten(X, y, rho):
    """Whiten AR(1) errors: row 1 times sqrt(1 - rho^2), row t - rho row t-1."""
    Xs, ys = X.copy(), y.copy()
    Xs[1:], ys[1:] = X[1:] - rho * X[:-1], y[1:] - rho * y[:-1]
    Xs[0], ys[0] = np.sqrt(1 - rho**2) * X[0], np.sqrt(1 - rho**2) * y[0]
    return Xs, ys


def ols(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ b
    cov = (e @ e / (len(y) - X.shape[1])) * np.linalg.inv(X.T @ X)
    return b, e, cov


b_ols, e, V_ols = ols(X, du)
# regress e_t on e_(t-1)
rho = np.sum(e[1:] * e[:-1]) / np.sum(e[:-1] ** 2)
# two-step Prais-Winsten
b_pw, e_pw, V_pw = ols(*prais_winsten(X, du, rho))
print(f"rho-hat {rho:.3f}")
print(f"OLS slope {b_ols[1]:.4f} (se {np.sqrt(V_ols[1, 1]):.4f}); "
      f"Prais-Winsten {b_pw[1]:.4f} (se {np.sqrt(V_pw[1, 1]):.4f})")
```

```{.python .run #cell-ar1-errors-co}
def cond_ss(b, r):
    """Conditional sum of squares S(beta, rho) over t = 2, ..., n."""
    u = (du[1:] - r * du[:-1]) - (X[1:] - r * X[:-1]) @ b
    return u @ u


b_co, r_co, path = b_ols, 0.0, []
for it in range(100):
    res = du - X @ b_co
    r_new = np.sum(res[1:] * res[:-1]) / np.sum(res[:-1] ** 2)   # rho step
    path.append(cond_ss(b_co, r_new))
    # beta step
    b_co, *_ = np.linalg.lstsq(X[1:] - r_new * X[:-1], du[1:] - r_new * du[:-1],
                               rcond=None)
    path.append(cond_ss(b_co, r_new))
    if abs(r_new - r_co) < 1e-10:
        break
    r_co = r_new
print(f"Cochrane-Orcutt after {it + 1} iterations: rho {r_co:.4f}, "
      f"slope {b_co[1]:.4f}")
```

## Efficiency and validity of least squares

How much is lost by ignoring AR(1) errors? There are two separate questions: how much less efficient least squares
is than generalized least squares, and how wrong its nominal standard error \( s\{(\X\T\X)^{-1}\}_{jj}^{1/2} \) is. The
answers depend on the regressor, and they run in opposite directions.

::: {#exm-het-ar1-efficiency}
[Three regressors]

Take \( n=100 \), an intercept, and one of three regressors: a linear trend, a fixed draw of white noise, and a
fixed draw of an AR(1) series with coefficient \( 0.8 \) (its sample lag-one autocorrelation is \( 0.87 \)). For AR(1)
errors, [Figure 21.6.1](#fig-het-ar1) shows (a) the ratio of the generalized to the ordinary least squares variance of
the slope, and (b) the ratio of the true least squares variance to the expected nominal one,
\( \E(s^2)\{(\X\T\X)^{-1}\}_{22} \), with \( \E s^2 \) from @thm-rv-quadform-mean. All are exact computations.

At \( \rho=0.6 \): for the trend, least squares is \( 0.970 \) efficient but its nominal variance is too small by a
factor of \( 4.01 \). For white noise it is only \( 0.462 \) efficient, yet its nominal variance is nearly right (factor
\( 1.11 \)). The autocorrelated regressor suffers both ways: efficiency \( 0.650 \), variance understated by a factor
of \( 3.19 \). At \( \rho=0.9 \) the factors for the trend become \( 0.842 \) and \( 19.59 \).

The trend repeats the lesson of the autocorrelated-errors example of
[Section 19.2](../ch19-theory-of-departures/02-wrong-covariance.html); the other two regressors are new.
:::

::: {when-format="html"}
![**Figure 21.6.1.** Least squares with AR(1) errors, \( n=100 \), for three regressors. (a) Efficiency: variance of the
generalized least squares slope divided by that of the ordinary one. (b) True variance of the least squares slope
divided by its expected nominal variance (log scale).](ar1_efficiency.svg){#fig-het-ar1 width=100%}
:::

::: {when-format="pdf"}
![Least squares with AR(1) errors, \( n=100 \), for three regressors. (a) Efficiency: variance of the
generalized least squares slope divided by that of the ordinary one. (b) True variance of the least squares slope
divided by its expected nominal variance (log scale).](ar1_efficiency.pdf){width=100%}
:::

The pattern has a simple explanation. The variance of \( \hat\beta_1 \) in a simple regression is
\( \sum_{s,t}(x_s-\bar x)(x_t-\bar x)\Cov(\varepsilon_s,\varepsilon_t)/S_{xx}^2 \). When the regressor is itself autocorrelated,
neighbouring terms add up, and for a stationary regressor with lag-one autocorrelation \( \phi \) the ratio of true to
nominal variance is approximately \( (1+\rho\phi)/(1-\rho\phi) \) (@exr-het-ar1-inflation). For the AR(0.8) regressor
this gives \( 3.18 \) at \( \rho=0.6 \), against the exact \( 3.19 \). A white-noise regressor has \( \phi\approx0 \), so the nominal
variance is nearly right; but generalized least squares, which quasi-differences \( x_t \), turns a rough regressor into
an even rougher one with more information, and least squares' efficiency tends to \( (1-\rho^2)/(1+\rho^2) \), which is
\( 0.471 \) at \( \rho=0.6 \) (@exr-het-ar1-white). For polynomial trends, Grenander (1954) showed that least squares is
asymptotically fully efficient, which the trend curve in panel (a) approaches; but its standard error is then at its
most misleading.

## Newey–West standard errors

The analogue of the sandwich of [Section 21.4](04-sandwich-estimators.html) for serially correlated errors keeps
least squares and estimates the middle matrix \( \X\T\Cov(\be)\X=\sum_{s,t}\Cov(\varepsilon_s,\varepsilon_t)\x_{(s)}\x_{(t)}\T \)
without a parametric model. The products \( \hat{\varepsilon}_s\hat{\varepsilon}_t \) can stand in for the covariances only
when \( s \) and \( t \) are close, since distant covariances are small but their estimates are not. With \( \mathbf{v}_t=\x_{(t)}\hat{\varepsilon}_t \)
and a truncation lag \( L \), Newey and West (1987) proposed
\[
\begin{gathered}
\hat{\boldsymbol{\Omega}}_L=\sum_{k=-L}^{L}\Bigl(1-\frac{\lvert k\rvert}{L+1}\Bigr)\hat{\boldsymbol{\Gamma}}_k,\\
\hat{\boldsymbol{\Gamma}}_k=\sum_{t=k+1}^n\mathbf{v}_t\mathbf{v}_{t-k}\T,\quad \hat{\boldsymbol{\Gamma}}_{-k}=\hat{\boldsymbol{\Gamma}}_k\T ,
\end{gathered}
\]{#eq-het-newey-west}

and the covariance estimate \( (\X\T\X)^{-1}\hat{\boldsymbol{\Omega}}_L(\X\T\X)^{-1} \). The declining **Bartlett weights**
\( 1-\lvert k\rvert/(L+1) \) are not a matter of taste.

::: {#prp-het-newey-west}
[Newey–West estimates are nonnegative definite]

For any vectors \( \mathbf{v}_1,\dots,\mathbf{v}_n \) and any \( L\ge0 \), \( \hat{\boldsymbol{\Omega}}_L \) in @eq-het-newey-west is nonnegative
definite. Without the weights, \( \sum_{\lvert k\rvert\le L}\hat{\boldsymbol{\Gamma}}_k \) need not be.
:::

::: {.proof}
Set \( \mathbf{v}_t=\bzero \) for \( t\le0 \) and \( t>n \), and let \( \mathbf{s}_j=\sum_{i=0}^{L}\mathbf{v}_{j-i} \) for every integer \( j \). Then
\[
\sum_j\mathbf{s}_j\mathbf{s}_j\T=\sum_j\sum_{i,l=0}^{L}\mathbf{v}_{j-i}\mathbf{v}_{j-l}\T=\sum_{k=-L}^{L}(L+1-\lvert k\rvert)\sum_t\mathbf{v}_t\mathbf{v}_{t-k}\T ,
\]
because for fixed \( t=j-i \) and lag \( k=l-i \) there are exactly \( L+1-\lvert k\rvert \) pairs \( (i,l) \) in \( \{0,\dots,L\}^2 \)
with \( l-i=k \). Dividing by \( L+1 \) shows that \( \hat{\boldsymbol{\Omega}}_L=(L+1)^{-1}\sum_j\mathbf{s}_j\mathbf{s}_j\T \), a sum of nonnegative
definite matrices. @exr-het-hac-negative gives an example for the unweighted sum.
:::

Newey and West (1987) proved that the estimator is consistent when \( L \) grows with \( n \) slowly enough, for errors
and regressors that are stationary with dependence dying out fast enough. We do not prove this. Its practical
weakness is the choice of \( L \), which trades bias from truncation against variance from estimating many
autocovariances; Andrews (1991) studied data-based choices. A common rule of thumb is
\( L=\lfloor4(n/100)^{2/9}\rfloor \). For the Okun regression, with \( L=4 \), the Newey–West standard error of the
slope is \( 0.0086 \), against the nominal \( 0.0050 \). With \( L=0 \) the estimator is HC0 (@exr-het-hac-hc0).

```{.python .run #cell-ar1-errors-nw}
def newey_west(X, e, L):
    """Bartlett-weighted HAC covariance of the OLS estimate, truncation lag L."""
    v = X * e[:, None]
    S = v.T @ v
    for k in range(1, L + 1):
        G = v[k:].T @ v[:-k]
        S += (1 - k / (L + 1)) * (G + G.T)
    B = np.linalg.inv(X.T @ X)
    return B @ S @ B


V_nw = newey_west(X, e, 4)
print(f"Newey-West (L = 4) se of slope {np.sqrt(V_nw[1, 1]):.4f}")
```

::: {#exm-het-ar1-coverage}
[Coverage with AR(1) errors]

Let the regressor be a fixed draw of an AR(1) series with coefficient \( 0.8 \), and the errors AR(1) with \( \rho=0.6 \)
and normal innovations. Over \( 10{,}000 \) simulated data sets the nominal 95% intervals for the slope covered as
follows. For \( n=50 \): least squares with the usual standard error \( 0.747 \), least squares with Newey–West
\( 0.784 \), two-step Prais–Winsten \( 0.913 \), and generalized least squares with the true \( \rho \) \( 0.946 \). For
\( n=200 \): \( 0.746 \), \( 0.871 \), \( 0.944 \) and \( 0.950 \). The Newey–West interval is valid only asymptotically, and
converges slowly; when the AR(1) model is right, using it pays.
:::

```{.python .run #cell-ar1-errors-coverage}
def ar1_draws(n, rho, reps, rng):
    u = rng.normal(size=(reps, n))
    x = np.empty((reps, n))
    x[:, 0] = u[:, 0] / np.sqrt(1 - rho**2)
    for t in range(1, n):
        x[:, t] = rho * x[:, t - 1] + u[:, t]
    return x


def coverage(n, reps, rng):
    """Coverage of 95% slope intervals: OLS, Newey-West, Prais-Winsten, GLS."""
    x = ar1_draws(n, 0.8, 1, rng)[0]
    Xc = np.column_stack([np.ones(n), x])
    L = int(4 * (n / 100) ** (2 / 9))
    q = stats.t.ppf(0.975, n - 2)
    hits = np.zeros(4)
    for y in ar1_draws(n, 0.6, reps, rng):          # true slope 0
        b, res, V = ols(Xc, y)
        r = np.sum(res[1:] * res[:-1]) / np.sum(res[:-1] ** 2)
        b_f, _, V_f = ols(*prais_winsten(Xc, y, r))
        b_g, _, V_g = ols(*prais_winsten(Xc, y, 0.6))
        ses = [V[1, 1], newey_west(Xc, res, L)[1, 1], V_f[1, 1], V_g[1, 1]]
        hits += [abs(bb[1]) <= q * np.sqrt(v)
                 for bb, v in zip((b, b, b_f, b_g), ses)]
    return hits / reps


rng = np.random.default_rng(2109)
for n_sim in (50, 200):
    print(n_sim, "OLS, Newey-West, Prais-Winsten, GLS:",
          np.round(coverage(n_sim, 500, rng), 3))
```

## Serial correlation as a symptom

An AR(1) error model is a restriction on a dynamic regression. Multiply the model at time \( t-1 \) by \( \rho \) and
subtract:
\[
y_t=\rho y_{t-1}+\x_{(t)}\T\bbeta-\rho\,\x_{(t-1)}\T\bbeta+u_t .
\]
This is a regression of \( y_t \) on its own lag and on current and lagged regressors, with the coefficients of the
lagged regressors constrained to equal \( -\rho \) times those of the current ones, the **common factor restriction**.
Residual autocorrelation is often the trace of a missing lag of \( y \) or of a regressor, and then the dynamic
regression without the restriction is the better model. Hendry and Mizon (1978) argued that the restriction should
be tested before the error model is adopted.

::: {#exm-het-okun-dynamic}
[Okun's law, dynamically]

Regressing \( \Delta u_t \) on an intercept, \( \Delta u_{t-1} \), \( g_t \) and \( g_{t-1} \) gives coefficients \( 0.2160 \), \( 0.3571 \)
(standard error \( 0.0574 \)), \( -0.0487 \) (\( 0.0043 \)) and \( -0.0162 \) (\( 0.0055 \)). The AR(1) model would require the last
coefficient to be \( -0.3571\times(-0.0487)=0.0174 \). The estimate has the opposite sign, and the delta-method test of
the restriction gives \( z=-8.1 \). The unrestricted dynamic regression has no remaining residual autocorrelation: its
Breusch–Godfrey statistic with four lags is \( 4.42 \), with p-value \( 0.35 \). Unemployment responds to growth with a lag
and with persistence; the AR(1) error model mistook this for correlated noise, and the slope estimates of
@exm-het-okun-pw, which move with \( \hat\rho \), were artefacts of that mistake.
:::

```{.python .run #cell-ar1-errors-dynamic}
# the AR(1) error model is the dynamic regression
#   du_t = a + rho du_(t-1) + b g_t + c g_(t-1) + u_t
#   with the restriction c = -rho b
Xd = np.column_stack([np.ones(n - 1), du[:-1], growth[1:], growth[:-1]])
bd, ed, Vd = ols(Xd, du[1:])
# zero under the restriction
restr = bd[3] + bd[1] * bd[2]
grad = np.array([0.0, bd[2], bd[1], 1.0])                      # delta method
z_cf = restr / np.sqrt(grad @ Vd @ grad)
print("dynamic regression: const, du(t-1), g(t), g(t-1) =", np.round(bd, 4))
print(f"common-factor restriction c + rho b = {restr:.4f}, z = {z_cf:.1f}")
```

An extreme case is the regression of one random walk on another, independent one, which gives large \( t \)
statistics and a small \( d \) far more often than the nominal level suggests (Granger and Newbold 1974). There the
residual autocorrelation signals a meaningless regression, and no correction of the standard errors repairs it.

## Exercises

### A. Check your understanding

::: {#exr-het-ar1-inverse}
[A1]

Compute \( \bL_\rho\T\bL_\rho \) and show that it is tridiagonal with diagonal \( (1,1+\rho^2,\dots,1+\rho^2,1) \) and off-diagonal
entries \( -\rho \). Check \( \V_\rho\bL_\rho\T\bL_\rho=\I \) for \( n=2 \).
:::

::: {.solution}
Column \( t \) of \( \bL_\rho \) has entry \( 1 \) (or \( \sqrt{1-\rho^2} \) for \( t=1 \)) in row \( t \) and \( -\rho \) in row \( t+1 \) (for \( t<n \)).
The inner product of columns \( t \) and \( t \) is \( 1+\rho^2 \) for \( 1<t<n \), \( (1-\rho^2)+\rho^2=1 \) for \( t=1 \) and \( 1 \) for \( t=n \);
columns \( t \) and \( t+1 \) share only row \( t+1 \), giving \( -\rho \); other pairs share nothing. For \( n=2 \),
\( \V_\rho=(1-\rho^2)^{-1}\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix} \) and
\( \bL_\rho\T\bL_\rho=\begin{pmatrix}1&-\rho\\-\rho&1\end{pmatrix} \), whose product is \( (1-\rho^2)^{-1}(1-\rho^2)\I=\I \).
:::

::: {#exr-het-hac-hc0}
[A2]

Show that the Newey–West estimator with \( L=0 \) is HC0, and that with \( \X=\bone \) and \( L\ge1 \) it estimates the long-run
variance \( \sum_k\gamma_k \) of a stationary series, scaled by \( n \).
:::

### B. Practice

::: {#exr-het-ar1-white}
[B1]

Let the regressor \( x_t \) be white noise with variance \( \tau^2 \), independent of AR(1) errors, and consider the slope in
a regression through the origin. Show heuristically that
\( n\Var(\hat\beta_{\text{OLS}})\to\sigma_u^2/\{(1-\rho^2)\tau^2\} \) and
\( n\Var(\hat\beta_{\text{GLS}})\to\sigma_u^2/\{(1+\rho^2)\tau^2\} \), so that the efficiency of least squares, the ratio bounded in @thm-dep-efficiency, tends to
\( (1-\rho^2)/(1+\rho^2) \).
:::

::: {#exr-het-ar1-inflation}
[B2]

Let the regressor be a stationary AR(1) series with coefficient \( \phi \), independent of AR(1) errors with coefficient
\( \rho \). Show heuristically that the ratio of the true variance of the least squares slope to its nominal variance
tends to \( (1+\rho\phi)/(1-\rho\phi) \).
:::

::: {.solution}
Conditionally on the regressor,
\( \Var(\hat\beta)=\sum_{s,t}\tilde x_s\tilde x_t\gamma_0\rho^{\lvert s-t\rvert}/S_{xx}^2 \), where \( \tilde x_t=x_t-\bar x \). For large
\( n \), \( \sum_t\tilde x_t\tilde x_{t-k}\approx S_{xx}\phi^{\lvert k\rvert} \), so the variance is about
\( (\gamma_0/S_{xx})\sum_k(\rho\phi)^{\lvert k\rvert}=(\gamma_0/S_{xx})(1+\rho\phi)/(1-\rho\phi) \). The nominal variance
\( \E(s^2)/S_{xx} \) is about \( \gamma_0/S_{xx} \), since \( s^2 \) estimates the marginal error variance \( \gamma_0 \) when
\( p/n \) is small.
:::

::: {#exr-het-co-first}
[B3]

In the model \( y_t=\alpha+\beta t+\varepsilon_t \), \( t=1,\dots,n \), with AR(1) errors, show that the Cochrane–Orcutt
regression has columns \( (1-\rho)\bone \) and \( (1-\rho)t+\rho \) on \( t=2,\dots,n \), and that the information it
carries about \( \beta \) is \( (1-\rho)^2\sum_{t\ge2}(t-\bar t)^2/\sigma_u^2 \), where \( \bar t \) is the mean of \( 2,\dots,n \). What
happens as \( \rho\to1 \) with \( n \) fixed, and why does the Prais–Winsten first row change the picture?
:::

### C. Going deeper

::: {#exr-het-hac-negative}
[C1]

Give vectors \( v_1,v_2,v_3 \) (with \( p=1 \)) for which the unweighted sum \( \sum_{\lvert k\rvert\le1}\hat{\Gamma}_k \) is negative.
:::

::: {.solution}
Take \( (v_1,v_2,v_3)=(1,-1,1) \). Then \( \hat\Gamma_0=3 \) and \( \hat\Gamma_1=\hat\Gamma_{-1}=v_2v_1+v_3v_2=-2 \), so the unweighted sum
is \( 3-4=-1 \). With Bartlett weights it is \( 3-2=1 \).
:::
