# The Box–Cox family

Choosing a rung of the ladder by eye works when the evidence is strong, but it gives no measure of
uncertainty and no way to compare rungs that fit almost equally well. Box and Cox (1964) turned the
choice into a parameter: embed the powers in a family indexed by a real number \( \lambda \), assume that
for some \( \lambda \) the normal linear model holds on the transformed scale, and estimate \( \lambda \) by
maximum likelihood along with everything else.

## The family

::: {#def-tr-box-cox}
[The Box–Cox family]

For \( y>0 \) and \( \lambda\in\Real \), the **Box–Cox transformation** of \( y \) is
\[
y^{(\lambda)}=\begin{cases}\dfrac{y^{\lambda}-1}{\lambda},&\lambda\ne0,\\[1.2ex] \log y,&\lambda=0.\end{cases}
\]{#eq-tr-box-cox}

For a sample \( y_1,\dots,y_n>0 \) with geometric mean \( \dot{y}=\exp\bigl(n^{-1}\sum_i\log y_i\bigr) \), the
**normalized** transformation is \( z_i^{(\lambda)}=y_i^{(\lambda)}/\dot{y}^{\,\lambda-1} \), and
\( \mathbf{z}^{(\lambda)}=(z_1^{(\lambda)},\dots,z_n^{(\lambda)})\T \).
:::

In a model with an intercept \( y^{(\lambda)} \) and \( y^{\lambda} \) give the same fit. Subtracting \( 1 \) and dividing by
\( \lambda \) make the family *continuous* in \( \lambda \), with the logarithm in its place, and increasing for every \( \lambda \).

::: {#prp-tr-family}
[Properties of the family]

Let \( y>0 \).

::: {.enumerate options="label=(\alph*)"}
1. \( y^{(\lambda)}\to\log y \) as \( \lambda\to0 \), so \( \lambda\mapsto y^{(\lambda)} \) is continuous.

2. For each \( \lambda \), \( y\mapsto y^{(\lambda)} \) is strictly increasing on \( (0,\infty) \), with derivative
   \( y^{\lambda-1} \). The normalized transformation has derivative \( (y/\dot{y})^{\lambda-1} \), which equals \( 1 \) at
   \( y=\dot{y} \) for every \( \lambda \).

3. If \( p<q \), then \( y^{(q)}=\varphi\bigl(y^{(p)}\bigr) \) for a strictly increasing and strictly convex function
   \( \varphi \). Equivalently, moving down the ladder from \( q \) to \( p \) applies a strictly concave increasing function.

4. As \( y \) ranges over \( (0,\infty) \), \( y^{(\lambda)} \) ranges over \( (-1/\lambda,\infty) \) if \( \lambda>0 \), over
   \( (-\infty,-1/\lambda) \) if \( \lambda<0 \), and over \( \Real \) if \( \lambda=0 \).
:::

:::

::: {.proof}
(a) \( y^{\lambda}=\exp(\lambda\log y)=1+\lambda\log y+O(\lambda^2) \), so \( (y^{\lambda}-1)/\lambda=\log y+O(\lambda) \).
(b) For \( \lambda\ne0 \) the derivative of \( (y^{\lambda}-1)/\lambda \) is \( y^{\lambda-1}>0 \), and for \( \lambda=0 \) the
derivative of \( \log y \) is \( y^{-1} \). Dividing by the constant \( \dot{y}^{\,\lambda-1} \) gives the second statement.
(c) Let \( u=y^{(p)} \). By (b) the map \( y\mapsto u \) has an increasing inverse \( y=\psi(u) \) with
\( \psi'(u)=y^{1-p} \). Then \( \varphi=(\cdot)^{(q)}\circ\psi \) has derivative \( \varphi'(u)=y^{q-1}y^{1-p}=\psi(u)^{q-p} \). Since
\( q-p>0 \) and \( \psi \) is strictly increasing, \( \varphi' \) is positive and strictly increasing, so \( \varphi \) is strictly
increasing and strictly convex. The inverse of a strictly increasing convex function is strictly increasing and
concave.
(d) \( y^{\lambda} \) takes every value in \( (0,\infty) \). For \( \lambda>0 \) subtracting \( 1 \) and dividing by \( \lambda \) gives
\( (-1/\lambda,\infty) \), and for \( \lambda<0 \) the division reverses the interval to \( (-\infty,-1/\lambda) \).
:::

Part (c) is the ladder's advice made precise: a plot of \( y^{(q)} \) that bends upwards is a convex function of
\( y^{(p)} \), and a lower power removes some of the bend. Part (d) is a warning. For \( \lambda\ne0 \) the transformed
response is bounded on one side, so it *cannot* be exactly normal. The Box–Cox model is always an
approximation, harmless when the bound is many standard deviations from the data.

## The likelihood

The working model is that for some \( \lambda \) the \( Y_i^{(\lambda)} \) are independent with
\( Y_i^{(\lambda)}\sim\Normal(\x_{(i)}\T\bbeta,\sigma^2) \), ignoring the small normal probability beyond the bound of
@prp-tr-family(d). The likelihood must be written for the *observed* \( y_i \), since otherwise changing \( \lambda \)
changes the data being explained. That is where the Jacobian enters, and the geometric mean with it.

::: {#thm-tr-box-cox}
[Profile likelihood of the Box–Cox parameter]

Let \( \X \) have rank \( r<n \), let \( y_1,\dots,y_n>0 \), and for each \( \lambda \) let
\( \text{SSE}_z(\lambda)=\norm{(\I-\M)\mathbf{z}^{(\lambda)}}^2 \), assumed positive. Under the model above:

::: {.enumerate options="label=(\alph*)"}
1. the log-likelihood of \( (\bbeta,\sigma^2,\lambda) \) is
   \[
\begin{aligned}
\ell(\bbeta,\sigma^2,\lambda)={}&-\frac n2\log(2\pi\sigma^2)-\frac{\norm{\y^{(\lambda)}-\X\bbeta}^2}{2\sigma^2}\\
&+(\lambda-1)\sum_{i=1}^n\log y_i ;
\end{aligned}
\]{#eq-tr-box-cox-loglik}

2. for fixed \( \lambda \) it is maximized by the least squares fit of \( \y^{(\lambda)} \) on \( \X \), and the profile
   log-likelihood of \( \lambda \) is
   \[
\ell_p(\lambda)=-\frac n2\Bigl\{\log\Bigl(\frac{2\pi\,\text{SSE}_z(\lambda)}{n}\Bigr)+1\Bigr\} ,
\]{#eq-tr-profile}

   so the maximum likelihood estimate \( \hat{\lambda} \) is the value that minimizes \( \text{SSE}_z(\lambda) \);

3. the likelihood ratio statistic for \( H:\lambda=\lambda_0 \) is
   \( 2\{\ell_p(\hat{\lambda})-\ell_p(\lambda_0)\}=n\log\{\text{SSE}_z(\lambda_0)/\text{SSE}_z(\hat{\lambda})\} \), and the
   likelihood ratio confidence set of level \( 1-\alpha \) is
   \[
\Bigl\{\lambda_0:\ n\log\frac{\text{SSE}_z(\lambda_0)}{\text{SSE}_z(\hat{\lambda})}\le\chi^2_{1,\alpha}\Bigr\},
\]{#eq-tr-lr-interval}

   where \( \chi^2_{1,\alpha} \) is the upper \( \alpha \) point of \( \chi^2(1) \).
:::

:::

::: {.proof}
(a) The map \( y\mapsto y^{(\lambda)} \) is strictly increasing and differentiable with derivative \( y^{\lambda-1} \) (@prp-tr-family(b)). By the change-of-variables formula the density of \( Y_i \) at \( y_i \) is the
\( \Normal(\x_{(i)}\T\bbeta,\sigma^2) \) density at \( y_i^{(\lambda)} \) times \( y_i^{\lambda-1} \). The \( Y_i \) are independent, so
the log-likelihood is the sum of the logs, which is @eq-tr-box-cox-loglik.
(b) The Jacobian term does not involve \( (\bbeta,\sigma^2) \), and the other two terms are the log-likelihood of the
normal linear model with response \( \y^{(\lambda)} \). By @thm-opt-mle its maximum over \( (\bbeta,\sigma^2) \) is
attained at any least squares fit and equals \( -\frac n2\{\log(2\pi\,\text{SSE}(\y^{(\lambda)})/n)+1\} \). Now
\( \y^{(\lambda)}=\dot{y}^{\,\lambda-1}\mathbf{z}^{(\lambda)} \), so
\( \text{SSE}(\y^{(\lambda)})=\dot{y}^{\,2(\lambda-1)}\text{SSE}_z(\lambda) \), and
\[
-\frac n2\log\dot{y}^{\,2(\lambda-1)}=-n(\lambda-1)\log\dot{y}=-(\lambda-1)\sum_{i=1}^n\log y_i .
\]
This cancels the Jacobian term exactly and leaves @eq-tr-profile, which is a decreasing function of
\( \text{SSE}_z(\lambda) \).
(c) Subtract two values of @eq-tr-profile.
:::

The theorem is algebraic. The *coverage* of @eq-tr-lr-interval is a separate, asymptotic matter: if the
model holds at \( \lambda_0 \) and the regularity conditions of Wilks's theorem are satisfied, the probability
that the set contains \( \lambda_0 \) tends to \( 1-\alpha \) (van der Vaart 1998, chapter 16). That result is not proved
here, and by @prp-tr-family(d) its hypothesis can hold exactly only at \( \lambda_0=0 \). Otherwise, as Hernandez and
Johnson (1980) showed, \( \hat{\lambda} \) converges to the value whose normal model is closest in Kullback–Leibler
divergence to the true distribution, and the interval describes the precision of that best approximation.

::: {.remark}
[Why normalize]

Residual sums of squares of \( \y^{(\lambda)} \) are in the units of \( y^{\lambda} \) and cannot be compared across
\( \lambda \). The geometric mean builds the Jacobian into the response: every \( \mathbf{z}^{(\lambda)} \) is in the units of
\( y \), with unit slope at \( \dot{y} \) (@prp-tr-family(b)), so \( s_z(\lambda)=\{\text{SSE}_z(\lambda)/(n-r)\}^{1/2} \) can be
compared across \( \lambda \), and minimizing it is maximizing the likelihood.
:::

The normalization also makes the analysis independent of the units in which \( y \) is recorded.

::: {#prp-tr-units}
[Units of the response]

Suppose \( \bone\in\C(\X) \). If every \( y_i \) is multiplied by the same \( c>0 \), then \( \text{SSE}_z(\lambda) \) is
multiplied by \( c^2 \) for every \( \lambda \). Consequently \( \ell_p \) changes by a constant, and \( \hat{\lambda} \), the
likelihood ratio statistics and the set @eq-tr-lr-interval are unchanged.
:::

::: {.proof}
The geometric mean of \( cy_1,\dots,cy_n \) is \( c\dot{y} \). For \( \lambda\ne0 \),
\[
\frac{(cy_i)^{\lambda}-1}{\lambda(c\dot{y})^{\lambda-1}}
=c\,\frac{y_i^{\lambda}-1}{\lambda\dot{y}^{\,\lambda-1}}+c\,\frac{1-c^{-\lambda}}{\lambda\dot{y}^{\,\lambda-1}},
\]
and for \( \lambda=0 \), \( c\dot{y}\log(cy_i)=c\,\dot{y}\log y_i+c\,\dot{y}\log c \). In both cases the new normalized vector
is \( c\,\mathbf{z}^{(\lambda)}+k\bone \) for a constant \( k \). Since \( (\I-\M)\bone=\bzero \), its residual vector is
\( c(\I-\M)\mathbf{z}^{(\lambda)} \), and \( \text{SSE}_z \) is multiplied by \( c^2 \). By @eq-tr-profile, \( \ell_p \) changes by
\( -n\log c \), a constant.
:::

::: {#exm-tr-engel-box-cox}
[Box–Cox for Engel's data]

First let transformed food expenditure be linear in income. The profile log-likelihood, panel (a) of [Figure 22.2.1](#fig-tr-box-cox-profile), peaks at
\( \hat{\lambda}=0.63 \), and the \( 95\% \) interval given by @eq-tr-lr-interval is \( [0.45,\ 0.80] \). Both the raw scale and the
logarithm are rejected: the likelihood ratio statistics are \( 18.8 \) for \( \lambda=1 \) and \( 37.2 \) for
\( \lambda=0 \), against the critical value \( 3.84 \).

Now let food expenditure be linear in *log* income. The profile, panel (b), peaks at \( \hat{\lambda}=0.025 \), with
interval \( [-0.117,\ 0.162] \). The logarithm is well inside, with statistic \( 0.12 \), while the square root
(\( 46.6 \)) and the raw scale (\( 181.4 \)) are firmly rejected. With log income as the regressor, the
normalized residual standard deviation is \( 78.6 \) on the log scale and \( 115.6 \) on the raw scale, both in the
currency of the survey (the geometric mean of food expenditure is \( 574.4 \)).

The answers differ because a response transformation is chosen to make a *particular* mean function
linear. [Section 22.3](03-transforming-regressors.html) estimates the two powers together.
:::

::: {when-format="html"}
![**Figure 22.2.1.** Box–Cox profile log-likelihoods for Engel's data, relative to their maxima. (a) Food
expenditure linear in income. (b) Linear in log income. The dashed line is at \( -\chi^2_{1,0.05}/2 \), and the
vertical lines mark the \( 95\% \) likelihood ratio intervals.](box_cox_profile.svg){#fig-tr-box-cox-profile width=100%}
:::

::: {when-format="pdf"}
![Box–Cox profile log-likelihoods for Engel's data, relative to their maxima. (a) Food
expenditure linear in income. (b) Linear in log income. The dashed line is at \( -\chi^2_{1,0.05}/2 \), and the
vertical lines mark the \( 95\% \) likelihood ratio intervals.](box_cox_profile.pdf){width=100%}
:::

```{.python .run #cell-box-cox-profile}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)


def box_cox_z(y, lam):
    """Box-Cox transform normalized by the geometric mean, so every lambda is in the units of y."""
    gm = np.exp(np.mean(np.log(y)))
    if abs(lam) < 1e-12:
        return gm * np.log(y)
    return np.expm1(lam * np.log(y)) / (lam * gm ** (lam - 1))   # (y^lam - 1)/(lam gm^(lam-1))


def sse(v, X):
    coef, *_ = np.linalg.lstsq(X, v, rcond=None)
    r = v - X @ coef
    return r @ r


def profile(lam, y, X):
    """Profile log-likelihood of lambda, up to an additive constant."""
    return -0.5 * len(y) * np.log(sse(box_cox_z(y, lam), X))


def box_cox_fit(y, X, level=0.95):
    """MLE of lambda and the likelihood-ratio interval {lambda: 2(l_max - l_p) <= chi2_1}."""
    lam_hat = optimize.minimize_scalar(lambda l: -profile(l, y, X), bounds=(-3, 3), method="bounded",
                                       options={"xatol": 1e-10}).x
    l_max = profile(lam_hat, y, X)
    cut = stats.chi2.ppf(level, 1) / 2
    g = lambda l: l_max - profile(l, y, X) - cut
    lo = optimize.brentq(g, lam_hat - 3, lam_hat)
    hi = optimize.brentq(g, lam_hat, lam_hat + 3)
    return lam_hat, (lo, hi), l_max


X_lin = np.column_stack([np.ones(n), income])            # food linear in income
X_log = np.column_stack([np.ones(n), np.log(income)])    # food linear in log income
for name, X in [("income", X_lin), ("log income", X_log)]:
    lam_hat, (lo, hi), l_max = box_cox_fit(food, X)
    lr = {lam: 2 * (l_max - profile(lam, food, X)) for lam in (0, 0.5, 1)}
    print(f"{name:10s} lambda_hat {lam_hat:.3f}  95% interval [{lo:.3f}, {hi:.3f}]  LR at 0, 1/2, 1:",
          " ".join(f"{v:.1f}" for v in lr.values()))
```

## Choosing a convenient value

Nobody thinks in units of food expenditure raised to the power \( 0.025 \). The usual practice is to pick a
simple value inside the interval, here \( \lambda=0 \), and to analyse on that scale as if it had been chosen in
advance. Simple values have interpretations: the logarithm gives percentages and elasticities, the
reciprocal turns a time into a rate.

Treating the chosen scale as known ignores that the data chose it. Bickel and Doksum (1981) showed that
when \( \lambda \) is estimated, the unconditional variance of \( \hbeta \) can be many times its fixed-\( \lambda \) value.
Box and Cox (1982) replied that the unconditional variance of a coefficient whose *units* change with
\( \lambda \) has no scientific meaning, and Hinkley and Runger (1984) argued for inference conditional on the
selected scale. Quantities whose meaning does not depend on \( \lambda \), such as a prediction of \( y \) or the hypothesis
that a regressor has no effect at all (@exr-tr-no-effect), avoid the difficulty, and resampling
([Chapter 23](../ch23-resampling-inference/index.html)) can include the choice of scale in their
uncertainty.

## The constructed variable

A single extra regression gives a quick check of a proposed \( \lambda_0 \) and shows which observations drive
the evidence. Let
\[
\bw^{(\lambda_0)}=\frac{\partial\mathbf{z}^{(\lambda)}}{\partial\lambda}\Big|_{\lambda=\lambda_0},
\]{#eq-tr-constructed}

computed entry by entry. Differentiating \( (y^{\lambda}-1)\dot{y}^{\,1-\lambda}/\lambda \) at \( \lambda_0=1 \) gives
\( w_i=y_i\{\log(y_i/\dot{y})-1\}+1+\log\dot{y} \), and at \( \lambda_0=0 \) it gives
\( w_i=\tfrac12\dot{y}\{\log(y_i/\dot{y})\}^2-\tfrac12\dot{y}(\log\dot{y})^2 \) (@exr-tr-constructed-zero). In a model with an
intercept the additive constants do not matter, so the **constructed variables** are
\( y_i\{\log(y_i/\dot{y})-1\} \) for checking the raw scale and \( \tfrac12\dot{y}\{\log(y_i/\dot{y})\}^2 \) for checking the
logarithm.

::: {#prp-tr-constructed}
[A constructed variable for the response power]

Let \( \he_0=(\I-\M)\mathbf{z}^{(\lambda_0)} \) and \( \bw=\bw^{(\lambda_0)} \), with \( \bw\notin\C(\X) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \dfrac{d}{d\lambda}\text{SSE}_z(\lambda)\Big|_{\lambda_0}=2\he_0\T\bw \), and so
   \( \ell_p'(\lambda_0)=-n\,\he_0\T\bw/\text{SSE}_z(\lambda_0) \).

2. In the least squares fit of \( \mathbf{z}^{(\lambda_0)} \) on \( [\X,\bw] \), the coefficient of \( \bw \) is
   \( \hat{\gamma}=\he_0\T\bw/\norm{(\I-\M)\bw}^2 \). Hence \( \hat{\gamma} \) and \( \ell_p'(\lambda_0) \) have opposite signs, and
   \( \hat{\gamma}=0 \) iff \( \lambda_0 \) is a stationary point of the profile log-likelihood. In particular
   \( \hat{\gamma}=0 \) at an interior maximum \( \lambda_0=\hat{\lambda} \).
:::

:::

::: {.proof}
(a) \( \text{SSE}_z(\lambda)=\mathbf{z}^{(\lambda)\top}(\I-\M)\mathbf{z}^{(\lambda)} \). Its derivative is
\( 2\,\mathbf{z}^{(\lambda)\top}(\I-\M)\,\partial\mathbf{z}^{(\lambda)}/\partial\lambda \), because \( \I-\M \) is symmetric. At
\( \lambda_0 \) this is \( 2\he_0\T\bw \). The chain rule applied to @eq-tr-profile gives the derivative of \( \ell_p \).
(b) By @thm-proj-fwl the coefficient of \( \bw \) is \( \tilde{\bw}\T\mathbf{z}^{(\lambda_0)}/\norm{\tilde{\bw}}^2 \) with
\( \tilde{\bw}=(\I-\M)\bw\ne\bzero \), and \( \tilde{\bw}\T\mathbf{z}^{(\lambda_0)}=\bw\T(\I-\M)\mathbf{z}^{(\lambda_0)}=\bw\T\he_0 \).
Comparing with (a) gives the rest.
:::

To first order \( \mathbf{z}^{(\lambda)}\approx\mathbf{z}^{(\lambda_0)}+(\lambda-\lambda_0)\bw \), so if the model holds at \( \lambda \),
then \( \mathbf{z}^{(\lambda_0)}\approx\X\bbeta-(\lambda-\lambda_0)\bw+\be \) and \( \lambda_0-\hat{\gamma} \) is a rough one-step estimate of
\( \lambda \). The added-variable plot of \( \bw \) (@thm-proj-fwl) shows which observations carry the evidence, the use
Atkinson (1985) made of it. The \( t \) statistic of \( \hat{\gamma} \) tests \( \lambda=\lambda_0 \), but not exactly: \( \bw \) is a
function of \( \y \) itself, not only of the fitted values, so @thm-cor-fitted-regressors does not apply.

::: {#exm-tr-engel-constructed}
[The constructed variable for Engel's data]

In the model with log income, the constructed variable for \( \lambda_0=1 \) has coefficient
\( \hat{\gamma}=0.831 \) and \( t=17.96 \): overwhelming evidence against the raw scale. The one-step estimate
\( 1-0.831=0.169 \) is far from \( \hat{\lambda}=0.025 \), as expected from so far away. For \( \lambda_0=0 \) the coefficient is
\( -0.031 \) with \( t=-0.38 \).

In \( 4000 \) data sets simulated from the fitted log–log model, with normal errors and Engel's incomes (so
\( \lambda=0 \) is true), the \( 95\% \) likelihood ratio interval covers \( 0 \) in a fraction \( 0.945 \), close to nominal,
while the nominal \( 5\% \) constructed-variable test rejects in a fraction \( 0.075 \). The constructed variable is a
diagnostic; the likelihood interval is for inference.

Deleting one household at a time, the largest change moves \( \hat{\lambda} \) from \( 0.025 \) to \( 0.081 \), when the
household with income \( 2823 \) and food expenditure \( 2033 \), the largest in the survey, is removed. It lies well
above the line in panels (a) and (c) of [Figure 22.1.1](01-transforming-the-response.html#fig-tr-engel-scales),
at a point of high leverage (the leverage of Engel's households is examined in
[Section 20.2](../ch20-residuals-leverage-influence/02-leverage.html)). The logarithm stays inside the interval whichever household is removed. For a
fixed scale [Chapter 20](../ch20-residuals-leverage-influence/index.html) gives deletion formulas (@thm-res-deletion, and @prp-cmp-loo for their computation); for \( \hat{\lambda} \) the refits are cheap enough to do directly.
:::

```{.python .run #cell-box-cox-constructed}
def constructed(y, lam0, h=1e-5):
    """w = d z(lambda)/d lambda at lam0 (numerical derivative of the normalized transform)."""
    return (box_cox_z(y, lam0 + h) - box_cox_z(y, lam0 - h)) / (2 * h)


def constructed_test(y, X, lam0):
    """Regress z(lam0) on [X, w]: the coefficient of w and its t statistic."""
    w = constructed(y, lam0)
    fit = sm.OLS(box_cox_z(y, lam0), np.column_stack([X, w])).fit()
    return fit.params[-1], fit.tvalues[-1]


for lam0 in (1, 0):
    g, t = constructed_test(food, X_log, lam0)
    print(f"lambda0 = {lam0}: coefficient of w {g:.3f}, t = {t:.2f}, one-step estimate {lam0 - g:.3f}")
```

In the browser the simulation uses \( 200 \) replications, so its fractions are rough.

```{.python .run #cell-box-cox-simulation}
def simulate(reps, seed):
    rng = np.random.default_rng(seed)
    fit = sm.OLS(np.log(food), X_log).fit()           # the log-log fit: data generated with lambda = 0
    mu, sigma = fit.fittedvalues, np.sqrt(fit.scale)
    cover = reject = 0
    for _ in range(reps):
        y_sim = np.exp(mu + sigma * rng.standard_normal(n))
        lam_s, (lo_s, hi_s), _ = box_cox_fit(y_sim, X_log)
        cover += lo_s <= 0 <= hi_s
        reject += abs(constructed_test(y_sim, X_log, 0)[1]) > stats.t.ppf(0.975, n - 3)
    return cover / reps, reject / reps


coverage, size = simulate(reps=200, seed=2201)
print(f"coverage of the 95% LR interval {coverage:.3f}, size of the 5% constructed-variable test {size:.3f}")
```

::: {.warning}
\( \hat{\lambda} \) responds to anything that spoils the normal linear model: curvature, a variance that changes with
the level, skewness, single outlying cases. A value that makes the residuals look better has not
necessarily corrected the defect that mattered. After transforming, look at the residual plots of
[Chapter 20](../ch20-residuals-leverage-influence/index.html) and
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) again, on the new scale.
:::

## Exercises

### A. Check your understanding

::: {#exr-tr-box-cox-values}
[A1]

Compute \( 4^{(\lambda)} \) for \( \lambda=-1,0,\tfrac12,1,2 \). Check that the values increase with \( \lambda \), and explain,
using the definition, why \( y^{(\lambda)} \) is increasing in \( \lambda \) for every fixed \( y>0 \).
:::

::: {#exr-tr-why-normalize}
[A2]

Show that if every \( y_i \) is multiplied by \( c>0 \), the residual sum of squares of the *unnormalized*
\( \y^{(\lambda)} \) in a model with an intercept is multiplied by \( c^{2\lambda} \). Conclude that choosing \( \lambda \) by
minimizing the unnormalized residual sum of squares would give an answer that depends on the units
of \( y \).
:::

::: {.solution}
For \( \lambda\ne0 \), \( ((cy)^{\lambda}-1)/\lambda=c^{\lambda}(y^{\lambda}-1)/\lambda+(c^{\lambda}-1)/\lambda \), a multiple \( c^{\lambda} \) of
\( y^{(\lambda)} \) plus a constant, which the intercept absorbs. So the residual sum of squares is multiplied by
\( c^{2\lambda} \), a factor that depends on \( \lambda \). Changing the units reweights the comparison between different
values of \( \lambda \): as \( c \) grows the comparison is pushed towards ever more negative \( \lambda \), and as \( c \)
shrinks towards ever larger \( \lambda \). By @prp-tr-units the
normalized version has the common factor \( c^2 \) and no such problem.
:::

### B. Practice

::: {#exr-tr-constructed-zero}
[B1]

Derive the constructed variable at \( \lambda_0=0 \). Write \( z^{(\lambda)}=f(\lambda)\,\dot{y}^{\,1-\lambda} \) with
\( f(\lambda)=(y^{\lambda}-1)/\lambda \), expand \( f(\lambda)=\log y+\tfrac12\lambda(\log y)^2+O(\lambda^2) \), and show that
\( \partial z^{(\lambda)}/\partial\lambda \) at \( 0 \) is \( \tfrac12\dot{y}\{\log(y/\dot{y})\}^2-\tfrac12\dot{y}(\log\dot{y})^2 \).
:::

::: {.solution}
The expansion follows from \( y^{\lambda}=1+\lambda\log y+\tfrac12\lambda^2(\log y)^2+O(\lambda^3) \), so \( f(0)=\log y \) and
\( f'(0)=\tfrac12(\log y)^2 \). Also \( \dot{y}^{\,1-\lambda}=\dot{y}\exp(-\lambda\log\dot{y}) \) has value \( \dot{y} \) and derivative
\( -\dot{y}\log\dot{y} \) at \( 0 \). By the product rule the derivative of \( z^{(\lambda)} \) at \( 0 \) is
\[
\dot{y}\bigl\{\tfrac12(\log y)^2-\log y\log\dot{y}\bigr\}=\tfrac12\dot{y}\bigl\{(\log y-\log\dot{y})^2-(\log\dot{y})^2\bigr\}.
\]
:::

::: {#exr-tr-lr-difference}
[B2]

In @exm-tr-engel-box-cox, with log income as regressor, the normalized residual standard deviations
are \( 78.6 \) at \( \lambda=0 \) and \( 115.6 \) at \( \lambda=1 \). Show from @eq-tr-profile that the difference between the
likelihood ratio statistics at \( \lambda=1 \) and \( \lambda=0 \) is \( 2n\log\{s_z(1)/s_z(0)\} \), and check it against the
values \( 181.4 \) and \( 0.12 \) of the example.
:::

::: {#exr-tr-shifted}
[B3]

For data that may include zeros, consider the shifted family \( (y+c)^{(\lambda)} \) with an unknown \( c \) such that
every \( y_i+c>0 \). Write the log-likelihood of \( (\bbeta,\sigma^2,\lambda,c) \), show that the profile log-likelihood of
\( (\lambda,c) \) is @eq-tr-profile with \( \mathbf{z}^{(\lambda)} \) computed from \( y_i+c \) and its geometric mean, and state the
likelihood ratio confidence region for \( (\lambda,c) \).
:::

::: {.solution}
The Jacobian of \( y\mapsto(y+c)^{(\lambda)} \) is \( (y+c)^{\lambda-1} \), so the log-likelihood is
@eq-tr-box-cox-loglik with \( y_i \) replaced by \( y_i+c \) throughout. The proof of @thm-tr-box-cox goes through word
for word for each fixed \( c \), with \( \dot{y} \) the geometric mean of the \( y_i+c \). The region is
\( \{(\lambda,c):n\log\{\text{SSE}_z(\lambda,c)/\text{SSE}_z(\hat{\lambda},\hat{c})\}\le\chi^2_{2,\alpha}\} \), two degrees of freedom because
two parameters are profiled. The maximization needs care near the boundary \( c\to-\min_iy_i \): there
\( \log(\min_iy_i+c)\to-\infty \), so for \( \lambda<1 \) the Jacobian term \( (\lambda-1)\sum_i\log(y_i+c) \) grows without
bound, and the likelihood can be unbounded, as it is for the three-parameter lognormal distribution.
:::

### C. Going deeper

::: {#exr-tr-no-effect}
[C1]

Suppose that the conditional distribution of \( Y \) given \( (x_1,\dots,x_k) \) does not depend on \( x_1 \). Show that the
same holds for \( g(Y) \) for every strictly increasing \( g \), so the hypothesis that \( x_1 \) has no effect is
the same hypothesis on every scale. Why does the converse direction, from a zero coefficient on
one scale to a zero coefficient on another, fail when the conditional distribution does depend on
\( x_1 \)? Give an example with two groups whose medians are equal but whose means are not.
:::

::: {#exr-tr-box-cox-influence}
[C2]

Let \( \hat{\lambda} \) be an interior maximizer of \( \ell_p \), and write \( s(\lambda)=\he(\lambda)\T\bw^{(\lambda)} \), which is zero at
\( \hat{\lambda} \) by @prp-tr-constructed. Suppose case \( i \) is deleted. Show, by one Newton step for the equation
\( s=0 \), that the change in \( \hat{\lambda} \) is approximately proportional to the contribution
\( \hat{\varepsilon}_{(i)}\tilde{w}_{(i)} \) of that case, where both factors are computed with the case deleted. Which
cases does this predict to be influential for \( \hat{\lambda} \), and how does leverage in the regression of
\( \bw \) on \( \X \) enter?
:::
