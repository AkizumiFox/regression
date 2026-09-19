# Maximum likelihood under normality

The Gauss–Markov theorem needs only two moments. That is its strength, and also the reason it can
say nothing about \( \sigma^2 \), nothing about nonlinear estimators, and nothing about the distribution
of the estimates. From now on we assume more: the errors are normal. The linear model becomes a
fully specified family of distributions,
\[
\Y\sim\Normal_n(\X\bbeta,\sigma^2\I),\qquad \bbeta\in\Real^p,\ \sigma^2>0,
\]{#eq-opt-normal-model}

the normal linear model of @def-lm-linear-model. The first general-purpose method of estimation in
such a family is maximum likelihood. This section shows that it gives back least squares for
\( \bbeta \), and gives a slightly different estimate of \( \sigma^2 \).

## The likelihood

By @thm-mvn-density with \( \bSigma=\sigma^2\I \), which is positive definite with
\( \det(\sigma^2\I)=\sigma^{2n} \), the density of \( \Y \) at the observed \( \y \) is
\[
L(\bbeta,\sigma^2)=(2\pi\sigma^2)^{-n/2}\exp\Bigl\{-\frac{\norm{\y-\X\bbeta}^2}{2\sigma^2}\Bigr\}.
\]
As a function of the parameters with \( \y \) held fixed, this is the **likelihood**. Its logarithm is
\[
\ell(\bbeta,\sigma^2)=-\frac n2\log(2\pi)-\frac n2\log\sigma^2-\frac{\norm{\y-\X\bbeta}^2}{2\sigma^2}.
\]{#eq-opt-loglik}

The data enter the log-likelihood only through the sum of squares \( \norm{\y-\X\bbeta}^2 \). That is the first
sign that normality and least squares belong together. The Gaussian density penalizes a residual by
its square, so the most likely parameter values are those with the smallest squared residuals.

One technical point matters. If \( \y\in\C(\X) \), the residual sum of squares at the least squares fit is
zero, and \( \ell\to\infty \) as \( \sigma^2\to0 \), so no maximum exists. This happens with probability
zero when \( n>r \). The subspace \( \C(\X) \) has dimension \( r<n \), so it has Lebesgue measure zero, and
\( \Y \) has a density. We therefore assume \( n>r \) throughout, and \( \text{SSE}>0 \) for the observed data.

## The maximum likelihood estimates

::: {#thm-opt-mle}
[Maximum likelihood in the normal linear model]

Assume @eq-opt-normal-model with \( n>r \), and let \( \y \) be such that
\( \text{SSE}=\norm{(\I-\M)\y}^2>0 \). Then \( (\tilde{\bbeta},\tilde{\sigma}^2) \) maximizes the likelihood iff
\( \tilde{\bbeta} \) is a least squares estimate and
\[
\tilde{\sigma}^2=\hat{\sigma}^2_{\text{ML}}=\frac{\text{SSE}}{n}.
\]
The maximized log-likelihood is
\( \ell_{\max}=-\frac n2\bigl\{\log(2\pi\,\text{SSE}/n)+1\bigr\} \). The maximum likelihood estimate of
\( \X\bbeta \) is \( \M\y \), and for estimable \( \blambda\T\bbeta \) it is \( \blambda\T\hbeta \). Both are unique.
:::

::: {.proof}
Fix \( \sigma^2>0 \). By @eq-opt-loglik, \( \ell(\bb,\sigma^2) \) is a decreasing function of
\( \norm{\y-\X\bb}^2 \). So it is maximized over \( \bb \) exactly at the least squares estimates, and there
\( \norm{\y-\X\bb}^2=\text{SSE} \) (@thm-proj-ls-projection). This gives
\[
\begin{aligned}
\ell(\bb,\sigma^2)&\le\ell(\hbeta,\sigma^2)=g(\sigma^2)\\
&:=-\frac n2\log(2\pi)-\frac n2\log\sigma^2-\frac{\text{SSE}}{2\sigma^2},
\end{aligned}
\]
with equality iff \( \bb \) is a least squares estimate. Now put \( v_0=\text{SSE}/n \) and
\( u=v_0/\sigma^2>0 \). Then \( \text{SSE}/(2\sigma^2)=nu/2 \) and
\( \log\sigma^2=\log v_0-\log u \), so
\[
g(v_0)-g(\sigma^2)=\frac n2\bigl(u-1-\log u\bigr)\ge0,
\]
because \( \log u\le u-1 \) for all \( u>0 \), with equality only at \( u=1 \). Hence
\( \ell(\bb,\sigma^2)\le g(\sigma^2)\le g(v_0) \), with equality throughout iff \( \bb \) is a least squares
estimate and \( \sigma^2=v_0 \). The value \( g(v_0) \) is the stated maximum. The least squares estimates
all give the fitted vector \( \M\y \) (@thm-proj-ls-projection), and all give the same value of
\( \blambda\T\hbeta \) for \( \blambda\in\C(\X\T) \) (@thm-proj-invariant-functions).
:::

The proof finds a global maximum without solving any equations. Calculus gives the same answer and
is the usual route in a first course: set the gradient
\( \partial\ell/\partial\bbeta=\X\T(\y-\X\bbeta)/\sigma^2 \) to zero, which gives the normal equations, and
\( \partial\ell/\partial\sigma^2=-n/(2\sigma^2)+\norm{\y-\X\bbeta}^2/(2\sigma^4) \) to zero. But the
stationarity argument still has to be completed by showing that the stationary point is a
maximum (@exr-opt-mle-calculus). The inequality \( \log u\le u-1 \) handles existence, maximality and uniqueness
together.

When \( \X \) is rank deficient, the likelihood is constant along the affine set of least squares
estimates \( \hbeta+\Null(\X) \), so the MLE of \( \bbeta \) is not unique. This is a general fact about
nonidentified parameters. If \( \X\bbeta_1=\X\bbeta_2 \), the two parameter values give the same distribution
of \( \Y \) and hence the same likelihood for every data set. Maximum likelihood cannot separate what the
model itself does not separate. What it does determine uniquely, \( \X\bbeta \) and the estimable functions,
is exactly what least squares determines uniquely.

::: {#exm-opt-stackloss-mle}
[Stack loss]

Brownlee's stack loss data record \( n=21 \) days of operation of a plant that oxidizes
ammonia to nitric acid. The response is the percentage of ammonia lost (times ten). The regressors
are air flow, the temperature of the cooling water, and the acid concentration. With an intercept,
\( p=4 \), and least squares gives
\[
\hbeta=(-39.9197,\ 0.7156,\ 1.2953,\ {-}0.1521)\T,\qquad \text{SSE}=178.830 .
\]
By @thm-opt-mle these are also the maximum likelihood estimates of \( \bbeta \). The MLE of \( \sigma^2 \) is
\( 178.830/21=8.516 \), while the unbiased estimate is
\( s^2=178.830/17=10.519 \). The maximized log-likelihood is
\( -52.288 \). A general-purpose numerical optimizer, started far from the answer and
knowing nothing about least squares, arrives at the same point.
:::

```{.python .run #cell-mle-profile-fit}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape

beta_hat, *_ = np.linalg.lstsq(X, y, rcond=None)
sse = np.sum((y - X @ beta_hat) ** 2)
sigma2_mle, s2 = sse / n, sse / (n - p)
loglik_max = -n / 2 * (np.log(2 * np.pi * sigma2_mle) + 1)
print("beta_hat   ", np.round(beta_hat, 4))
print(f"SSE = {sse:.3f}   MLE of sigma^2 = {sigma2_mle:.3f}   s^2 = {s2:.3f}")
print(f"maximized log-likelihood = {loglik_max:.3f}")
```

```{.python .run #cell-mle-profile-numerical}
def negloglik(theta):
    """Minus the normal log-likelihood; theta = (beta, log sigma^2)."""
    b, log_s2 = theta[:p], theta[p]
    r = y - X @ b
    return n / 2 * (np.log(2 * np.pi) + log_s2) + r @ r / (2 * np.exp(log_s2))

start = np.r_[np.zeros(p), np.log(np.var(y))]
opt = optimize.minimize(negloglik, start, method="BFGS", options={"gtol": 1e-8})
print("numerical MLE of beta  ", np.round(opt.x[:p], 4))
print("numerical MLE sigma^2  ", round(float(np.exp(opt.x[p])), 3))
```

## The bias of the variance estimate

Unlike least squares for \( \bbeta \), the MLE of \( \sigma^2 \) is biased. Since
\( \E(\text{SSE})=\sigma^2(n-r) \) (@thm-lm-sigma2 in full rank, and [Section 6.4](../ch06-projections/04-least-squares.html)
in general),
\[
\E\hat{\sigma}^2_{\text{ML}}=\frac{n-r}{n}\,\sigma^2 .
\]
The bias is downward, by the fraction \( r/n \). The residual vector lives in the \( (n-r) \)-dimensional space
\( \C(\X)\perpc \), and the MLE divides its squared length by \( n \) as if it had \( n \) dimensions. The \( r \)
dimensions spent on fitting the mean are the dimensions in which the fit has pulled the data toward
itself. With a few regressors and many observations the bias is negligible. It is not negligible when the
number of parameters grows with the sample size.

::: {#exm-opt-neyman-scott}
[Many means, one variance]

Suppose \( m \) laboratories each make two measurements, \( Y_{i1},Y_{i2}\iid\Normal(\mu_i,\sigma^2) \), with
a different \( \mu_i \) for each laboratory and a common measurement variance. This is a linear model with
\( n=2m \) and \( r=m \), one indicator column per laboratory. The fitted value for each pair is its mean,
and
\[
\begin{aligned}
\hat{\sigma}^2_{\text{ML}}&=\frac{1}{2m}\sum_{i=1}^m\sum_{j=1}^2(Y_{ij}-\bar{Y}_{i\cdot})^2
=\frac1{4m}\sum_{i=1}^m(Y_{i1}-Y_{i2})^2,\\
\E\hat{\sigma}^2_{\text{ML}}&=\frac{\sigma^2}{2}.
\end{aligned}
\]
The bias does not disappear as \( m\to\infty \). By the law of large numbers, \( \hat{\sigma}^2_{\text{ML}}\to\sigma^2/2 \)
with probability one, so the MLE is *inconsistent*. The unbiased \( s^2=\text{SSE}/(n-r) \) is exactly twice as large
and is consistent. This is the example of Neyman and Scott (1948). It shows that the good large-sample
behaviour of maximum likelihood can fail when the number of nuisance parameters grows with the data.
Restricted maximum likelihood, which removes the mean before estimating the variance
(@exr-opt-reml), repairs it, and Chapter 32 builds variance-component estimation on that idea.
:::

Unbiasedness is not the only possible criterion for choosing a divisor. Under normality we can compare
estimators of the form \( \text{SSE}/d \) by mean squared error.

::: {#prp-opt-divisor}
[The best multiple of SSE]

Assume @eq-opt-normal-model with \( k=n-r\ge1 \). Among estimators \( \text{SSE}/d \) with \( d>0 \),
\[
\frac{\E\bigl(\text{SSE}/d-\sigma^2\bigr)^2}{\sigma^4}=\frac{2k}{d^2}+\Bigl(\frac kd-1\Bigr)^2,
\]
which is minimized at \( d=k+2 \), with minimum value \( 2/(k+2) \). The unbiased choice \( d=k \) gives \( 2/k \),
and the maximum likelihood choice \( d=n \) gives \( (2k+r^2)/n^2 \).
:::

::: {.proof}
Write \( \text{SSE}=\Y\T(\I-\M)\Y \). By @thm-qf-mean-var with \( \A=\I-\M \), \( \bSigma=\sigma^2\I \) and
\( \bmu=\X\bbeta \), for which \( \A\bmu=\bzero \),
\[
\E(\text{SSE})=\sigma^2\tr(\I-\M)=\sigma^2k,\qquad
\Var(\text{SSE})=2\sigma^4\tr\bigl((\I-\M)^2\bigr)=2\sigma^4k .
\]
The mean squared error is variance plus squared bias:
\( 2\sigma^4k/d^2+(k/d-1)^2\sigma^4 \). As a function of \( t=1/d \), this is
\( \sigma^4\{(2k+k^2)t^2-2kt+1\} \), a convex quadratic minimized at \( t=k/(2k+k^2)=1/(k+2) \). Substituting
\( d=k+2 \), \( k \) and \( n=k+r \) gives the three values. For the last,
\( 2k/n^2+(k-n)^2/n^2=(2k+r^2)/n^2 \).
:::

For the stack loss data, \( k=17 \), and the three relative mean squared errors are
\( 0.1176 \) for \( s^2 \), \( 0.1134 \) for the MLE and
\( 0.1053 \) for \( \text{SSE}/(k+2) \). The unbiased estimator is the worst of the three by
this criterion. Comparing \( (2k+r^2)/n^2 \) with \( 2/k \) shows that the MLE beats \( s^2 \) iff
\( (r-4)k<2r \). So it always does when \( r\le4 \), and it loses once many parameters are fitted relative
to the residual degrees of freedom, as in @exm-opt-neyman-scott. In practice \( s^2 \) is used anyway, because it is what makes the
\( t \) and \( F \) statistics of [Section 7.5](05-sampling-distributions.html) exact. The choice of divisor
matters for point estimation of \( \sigma^2 \) and for nothing else. Any fixed multiple of SSE leads to the same
tests and intervals once the matching distribution is used.

## Profile likelihood

In most problems only some parameters are of interest. The others, called **nuisance parameters**, still
appear in the likelihood. A simple and general way to eliminate them is to maximize them out.

::: {#def-opt-profile}
[Profile likelihood]

Let the parameter be \( \boldsymbol{\uptheta}=(\psi,\boldsymbol{\upeta}) \), where \( \psi \) is of interest and \( \boldsymbol{\upeta} \) is a nuisance
parameter. The **profile log-likelihood** of \( \psi \) is
\( \ell_p(\psi)=\sup_{\boldsymbol{\upeta}}\ell(\psi,\boldsymbol{\upeta}) \).
:::

The profile log-likelihood is maximized at the MLE of \( \psi \), with maximum \( \ell_{\max} \). The drop
\( \ell_{\max}-\ell_p(\psi) \) measures how much less plausible the value \( \psi \) is than the best value.
In the linear model the profiles have closed forms. The proof of @thm-opt-mle computed two of them.
With \( \sigma^2 \) as nuisance, the profile of \( \bbeta \) is
\[
\ell_p(\bbeta)=-\frac n2\Bigl\{\log\Bigl(\frac{2\pi\norm{\y-\X\bbeta}^2}{n}\Bigr)+1\Bigr\},
\]
a decreasing function of \( \norm{\y-\X\bbeta}^2 \). Its level sets are the level sets of the residual sum of
squares. When \( \X \) has full rank they are ellipsoids centred at \( \hbeta \), and they are the basis of the
confidence ellipsoids of [Chapter 12](../ch12-intervals-and-bands/index.html) (@thm-ci-ellipsoid). With \( \bbeta \) as nuisance, the profile of \( \sigma^2 \) is the function \( g \) in the proof.

The most useful case is a single estimable function. Here the profile likelihood turns out to be a
function of the familiar \( t \) statistic.

::: {#prp-opt-profile-t}
[Profile likelihood of an estimable function]

Assume @eq-opt-normal-model with \( n>r \), and let \( \blambda\in\C(\X\T) \), \( \blambda\ne\bzero \). For
\( \psi\in\Real \), put
\[
t(\psi)=\frac{\blambda\T\hbeta-\psi}{s\sqrt{\blambda\T(\X\T\X)\ginv\blambda}},\qquad s^2=\frac{\text{SSE}}{n-r}.
\]
Then the smallest residual sum of squares subject to \( \blambda\T\bbeta=\psi \) is
\[
\text{SSE}(\psi)=\text{SSE}+\frac{(\blambda\T\hbeta-\psi)^2}{\blambda\T(\X\T\X)\ginv\blambda},
\]{#eq-opt-constrained-sse}

and the profile log-likelihood of \( \psi=\blambda\T\bbeta \), with \( (\bbeta,\sigma^2) \) otherwise free, satisfies
\[
2\{\ell_{\max}-\ell_p(\psi)\}=n\log\Bigl(1+\frac{t(\psi)^2}{n-r}\Bigr).
\]{#eq-opt-profile-t}

:::

::: {.proof}
Write \( \blambda=\X\T\boldsymbol{\uprho} \) and \( \mathbf{c}=\M\boldsymbol{\uprho} \). Then \( \mathbf{c}\ne\bzero \), since \( \X\T\mathbf{c}=\X\T\boldsymbol{\uprho}=\blambda\ne\bzero \),
and \( \norm{\mathbf{c}}^2=\blambda\T(\X\T\X)\ginv\blambda \) as in the proof of @thm-opt-gauss-markov(a). For any
\( \bb \), @eq-proj-distance-split gives
\( \norm{\y-\X\bb}^2=\text{SSE}+\norm{\X\hbeta-\X\bb}^2 \). Put \( \bu=\X(\hbeta-\bb)\in\C(\X) \). The constraint
\( \blambda\T\bb=\psi \) reads \( \boldsymbol{\uprho}\T\X(\hbeta-\bb)=\blambda\T\hbeta-\psi=:\delta \), that is,
\( \boldsymbol{\uprho}\T\bu=\delta \). Since \( \bu\in\C(\X) \), \( \boldsymbol{\uprho}\T\bu=\boldsymbol{\uprho}\T\M\bu=\mathbf{c}\T\bu \). By the Cauchy–Schwarz
inequality, \( \lvert\delta\rvert=\lvert\mathbf{c}\T\bu\rvert\le\norm{\mathbf{c}}\norm{\bu} \), so
\( \norm{\bu}^2\ge\delta^2/\norm{\mathbf{c}}^2 \). Equality holds for \( \bu=\delta\mathbf{c}/\norm{\mathbf{c}}^2 \), which lies in
\( \C(\X) \) and so equals \( \X(\hbeta-\bb) \) for some \( \bb \). That \( \bb \) satisfies the constraint. This proves @eq-opt-constrained-sse.
Maximizing over \( \sigma^2 \) as in the proof of @thm-opt-mle, with SSE replaced by
\( \text{SSE}(\psi) \), gives \( \ell_p(\psi)=-\frac n2\{\log(2\pi\,\text{SSE}(\psi)/n)+1\} \). Hence
\[
\begin{aligned}
2\{\ell_{\max}-\ell_p(\psi)\}&=n\log\frac{\text{SSE}(\psi)}{\text{SSE}}
=n\log\Bigl(1+\frac{\delta^2}{\norm{\mathbf{c}}^2\,\text{SSE}}\Bigr)\\
&=n\log\Bigl(1+\frac{t(\psi)^2}{n-r}\Bigr),
\end{aligned}
\]
since \( \delta^2/(\norm{\mathbf{c}}^2\,\text{SSE})=t(\psi)^2/(n-r) \).
:::

Two consequences show up immediately. First, any procedure based on the profile likelihood of
\( \blambda\T\bbeta \) is a procedure based on \( \lvert t\rvert \), because @eq-opt-profile-t is increasing in \( t^2 \).
The likelihood ratio test of \( \blambda\T\bbeta=\psi_0 \) is therefore the \( t \) test, a fact [Chapter 11](../ch11-general-linear-hypothesis/index.html) extends
to general hypotheses (@thm-glh-lrt). Second, the large-sample calibration of the likelihood ratio, "reject when
\( 2\{\ell_{\max}-\ell_p\} \) exceeds a \( \chi^2(1) \) quantile", is only approximate. The exact calibration uses the \( t(n-r) \)
distribution of [Section 7.5](05-sampling-distributions.html).

::: {#exm-opt-stackloss-profile}
[Profile likelihood for the air-flow coefficient]

For the stack loss model, [Figure 7.3.1](#fig-opt-profile) plots
\( \ell_p(\beta_1)-\ell_{\max} \) for the air-flow coefficient \( \beta_1 \), computed by refitting the model with
\( \beta_1 \) fixed on a grid. The listing confirms @eq-opt-profile-t to rounding error. The standard error of
\( \hat{\beta}_1 \) is \( 0.1349 \). The set of \( \beta_1 \) with
\( 2\{\ell_{\max}-\ell_p\}\le3.841 \), the \( 95\% \) point of \( \chi^2(1) \), is the likelihood interval
\( (0.467,\ 0.965) \). By @eq-opt-profile-t it consists of the values with
\( \lvert t\rvert\le1.847 \). The exact \( 95\% \) interval uses the \( t(17) \) quantile
\( 2.110 \) instead, and is \( (0.431,\ 1.000) \). With only
\( 17 \) residual degrees of freedom, the \( \chi^2 \) calibration makes the interval about \( 12\% \)
too short.
:::

```{.python .run #cell-mle-profile-profile}
def profile_loglik(b1):
    """Profile log-likelihood of the air-flow coefficient: maximize over everything else."""
    others = np.delete(X, 1, axis=1)
    z = y - b1 * X[:, 1]                                 # move the fixed term to the left side
    coef, *_ = np.linalg.lstsq(others, z, rcond=None)
    sse_b = np.sum((z - others @ coef) ** 2)             # SSE(b1): smallest SSE with beta_1 = b1
    return -n / 2 * (np.log(2 * np.pi * sse_b / n) + 1)

se1 = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
grid = np.linspace(beta_hat[1] - 4 * se1, beta_hat[1] + 4 * se1, 201)
lp = np.array([profile_loglik(b) for b in grid])
t = (beta_hat[1] - grid) / se1
print("max |2(l_max - l_p) - n log(1 + t^2/(n-p))| =",
      np.abs(2 * (loglik_max - lp) - n * np.log1p(t ** 2 / (n - p))).max())
```

::: {when-format="html"}
![**Figure 7.3.1.** Profile log-likelihood of the air-flow coefficient in the stack loss regression,
relative to its maximum. The ticks on the dashed cutoff line mark the \( \chi^2 \)-calibrated likelihood interval,
and the bar below it marks the exact \( 95\% \) \( t \) interval, which is wider.](profile_likelihood.svg){#fig-opt-profile width=70%}
:::

::: {when-format="pdf"}
![Profile log-likelihood of the air-flow coefficient in the stack loss regression,
relative to its maximum. The ticks on the dashed cutoff line mark the \( \chi^2 \)-calibrated likelihood interval,
and the bar below it marks the exact \( 95\% \) \( t \) interval, which is wider.](profile_likelihood.pdf){width=70%}
:::

::: {.remark}
[Other error distributions]

The agreement between maximum likelihood and least squares is special to normal errors. If the errors
have the Laplace density \( (2\tau)^{-1}e^{-\lvert\varepsilon\rvert/\tau} \), the log-likelihood is
\( -n\log(2\tau)-\sum_i\lvert y_i-\x_{(i)}\T\bbeta\rvert/\tau \), and maximum likelihood minimizes the sum of *absolute*
residuals; @exr-opt-laplace-mle works out the location case. For errors with a \( t \) distribution it downweights large residuals. For each error
distribution the likelihood picks its own loss function. Least squares is the one that goes with the
Gaussian. [Section 20.6](../ch20-residuals-leverage-influence/06-robust.html) uses robust fitting as a diagnostic, and Chapter 45 treats median and quantile regression.
:::

## Exercises

### A. Check your understanding

::: {#exr-opt-location-mle}
[A1]

For \( Y_1,\dot{s},Y_n\iid\Normal(\mu,\sigma^2) \), write the model in the form @eq-opt-normal-model and read off
the MLEs from @thm-opt-mle. Show that \( \hat{\sigma}^2_{\text{ML}}=(n-1)s^2/n \), where \( s^2 \) is the sample
variance.
:::

::: {#exr-opt-lr-nested}
[A2]

Let \( \C(\X_0)\subseteq\C(\X) \) with residual sums of squares \( \text{SSE}_0\ge\text{SSE} \). Show that the ratio of
maximized likelihoods is \( L_{0,\max}/L_{\max}=(\text{SSE}/\text{SSE}_0)^{n/2} \), and write it as a decreasing
function of \( (\text{SSE}_0-\text{SSE})/\text{SSE} \).
:::

::: {.solution}
By @thm-opt-mle, \( L_{\max}=(2\pi\,\text{SSE}/n)^{-n/2}e^{-n/2} \), and similarly for the smaller model. The
ratio is \( (\text{SSE}_0/\text{SSE})^{-n/2}=\{1+(\text{SSE}_0-\text{SSE})/\text{SSE}\}^{-n/2} \), which decreases in
\( (\text{SSE}_0-\text{SSE})/\text{SSE} \). That quantity is, up to degrees of freedom, the \( F \) statistic of @thm-qf-nested-f.
So the likelihood ratio test for nested normal linear models is the \( F \) test.
[Chapter 11](../ch11-general-linear-hypothesis/index.html) develops this (@thm-glh-lrt).
:::

### B. Practice

::: {#exr-opt-mle-calculus}
[B1]

Assume \( \rank(\X)=p \). Compute the gradient and the Hessian of \( \ell(\bbeta,\sigma^2) \), show that the only
stationary point is \( (\hbeta,\text{SSE}/n) \), and show that the Hessian there is negative definite. Why
does this not by itself prove that the stationary point is a global maximum? What additional
argument completes the proof?
:::

::: {#exr-opt-neyman-scott}
[B2]

In @exm-opt-neyman-scott, verify the identity
\( \sum_j(Y_{ij}-\bar{Y}_{i\cdot})^2=(Y_{i1}-Y_{i2})^2/2 \), compute \( \E\hat{\sigma}^2_{\text{ML}} \) and \( \Var\hat{\sigma}^2_{\text{ML}} \),
and prove that \( \hat{\sigma}^2_{\text{ML}}\to\sigma^2/2 \) in probability as \( m\to\infty \).
:::

::: {.solution}
\( Y_{i1}-\bar{Y}_{i\cdot}=(Y_{i1}-Y_{i2})/2=-(Y_{i2}-\bar{Y}_{i\cdot}) \), so the sum of squares is
\( 2\cdot(Y_{i1}-Y_{i2})^2/4 \). Now \( D_i=Y_{i1}-Y_{i2}\sim\Normal(0,2\sigma^2) \) independently, so
\( D_i^2/(2\sigma^2)\sim\chi^2(1) \), with mean \( 1 \) and variance \( 2 \). Hence
\( \hat{\sigma}^2_{\text{ML}}=\frac1{4m}\sum_iD_i^2 \) has mean \( 2\sigma^2/4=\sigma^2/2 \) and variance
\( m\cdot4\sigma^4\cdot2/(16m^2)=\sigma^4/(2m) \). The variance tends to zero, so by Chebyshev's inequality
\( \hat{\sigma}^2_{\text{ML}}\to\sigma^2/2 \) in probability.
:::

::: {#exr-opt-gls-mle}
[B3]

Suppose \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\V) \) with \( \V \) known and positive definite and \( \X \) of full rank. Find the
MLEs of \( \bbeta \) and \( \sigma^2 \).
:::

### C. Going deeper

::: {#exr-opt-laplace-mle}
[C1]

In the location model \( Y_i=\mu+\varepsilon_i \), \( i=1,\dot{s},n \), let the errors be independent with the Laplace
density \( (2\tau)^{-1}e^{-\lvert\varepsilon\rvert/\tau} \), so \( \Var(\varepsilon_i)=2\tau^2 \).

::: {.enumerate options="label=(\alph*)"}
1. Show that every sample median maximizes the likelihood in \( \mu \), whatever \( \tau \) is, and that the MLE of
   \( \tau \) is the mean absolute deviation \( n^{-1}\sum_i\lvert Y_i-\tilde{\mu}\rvert \) about a median \( \tilde{\mu} \).

2. The sample median of \( n \) observations from a density \( f \) that is positive and continuous at its median
   \( \mu \) has asymptotic variance \( 1/\{4nf(\mu)^2\} \). Compute this for the Laplace density and compare it with
   \( \Var(\bar{Y}) \). Which is larger, and by what factor?

3. \( \bar{Y} \) is the BLUE of \( \mu \) (@cor-opt-gm-consequences). Explain why (b) does not contradict
   @thm-opt-gauss-markov, and contrast the answer with @exm-opt-midrange, where the winning nonlinear
   estimator used the extremes rather than the middle of the sample.
:::
:::

::: {.solution}
(a) The log-likelihood is \( -n\log(2\tau)-\tau^{-1}\sum_i\lvert Y_i-\mu\rvert \). For fixed \( \tau \) it is maximized
by minimizing \( g(\mu)=\sum_i\lvert Y_i-\mu\rvert \). The function \( g \) is convex and piecewise linear, with slope
\( \#\{i:Y_i<\mu\}-\#\{i:Y_i>\mu\} \) between data points. The slope is \( \le0 \) to the left of every median and
\( \ge0 \) to the right, so the minimizers are exactly the sample medians. Substituting \( \tilde{\mu} \) and maximizing
\( -n\log\tau-g(\tilde{\mu})/\tau \) over \( \tau>0 \) gives \( \hat{\tau}=g(\tilde{\mu})/n \).
(b) \( f(\mu)=1/(2\tau) \), so the median has asymptotic variance \( \tau^2/n \), while
\( \Var(\bar{Y})=2\tau^2/n \). The mean has twice the variance of the median.
(c) The median is not a linear function of \( \Y \): its weights depend on the ordering of the sample. The
Gauss–Markov theorem compares \( \bar{Y} \) only with linear unbiased estimators. With Laplace errors, which have
heavier tails than the normal, the extremes are unreliable and the median gains by ignoring them; with uniform
errors the extremes are the most informative observations and the midrange gains by using only them.
:::

::: {#exr-opt-reml}
[C2]

Let \( \Q_2 \) be an \( n\times(n-r) \) matrix whose columns are an orthonormal basis of \( \C(\X)\perpc \), and let
\( \Z=\Q_2\T\Y \) ("error contrasts"). Show that \( \Z\sim\Normal_{n-r}(\bzero,\sigma^2\I) \), whatever \( \bbeta \) is, and that
\( \norm{\Z}^2=\text{SSE} \). Show that the MLE of \( \sigma^2 \) based on \( \Z \) alone is \( \text{SSE}/(n-r) \). This is
**restricted** (or residual) maximum likelihood.
:::

::: {.solution}
\( \Q_2\T\X=\mathbf{0} \), so \( \Z=\Q_2\T\be\sim\Normal_{n-r}(\bzero,\sigma^2\Q_2\T\Q_2)=\Normal_{n-r}(\bzero,\sigma^2\I) \)
by @thm-mvn-linear. Also \( \Q_2\Q_2\T=\I-\M \) (@prp-proj-orthonormal-formula), so
\( \norm{\Z}^2=\Y\T\Q_2\Q_2\T\Y=\Y\T(\I-\M)\Y=\text{SSE} \). The log-likelihood of \( \Z \) is
\( -\frac{n-r}2\log(2\pi\sigma^2)-\norm{\Z}^2/(2\sigma^2) \), which is the function \( g \) of the proof of
@thm-opt-mle with \( n \) replaced by \( n-r \). It is maximized at \( \norm{\Z}^2/(n-r)=\text{SSE}/(n-r) \). In
@exm-opt-neyman-scott this gives \( s^2 \), which is consistent.
:::
