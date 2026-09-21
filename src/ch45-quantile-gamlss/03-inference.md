# Inference for quantile regression

The estimator of [Section 45.2](02-quantile-regression.html) minimizes a criterion
with no derivative at the kink and no closed form, yet its large-sample behaviour is
as clean as that of least squares, with one new ingredient: the asymptotic variance
involves the *density* of the response at the fitted quantile, which is hard to
estimate well. Everything practical below is an attempt either to estimate that
density or to avoid it.

## An identity of Knight's

The analysis rests on an exact decomposition of the increment of the check loss into
a linear part and a controllable remainder.

::: {#lem-qnt-knight}
[Knight's identity]

For all real \( u \) and \( v \), with \( \psi_{\tau}(u)=\tau-1\{u<0\} \),
\[
\rho_{\tau}(u-v)-\rho_{\tau}(u)
=-v\,\psi_{\tau}(u)+\int_{0}^{v}\bigl\{1\{u\le s\}-1\{u\le0\}\bigr\}\,ds .
\]{#eq-qnt-knight}

:::

::: {.proof}
Both sides vanish at \( v=0 \), so it suffices to check four cases. Suppose first
\( u>0 \). If \( 0<v\le u \) the left side is \( \tau(u-v)-\tau u=-\tau v \) and the
integrand vanishes for \( s\le v<u \), so the right side is \( -\tau v \) too. If
\( v>u \) the left side is \( (\tau-1)(u-v)-\tau u=-u-\tau v+v \), while the integral
contributes \( v-u \), giving \( -\tau v+(v-u) \): the same. If \( v<0 \) the integral
runs backwards over \( s\in(v,0) \), where \( s<0<u \) makes the integrand zero, and
the left side is \( -\tau v \) again.

Now suppose \( u\le 0 \), so \( \psi_{\tau}(u)=\tau-1 \) and
\( 1\{u\le0\}=1 \). If \( v>0 \) then for \( s\in(0,v) \) we have \( u\le0<s \), the
integrand is zero, and both sides equal \( -v(\tau-1) \). If \( u\le v<0 \) the same
computation holds with the integral over \( (v,0) \), where \( u\le v<s \) again makes
the integrand vanish. If \( v<u\le0 \), the left side is
\( \tau(u-v)-u(\tau-1)=u-\tau v \); the integrand equals \( -1 \) on \( (v,u) \) and
\( 0 \) on \( (u,0) \), so the backwards integral contributes \( u-v \), and the right
side is \( -v(\tau-1)+(u-v)=u-\tau v \).
:::

## The limiting distribution

::: {#thm-qnt-asymptotics}
[Asymptotic normality of a regression quantile]

Fix \( \tau \) and suppose the rows \( \x_{(i)} \) are fixed and the errors
\( u_i=y_i-\x_{(i)}\T\bbeta_{\tau} \) are independent with distribution functions
\( F_i \) satisfying \( F_i(0)=\tau \). Assume

::: {.enumerate options="label=(C\arabic*)"}
1. each \( F_i \) has a density \( f_i \) in a neighbourhood of zero with
   \( f_i(0)>0 \); the \( f_i \) are equicontinuous at zero, that is
   \( \sup_i\sup_{\lvert s\rvert\le t}\lvert f_i(s)-f_i(0)\rvert\to0 \) as
   \( t\downarrow0 \); and
   \( \sup_i\sup_{\lvert s\rvert<\epsilon}f_i(s)<\infty \) for some \( \epsilon>0 \);

2. \( n^{-1}\sum_i\x_{(i)}\x_{(i)}\T\to\bD_0 \) and
   \( n^{-1}\sum_if_i(0)\,\x_{(i)}\x_{(i)}\T\to\bD_1 \), both positive definite;

3. \( n^{-1/2}\max_i\norm{\x_{(i)}}\to0 \).
:::

Then
\[
\sqrt{n}\,(\hbeta_{\tau}-\bbeta_{\tau})
\;\xrightarrow{d}\;\Normal_p\bigl(\bzero,\;\tau(1-\tau)\,\bD_1^{-1}\bD_0\bD_1^{-1}\bigr).
\]{#eq-qnt-sandwich}

If the errors are identically distributed with density \( f \), then
\( \bD_1=f(0)\bD_0 \) and the covariance collapses to
\[
\frac{\tau(1-\tau)}{f(0)^2}\,\bD_0^{-1}
=\tau(1-\tau)\,s(\tau)^2\,\bD_0^{-1},
\qquad s(\tau)=\frac{1}{f\{F^{-1}(\tau)\}},
\]{#eq-qnt-iid}

where \( s \) is the **sparsity function**, the derivative of the quantile function.
:::

::: {.proof}
Write \( \boldsymbol{\updelta}\in\Real^p \), \( c_i=\x_{(i)}\T\boldsymbol{\updelta}/\sqrt{n} \) and
\[
Z_n(\boldsymbol{\updelta})=\sum_{i=1}^n\bigl\{\rho_{\tau}(u_i-c_i)-\rho_{\tau}(u_i)\bigr\},
\]
a convex function of \( \boldsymbol{\updelta} \) minimized at
\( \hat{\boldsymbol{\updelta}}_n=\sqrt{n}(\hbeta_{\tau}-\bbeta_{\tau}) \). By @eq-qnt-knight,
\[
\begin{gathered}
Z_n(\boldsymbol{\updelta})=-\boldsymbol{\updelta}\T\bw_n+R_n(\boldsymbol{\updelta}),
\qquad \bw_n=\frac{1}{\sqrt{n}}\sum_i\x_{(i)}\psi_{\tau}(u_i),\\
R_n(\boldsymbol{\updelta})=\sum_i\int_0^{c_i}
   \bigl\{1\{u_i\le s\}-1\{u_i\le0\}\bigr\}\,ds .
\end{gathered}
\]
Because \( F_i \) is continuous at zero and \( F_i(0)=\tau \), the variables
\( \psi_{\tau}(u_i) \) are independent with mean zero and variance \( \tau(1-\tau) \),
and they are bounded. Condition (C3) makes the Lindeberg condition hold for
\( \bw_n \), which by (C2) and the Cramér–Wold device gives
\( \bw_n\xrightarrow{d}\Normal_p\{\bzero,\tau(1-\tau)\bD_0\} \).

For the remainder, \( \E R_n(\boldsymbol{\updelta})=\sum_i\int_0^{c_i}\{F_i(s)-F_i(0)\}\,ds \).
By the equicontinuity in (C1), \( F_i(s)-F_i(0)=f_i(0)s+o(s) \) uniformly in
\( i \), and
\( \max_i\lvert c_i\rvert\le n^{-1/2}\norm{\boldsymbol{\updelta}}\max_i\norm{\x_{(i)}}\to0 \) by
(C3), so
\[
\begin{gathered}
\E R_n(\boldsymbol{\updelta})=\sum_i\Bigl\{\tfrac12f_i(0)c_i^2+o(c_i^2)\Bigr\}\\
=\tfrac12\boldsymbol{\updelta}\T\Bigl\{\tfrac1n\sum_if_i(0)\x_{(i)}\x_{(i)}\T\Bigr\}
   \boldsymbol{\updelta}+o(1)
\longrightarrow \tfrac12\boldsymbol{\updelta}\T\bD_1\boldsymbol{\updelta} .
\end{gathered}
\]
Each summand of \( R_n \) is bounded in absolute value by
\( \lvert c_i\rvert\,1\{\lvert u_i\rvert\le\lvert c_i\rvert\} \), whose second
moment is at most \( c_i^2\cdot 2\lvert c_i\rvert\sup f_i \) once
\( \lvert c_i\rvert<\epsilon \). Summing the independent terms,
\[
\Var R_n(\boldsymbol{\updelta})\le 2\bigl(\sup_i\sup_{\lvert s\rvert<\epsilon}f_i(s)\bigr)
   \max_i\lvert c_i\rvert\sum_ic_i^2\longrightarrow0,
\]
because \( \sum_ic_i^2\to\boldsymbol{\updelta}\T\bD_0\boldsymbol{\updelta} \) and \( \max_i\lvert c_i\rvert\to0 \).
Hence \( R_n(\boldsymbol{\updelta})\to\frac12\boldsymbol{\updelta}\T\bD_1\boldsymbol{\updelta} \) in probability, and for each
fixed \( \boldsymbol{\updelta} \),
\[
Z_n(\boldsymbol{\updelta})\;\xrightarrow{d}\;Z(\boldsymbol{\updelta})=-\boldsymbol{\updelta}\T\bw+\tfrac12\boldsymbol{\updelta}\T\bD_1\boldsymbol{\updelta},
\qquad \bw\sim\Normal_p\{\bzero,\tau(1-\tau)\bD_0\},
\]
jointly over finitely many \( \boldsymbol{\updelta} \), since \( \bw_n \) is the only random part in
the limit. The limit \( Z \) is a strictly convex quadratic with unique minimizer
\( \bD_1^{-1}\bw \).

The last step is the one this book takes on trust. Each \( Z_n \) is convex, so
pointwise convergence in distribution to a limit with a unique minimizer implies
convergence in distribution of the minimizers: this is the convexity lemma of
Pollard (1991), proved there from the fact that convex functions converging
pointwise on a dense set converge uniformly on compacta. Applying it,
\( \hat{\boldsymbol{\updelta}}_n\xrightarrow{d}\bD_1^{-1}\bw \), whose covariance
is @eq-qnt-sandwich. The identically distributed case follows by substituting
\( f_i=f \) and noting that \( s(\tau)=1/f\{F^{-1}(\tau)\} \) is the derivative of
\( F^{-1} \) at \( \tau \).
:::

::: {.warning}
[What the sparsity function does to precision]

In @eq-qnt-iid the factor \( \tau(1-\tau) \) shrinks as \( \tau \) approaches
\( 0 \) or \( 1 \), suggesting that extreme quantiles are estimated *more*
precisely. They are not: \( s(\tau)^2 \) blows up faster, because the tail density is
small. For a standard normal \( \tau(1-\tau)s(\tau)^2 \) is
\( 0.25/\varphi(0)^2\approx1.571 \) at the median and about \( 13.9 \) at
\( \tau=0.99 \). Extreme conditional quantiles need far more data, or a model that
borrows strength across levels.
:::

## Estimating the sparsity function

The iid formula @eq-qnt-iid requires \( f\{F^{-1}(\tau)\} \), a density at a point.
The standard device replaces it by a difference quotient of the empirical quantile
function of the residuals,
\[
\hat s(\tau)=\frac{\hat F_n^{-1}(\tau+h_n)-\hat F_n^{-1}(\tau-h_n)}{2h_n},
\]{#eq-qnt-sparsity}

with a bandwidth \( h_n\to0 \) slowly; Hall and Sheather (1988) derived the choice
\( h_n=n^{-1/3}z_{\alpha}^{2/3}\{1.5\varphi(z_{\tau})^2/(2z_{\tau}^2+1)\}^{1/3} \)
from an Edgeworth expansion for studentized quantiles. The non-identically-distributed
case needs \( \bD_1 \), estimated by replacing \( f_i(0) \) with a kernel estimate of
the density of the \( i \)th residual at zero, in the form Koenker (2005, chapter 4)
attributes to Powell. Both involve a bandwidth and both converge slowly, which is the
reason for the alternatives below.

::: {#exm-qnt-three-ses}
[Three standard errors for one slope]

At \( \tau=0.5 \) the fitted income slope is \( 0.5602 \), and the three routes
give very different precisions:

| method | \( \tau=0.5 \) | \( \tau=0.9 \) |
|---|---|---|
| iid formula @eq-qnt-iid | 0.0121 | 0.0143 |
| sandwich @eq-qnt-sandwich | 0.0308 | 0.0300 |
| pairs bootstrap, \( B=999 \) | 0.0344 | 0.0264 |

At the median the iid formula is optimistic by a factor of \( 2.5 \) against the
sandwich and rather more against the bootstrap, which is no surprise:
@exm-qnt-engel-bands showed that these errors are anything but identically
distributed. The least squares slope has standard error \( 0.0144 \) under the same
false assumption. The estimated sparsity is \( 192.2 \) francs at \( \tau=0.5 \) and
\( 378.9 \) at \( \tau=0.9 \).
:::

```{.python .run #cell-qr-inference-sparsity}
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import linprog
data = sm.datasets.engel.load_pandas().data

def qreg(X, y, tau):
    """Quantile regression by linear programming (as in Section 45.2)."""
    n, p = X.shape
    cost = np.concatenate([np.zeros(2 * p), tau * np.ones(n), (1 - tau) * np.ones(n)])
    A = np.hstack([X, -X, np.eye(n), -np.eye(n)])
    out = linprog(cost, A_eq=A, b_eq=y, bounds=(0, None), method="highs")
    return out.x[:p] - out.x[p:2 * p]

income = data["income"].to_numpy()
food = data["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])

def hall_sheather(n, tau, alpha=0.05):
    """Bandwidth in the tau direction for a difference quotient of the quantile function."""
    z, za = stats.norm.ppf(tau), stats.norm.ppf(1 - alpha / 2)
    return n ** (-1 / 3) * za ** (2 / 3) * (1.5 * stats.norm.pdf(z) ** 2 / (2 * z ** 2 + 1)) ** (1 / 3)

def sparsity(resid, tau, n):
    """Siddiqui's estimate of s(tau) = 1 / f(F^{-1}(tau)): a difference quotient."""
    h = min(hall_sheather(n, tau), tau, 1 - tau)
    hi, lo = np.quantile(resid, [min(tau + h, 1.0), max(tau - h, 0.0)])
    return (hi - lo) / (2 * h)

def se_iid(X, y, tau):
    """Standard errors under independent and identically distributed errors."""
    beta = qreg(X, y, tau)
    s = sparsity(y - X @ beta, tau, len(y))
    return beta, np.sqrt(tau * (1 - tau) * s ** 2 * np.diag(np.linalg.inv(X.T @ X)))

def se_sandwich(X, y, tau):
    """Powell's sandwich: a kernel estimate of the matrix sum_i f_i(0) x_i x_i'."""
    beta = qreg(X, y, tau)
    n = len(y)
    r = y - X @ beta
    h = sparsity(r, tau, n) * hall_sheather(n, tau)    # a window on the residual scale
    dens = np.exp(-0.5 * (r / h) ** 2) / (h * np.sqrt(2 * np.pi))
    D1 = np.linalg.inv(X.T @ (dens[:, None] * X))
    return beta, np.sqrt(tau * (1 - tau) * np.diag(D1 @ (X.T @ X) @ D1))

def se_bootstrap(X, y, tau, B, rng):
    """The pairs bootstrap of chapter 23: resample (y_i, x_i) with replacement and refit."""
    draws = np.array([qreg(X[i], y[i], tau)
                      for i in rng.integers(0, len(y), (B, len(y)))])
    return draws.std(axis=0, ddof=1)

beta50, se50_iid = se_iid(X, food, 0.50)
_, se50_sand = se_sandwich(X, food, 0.50)
se50_small = se_bootstrap(X, food, 0.50, 199, np.random.default_rng(4545))
print(f"slope at tau = 0.5: {beta50[1]:.4f}")
print(f"  se: iid errors {se50_iid[1]:.4f}, sandwich {se50_sand[1]:.4f},"
      f" bootstrap with B = 199 {se50_small[1]:.4f}")
```

## Avoiding the density

Two families of methods sidestep @eq-qnt-sparsity altogether.

**Rank-score tests.** The dual variables of @eq-qnt-lp, described
in @exr-qnt-lp-dual, are a regression analogue of ranks: they depend on the data only
through each observation's position relative to the fitted hyperplanes. Gutenbrunner
and Jurečková (1992) showed that tests built from them are distribution-free in the
sense classical rank tests are, and intervals follow by inverting the test over a
grid. No density is estimated, which is why Koenker (2005, chapter 3) makes them the
default; the cost is a grid search per interval.

**The bootstrap.** Resampling \( (y_i,\x_{(i)}) \) pairs — scheme (a)
of @def-bs-schemes — needs only independence of the pairs, and Hahn (1995) proved that
the bootstrap distribution of \( \sqrt{n}(\hbeta_{\tau}^{*}-\hbeta_{\tau}) \) has
the same limit as that of \( \sqrt{n}(\hbeta_{\tau}-\bbeta_{\tau}) \). Resampling
*residuals* is not valid here, for the reason
[Section 23.1](../ch23-resampling-inference/01-cases-and-residuals.html) gave in the
heteroscedastic linear model: it imposes the identical-distribution structure quantile
regression was chosen to avoid.

::: {#exm-qnt-coverage}
[Which interval keeps its promise]

Generate \( n=150 \) observations from \( y=1+x+(0.2+1.2x)\varepsilon \) with
\( x \) uniform on \( (0,2) \) and \( \varepsilon \) standard normal, so the true
slope of the \( \tau \)-quantile is \( 1+1.2z_{\tau} \): \( 1.000 \) at the median,
\( 2.538 \) at \( \tau=0.9 \). Over \( 800 \) replicates (\( 150 \) for the
bootstrap, \( B=99 \)), the nominal \( 95\% \) intervals cover:

| method | \( \tau=0.5 \) | \( \tau=0.9 \) |
|---|---|---|
| iid formula | 0.899 | 0.911 |
| sandwich | 0.965 | 0.911 |
| pairs bootstrap | 0.933 | 0.940 |

The iid formula undercovers at both levels, as it must. The sandwich is conservative
at the median (mean width \( 0.989 \) against \( 0.743 \)) and no better than the
iid formula at \( \tau=0.9 \), where its kernel estimate of \( \bD_1 \) is fed by few
residuals near zero. The bootstrap is the steadiest, at \( B \) extra linear
programs.
:::

## Quantile curves that cross

Nothing in separate fitting forces the fitted curves to be ordered, yet the objects
they estimate certainly are.

::: {#prp-qnt-crossing}
[Ordering, and how to restore it]

::: {.enumerate options="label=(\alph*)"}
1. Let \( \bar{\x}=n^{-1}\sum_i\x_{(i)} \) and let \( \hbeta_{\tau} \)
   solve @eq-qnt-sample. Then \( \tau\mapsto\bar{\x}\T\hbeta_{\tau} \) is nondecreasing: at
   the average design point the fitted quantiles never cross.

2. Fix \( \x \) and levels \( \tau_1<\dots<\tau_m \), let
   \( q_j=\x\T\hbeta_{\tau_j} \) be the fitted values and let \( q_{(1)}\le\dots\le
   q_{(m)} \) be the same numbers sorted. For any nondecreasing sequence
   \( g_1\le\dots\le g_m \) and any \( p\ge1 \),
   \[
   \sum_{j=1}^m\bigl\lvert q_{(j)}-g_j\bigr\rvert^{p}
   \;\le\;\sum_{j=1}^m\bigl\lvert q_j-g_j\bigr\rvert^{p}.
   \]{#eq-qnt-rearrange}

   Sorting the fitted quantiles therefore moves them no further from the true
   quantile curve, in every \( L^p \) sense at once.
:::

:::

::: {.proof}
(a) Split the check loss as \( \rho_{\tau}(u)=\tau u+u_{-} \), where
\( u_{-}=(-u)_{+} \); this is an identity, since for \( u\ge0 \) both sides are
\( \tau u \) and for \( u<0 \) both are \( \tau u-u \). Summing over the residuals at
\( \bb \),
\[
R_{\tau}(\bb)=\tau\sum_i\bigl(y_i-\x_{(i)}\T\bb\bigr)+S(\bb)
=\tau n\bigl(\bar y-\bar{\x}\T\bb\bigr)+S(\bb),
\qquad S(\bb)=\sum_i\bigl(\x_{(i)}\T\bb-y_i\bigr)_{+}.
\]
Let \( \tau_1<\tau_2 \) and abbreviate \( \hbeta_j=\hbeta_{\tau_j} \). Optimality gives
\( R_{\tau_2}(\hbeta_2)\le R_{\tau_2}(\hbeta_1) \) and
\( R_{\tau_1}(\hbeta_1)\le R_{\tau_1}(\hbeta_2) \). Adding the two inequalities, the
terms \( S(\hbeta_1) \), \( S(\hbeta_2) \) and \( \bar y \) cancel and there remains
\[
-\tau_2\bar{\x}\T\hbeta_2-\tau_1\bar{\x}\T\hbeta_1
\le-\tau_2\bar{\x}\T\hbeta_1-\tau_1\bar{\x}\T\hbeta_2,
\]
that is \( (\tau_2-\tau_1)\bigl(\bar{\x}\T\hbeta_1-\bar{\x}\T\hbeta_2\bigr)\le0 \).
Since \( \tau_2>\tau_1 \), \( \bar{\x}\T\hbeta_1\le\bar{\x}\T\hbeta_2 \).

(b) Sorting a finite sequence is a finite composition of transpositions of adjacent
inverted pairs, so it suffices to show that exchanging one inverted pair does not
increase the sum. Suppose \( i<j \) and \( q_i>q_j \), and write \( \phi(t)=\lvert
t\rvert^{p} \), which is convex. Put \( a=q_j-g_j \), \( b=q_i-g_i \),
\( a'=q_i-g_j \), \( b'=q_j-g_i \). Then \( a'+b'=a+b \), and since
\( g_i\le g_j \) and \( q_j<q_i \) both \( a' \) and \( b' \) lie between
\( \min(a,b) \) and \( \max(a,b) \). Two numbers with the same sum lying inside the
interval spanned by another two can be written as
\( a'=\lambda a+(1-\lambda)b \) and \( b'=(1-\lambda)a+\lambda b \) for some
\( \lambda\in[0,1] \), so convexity gives
\( \phi(a')+\phi(b')\le\phi(a)+\phi(b) \). Exchanging the pair replaces
\( \phi(b)+\phi(a) \) by \( \phi(a')+\phi(b') \) and leaves every other term alone.
:::

Part (b) is the **rearrangement** repair of Chernozhukov, Fernández-Val and Galichon
(2010): fit the levels separately, then sort the fitted values at each \( \x \). It
costs nothing, can only help, and produces a genuine quantile function. Other repairs
impose the ordering during fitting — He (1997) through a location-scale model,
Bondell, Reich and Wang (2010) through linear constraints — at the price of a stronger
model or a larger optimization.

::: {#exm-qnt-crossing-engel}
[Where Engel's lines cross]

The seven fitted lines of @exm-qnt-lp are not parallel, so every neighbouring pair
crosses somewhere: at incomes of \( 252 \), \( 203 \), \( 163 \), \( 228 \),
\( -117 \) and \( 143 \) francs. The smallest income in the data is \( 377 \), so
every crossing lies outside the data, where the model was never
tested ([Figure 45.3.1](03-inference.html#fig-qnt-inference)(b)). Crossing reaches
inside when the predictor is more flexible: adding a quadratic term in income leaves
the curves out of order on \( 0.035 \) of a grid spanning the lowest \( 99\% \) of
incomes, beginning at the very edge. Rearranging repairs the ordering and moves the
curves by at most \( 4.0 \) francs.
:::

::: {when-format="html"}
![**Figure 45.3.1.** (a) The fitted income slope against \( \tau \), with pointwise
\( 95\% \) pairs-bootstrap bands and the least squares slope. (b) The seven fitted
lines extended left of the smallest observed income: the fan closes and
crosses.](qr_inference.svg){#fig-qnt-inference width=100%}
:::

::: {when-format="pdf"}
![(a) The fitted income slope against \( \tau \), with pointwise \( 95\% \)
pairs-bootstrap bands and the least squares slope. (b) The seven fitted lines
extended left of the smallest observed income: the fan closes and
crosses.](qr_inference.pdf){width=100%}
:::

```{.python .run #cell-qr-inference-crossing}
TAUS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
beta = {tau: qreg(X, food, tau) for tau in TAUS}
crossings = {}
for lo, hi in zip(TAUS[:-1], TAUS[1:]):
    d0, d1 = beta[hi] - beta[lo]
    crossings[(lo, hi)] = -d0 / d1                 # where two fitted lines meet
print(f"incomes in the data run from {income.min():.0f} to {income.max():.0f}")
print("neighbouring fitted lines meet at income",
      " ".join(f"{v:.0f}" for v in crossings.values()))
```

## Two extensions, briefly

**A Bayesian reading.** Minimizing @eq-qnt-sample maximizes the likelihood of an
asymmetric Laplace density with location \( \x_{(i)}\T\bbeta_{\tau} \) and asymmetry
\( \tau \), so a flat prior makes the regression quantile a posterior mode and
[Chapter 39](../ch39-glms-in-practice-bayes/index.html) (@def-prc-posterior,
@def-prc-mcmc) applies (Yu and Moyeed 2001; Kozumi and Kobayashi 2011 give a scale
mixture yielding conditionally normal Gibbs steps). Only one \( \tau \) can have the
assumed error law, so the credible intervals deserve the scepticism due to any
misspecified likelihood.

**Nonlinear predictors.** A roughness penalty on the check loss gives additive
quantile regression — a linear program under an \( L^1 \) penalty (Koenker, Ng and
Portnoy 1994), a quadratic program under the penalties of
[Chapter 43](../ch43-smoothing/index.html) — and structured additive
predictors (@def-add-star) carry over. The difficulty is the smoothing
parameter: @thm-sel-loocv needs a quadratic criterion, so cross-validation must
refit, or the problem must be handed to [Section 45.5](05-gamlss.html).

## Exercises

### A. Check your understanding

::: {#exr-qnt-normal-median}
[A1]

Use @eq-qnt-iid to show that for iid normal errors the asymptotic variance of a
regression median exceeds that of least squares by the factor \( \pi/2 \). Which
assumption of @thm-qnt-asymptotics is doing the work?
:::

::: {.solution}
With \( f(0)=1/(\sigma\sqrt{2\pi}) \) and \( \tau=1/2 \), @eq-qnt-iid gives
\( (1/4)\,2\pi\sigma^2\bD_0^{-1}=(\pi/2)\sigma^2\bD_0^{-1} \), against
\( \sigma^2\bD_0^{-1} \) for least squares. The identical-distribution assumption is
what makes the comparison meaningful; under heteroscedasticity the comparison can go
the other way, as @exm-qnt-three-ses shows in reverse.
:::

### B. Practice

::: {#exr-qnt-extreme}
[B1]

For a standard normal error law, compute \( \tau(1-\tau)s(\tau)^2 \) at
\( \tau=0.5,0.9,0.99 \) and confirm the claim in the warning above that precision
deteriorates in the tail. At which \( \tau \) is the factor smallest?
:::

::: {.solution}
\( s(\tau)=1/\varphi\{\Phi^{-1}(\tau)\} \), so the factor is
\( \tau(1-\tau)/\varphi\{\Phi^{-1}(\tau)\}^2 \). At \( \tau=0.5 \) it is
\( 0.25/(0.3989)^2\approx1.571 \); at \( \tau=0.9 \), \( 0.09/(0.1755)^2\approx2.92 \);
at \( \tau=0.99 \), \( 0.0099/(0.0267)^2\approx13.9 \). Differentiating shows the
minimum is at \( \tau=1/2 \) by symmetry.
:::

::: {#exr-qnt-crossing-two-lines}
[B2]

Two fitted quantile lines \( a_1+b_1x \) and \( a_2+b_2x \) with \( b_1<b_2 \) cross at
\( x^{*}=-(a_2-a_1)/(b_2-b_1) \). Using the numbers of @exm-qnt-lp, verify the
crossing income of the \( \tau=0.05 \) and \( \tau=0.10 \) lines, and explain
why @prp-qnt-crossing(a) guarantees that \( x^{*} \) cannot lie at the mean income.
:::

### C. Going deeper

::: {#exr-qnt-knight-mean}
[C1]

Use @eq-qnt-knight to show that the estimating equation for \( \bbeta_{\tau} \) in the
population is \( \E\{\x\,\psi_{\tau}(Y-\x\T\bbeta_{\tau})\}=\bzero \), and that this is
the exact analogue of the normal equations with \( \psi_{\tau} \) in place of the
identity. Compare with the estimating-function framework of
[Chapter 38](../ch38-quasi-likelihood/index.html) (@def-ql-quasi): is \( \psi_{\tau} \)
a quasi-score?
:::

