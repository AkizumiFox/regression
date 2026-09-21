# Expectile regression

The check loss buys freedom from distributional assumptions at the price of a kink,
and the kink makes the estimator awkward: no closed form, no hat matrix, no leverage
shortcut for cross-validation (@thm-sel-loocv), no easy way to attach the quadratic
penalties that [Chapter 30](../ch30-regularization-boosting/index.html) and
[Chapter 43](../ch43-smoothing/index.html) are built on. Newey and Powell (1987)
replaced the absolute deviation by a squared one, keeping the asymmetry. The answer
is a family of location measures called expectiles, fitted by weighted least
squares.

## Asymmetric squared loss

::: {#def-qnt-expectile}
[Expectiles and expectile regression]

For \( \tau\in(0,1) \) define the **asymmetric squared loss**
\[
\rho^{E}_{\tau}(u)=\bigl\lvert\tau-1\{u<0\}\bigr\rvert\,u^{2}
=\begin{cases}\tau u^{2},&u\ge0,\\(1-\tau)u^{2},&u<0.\end{cases}
\]{#eq-qnt-expectile-loss}

For a random variable \( Y \) with \( \E Y^{2}<\infty \), the \( \tau \)**-expectile**
\( \mu_{\tau} \) is a minimizer of \( S(\mu)=\E\rho^{E}_{\tau}(Y-\mu) \). The
**expectile regression** model at level \( \tau \) asserts that the \( \tau \)-expectile
of the conditional law of \( Y \) is \( \x\T\bbeta_{\tau} \), and the estimate is
\[
\hbeta^{E}_{\tau}\in\argmin_{\bb}\;\sum_{i=1}^n
   \rho^{E}_{\tau}\bigl(y_i-\x_{(i)}\T\bb\bigr).
\]{#eq-qnt-expectile-sample}

:::

The weight \( \tau \) applies above the fit and \( 1-\tau \) below, as in the check
loss; only the power has changed, and at \( \tau=1/2 \)
problem @eq-qnt-expectile-sample is least squares.

## What an expectile is

::: {#prp-qnt-expectile}
[Properties of expectiles]

Let \( \E Y^{2}<\infty \) and let \( Y \) be non-degenerate.

::: {.enumerate options="label=(\alph*)"}
1. **Existence, uniqueness, monotonicity.** \( S \) is finite, differentiable and
   strictly convex, and \( \mu_{\tau} \) is its unique stationary point, characterized
   by
   \[
   \tau\,\E(Y-\mu_{\tau})_{+}=(1-\tau)\,\E(\mu_{\tau}-Y)_{+}.
   \]{#eq-qnt-expectile-equation}

   The map \( \tau\mapsto\mu_{\tau} \) is strictly increasing on \( (0,1) \).

2. **The mean is an expectile.** \( \mu_{1/2}=\E Y \).

3. **The level of a point.** For every \( \mu \) with \( P(Y<\mu)>0 \) and
   \( P(Y>\mu)>0 \) — equivalently, every \( \mu \) interior to the convex hull of
   the support of \( Y \) — there is exactly one \( \tau \) with
   \( \mu_{\tau}=\mu \), namely
   \[
   \tau=\frac{\E(\mu-Y)_{+}}{\E\lvert Y-\mu\rvert}.
   \]{#eq-qnt-expectile-level}

4. **Equivariance.** For \( a>0 \) and \( b\in\Real \),
   \( \mu_{\tau}(aY+b)=a\mu_{\tau}(Y)+b \), and \( \mu_{\tau}(-Y)=-\mu_{1-\tau}(Y) \).
   Expectiles are *not* equivariant under nonlinear increasing transformations.
:::

:::

::: {.proof}
(a) Write \( S(\mu)=\tau\E(Y-\mu)_{+}^{2}+(1-\tau)\E(\mu-Y)_{+}^{2} \), finite
because \( (y-\mu)^2\le2y^2+2\mu^2 \). The loss is continuously differentiable with
\( \rho^{E\prime}_{\tau}(u)=2u\lvert\tau-1\{u<0\}\rvert \), continuous at the
origin because both one-sided derivatives vanish, and nondecreasing, so
\( \rho^{E}_{\tau} \) is convex. Dominated convergence permits differentiation under
the expectation, giving
\[
\tfrac12 S'(\mu)=-\tau\E(Y-\mu)_{+}+(1-\tau)\E(\mu-Y)_{+}=:-H(\mu).
\]
Now \( H \) is continuous, and for \( \mu<\mu' \),
\[
\begin{gathered}
H(\mu)-H(\mu')=\tau\,\E\bigl\{(Y-\mu)_{+}-(Y-\mu')_{+}\bigr\}\\
   +\;(1-\tau)\,\E\bigl\{(\mu'-Y)_{+}-(\mu-Y)_{+}\bigr\}\ge0,
\end{gathered}
\]
each expectation being nonnegative because \( t\mapsto(y-t)_{+} \) is nonincreasing
and \( t\mapsto(t-y)_{+} \) nondecreasing, and strict unless \( Y \) puts no mass
outside \( [\mu,\mu'] \) on either side, which a non-degenerate \( Y \) cannot do for
every pair. Moreover, as \( \mu\to-\infty \) monotone convergence
gives \( \E(\mu-Y)_{+}\downarrow0 \), while
\( \E(Y-\mu)_{+}\ge\E Y-\mu\to+\infty \), so \( H(\mu)\to+\infty \); the same
argument reflected gives \( H(\mu)\to-\infty \) as \( \mu\to+\infty \). So
\( H \) has a unique zero, \( S \) is strictly convex there,
and @eq-qnt-expectile-equation is exactly \( H(\mu_{\tau})=0 \). For monotonicity, let
\( \tau<\tau' \) and write \( H_{\tau} \) for the function above. Then
\[
\begin{gathered}
H_{\tau'}(\mu_{\tau})=H_{\tau}(\mu_{\tau})
   +(\tau'-\tau)\,\E\lvert Y-\mu_{\tau}\rvert\\
=(\tau'-\tau)\,\E\lvert Y-\mu_{\tau}\rvert>0,
\end{gathered}
\]
using \( \E(Y-\mu)_{+}+\E(\mu-Y)_{+}=\E\lvert Y-\mu\rvert \); since \( H_{\tau'} \) is
decreasing with a unique zero, that zero \( \mu_{\tau'} \) lies to the right of
\( \mu_{\tau} \).

(b) At \( \tau=1/2 \), @eq-qnt-expectile-equation reads
\( \E(Y-\mu)_{+}=\E(\mu-Y)_{+} \), that is \( \E(Y-\mu)=0 \).

(c) Put \( N=\E(Y-\mu)_{+} \) and \( P=\E(\mu-Y)_{+} \); both are strictly positive
under the stated condition on \( \mu \), so \( N+P=\E\lvert Y-\mu\rvert>0 \).
Equation @eq-qnt-expectile-equation says \( \tau N=(1-\tau)P \), hence
\( \tau(N+P)=P \), and @eq-qnt-expectile-level follows. Conversely, given such a
\( \mu \), the ratio \( P/(N+P) \) lies strictly between \( 0 \) and \( 1 \) and
solves the equation, and by (a) the solution for that \( \tau \) is unique.

(d) For \( a>0 \), \( \rho^{E}_{\tau}(au)=a^{2}\rho^{E}_{\tau}(u) \), so
\( \E\rho^{E}_{\tau}(aY+b-\mu)=a^{2}\E\rho^{E}_{\tau}\{Y-(\mu-b)/a\} \) and the
minimizers correspond as claimed; \( \rho^{E}_{\tau}(-u)=\rho^{E}_{1-\tau}(u) \) gives
the reflection. For the failure under nonlinear maps it is enough to look at
\( \tau=1/2 \): by (b) the \( 1/2 \)-expectile is the mean, and for strictly convex
\( h \) Jensen's inequality gives \( \E h(Y)>h(\E Y) \), so
\( \mu_{1/2}\{h(Y)\}\ne h\{\mu_{1/2}(Y)\} \).
:::

Part (c) says what level an expectile sits at. Compare the two characterizations:
\[
\tau=\frac{\E\,1\{Y\le q_{\tau}\}}{\E\,1}\quad\text{(quantile)},
\qquad
\tau=\frac{\E(\mu_{\tau}-Y)_{+}}{\E\lvert Y-\mu_{\tau}\rvert}\quad\text{(expectile)} .
\]
A quantile counts the mass below the point; an expectile weighs it, by distance. So an
expectile depends on the *entire* distribution and has none of the robustness
of @def-qnt-quantile: one response moved to \( +10^{6} \) moves every expectile. In
exchange it is smooth.

::: {.remark}
[Expectiles as quantiles of another law]

Since @eq-qnt-expectile-level determines \( \tau \) from \( \mu \) and \( F \), the
expectiles of \( F \) are the quantiles of another distribution built from \( F \) by
partial moments — the observation of Jones (1994). The catch is that the relabelling
depends on \( F \): the \( \tau \)-quantile is at level \( \tau \) for every
distribution, while the \( \tau \)-expectile is at a level only the distribution can
tell you.
:::

::: {#exm-qnt-expectile-levels}
[Which quantile is the nine-tenths expectile?]

For a standard normal the \( 0.9 \)-expectile sits at the \( 0.806 \)-quantile and
the \( 0.5 \)-expectile at the \( 0.500 \)-quantile, as symmetry requires. For a
lognormal with log-scale \( 0.8 \) the \( 0.9 \)-expectile sits at the
\( 0.892 \)-quantile: the heavy right tail drags it further out. For the residuals of
the least squares fit to Engel's budgets the \( 0.9 \)-expectile sits at the
\( 0.851 \)-quantile and the \( 0.1 \)-expectile at the \( 0.170 \)-quantile
([Figure 45.4.1](04-expectiles.html#fig-qnt-engel-expectiles)(b)). Three
distributions, three answers to the same question.
:::

## Fitting by weighted least squares

::: {#prp-qnt-iwls}
[The estimating equation and its iteration]

Let \( \rank(\X)=p \) and \( Q(\bb)=\sum_i\rho^{E}_{\tau}(y_i-\x_{(i)}\T\bb) \).

::: {.enumerate options="label=(\alph*)"}
1. \( Q \) is continuously differentiable and strictly convex, so it has a unique
   minimizer, characterized by
   \[
   \X\T\W(\bb)\,(\y-\X\bb)=\bzero,\qquad
   \W(\bb)=\diag\bigl\{w_i(\bb)\bigr\},\quad
   w_i(\bb)=\bigl\lvert\tau-1\{y_i<\x_{(i)}\T\bb\}\bigr\rvert .
   \]{#eq-qnt-expectile-normal}

2. The **iterated weighted least squares** step
   \( \bb^{(t+1)}=\{\X\T\W(\bb^{(t)})\X\}^{-1}\X\T\W(\bb^{(t)})\y \) has the property
   that if \( \W(\bb^{(t+1)})=\W(\bb^{(t)}) \), then \( \bb^{(t+1)} \) is the exact
   global minimizer of \( Q \).
:::

:::

::: {.proof}
(a) The scalar loss is continuously differentiable with derivative
\( 2u\lvert\tau-1\{u<0\}\rvert \), which is strictly increasing, so
\( \rho^{E}_{\tau} \) is strictly convex; composing with the affine map
\( \bb\mapsto y_i-\x_{(i)}\T\bb \) and summing gives a convex \( Q \), strictly convex
because \( \rank(\X)=p \) makes \( \bb\mapsto\X\bb \) injective. It is also
coercive: \( \rho^{E}_{\tau}(u)\ge\min(\tau,1-\tau)u^{2} \) gives
\( Q(\bb)\ge\min(\tau,1-\tau)\norm{\y-\X\bb}^{2} \), which tends to infinity as
\( \norm{\bb}\to\infty \) because \( \X \) has full column rank, so a minimizer
exists, and strict convexity makes it unique. The chain rule gives
\( \nabla Q(\bb)=-2\X\T\W(\bb)(\y-\X\bb) \), and for a differentiable convex function
a stationary point is a global minimum.

(b) The weighted least squares step solves \( \X\T\W_t\X\bb=\X\T\W_t\y \) with
\( \W_t=\W(\bb^{(t)}) \), so \( \X\T\W_t(\y-\X\bb^{(t+1)})=\bzero \). If
\( \W(\bb^{(t+1)})=\W_t \), this is @eq-qnt-expectile-normal at \( \bb^{(t+1)} \), so
\( \nabla Q(\bb^{(t+1)})=\bzero \) and by (a) the point is the global minimizer.
:::

::: {.warning}
[What is not proved here]

Part (b) says the iteration is *correct when it stops*, not that it stops. The
weights depend on \( \bb \) only through the signs of the residuals, so the step is a
Newton step for a piecewise linear system with finitely many weight patterns; global
convergence is not proved in this book. In practice it is fast: for the seven fits
below none needed more than \( 6 \) reweightings from a least squares start.
Compare @prp-res-irls, whose majorization argument does give a descent guarantee but
requires \( \rho(u)=g(u^2) \) with \( g \) concave, which @eq-qnt-expectile-loss is
not.
:::

::: {#exm-qnt-engel-expectiles}
[Expectile lines for Engel's budgets]

The same \( 235 \) households, fitted at the same seven levels
by @eq-qnt-expectile-sample:

| \( \tau \) | 0.05 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 0.95 |
|---|---|---|---|---|---|---|---|
| slope | 0.3617 | 0.3830 | 0.4293 | 0.4852 | 0.5454 | 0.6017 | 0.6322 |
| intercept | 157.9 | 162.6 | 159.0 | 147.5 | 128.8 | 109.0 | 99.8 |

The level \( \tau=1/2 \) reproduces least squares exactly, slope \( 0.4852 \). The
fan opens in the same direction as the quantile fan of @exm-qnt-lp but less widely,
from \( 0.3617 \) to \( 0.6322 \) against \( 0.3434 \) to \( 0.7091 \) — not a
discrepancy, since by @exm-qnt-expectile-levels the level \( 0.9 \) means a more
central position for an expectile.
:::

::: {when-format="html"}
![**Figure 45.4.1.** (a) Seven fitted expectile lines for Engel's budgets. (b) The
quantile level at which the \( \tau \)-expectile sits, for three laws; the diagonal
is where the levels agree.](engel_expectiles.svg){#fig-qnt-engel-expectiles width=100%}
:::

::: {when-format="pdf"}
![(a) Seven fitted expectile lines for Engel's budgets. (b) The quantile level at
which the \( \tau \)-expectile sits, for three laws; the diagonal is where the levels
agree.](engel_expectiles.pdf){width=100%}
:::

```{.python .run #cell-expectiles-iwls}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.engel.load_pandas().data
TAUS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

income = data["income"].to_numpy()
food = data["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])

def expectile_reg(X, y, tau, tol=1e-12, max_iter=100):
    """Expectile regression by iterated weighted least squares.

    The weight of observation i is tau if it sits above the current fit and
    1 - tau if it sits below; each step is a weighted least squares fit.
    """
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    for step in range(1, max_iter + 1):
        w = np.where(y >= X @ beta, tau, 1 - tau)
        new, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
        if np.max(np.abs(new - beta)) <= tol * (1 + np.max(np.abs(beta))):
            return new, step
        beta = new
    raise RuntimeError("iterated weighted least squares did not converge")

fits = {tau: expectile_reg(X, food, tau) for tau in TAUS}
for tau in (0.10, 0.50, 0.90):
    beta, steps = fits[tau]
    print(f"tau = {tau:.2f}: intercept {beta[0]:7.2f}  slope {beta[1]:.4f}  ({steps} steps)")
```

```{.python .run #cell-expectiles-level}
def expectile(y, tau):
    """The sample tau-expectile of a vector: the scalar case of the same iteration."""
    m = y.mean()
    for _ in range(200):
        w = np.where(y >= m, tau, 1 - tau)
        new = (w * y).sum() / w.sum()
        if abs(new - m) < 1e-12 * (1 + abs(m)):
            return new
        m = new
    return m

def level(y, mu):
    """tau = E(mu - Y)_+ / E|Y - mu|: the expectile level of the point mu."""
    return np.mean(np.maximum(mu - y, 0)) / np.mean(np.abs(y - mu))

ols, _ = np.linalg.lstsq(X, food, rcond=None)[:2]
resid = food - X @ ols
for tau in (0.1, 0.5, 0.9):
    assert abs(level(resid, expectile(resid, tau)) - tau) < 1e-9
print("the identity tau = E(mu - Y)_+ / E|Y - mu| holds at every level tried")
```

## When to prefer an expectile

The case for expectiles is about what a differentiable loss makes possible.

**Penalties and smoothing.** A penalized spline turns @eq-qnt-expectile-sample into
penalized weighted least squares, which is @def-reg-penalized with weights. Every
step is a linear smoother, so the effective degrees of freedom is a trace and
generalized cross-validation @def-sel-gcv can choose the smoothing parameter through
the leverage shortcut @thm-sel-loocv without refitting. Schnabel and Eilers (2009)
built expectile smoothing on this, and Sobotka and Kneib (2012) extended it to
structured additive predictors, where the smoothing parameters go to the mixed-model
representation (@thm-smo-mixed, using @def-mix-reml). None of this is available for
the check loss without approximation.

**Efficiency and risk.** The expectile estimator uses the magnitude of every
residual, the quantile estimator only its sign, whose asymptotic variance carries the
sparsity factor of @eq-qnt-iid. Which is the more precise depends on the error law,
and the comparison is meaningful only at matched levels, since at a common
\( \tau \) the two estimate different functionals; Newey and Powell (1987, §3) give
the asymptotic variance from which such comparisons are made. A functional is
*elicitable* if it minimizes an expected loss — the mean by squared error, quantiles
by the check loss (@thm-qnt-check), expectiles by @eq-qnt-expectile-loss — and
Ziegel (2016) showed that the \( \tau \)-expectiles with \( \tau\ge1/2 \) are the
only law-invariant coherent risk measures that are also elicitable, which is what put
them into financial risk management.

The case against is short and serious: nobody can say in one sentence what the
\( 0.9 \)-expectile is. "Ninety per cent of households fall below this line" a reader
understands; "the lower partial moment at this line is nine times the upper one" is
not.

## Exercises

### A. Check your understanding

::: {#exr-qnt-expectile-two-point}
[A1]

Let \( Y \) take the values \( 0 \) and \( 1 \) with probability \( 1/2 \) each.
Compute \( \mu_{\tau} \) from @eq-qnt-expectile-equation as a function of \( \tau \),
and compare with the \( \tau \)-quantile.
:::

::: {.solution}
For \( \mu\in(0,1) \), \( \E(Y-\mu)_{+}=(1-\mu)/2 \) and \( \E(\mu-Y)_{+}=\mu/2 \),
so @eq-qnt-expectile-equation gives \( \tau(1-\mu)=(1-\tau)\mu \), that is
\( \mu_{\tau}=\tau \). The expectile moves continuously from \( 0 \) to \( 1 \), while
the \( \tau \)-quantile jumps from \( 0 \) to \( 1 \) at \( \tau=1/2 \). Expectiles
smooth over the gaps of a discrete law.
:::

### B. Practice

::: {#exr-qnt-expectile-normal}
[B1]

For \( Y\sim\Normal(0,1) \), show that @eq-qnt-expectile-equation becomes
\( \tau\{\varphi(\mu)-\mu(1-\Phi(\mu))\}=(1-\tau)\{\varphi(\mu)+\mu\Phi(\mu)\} \), and
verify numerically that \( \mu_{0.9}=0.8616 \) lies at the
\( 0.806 \)-quantile, as reported in @exm-qnt-expectile-levels.
:::

::: {.solution}
For the standard normal, \( \E(Y-\mu)_{+}=\varphi(\mu)-\mu\{1-\Phi(\mu)\} \) and
\( \E(\mu-Y)_{+}=\varphi(\mu)+\mu\Phi(\mu) \), both by integration by parts.
Substituting into @eq-qnt-expectile-equation gives the stated identity; solving it
numerically at \( \tau=0.9 \) gives \( \mu_{0.9}=0.8616 \), and
\( \Phi(0.8616)=0.806 \).
:::

::: {#exr-qnt-expectile-symmetry}
[B2]

Show that if the law of \( Y \) is symmetric about \( m \), then
\( \mu_{\tau}+\mu_{1-\tau}=2m \), and that the \( \tau \)-expectile lies at the
\( \alpha(\tau) \)-quantile with \( \alpha(1/2)=1/2 \). Does symmetry force
\( \alpha(\tau)=\tau \)?
:::

### C. Going deeper

::: {#exr-qnt-expectile-heavy}
[C1]

Let \( Y \) have a Cauchy distribution. Which of the objects in this chapter still
exist: the conditional mean, the \( \tau \)-quantiles, the \( \tau \)-expectiles? What
does this say about the price of the smooth loss?
:::

::: {.solution}
Quantiles exist and are finite for every \( \tau\in(0,1) \),
since @thm-qnt-check(a) needs only a distribution function (and the remark after it
removes even the moment condition). The mean does not exist, so neither does the
\( 1/2 \)-expectile; indeed @def-qnt-expectile assumes \( \E Y^{2}<\infty \), and
\( \E\lvert Y\rvert=\infty \) already destroys @eq-qnt-expectile-equation, whose two
sides are both infinite. The smooth loss buys computation at the cost of moments.
:::
