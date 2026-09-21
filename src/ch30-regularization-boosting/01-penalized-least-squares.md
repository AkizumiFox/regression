# Penalized least squares in general

Ridge regression and the lasso of [Chapter 27](../ch27-shrinkage/index.html) raise the same questions: is
there a minimizer, is it unique, and how do we recognize one? This section answers them for a general, mainly convex,
penalty, reads a penalty as a prior, and ends with penalties that give up convexity on purpose. Throughout, \( \X \) is \( n\times p \) of any rank, and \( p \) may
exceed \( n \).

## The criterion

::: {#def-reg-penalized}
[Penalized least squares]

A **penalty** is a function \( P:\Real^p\to[0,\infty) \). For \( \lambda\ge0 \), the **penalized least
squares criterion** is
\[
Q_\lambda(\bb)=\tfrac12\norm{\y-\X\bb}^2+\lambda P(\bb),
\]{#eq-reg-pls}

and a **penalized least squares estimate** is any minimizer \( \hbeta_\lambda \) of \( Q_\lambda \). The penalty
is **convex** if \( P(\alpha\bb+(1-\alpha)\mathbf{c})\le\alpha P(\bb)+(1-\alpha)P(\mathbf{c}) \) for all
\( \bb,\mathbf{c} \) and \( \alpha\in(0,1) \), and **strictly convex** if the inequality is strict whenever
\( \bb\ne\mathbf{c} \). The penalization is **partial** if \( \bb=(\bb_1\T,\bb_2\T)\T \) and \( P \) depends on
\( \bb_2 \) alone; the coefficients in \( \bb_1 \) are then **unpenalized**.
:::

Other scalings (no \( \tfrac12 \), or a divisor \( n \)) multiply \( \lambda \) by a constant; results are quoted here in
the scaling of @eq-reg-pls. Ridge (\( P=\tfrac12\norm{\bb}^2 \), @thm-shr-ridge) and the lasso (\( P=\norm{\bb}_1 \), @def-shr-lasso)
are convex; the best-subset penalty \( \#\{j:b_j\neq0\} \) of [Chapter 29](../ch29-model-selection/index.html) and the bridge
penalties \( \sum_j\lvert b_j\rvert^q \), \( q<1 \), are not. Sections [30.2](02-elastic-net.html) and [30.3](03-group-penalties.html)
add penalties with structure.

**Intercept and scaling.** An intercept is not penalized. Minimizing over it first gives \( b_0=\bar y-\bar{\x}\T\bb \)
and leaves @eq-reg-pls with \( \y \) and the columns of \( \X \) centred (@exr-reg-intercept), the partialling-out of
@thm-proj-fwl, so we centre and drop it. A penalty that treats coordinates alike is not equivariant under rescaling a column,
so the columns are also scaled to a common length, a choice that is part of the model.

## Existence and uniqueness

For least squares a minimizer always exists and its fit is unique (@thm-proj-ls-projection). For penalized criteria:

::: {#thm-reg-existence}
[Existence and uniqueness of penalized estimates]

Let \( P \) be continuous and \( \lambda\ge0 \).

::: {.enumerate options="label=(\alph*)"}
1. If \( Q_\lambda(\bb)\to\infty \) as \( \norm{\bb}\to\infty \), a minimizer exists. This holds in each of three
   cases: (i) \( \X \) has full column rank; (ii) \( \lambda>0 \) and \( P(\bb)\ge c\norm{\bb} \) for
   some \( c>0 \); (iii) the penalization is partial, \( \X=[\X_1,\X_2] \) with \( \X_1 \) of full column rank,
   \( \lambda>0 \), and \( P(\bb)\ge c\norm{\bb_2} \) for some \( c>0 \).

2. If \( P \) is convex, the set of minimizers is convex, and all minimizers have the same fitted
   vector \( \X\hbeta_\lambda \). If in addition \( \lambda>0 \), they all have the same penalty \( P(\hbeta_\lambda) \).

3. If \( P \) is convex and either \( \X \) has full column rank, or \( \lambda>0 \) and \( P \) is strictly convex, the
   minimizer is unique.
:::

:::

::: {.proof}
(a) Let \( \mathcal S=\{\bb:Q_\lambda(\bb)\le Q_\lambda(\bzero)\} \). It is closed because \( Q_\lambda \) is continuous,
and bounded because \( Q_\lambda\to\infty \), so \( Q_\lambda \) attains its minimum on \( \mathcal S \). Any point
outside \( \mathcal S \) has a larger value than \( \bzero \), so this minimum is a global minimizer.
In case (i), the smallest singular value \( s \) of \( \X \) is positive (@thm-mat-svd), so
\( \norm{\y-\X\bb}\ge s\norm{\bb}-\norm{\y} \), which tends to infinity. In case (ii), \( Q_\lambda(\bb)\ge\lambda c\norm{\bb} \).
In case (iii), suppose \( \norm{\bb^{(k)}}\to\infty \) along a sequence on which \( Q_\lambda \) stays bounded.
Then \( \lambda c\norm{\bb_2^{(k)}} \) is bounded, so \( \norm{\bb_1^{(k)}}\to\infty \). With \( s_1>0 \) the smallest
singular value of \( \X_1 \),
\( \norm{\y-\X\bb^{(k)}}\ge s_1\norm{\bb_1^{(k)}}-\norm{\y}-\norm{\X_2\bb_2^{(k)}}\to\infty \), a contradiction.

(b) If \( P \) is convex, so is \( Q_\lambda \), as a sum of convex functions, and the set of minimizers of a convex
function is convex. Let \( \bb \) and \( \mathbf{c} \) be minimizers with value \( q \), and put \( \mathbf{m}=\tfrac12(\bb+\mathbf{c}) \).
With \( \mathbf{u}=\y-\X\bb \) and \( \mathbf{v}=\y-\X\mathbf{c} \), the parallelogram identity gives
\[
\norm{\y-\X\mathbf{m}}^2=\bigl\lVert\tfrac12(\mathbf{u}+\mathbf{v})\bigr\rVert^2
=\tfrac12\norm{\mathbf{u}}^2+\tfrac12\norm{\mathbf{v}}^2-\tfrac14\norm{\X(\bb-\mathbf{c})}^2 .
\]
Together with \( P(\mathbf{m})\le\tfrac12P(\bb)+\tfrac12P(\mathbf{c}) \) this yields
\( q\le Q_\lambda(\mathbf{m})\le q-\tfrac18\norm{\X(\bb-\mathbf{c})}^2 \). Hence \( \X\bb=\X\mathbf{c} \). The two minimizers
then have equal residual sums of squares, so when \( \lambda>0 \) their penalties are equal too.

(c) If \( \X \) has full column rank, \( \X\bb=\X\mathbf{c} \) forces \( \bb=\mathbf{c} \). If \( P \) is strictly convex and
\( \bb\neq\mathbf{c} \), the inequality for \( P(\mathbf{m}) \) in (b) is strict and \( Q_\lambda(\mathbf{m})<q \), which is
impossible.
:::

Norms satisfy (a)(ii), so the lasso, group lasso and elastic net always have minimizers. By (b), when \( p>n \) the lasso's
minimizers, possibly many, share their fit and \( \ell_1 \) norm (they are unique in general position; Tibshirani, 2013).

## Optimality conditions

The corner of \( \lvert b\rvert \) at zero, which makes the lasso set coefficients to zero, calls for a substitute
for the gradient.

::: {#def-reg-subgradient}
[Subgradient]

Let \( P \) be convex. A vector \( \mathbf{s} \) is a **subgradient** of \( P \) at \( \bb \) if
\( P(\mathbf{c})\ge P(\bb)+\mathbf{s}\T(\mathbf{c}-\bb) \) for every \( \mathbf{c}\in\Real^p \). The set of all subgradients
at \( \bb \) is the **subdifferential** \( \partial P(\bb) \).
:::

Where \( P \) is differentiable, \( \partial P(\bb)=\{\nabla P(\bb)\} \). For \( P(b)=\lvert b\rvert \) on \( \Real \), the subdifferential is
\( \{\operatorname{sign}(b)\} \) if \( b\neq0 \) and \( [-1,1] \) at \( 0 \).

::: {#lem-reg-optimality}
[Optimality condition]

Let \( P \) be convex and \( \lambda>0 \). A vector \( \hbeta \) minimizes \( Q_\lambda \) iff
\[
\X\T(\y-\X\hbeta)\in\lambda\,\partial P(\hbeta).
\]{#eq-reg-kkt}

In that case \( Q_\lambda(\bb)-Q_\lambda(\hbeta)\ge\tfrac12\norm{\X(\bb-\hbeta)}^2 \) for every \( \bb \).
:::

::: {.proof}
Write \( \br=\y-\X\hbeta \) and \( \mathbf{d}=\bb-\hbeta \). Expanding the square,
\[
Q_\lambda(\bb)-Q_\lambda(\hbeta)=-\mathbf{d}\T\X\T\br+\tfrac12\norm{\X\mathbf{d}}^2+\lambda\bigl(P(\bb)-P(\hbeta)\bigr).
\]{#eq-reg-expansion}

If \( \X\T\br=\lambda\mathbf{s} \) with \( \mathbf{s}\in\partial P(\hbeta) \), then
\( \lambda(P(\bb)-P(\hbeta))\ge\lambda\mathbf{s}\T\mathbf{d}=\mathbf{d}\T\X\T\br \), and @eq-reg-expansion gives the stated lower bound,
which is nonnegative. Conversely, let \( \hbeta \) be a minimizer and fix \( \bb \). For \( 0<\tau<1 \), convexity gives
\( P(\hbeta+\tau\mathbf{d})-P(\hbeta)\le\tau(P(\bb)-P(\hbeta)) \). Apply @eq-reg-expansion to \( \hbeta+\tau\mathbf{d} \) in
place of \( \bb \):
\[
\begin{aligned}
0&\le Q_\lambda(\hbeta+\tau\mathbf{d})-Q_\lambda(\hbeta)\\
&\le-\tau\mathbf{d}\T\X\T\br+\tfrac{\tau^2}2\norm{\X\mathbf{d}}^2+\lambda\tau\bigl(P(\bb)-P(\hbeta)\bigr).
\end{aligned}
\]
Divide by \( \tau \) and let \( \tau\to0 \). The result is \( \lambda(P(\bb)-P(\hbeta))\ge\mathbf{d}\T\X\T\br \) for every \( \bb \),
which says that \( \X\T\br/\lambda\in\partial P(\hbeta) \).
:::

Condition @eq-reg-kkt is the **Karush–Kuhn–Tucker** (KKT) condition; with \( \lambda=0 \) it is the normal
equations. Most penalties below are sums of Euclidean norms of blocks.

::: {#lem-reg-norm-subgradient}
[Subdifferential of a sum of block norms]

Let \( G_1,\dots,G_K \) partition \( \{1,\dots,p\} \), let \( w_1,\dots,w_K>0 \), and let
\( P(\bb)=\sum_kw_k\norm{\bb_{G_k}} \), where \( \bb_{G_k} \) is the subvector of \( \bb \) with coordinates in
\( G_k \). Then \( \mathbf{s}\in\partial P(\bb) \) iff, for every \( k \),
\[
\mathbf{s}_{G_k}=w_k\,\frac{\bb_{G_k}}{\norm{\bb_{G_k}}}\quad\text{if }\bb_{G_k}\neq\bzero,
\qquad \norm{\mathbf{s}_{G_k}}\le w_k\quad\text{if }\bb_{G_k}=\bzero .
\]
With singleton blocks and unit weights, \( P=\norm{\cdot}_1 \), and \( \mathbf{s}\in\partial\norm{\bb}_1 \) iff
\( s_j=\operatorname{sign}(b_j) \) when \( b_j\neq0 \) and \( \lvert s_j\rvert\le1 \) when \( b_j=0 \).
:::

::: {.proof}
First, \( \mathbf{s}\in\partial P(\bb) \) iff \( \mathbf{s}_{G_k}\in\partial(w_k\norm{\cdot})(\bb_{G_k}) \) for every \( k \). Summing the
blockwise inequalities gives one direction. For the other, apply the subgradient inequality to a
\( \mathbf{c} \) that differs from \( \bb \) only in block \( k \). So it suffices to describe the subdifferential of
\( f(\bu)=w\norm{\bu} \) on \( \Real^d \). We show that \( \mathbf{s}\in\partial f(\bu) \) iff \( \norm{\mathbf{s}}\le w \) and
\( \mathbf{s}\T\bu=w\norm{\bu} \).

If these two conditions hold, then for every \( \bv \), by Cauchy–Schwarz (@prp-mat-cauchy-schwarz),
\( f(\bv)=w\norm{\bv}\ge\mathbf{s}\T\bv=f(\bu)+\mathbf{s}\T(\bv-\bu) \). Conversely, let \( \mathbf{s}\in\partial f(\bu) \).
Taking \( \bv=\bzero \) and \( \bv=2\bu \) gives \( \mathbf{s}\T\bu\ge w\norm{\bu} \) and \( \mathbf{s}\T\bu\le w\norm{\bu} \). Then
\( w\norm{\bu+\tau\bv}\ge w\norm{\bu}+\tau\mathbf{s}\T\bv \) for all \( \tau>0 \). The triangle inequality bounds the
left side by \( w\norm{\bu}+\tau w\norm{\bv} \), so \( \mathbf{s}\T\bv\le w\norm{\bv} \) for every \( \bv \). Taking \( \bv=\mathbf{s} \) gives
\( \norm{\mathbf{s}}\le w \).

If \( \bu=\bzero \), the two conditions say \( \norm{\mathbf{s}}\le w \). If \( \bu\neq\bzero \), then \( \mathbf{s}\T\bu=w\norm{\bu}\ge\norm{\mathbf{s}}\norm{\bu} \)
is equality in Cauchy–Schwarz, so \( \mathbf{s} \) is a nonnegative multiple of \( \bu \) with norm \( w \), that is,
\( \mathbf{s}=w\bu/\norm{\bu} \).
:::

Combining the two lemmas, \( \hbeta \) is a lasso estimate iff \( \x_j\T(\y-\X\hbeta)=\lambda\operatorname{sign}(\hat\beta_j) \)
for every \( j \) with \( \hat\beta_j\neq0 \) and \( \lvert\x_j\T(\y-\X\hbeta)\rvert\le\lambda \) for every \( j \) with
\( \hat\beta_j=0 \), as in @def-shr-lasso: active regressors share the absolute correlation \( \lambda \) with the residual,
and no inactive one exceeds it.

## The Bayesian reading

Ridge regression is a posterior mean (@prp-opt-shrinkage(d)). For other penalties the matching statement is about
the posterior *mode*; part (b) below restates, for completeness, the reading of the lasso found in @exr-shr-laplace
([Section 27.4](../ch27-shrinkage/04-lasso.html)).

::: {#prp-reg-map}
[Penalized estimates as posterior modes]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \sigma^2 \) known, and give \( \bbeta \) a prior density
\( \pi(\bb)\propto\exp\{-(\lambda/\sigma^2)P(\bb)\} \), assumed integrable. Then the posterior density of
\( \bbeta \) is proportional to \( \exp\{-Q_\lambda(\bb)/\sigma^2\} \). So the penalized least squares estimates are
exactly the posterior modes. In particular:

::: {.enumerate options="label=(\alph*)"}
1. a normal prior \( \bbeta\sim\Normal_p(\bzero,\tau^2\I) \) gives ridge regression (@thm-shr-ridge) with \( \lambda=\sigma^2/\tau^2 \),
   and the posterior mode equals the posterior mean;

2. independent Laplace priors with density \( (a/2)e^{-a\lvert b_j\rvert} \) give the lasso with \( \lambda=a\sigma^2 \).
:::

:::

::: {.proof}
By Bayes' theorem the posterior density is proportional to the likelihood times the prior, that is, to
\( \exp\{-\norm{\y-\X\bb}^2/(2\sigma^2)\}\exp\{-(\lambda/\sigma^2)P(\bb)\}=\exp\{-Q_\lambda(\bb)/\sigma^2\} \), and
\( t\mapsto e^{-t/\sigma^2} \) is decreasing. For (a), \( \tau^{-2}\norm{\bb}^2/2=(\lambda/\sigma^2)\norm{\bb}^2/2 \) with
\( \lambda=\sigma^2/\tau^2 \). The posterior is normal, so its mode and mean coincide (@thm-opt-bayes-conjugate with \( \sigma^2 \) known, or @prp-opt-shrinkage(d)). For (b), as in @exr-shr-laplace,
\( a\norm{\bb}_1=(\lambda/\sigma^2)\norm{\bb}_1 \) with \( \lambda=a\sigma^2 \).
:::

Under squared error loss, though, a Bayesian reports the posterior mean, which under a Laplace prior is never sparse.

::: {#exm-reg-laplace-mean}
[Posterior mode and posterior mean under a Laplace prior]

Take one observation \( z\sim\Normal(\theta,1) \) and a Laplace prior on \( \theta \) with rate
\( a=1.5 \). By @prp-reg-map(b) the posterior mode is the lasso estimate with \( \lambda=1.5 \),
which is soft thresholding at \( 1.5 \) (@thm-shr-lasso-orthonormal). The posterior mean, computed by numerical
integration in the listing below, behaves differently:

| \( z \) | 0.5 | 1.0 | 2.0 | 4.0 |
|---|---|---|---|---|
| posterior mode | 0 | 0 | 0.500 | 2.500 |
| posterior mean | 0.174 | 0.364 | 0.859 | 2.509 |

The mode is zero for \( \lvert z\rvert\le1.5 \); the mean is zero only at \( z=0 \), and approaches \( z-1.5 \) for large \( z \)
([Figure 30.1.1](#fig-reg-thresholding)(c)). Exact zeros come from the mode, not from the prior; a prior that believes in
exact zeros needs a point mass, as in the spike-and-slab priors of George and McCulloch (1993).
[Chapter 39](../ch39-glms-in-practice-bayes/index.html) returns to priors as regularizers (@prp-prc-priors).
:::

```{.python .run #cell-penalties-laplace}
import numpy as np
from scipy import integrate

def soft(z, lam):
    return np.sign(z) * np.maximum(np.abs(z) - lam, 0.0)

def laplace_posterior_mean(z, rate):
    """E(theta | z) when z ~ N(theta, 1) and theta has density (rate/2) exp(-rate |theta|)."""
    w = lambda t: np.exp(-0.5 * (z - t) ** 2 - rate * abs(t))
    num = integrate.quad(lambda t: t * w(t), -30, 30, points=[0.0, z])[0]
    den = integrate.quad(w, -30, 30, points=[0.0, z])[0]
    return num / den

rate = 1.5
for z in [0.5, 1.0, 2.0, 4.0]:
    print(f"z = {z}: posterior mode {float(soft(z, rate)):.3f}, "
          f"posterior mean {laplace_posterior_mean(z, rate):.3f}")
```

## Penalties that are not convex

Under an orthonormal design the lasso subtracts \( \lambda \) from every surviving coefficient, however large (@thm-shr-lasso-orthonormal).
Fan and Li (2001) asked for a rule that is **sparse**, **continuous** in the data, and
**nearly unbiased** (large estimates left alone). Soft thresholding fails the third, and hard thresholding, the best-subset
rule (@thm-shr-lasso-orthonormal(c)), the second. No convex penalty has all three, since one that is flat far from zero on
both sides is constant. Two penalties \( P(\bb)=\sum_j\rho(b_j) \), with \( \rho \) concave on \( [0,\infty) \), are in common use.

**SCAD** (the smoothly clipped absolute deviation of Fan and Li, 2001), with a parameter \( a>2 \), is defined through
its derivative on \( t>0 \):
\[
\rho'(t)=\lambda\Bigl\{\mathbf 1(t\le\lambda)+\frac{(a\lambda-t)_+}{(a-1)\lambda}\mathbf 1(t>\lambda)\Bigr\},
\qquad \rho(0)=0,\quad \rho(-t)=\rho(t).
\]
It is the lasso penalty up to \( \lambda \); its slope then falls linearly to zero at \( a\lambda \). Here \( \lambda \) sits
inside \( \rho \), and the criterion is \( \tfrac12\norm{\y-\X\bb}^2+\sum_j\rho(b_j) \). Fan and Li suggest \( a=3.7 \).

**MCP** (the minimax concave penalty of Zhang, 2010), with a parameter \( \gamma>1 \), is
\[
\rho(t)=\begin{cases}\lambda\lvert t\rvert-t^2/(2\gamma), & \lvert t\rvert\le\gamma\lambda,\\ \gamma\lambda^2/2, & \lvert t\rvert>\gamma\lambda.\end{cases}
\]
Its slope \( \lambda-\lvert t\rvert/\gamma \) starts to fall at once, rather than after \( \lambda \) as for SCAD.

::: {#prp-reg-mcp}
[MCP: convexity and firm thresholding]

Let \( Q(\bb)=\tfrac12\norm{\y-\X\bb}^2+\sum_j\rho(b_j) \) with \( \rho \) the MCP penalty.

::: {.enumerate options="label=(\alph*)"}
1. If the smallest eigenvalue of \( \X\T\X \) is at least \( 1/\gamma \), then \( Q \) is convex. If it exceeds
   \( 1/\gamma \), then \( Q \) is strictly convex and has a unique minimizer.

2. If \( \X\T\X=\I \), the unique minimizer is \( \hat\beta_j=\phi(z_j) \) with \( \mathbf{z}=\X\T\y \), where \( \phi \) is **firm
   thresholding**:
   \[
   \phi(z)=\begin{cases}0,&\lvert z\rvert\le\lambda,\\ \operatorname{sign}(z)\dfrac{\lvert z\rvert-\lambda}{1-1/\gamma},&\lambda<\lvert z\rvert\le\gamma\lambda,\\ z,&\lvert z\rvert>\gamma\lambda.\end{cases}
   \]
:::

:::

::: {.proof}
(a) Write
\[
\begin{aligned}
Q(\bb)={}&\Bigl\{\tfrac12\bb\T\bigl(\X\T\X-\gamma^{-1}\I\bigr)\bb-\y\T\X\bb+\tfrac12\y\T\y\Bigr\}\\
&+\sum_j\Bigl\{\rho(b_j)+\frac{b_j^2}{2\gamma}\Bigr\}.
\end{aligned}
\]
The first brace is a quadratic with nonnegative definite matrix, hence convex, and strictly convex if the
matrix is positive definite. Each term of the second sum is \( \lambda\lvert t\rvert \) for \( \lvert t\rvert\le\gamma\lambda \) and
\( \gamma\lambda^2/2+t^2/(2\gamma) \) beyond. On \( t>0 \) its derivative is \( \lambda \) up to \( \gamma\lambda \) and \( t/\gamma \) after. This derivative
is continuous at \( \gamma\lambda \) and nondecreasing, so the function is convex (it is even, with a corner at 0). A sum of convex
functions is convex, and adding a strictly convex function makes the sum strictly convex. A strictly convex function that
tends to infinity has exactly one minimizer, by the argument of @thm-reg-existence(a) and (c).

(b) With \( \X\T\X=\I \), part (a) applies because \( \gamma>1 \). Moreover
\( Q(\bb)=\tfrac12\norm{\mathbf{z}-\bb}^2+\sum_j\rho(b_j)+\text{const} \), which separates into one problem per coordinate:
minimize \( g(\theta)=\tfrac12(z-\theta)^2+\rho(\theta) \). Take \( z\ge0 \). For \( \theta>0 \),
\( g(-\theta)-g(\theta)=2z\theta\ge0 \), so we may restrict to \( \theta\ge0 \). There
\[
g'(\theta)=\begin{cases}(1-1/\gamma)\theta-(z-\lambda),&0<\theta\le\gamma\lambda,\\ \theta-z,&\theta>\gamma\lambda,\end{cases}
\]
a continuous nondecreasing function. If \( z\le\lambda \), then \( g'\ge0 \) on \( (0,\infty) \), and the minimizer is
\( 0 \). If \( \lambda<z\le\gamma\lambda \), then \( g' \) vanishes at \( \theta=(z-\lambda)/(1-1/\gamma) \), which is at most
\( \gamma\lambda \). If \( z>\gamma\lambda \), then \( g'(\gamma\lambda)=\gamma\lambda-z<0 \), and \( g' \) vanishes at \( \theta=z \).
:::

Firm thresholding has all three properties; it tends to soft thresholding as \( \gamma\to\infty \) and to hard
thresholding as \( \gamma\downarrow1 \). The SCAD rule under an orthonormal design is (Fan and Li, 2001)
\[
\phi_{\mathrm{SCAD}}(z)=\begin{cases}\operatorname{sign}(z)(\lvert z\rvert-\lambda)_+,&\lvert z\rvert\le2\lambda,\\ \dfrac{(a-1)z-\operatorname{sign}(z)a\lambda}{a-2},&2\lambda<\lvert z\rvert\le a\lambda,\\ z,&\lvert z\rvert>a\lambda,\end{cases}
\]
which @exr-reg-scad derives. The script checks each rule of [Figure 30.1.1](#fig-reg-thresholding) against brute-force
minimization.

::: {when-format="html"}
![**Figure 30.1.1.** (a) The lasso, SCAD (\( a=3.7 \)), MCP (\( \gamma=3 \)) and \( \ell_0 \) penalties with
\( \lambda=1 \). (b) The thresholding rules they induce under an orthonormal design (soft, SCAD, firm, and
hard thresholding dashed); the grey line is the identity. (c) Posterior mode and posterior mean of \( \theta \) given \( z\sim\Normal(\theta,1) \) under a Laplace prior
with rate \( 1.5 \) (@exm-reg-laplace-mean).](thresholding.svg){#fig-reg-thresholding width=100%}
:::

::: {when-format="pdf"}
![(a) The lasso, SCAD (\( a=3.7 \)), MCP (\( \gamma=3 \)) and \( \ell_0 \) penalties with
\( \lambda=1 \). (b) The thresholding rules they induce under an orthonormal design (soft, SCAD, firm, and
hard thresholding dashed); the grey line is the identity. (c) Posterior mode and posterior mean of \( \theta \) given \( z\sim\Normal(\theta,1) \) under a Laplace prior
with rate \( 1.5 \) (@exm-reg-laplace-mean).](thresholding.pdf){width=100%}
:::

```{.python .run #cell-penalties-rules}
import numpy as np

def soft(z, lam):
    return np.sign(z) * np.maximum(np.abs(z) - lam, 0.0)

def mcp_rule(z, lam, gam):                    # firm thresholding
    az = np.abs(z)
    mid = np.sign(z) * (az - lam) / (1 - 1 / gam)
    return np.where(az <= lam, 0.0, np.where(az <= gam * lam, mid, z))

def scad_rule(z, lam, a):
    az = np.abs(z)
    mid = ((a - 1) * z - np.sign(z) * a * lam) / (a - 2)
    return np.where(az <= 2 * lam, soft(z, lam), np.where(az <= a * lam, mid, z))

def hard_rule(z, lam):                        # penalty (lam^2 / 2) * 1{theta != 0}
    return np.where(np.abs(z) > lam, z, 0.0)
```

Without the eigenvalue condition there can be several local minima (@exr-reg-mcp-nonconvex). Replacing \( \rho \) by its
tangent at a preliminary estimate gives a weighted lasso (Zou and Li, 2008), as does the adaptive lasso (Zou, 2006). Fan and Li (2001) proved an
*oracle property*, which we state without proof: with \( p \) fixed, \( \lambda_n/n\to0 \), \( \lambda_n/\sqrt n\to\infty \) (in the
scaling of @eq-reg-pls; Fan and Li write the penalty as \( n\sum_jp_\lambda(\lvert b_j\rvert) \), where the conditions read
\( \lambda_n\to0 \), \( \sqrt n\lambda_n\to\infty \)) and regularity conditions, some local SCAD minimizer finds the zero coefficients with probability tending to one and estimates the rest as well as
least squares on the true support. The convergence is not uniform in \( \bbeta \) (Leeb and Pötscher, 2005), and it does not justify
oracle standard errors ([Chapter 29](../ch29-model-selection/index.html)).

## Exercises

### A. Check your understanding

::: {#exr-reg-rules}
[A1]

With \( \lambda=1 \), compute the lasso, ridge (penalty \( \tfrac12b^2 \)), hard, firm (\( \gamma=3 \)) and SCAD
(\( a=3.7 \)) estimates of \( \theta \) for a single observation \( z=2.5 \) under an orthonormal design. Which rules shrink it,
and by how much?
:::

::: {.solution}
Lasso: \( 2.5-1=1.5 \). Ridge: \( 2.5/2=1.25 \). Hard: \( 2.5 \), since \( \lvert z\rvert>1 \). Firm: \( \lambda<z\le\gamma\lambda=3 \), so
\( (2.5-1)/(1-1/3)=2.25 \). SCAD: \( 2<z\le3.7 \), so \( (2.7\cdot2.5-3.7)/1.7=3.05/1.7=1.79 \). Only hard
thresholding leaves \( z \) alone. At \( z=4 \) the firm and SCAD rules would leave it alone as well, since \( 4 \) exceeds both
\( \gamma\lambda=3 \) and \( a\lambda=3.7 \).
:::

::: {#exr-reg-firm-continuous}
[A2]

Show that the firm thresholding rule \( \phi \) of @prp-reg-mcp(b) is continuous at \( \lvert z\rvert=\lambda \) and at
\( \lvert z\rvert=\gamma\lambda \), whereas the hard thresholding rule \( z\mathbf 1(\lvert z\rvert>\lambda) \) jumps by
\( \lambda \) at \( \lvert z\rvert=\lambda \). What is the largest slope of \( \phi \), and what happens to it as \( \gamma\downarrow1 \)?
:::

::: {.solution}
Take \( z\ge0 \). As \( z\downarrow\lambda \) the middle branch \( (z-\lambda)/(1-1/\gamma) \) tends to \( 0 \), the value below
\( \lambda \). At \( z=\gamma\lambda \) it equals \( (\gamma-1)\lambda/\{(\gamma-1)/\gamma\}=\gamma\lambda \), the value of the upper
branch. Hard thresholding is \( 0 \) at \( z=\lambda \) and tends to \( \lambda \) as \( z\downarrow\lambda \). The slope of \( \phi \) is
\( 0 \), then \( \gamma/(\gamma-1) \), then \( 1 \); the largest is \( \gamma/(\gamma-1) \), which tends to infinity as
\( \gamma\downarrow1 \), where \( \phi \) approaches the discontinuous hard rule.
:::

### B. Practice

::: {#exr-reg-scad}
[B1]

Derive the SCAD thresholding rule displayed above. Minimize \( g(\theta)=\tfrac12(z-\theta)^2+\rho(\theta) \) for
\( z\ge0 \) on the three pieces \( [0,\lambda] \), \( [\lambda,a\lambda] \) and \( [a\lambda,\infty) \), and check that \( g \) is convex on
\( [0,\infty) \) when \( a>2 \).
:::

::: {.solution}
On \( (0,\infty) \), \( g'(\theta)=\theta-z+\rho'(\theta) \). This is \( \theta-z+\lambda \) on \( (0,\lambda] \),
\( \theta-z+(a\lambda-\theta)/(a-1) \) on \( (\lambda,a\lambda] \), and \( \theta-z \) beyond. The middle piece has slope
\( 1-1/(a-1)=(a-2)/(a-1)>0 \) when \( a>2 \), and the pieces agree at \( \lambda \) and at \( a\lambda \). So \( g' \) is continuous
and increasing, and \( g \) is strictly convex on \( [0,\infty) \). As in the proof of @prp-reg-mcp, a negative
\( \theta \) is never better when \( z\ge0 \). If \( z\le\lambda \), then \( g'(0+)=\lambda-z\ge0 \) and \( \hat\theta=0 \). If \( \lambda<z\le2\lambda \),
the root \( z-\lambda \) lies in \( (0,\lambda] \). If \( 2\lambda<z\le a\lambda \), the middle piece vanishes at
\( \theta=\{(a-1)z-a\lambda\}/(a-2) \), which lies in \( (\lambda,a\lambda] \) exactly for these \( z \). If \( z>a\lambda \), then
\( \hat\theta=z \).
:::

::: {#exr-reg-intercept}
[B2]

Let \( P \) be any penalty on the slopes, and let the model have an unpenalized intercept. By @prp-shr-centring,
\( (\hat b_0,\hbeta) \) minimizes \( \tfrac12\sum_i(y_i-b_0-\x_{(i)}\T\bb)^2+\lambda P(\bb) \) iff \( \hbeta \) minimizes the centred
criterion and \( \hat b_0=\bar y-\bar{\x}\T\hbeta \). Deduce that the fitted values of penalized regression with an intercept
are unchanged when a constant is added to \( \y \), and that the slopes are unchanged when a constant is added to a column of
\( \X \). Does the second statement survive standardizing the columns before fitting?
:::

### C. Going deeper

::: {#exr-reg-mcp-nonconvex}
[C1]

Let \( n=p=2 \), \( \X\T\X=\begin{pmatrix}1&r\\r&1\end{pmatrix} \) with \( r\in(0,1) \), and use MCP with \( \gamma>1 \). For
which \( r \) does @prp-reg-mcp(a) guarantee convexity? Find a response \( \y \) for which the MCP criterion has two
distinct local minima when \( r \) exceeds that bound.
:::

