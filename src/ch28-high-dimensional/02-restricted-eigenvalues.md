# Restricted eigenvalue conditions

In Parts II and III, the precision of least squares was governed by the smallest eigenvalue of \( \X\T\X \): the
worst variance of \( \mathbf{a}\T\hbeta \) over unit vectors \( \mathbf{a} \) is \( \sigma^2/\lambda_{\min}(\X\T\X) \) (@prp-dep-collinear-directions).
When \( p>n \) that eigenvalue is zero. This section explains which weaker
condition takes its place once \( \bbeta \) is assumed sparse, and shows with numbers how demanding it is.

## Notation

Throughout Sections 28.2–28.5 the columns of \( \X \) are scaled so that \( \norm{\x_j}^2=n \) for every
\( j \), and \( \hat{\bSigma}=\X\T\X/n \) is the sample Gram matrix, with unit diagonal. For a vector
\( \bv\in\Real^p \) and a set \( S\subseteq\{1,\dots,p\} \), \( \bv_S \) is the vector that agrees with \( \bv \) on \( S \)
and is zero off \( S \), so that \( \bv=\bv_S+\bv_{S^c} \). The \( \ell_1 \) and \( \ell_\infty \) norms are
\( \norm{\bv}_1=\sum_j|v_j| \) and \( \norm{\bv}_\infty=\max_j|v_j| \), and \( \norm{\bv} \) is still the Euclidean
norm. The **support** of \( \bbeta \) is \( S=\{j:\beta_j\ne0\} \), and \( s=|S| \). A vector with at most \( s \) nonzero
entries is **\( s \)-sparse**. By Cauchy–Schwarz (@prp-mat-cauchy-schwarz), for every \( \bv \),
\[
\norm{\bv_S}_1\le\sqrt{s}\,\norm{\bv_S}\le\sqrt s\,\norm{\bv}.
\]{#eq-hd-l1-l2}

## Sparsity restores identifiability

If \( \bbeta \) is only known to lie in \( \Real^p \), then \( \bbeta \) and \( \bbeta+\bv \) with \( \bv\in\Null(\X) \) give the same
mean, and nothing can separate them. If \( \bbeta \) is known to be sparse, most of those alternatives are
excluded.

::: {#prp-hd-sparse-identifiable}
[Identifiability of sparse vectors]

Distinct \( s \)-sparse vectors \( \bbeta\ne\bbeta' \) always have \( \X\bbeta\ne\X\bbeta' \) iff every set of \( 2s \)
columns of \( \X \) is linearly independent. In particular this requires \( 2s\le n \).
:::

::: {.proof}
If \( \bbeta,\bbeta' \) are \( s \)-sparse and \( \X\bbeta=\X\bbeta' \), then \( \bv=\bbeta-\bbeta' \) is a vector in \( \Null(\X) \)
with at most \( 2s \) nonzero entries. It is zero if every \( 2s \) columns are independent. Conversely, if some
\( 2s \) columns are dependent, there is a nonzero \( \bv\in\Null(\X) \) supported on them. Split its support into two
sets of at most \( s \) indices, \( T_1 \) and \( T_2 \), and put \( \bbeta=\bv_{T_1} \), \( \bbeta'=-\bv_{T_2} \). These are distinct
\( s \)-sparse vectors with \( \X\bbeta-\X\bbeta'=\X\bv=\bzero \). Finally, \( 2s \) vectors in \( \Real^n \) can be independent
only if \( 2s\le n \).
:::

Identifiability is a statement about exactly sparse vectors and noiseless data. For
estimation we need a quantitative version that allows for noise and for estimates that are only
approximately sparse. The key observation, proved in [Section 28.3](03-lasso-bounds.html), is that
the error of the lasso, \( \boldsymbol{\Delta}=\hat{\bbeta}_\lambda-\bbeta \), lies with high probability in a
particular cone. Its \( \ell_1 \) mass off the support is at most three times its mass on the support.
Curvature of the least squares criterion is needed only in the directions of that cone.

## Two conditions on a cone

::: {#def-hd-re}
[Compatibility and restricted eigenvalue conditions]

Let \( S\subseteq\{1,\dots,p\} \) have \( s\ge1 \) elements, and let
\[
\mathcal{K}(S)=\bigl\{\bv\in\Real^p:\ \norm{\bv_{S^c}}_1\le3\norm{\bv_S}_1\bigr\}.
\]{#eq-hd-cone}

The **compatibility constant** and the **restricted eigenvalue** of the design for \( S \) are
\[
\phi^2(S)=\min_{\bv\in\mathcal{K}(S),\,\bv\ne\bzero}\ \frac{s\,\bv\T\hat{\bSigma}\bv}{\norm{\bv_S}_1^2},\qquad
\kappa^2(S)=\min_{\bv\in\mathcal{K}(S),\,\bv\ne\bzero}\ \frac{\bv\T\hat{\bSigma}\bv}{\norm{\bv}^2}.
\]
The design satisfies the **compatibility condition** for \( S \) if \( \phi^2(S)>0 \), and the **restricted
eigenvalue condition** for \( S \) if \( \kappa^2(S)>0 \). The same definitions with a nonnegative definite
matrix \( \bSigma \) in place of \( \hat{\bSigma} \) give \( \phi^2_{\bSigma}(S) \) and \( \kappa^2_{\bSigma}(S) \).
:::

A nonzero \( \bv\in\mathcal{K}(S) \) has \( \bv_S\ne\bzero \), since otherwise \( \norm{\bv_{S^c}}_1\le0 \). So both ratios are
defined. Both are unchanged when \( \bv \) is rescaled, and the sets \( \{\bv\in\mathcal{K}(S):\norm{\bv}=1\} \) and
\( \{\bv\in\mathcal{K}(S):\norm{\bv_S}_1=1\} \) are closed and bounded. So the minima are attained, which
justifies writing \( \min \). The constant \( 3 \) comes from the choice of penalty in [Section 28.3](03-lasso-bounds.html).
Other choices give other constants and the same theory.

The cone contains the \( s \)-dimensional subspace of vectors supported on \( S \) but no larger subspace (@exr-hd-cone-subspaces).
So the null space of \( \X \), of dimension at least \( p-n \), can still avoid it.

::: {#prp-hd-re-basic}
[Basic facts]

::: {.enumerate options="label=(\alph*)"}
1. \( \lambda_{\min}(\hat{\bSigma})\le\kappa^2(S)\le\phi^2(S) \). If \( p>n \), then \( \lambda_{\min}(\hat{\bSigma})=0 \).

2. \( \kappa^2(S)>0 \) iff \( \phi^2(S)>0 \) iff \( \Null(\X)\cap\mathcal{K}(S)=\{\bzero\} \).

3. *(Transfer.)* If \( \bSigma \) is symmetric and \( \max_{j,k}|\hat{\Sigma}_{jk}-\Sigma_{jk}|\le\delta \), then
   \( \kappa^2(S)\ge\kappa_{\bSigma}^2(S)-16s\delta \) and \( \phi^2(S)\ge\phi_{\bSigma}^2(S)-16s\delta \).

4. *(Coherence.)* If \( \mu=\max_{j\ne k}|\hat{\Sigma}_{jk}| \) is the **mutual coherence** of the scaled design,
   then \( \kappa^2(S)\ge1-16s\mu \) and \( \phi^2(S)\ge1-16s\mu \).
:::

:::

::: {.proof}
(a) The first inequality holds because \( \kappa^2(S) \) minimizes the Rayleigh quotient of \( \hat{\bSigma} \) over a subset of
nonzero vectors, while \( \lambda_{\min}(\hat{\bSigma}) \) minimizes it over all of them (@thm-mat-extremal-rayleigh).
For the second, @eq-hd-l1-l2 gives \( \norm{\bv_S}_1^2\le s\norm{\bv}^2 \), so
\( s\,\bv\T\hat{\bSigma}\bv/\norm{\bv_S}_1^2\ge\bv\T\hat{\bSigma}\bv/\norm{\bv}^2\ge\kappa^2(S) \) for every \( \bv\in\mathcal{K}(S) \).
If \( p>n \), then \( \rank(\hat{\bSigma})\le n<p \), so \( \hat{\bSigma} \) is singular.

(b) \( \bv\T\hat{\bSigma}\bv=\norm{\X\bv}^2/n \) vanishes iff \( \X\bv=\bzero \). Since the minima are attained, each
constant is zero iff some nonzero \( \bv\in\mathcal{K}(S) \) has \( \X\bv=\bzero \).

(c) For \( \bv\in\mathcal{K}(S) \), \( \norm{\bv}_1=\norm{\bv_S}_1+\norm{\bv_{S^c}}_1\le4\norm{\bv_S}_1 \). Hence
\[
\bigl|\bv\T(\hat{\bSigma}-\bSigma)\bv\bigr|\le\sum_{j,k}|v_j||v_k|\,\delta=\delta\norm{\bv}_1^2
\le16\delta\norm{\bv_S}_1^2\le16s\delta\norm{\bv}^2,
\]
the last step by @eq-hd-l1-l2. Dividing by \( \norm{\bv}^2 \) gives
\( \bv\T\hat{\bSigma}\bv/\norm{\bv}^2\ge\kappa^2_{\bSigma}(S)-16s\delta \). Dividing the middle bound by
\( \norm{\bv_S}_1^2/s \) gives \( s\,\bv\T\hat{\bSigma}\bv/\norm{\bv_S}_1^2\ge\phi^2_{\bSigma}(S)-16s\delta \).

(d) Take \( \bSigma=\I \) in (c), with \( \delta=\mu \) (the diagonals agree). Then \( \kappa^2_{\I}(S)=1 \), and
\( \phi^2_{\I}(S)\ge\kappa^2_{\I}(S)=1 \) by (a) applied to \( \I \).
:::

By (b) the two conditions agree about *whether* the constants vanish, but by (a) their sizes, which enter the
bounds, can differ a great deal: a lower bound on \( \kappa^2(S) \) is the stronger assumption. The fast prediction
rate of [Section 28.3](03-lasso-bounds.html) uses \( \phi^2(S) \), the \( \ell_2 \) bound \( \kappa^2(S) \). Part (c) carries a
population condition to the sample when \( \hat{\bSigma} \) is uniformly close to \( \bSigma \), and (d) is the case \( \bSigma=\I \).

## When the conditions hold

**Populations.** If \( \lambda_{\min}(\bSigma)>0 \), part (a) applied to \( \bSigma \) gives
\( \phi^2_{\bSigma}(S)\ge\kappa^2_{\bSigma}(S)\ge\lambda_{\min}(\bSigma) \). For example, the equicorrelation matrix
\( \bSigma=(1-\rho)\I+\rho\bone\bone\T \), \( 0\le\rho<1 \), has
\( \bv\T\bSigma\bv=(1-\rho)\norm{\bv}^2+\rho(\bone\T\bv)^2\ge(1-\rho)\norm{\bv}^2 \), so both constants are at least \( 1-\rho \)
for every \( p \) and \( S \). The difficulty lies in going from \( \bSigma \) to \( \hat{\bSigma} \) when \( p>n \).

**Random designs.** Suppose the rows of \( \X \) are independent \( \Normal_p(\bzero,\bSigma) \) vectors with unit
variances. A union bound over the \( p^2 \) entries of \( \hat{\bSigma}-\bSigma \), each an average of \( n \) products of
normal variables, gives \( \max_{jk}|\hat{\Sigma}_{jk}-\Sigma_{jk}|\le C\sqrt{\log p/n} \) with high probability once
\( n\ge C'\log p \). The tail bound for each entry is not the normal bound of the next section, since
products of normal variables have heavier tails; it is a Bernstein-type inequality, not proved here (see
Wainwright 2019, chapter 2). Part (c) then gives a positive restricted eigenvalue once \( n \) exceeds a multiple of \( s^2\log p \). This is wasteful,
because it bounds \( \bv\T(\hat{\bSigma}-\bSigma)\bv \) by the worst entry. Raskutti, Wainwright and Yu (2010) proved, by
an argument we do not reproduce, that with probability tending to one
\[
\frac{\norm{\X\bv}}{\sqrt n}\ \ge\ c_1\,\norm{\bSigma^{1/2}\bv}-c_2\sqrt{\frac{\log p}{n}}\,\norm{\bv}_1
\qquad\text{for all }\bv\in\Real^p ,
\]
for constants \( c_1,c_2>0 \). On the cone, \( \norm{\bv}_1\le4\sqrt s\norm{\bv} \), so this gives
\( \kappa(S)\ge c_1\lambda_{\min}(\bSigma)^{1/2}-4c_2\sqrt{s\log p/n} \), which is positive once \( n \) exceeds a multiple of
\( s\log p \), possibly far smaller than \( p \).

## When they fail

The conditions fail when a vector of the cone is (nearly) invisible to \( \X \).

*Duplicated columns.* If \( \x_k=\x_j \) with \( j\in S \) and \( k\notin S \), then \( \bv=\mathbf{e}_j-\mathbf{e}_k \) has
\( \norm{\bv_{S^c}}_1=\norm{\bv_S}_1=1 \), so it lies in the cone, and \( \X\bv=\bzero \). Both constants are zero.
Nearly equal columns give nearly zero constants (@exr-hd-near-duplicates). No estimator can say whether
the effect belongs to \( \x_j \) or to its copy.

*Short combinations of inactive columns.* More generally, suppose an active column is a combination
of inactive ones, \( \x_j=\X_{S^c}\mathbf{a} \) with \( \norm{\mathbf{a}}_1\le3 \). Then \( \bv=\mathbf{e}_j-\mathbf{a} \) (with \( \mathbf{a} \) placed in
the coordinates of \( S^c \)) lies in the cone and \( \X\bv=\bzero \). When \( p>n \) the inactive columns usually span
\( \Real^n \), so the conditions require every such combination to have a large \( \ell_1 \) norm.

*Too many active variables for the sample size.* Even with independent regressors and a population
compatibility constant of one, the sample constant can be close to zero if \( s \) is large relative to
\( n/\log p \). A combination of hundreds of inactive columns with total \( \ell_1 \) weight \( 3 \) can then imitate
\( \X_S\bv_S \) closely. The next example measures this.

::: {#exm-hd-compatibility-values}
[Compatibility constants of Gaussian designs]

Take \( p=1000 \) regressors, \( S=\{1,\dots,s\} \), and three random designs, with columns scaled so that
\( \norm{\x_j}^2=n \): independent standard normal entries with \( n=200 \) (the design used in
[Section 28.3](03-lasso-bounds.html)) and with \( n=400 \), and equicorrelated normal regressors with
\( \rho=0.5 \) and \( n=400 \). The compatibility constant is a minimum of a convex quadratic over the cone.
For each sign pattern of \( \bv_S \) the problem is convex, and a duality gap certifies the computed
value. The table gives \( \phi^2(S) \) and the mutual coherence \( \mu \); its rows are labelled by the common correlation \( \rho \).

| Design | \( \mu \) | \( s=1 \) | \( s=2 \) | \( s=5 \) | Population |
|---|---|---|---|---|---|
| \( \rho=0 \), \( n=200 \) | \( 0.350 \) | \( 0.328 \) | \( 0.146 \) | \( 0.0045 \) | \( 1 \) |
| \( \rho=0 \), \( n=400 \) | \( 0.238 \) | \( 0.479 \) | \( 0.384 \) | \( 0.177 \) | \( 1 \) |
| \( \rho=0.5 \), \( n=400 \) | \( 0.641 \) | \( 0.266 \) | \( 0.192 \) | \( 0.089 \) | \( \ge0.5 \) |

The sample constants lie far below the population values, even for \( s=1 \), and fall quickly with \( s \).
At \( n=200 \) and \( s=5 \) the constant \( 0.0045 \) makes the fast-rate bound of the next section
numerically useless for that design; doubling \( n \) helps a great deal. The coherence bound of
@prp-hd-re-basic(d) is useless throughout: with \( \mu=0.350 \), \( 1-16s\mu \) is \( -4.60 \) already at \( s=1 \).
:::

```{.python .run #cell-compatibility-small}
import numpy as np

def proj_simplex(u):
    """Euclidean projection onto {u >= 0, sum(u) = 1}."""
    v = np.sort(u)[::-1]
    css = np.cumsum(v) - 1
    k = np.flatnonzero(v - css / np.arange(1, len(u) + 1) > 0)[-1]
    return np.maximum(u - css[k] / (k + 1), 0)

def proj_l1(w, radius):
    """Euclidean projection onto the l1 ball of the given radius."""
    if np.abs(w).sum() <= radius:
        return w
    return np.sign(w) * proj_simplex(np.abs(w) / radius) * radius

def compatibility(X, S, radius=3.0, iters=4000):
    """Return phi^2(S) and a certified lower bound (value minus duality gap)."""
    n, p = X.shape
    s = len(S)
    Sc = np.setdiff1d(np.arange(p), S)
    A, C = X[:, S], X[:, Sc]
    L = 2 * s / n * np.linalg.norm(X, 2) ** 2         # Lipschitz constant of the gradient
    best, lower = np.inf, np.inf
    for code in range(2 ** (s - 1)):                  # sign patterns of v_S up to a global sign
        sg = np.array([1.0] + [-1.0 if (code >> k) & 1 else 1.0 for k in range(s - 1)])
        u, w = np.full(s, 1 / s), np.zeros(len(Sc))
        yu, yw, t = u.copy(), w.copy(), 1.0
        for _ in range(iters):                        # accelerated projected gradient
            r = A @ (sg * yu) + C @ yw
            gu, gw = 2 * s / n * sg * (A.T @ r), 2 * s / n * (C.T @ r)
            u_new, w_new = proj_simplex(yu - gu / L), proj_l1(yw - gw / L, radius)
            t_new = (1 + np.sqrt(1 + 4 * t * t)) / 2
            yu = u_new + (t - 1) / t_new * (u_new - u)
            yw = w_new + (t - 1) / t_new * (w_new - w)
            u, w, t = u_new, w_new, t_new
        r = A @ (sg * u) + C @ w
        value = s / n * r @ r
        gu, gw = 2 * s / n * sg * (A.T @ r), 2 * s / n * (C.T @ r)
        gap = gu @ u - gu.min() + gw @ w + radius * np.abs(gw).max()
        best, lower = min(best, value), min(lower, value - gap)
    return best, lower

rng = np.random.default_rng(2820)
n, p = 100, 300
X = rng.normal(size=(n, p))
X /= np.sqrt(np.sum(X ** 2, axis=0) / n)            # ||x_j||^2 = n
print("smallest eigenvalue of X^T X / n:", np.linalg.eigvalsh(X.T @ X / n)[0])
for s in [1, 2, 3]:
    value, lower = compatibility(X, np.arange(s), iters=1500)
    print(f"s = {s}: phi^2 = {value:.3f}")
```

::: {.remark}
[Checking the conditions in practice]

The constants depend on the unknown \( S \) and cannot be checked from data. They describe the designs for which
the theory is informative: those in which no short combination of columns imitates another. With groups of
nearly equivalent regressors, common in genomics, the lasso may still predict well (by the slow rate of
[Section 28.3](03-lasso-bounds.html)), but its choice among them is arbitrary.
:::

## Exercises

### A. Check your understanding

::: {#exr-hd-cone-nonconvex}
[A1]

Let \( p=2 \) and \( S=\{1\} \). Sketch \( \mathcal{K}(S) \) and show that it is not convex. Show that for general
\( S \) it is the union of \( 2^s \) convex cones, one for each sign pattern of \( \bv_S \). This is why the example
above solves one convex problem per sign pattern.
:::

### B. Practice

::: {#exr-hd-near-duplicates}
[B1]

Let \( j\in S \) and \( k\notin S \). Show that \( \phi^2(S)\le s\norm{\x_j-\x_k}^2/n \) and
\( \kappa^2(S)\le\norm{\x_j-\x_k}^2/(2n) \). If \( \x_j\T\x_k/n=r \) (the cosine of the angle between the two
columns, which is their sample correlation when the columns are centred), conclude that \( \phi^2(S)\le2s(1-r) \).
:::

::: {.solution}
Take \( \bv=\mathbf{e}_j-\mathbf{e}_k\in\mathcal{K}(S) \). Then \( \norm{\bv_S}_1=1 \), \( \norm{\bv}^2=2 \) and
\( \bv\T\hat{\bSigma}\bv=\norm{\x_j-\x_k}^2/n \), which gives both bounds. With \( \norm{\x_j}^2=\norm{\x_k}^2=n \),
\( \norm{\x_j-\x_k}^2/n=2-2\x_j\T\x_k/n=2(1-r) \).
:::

::: {#exr-hd-equicorrelation-re}
[B2]

For \( \bSigma=(1-\rho)\I+\rho\bone\bone\T \) with \( 0\le\rho<1 \) and \( p\ge2 \), show that
\( \kappa_{\bSigma}^2(S)=1-\rho \) for every \( S \). Is the equicorrelation structure a hard case for the
restricted eigenvalue condition?
:::

::: {.solution}
The lower bound \( 1-\rho \) was shown in the text. For equality, it is enough to find a nonzero \( \bv\in\mathcal{K}(S) \)
with \( \bone\T\bv=0 \), since then \( \bv\T\bSigma\bv=(1-\rho)\norm{\bv}^2 \). If \( s\ge2 \), take \( \mathbf{e}_j-\mathbf{e}_k \) with
\( j,k\in S \). If \( s=1 \), take \( \mathbf{e}_j-\mathbf{e}_k \) with \( S=\{j\} \) and any \( k\ne j \); it lies in the cone because
\( 1\le3 \). The population condition is therefore as good as the smallest eigenvalue, and equicorrelation
is not a hard case in the population. It is harder in the sample (the table above), and harder still for
exact support recovery, as [Section 28.4](04-support-recovery.html) shows for a related design.
:::

### C. Going deeper

::: {#exr-hd-cone-subspaces}
[C1]

Show that every linear subspace contained in \( \mathcal{K}(S) \) has dimension at most \( s \), and that the bound is
attained. Deduce that if \( \dim\Null(\X)>p-s \) then the conditions fail for every \( S \) of size \( s \),
whereas \( \dim\Null(\X)=p-n \) with \( n\ge s \) is not excluded by this dimension count.
:::

::: {#exr-hd-null-space-property}
[C2]

*(Exact recovery by \( \ell_1 \) minimization.)* Say \( \X \) has the **null space property** for \( S \) if every nonzero
\( \bv\in\Null(\X) \) has \( \norm{\bv_{S^c}}_1>\norm{\bv_S}_1 \). Show that this holds iff, for every \( \bbeta \) supported
on \( S \), \( \bbeta \) is the unique solution of: minimize \( \norm{\bb}_1 \) subject to \( \X\bb=\X\bbeta \). Relate the
property to @prp-hd-re-basic(b) with the constant \( 3 \) replaced by \( 1 \).
:::

::: {.solution}
Suppose the property holds, \( \bbeta \) is supported on \( S \), and \( \X\bb=\X\bbeta \) with \( \bb\ne\bbeta \). Then
\( \bv=\bb-\bbeta\in\Null(\X) \) is nonzero, and \( \bb_S=\bbeta+\bv_S \), \( \bb_{S^c}=\bv_{S^c} \), so
\( \norm{\bb}_1\ge\norm{\bbeta}_1-\norm{\bv_S}_1+\norm{\bv_{S^c}}_1>\norm{\bbeta}_1 \). So \( \bbeta \) is the unique minimizer.
Conversely, if some nonzero \( \bv\in\Null(\X) \) has \( \norm{\bv_{S^c}}_1\le\norm{\bv_S}_1 \), put \( \bbeta=\bv_S \) and
\( \bb=-\bv_{S^c} \). Then \( \X\bb=\X\bbeta-\X\bv=\X\bbeta \), \( \bb\ne\bbeta \) and \( \norm{\bb}_1\le\norm{\bbeta}_1 \), so \( \bbeta \) is not
the unique minimizer. The property says exactly that \( \Null(\X) \) meets the cone
\( \{\norm{\bv_{S^c}}_1\le\norm{\bv_S}_1\} \) only at zero. This is part (b) with the constant \( 1 \) in place of \( 3 \). The
noisy problem needs the larger cone because the penalty must also absorb the noise.
:::
