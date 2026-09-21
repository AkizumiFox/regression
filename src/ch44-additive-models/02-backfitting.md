# Backfitting and concurvity

[Section 44.1](01-additive-models.html) reduced the additive model to one penalized
regression with \( d=p+\sum_jd_j \) columns. When \( d \) is a few dozen, solving
@eq-add-normal directly is right and this section is only a change of view. When \( d \)
runs into the thousands — one coefficient per district, one per subject — or when a term's
smoother is not given by a basis at all, an older idea takes over: fit each term in turn to
what the others have left behind, and go round again.

## The algorithm

Write \( \bS_j \) for the smoother of the \( j \)th term, the linear operator that fits
that term alone to a response vector:
\[
\bS_j=\Z_j\bigl(\Z_j\T\Z_j+\lambda_j\bP_j\bigr)^{-1}\Z_j\T .
\]{#eq-add-blocksmoother}

The **partial residual** for term \( j \) is what the response has left once every other
part of the model is subtracted,
\( \br_j=\y-\X\hbeta-\sum_{k\ne j}\hat{\mathbf{f}}_k \); backfitting replaces
\( \hat{\mathbf{f}}_j \) by \( \bS_j\br_j \) and moves on.

::: {.algorithm}
**Backfitting.**

1. Initialize \( \hbeta \) by least squares of \( \y \) on \( \X \), and
   \( \hat{\mathbf{f}}_j=\bzero \) for \( j=1,\dots,q \).
2. For \( j=1,\dots,q \) in turn, set
   \( \br_j=\y-\X\hbeta-\sum_{k\ne j}\hat{\mathbf{f}}_k \) and
   \( \hat{\mathbf{f}}_j\leftarrow\bS_j\br_j \).
3. Update \( \hbeta\leftarrow(\X\T\X)^{-1}\X\T\bigl(\y-\sum_k\hat{\mathbf{f}}_k\bigr) \),
   and return to step 2 until the largest change in any \( \hat{\mathbf{f}}_j \) is below a
   tolerance.
:::

Each step is a one-dimensional smoothing problem of the kind
[Chapter 43](../ch43-smoothing/index.html) solved; the terms are never assembled into one
large matrix; and a smoother *not* of the form @eq-add-blocksmoother — a kernel estimate, a
monotone fit — drops in unchanged. What is not obvious is that the cycling converges, and
to what.

## Backfitting is Gauss–Seidel

::: {#thm-add-backfitting}
[Backfitting]

Consider the criterion @eq-add-criterion with blocks \( \Z_0=\X \),
\( \bP_0=\bzero \) and \( \Z_1,\dots,\Z_q \), and write
\( \A_{jj}=\Z_j\T\Z_j+\lambda_j\bP_j \), assumed nonsingular for every \( j \). Let
\( \A=\Z\T\Z+\bP \) as in @eq-add-normal, so that \( \A \) is nonnegative definite
with diagonal blocks \( \A_{jj} \).

::: {.enumerate options="label=(\alph*)"}
1. *(Coordinate descent.)* The backfitting update of the \( j \)th block is the exact
   minimizer of \( Q \) over \( \bgamma_j \) with all other blocks held fixed. In
   particular \( Q \) never increases along the iteration.

2. *(Gauss–Seidel.)* In terms of coefficients, the update is
   \( \bgamma_j\leftarrow\A_{jj}^{-1}\Z_j\T\bigl(\y-\sum_{k\ne j}\Z_k\bgamma_k\bigr) \),
   which is one block Gauss–Seidel sweep for the linear system
   \( \A\bgamma=\Z\T\y \). Fix any order of the blocks — the algorithm above takes them in
   the order \( 1,\dots,q,0 \) — and write \( \A=\bL+\bD+\bL\T \) with \( \bD \) the block
   diagonal and \( \bL \) the strictly block lower triangular part in that order. One full
   cycle sends
   \( \bgamma \) to \( \mathbf{T}\bgamma+(\bD+\bL)^{-1}\Z\T\y \) with
   \[
   \mathbf{T}=\I-(\bD+\bL)^{-1}\A .
   \]{#eq-add-gauss-seidel}

3. *(Convergence.)* Let \( \bgamma^{(t)} \) be the iterates from any starting point and
   let \( \hat{\bgamma} \) be any solution of \( \A\bgamma=\Z\T\y \). Then
   \( \A^{1/2}(\bgamma^{(t)}-\hat{\bgamma})\to\bzero \) geometrically. Consequently the
   fitted vector \( \Z\bgamma^{(t)} \) converges to \( \Z\hat{\bgamma} \), the criterion
   \( Q(\bgamma^{(t)}) \) converges to its minimum, and \( \bgamma^{(t)} \) itself
   converges to a minimizer of \( Q \), which may depend on the starting point.

4. *(Uniqueness.)* If \( \A \) is positive definite, the limit is the unique solution
   \( \A^{-1}\Z\T\y \) of @eq-add-normal, and the convergence is geometric in every norm.
:::

:::

::: {.proof}
(a) As a function of \( \bgamma_j \) alone, with \( \br_j=\y-\sum_{k\ne j}\Z_k\bgamma_k \),
\[
Q=\norm{\br_j-\Z_j\bgamma_j}^2+\lambda_j\bgamma_j\T\bP_j\bgamma_j+\text{const},
\]
strictly convex because \( \A_{jj} \) is nonsingular; its gradient vanishes at
\( \A_{jj}\bgamma_j=\Z_j\T\br_j \), the stated update, and
\( \Z_j\bgamma_j=\bS_j\br_j \) by @eq-add-blocksmoother. Each step minimizes \( Q \) in
one block, so \( Q \) is nonincreasing.

(b) The \( j \)th block row of \( \A\bgamma=\Z\T\y \) reads
\( \A_{jj}\bgamma_j+\sum_{k\ne j}\Z_j\T\Z_k\bgamma_k=\Z_j\T\y \), and solving it for
\( \bgamma_j \) with the other blocks at their current values is the update in (a); taking
the *new* values for \( k<j \) and the old for \( k>j \) is block Gauss–Seidel. Over one
full cycle the new iterate satisfies
\( (\bD+\bL)\bgamma^{\text{new}}=-\bL\T\bgamma^{\text{old}}+\Z\T\y \), which
rearranges to @eq-add-gauss-seidel.

(c) Put \( \mathbf{E}=\bD+\bL \), so that \( \mathbf{E}+\mathbf{E}\T-\A=\bD \). The error
\( \mathbf{e}^{(t)}=\bgamma^{(t)}-\hat{\bgamma} \) obeys \( \mathbf{e}^{(t+1)}=\mathbf{T}\mathbf{e}^{(t)} \).
For any \( \mathbf{e} \), with \( \bv=\mathbf{E}^{-1}\A\mathbf{e} \) so that
\( \mathbf{T}\mathbf{e}=\mathbf{e}-\bv \),

\[
\mathbf{e}\T\A\mathbf{e}-(\mathbf{T}\mathbf{e})\T\A(\mathbf{T}\mathbf{e})
=2\bv\T\A\mathbf{e}-\bv\T\A\bv
=2\bv\T\mathbf{E}\bv-\bv\T\A\bv
=\bv\T(\mathbf{E}+\mathbf{E}\T-\A)\bv=\bv\T\bD\bv ,
\]{#eq-add-energy}

using \( \A\mathbf{e}=\mathbf{E}\bv \) and the symmetry of \( \A \). The matrix \( \bD \)
is positive definite, so the \( \A \)-seminorm of the error never increases, and it
decreases strictly unless \( \bv=\bzero \), that is unless \( \A\mathbf{e}=\bzero \).

For geometric decay, set \( \bu=\A^{1/2}\mathbf{e} \) and
\( \mathbf{N}=\I-\A^{1/2}\mathbf{E}^{-1}\A^{1/2} \), so that
\( \A^{1/2}\mathbf{T}\mathbf{e}=\mathbf{N}\bu \). In terms of \( \bu \),
@eq-add-energy reads, for *every* \( \bu\in\Real^{d} \),
\[
\norm{\bu}^2-\norm{\mathbf{N}\bu}^2=\bv\T\bD\bv,\qquad \bv=\mathbf{E}^{-1}\A^{1/2}\bu .
\]{#eq-add-contraction}

Now \( \C(\A^{1/2}) \) is invariant under \( \mathbf{N} \), since
\( \mathbf{N}\bu=\bu-\A^{1/2}(\mathbf{E}^{-1}\A^{1/2}\bu) \), and on that subspace
\( \bv=\bzero \) forces \( \A^{1/2}\bu=\bzero \), hence \( \bu=\bzero \). So
@eq-add-contraction gives \( \norm{\mathbf{N}\bu}<\norm{\bu} \) for every nonzero
\( \bu\in\C(\A^{1/2}) \), and its unit sphere is compact, so
\( \rho:=\max\{\norm{\mathbf{N}\bu}:\bu\in\C(\A^{1/2}),\norm{\bu}=1\}<1 \). As
\( \bu^{(0)}\in\C(\A^{1/2}) \) and the subspace is invariant,
\( \norm{\A^{1/2}\mathbf{e}^{(t)}}\le\rho^{t}\norm{\A^{1/2}\mathbf{e}^{(0)}} \).

Finally, \( \norm{\A^{1/2}\mathbf{e}}^2=\mathbf{e}\T\A\mathbf{e}
=\norm{\Z\mathbf{e}}^2+\mathbf{e}\T\bP\mathbf{e} \), so both terms tend to zero,
which gives \( \Z\bgamma^{(t)}\to\Z\hat{\bgamma} \) and
\( Q(\bgamma^{(t)})\to Q(\hat{\bgamma}) \). For the coefficients themselves,
\[
\bgamma^{(t+1)}-\bgamma^{(t)}=\mathbf{T}\mathbf{e}^{(t)}-\mathbf{e}^{(t)}
=-\mathbf{E}^{-1}\A^{1/2}\bigl(\A^{1/2}\mathbf{e}^{(t)}\bigr),
\]
whose norm is at most \( \norm{\mathbf{E}^{-1}\A^{1/2}}\rho^{t}\norm{\A^{1/2}\mathbf{e}^{(0)}} \).
These increments are summable, so \( \bgamma^{(t)} \) converges to a point
\( \bgamma^{\infty} \) with \( \A\bgamma^{\infty}=\Z\T\y \): a minimizer of
\( Q \).

(d) If \( \A \) is positive definite then \( \A^{1/2} \) is nonsingular,
\( \C(\A^{1/2})=\Real^{d} \), and the argument above gives
\( \norm{\mathbf{e}^{(t)}}\le c\,\rho^{t}\norm{\mathbf{e}^{(0)}} \) for a constant \( c \)
depending only on the condition number of \( \A \).
:::

The unpenalized case is the picture [Chapter 6](../ch06-projections/index.html) draws.

::: {#cor-add-projection}
[Backfitting as alternating projection]

If \( \lambda_j=0 \) for every \( j \) and each \( \Z_j \) has full column rank, then each
\( \bS_j \) is the orthogonal projection onto
\( \mathcal{V}_j=\C(\Z_j) \) by @thm-proj-sym-idem, and backfitting converges to the
orthogonal projection of \( \y \) onto
\( \mathcal{V}_0+\mathcal{V}_1+\dots+\mathcal{V}_q \). The fitted components are
unique iff the only way to write \( \bzero \) as \( \sum_j\bv_j \) with
\( \bv_j\in\mathcal{V}_j \) is with every \( \bv_j=\bzero \).
:::

::: {.proof}
With \( \bP=\bzero \) the criterion is \( \norm{\y-\Z\bgamma}^2 \) and
\( \A=\Z\T\Z \), nonnegative definite with nonsingular diagonal blocks, which is what
@thm-add-backfitting(c) needs; its conclusion is \( \Z\bgamma^{(t)}\to\Z\hat{\bgamma} \),
and \( \Z\hat{\bgamma}=\M\y \) is the projection onto
\( \C(\Z)=\sum_j\mathcal{V}_j \) by @thm-proj-ls-projection. The components are the
coordinates of that projection in the sum, unique precisely when the sum is direct, which is
the stated condition.
:::

## What it costs, and when to use it

::: {#exm-add-backfit}
[Backfitting on two P-spline terms]

Two covariates on \( [0,1] \) come from a bivariate normal pair squashed through a
hyperbolic tangent, so one correlation \( \rho \) controls their dependence; the mean is
\( 3+f_1(z_1)+f_2(z_2) \), cubic-like in the first component and a bump in the second, with
normal errors of standard deviation \( 0.4 \). With \( n=400 \), \( 15 \)
cubic B-splines per term, second-difference penalties and
\( \lambda_1=\lambda_2=8.0 \), backfitting and one direct solve of @eq-add-normal agree to
\( 9.4\times 10^{-14} \) in the fitted values, as they must, using
\( 11.64 \) effective degrees of freedom: \( 1 \) for the intercept,
\( 5.34 \) and \( 5.30 \) for the two components.

The work is not shared. At \( \rho=0.0 \) the backfitting change falls below
\( 10^{-12} \) after \( 9 \) cycles, shrinking by a factor of about
\( 0.015 \) per cycle; at \( \rho=0.5 \) it takes \( 21 \)
cycles and the factor is \( 0.240 \); at \( \rho=0.9 \) it takes
\( 121 \) cycles and the factor is \( 0.674 \).
[Figure 44.2.1](02-backfitting.html#fig-add-backfitting) shows the fitted components and the
three convergence curves, straight on a logarithmic scale as @thm-add-backfitting(c)
predicts.
:::

::: {when-format="html"}
![**Figure 44.2.1.** Backfitting on two P-spline terms. (a) and (b) the fitted components against
the truth. (c) the largest change per cycle at three levels of dependence: geometric decay,
at a rate that worsens as the terms grow harder to tell apart.](backfitting.svg){#fig-add-backfitting width=100%}
:::

::: {when-format="pdf"}
![Backfitting on two P-spline terms. (a) and (b) the fitted components against
the truth. (c) the largest change per cycle at three levels of dependence: geometric decay,
at a rate that worsens as the terms grow harder to tell apart.](backfitting.pdf){width=100%}
:::

The blocks used throughout the chapter are built once and reused. A P-spline term is a
basis, a difference penalty and the constraint of @prp-add-identifiability(c) absorbed into
it:

```{.python .run #cell-psplines-basis}
import numpy as np
from scipy.interpolate import BSpline

def bspline_basis(x, n_basis=15, degree=3, lo=None, hi=None):
    """n_basis B-splines of the given degree on equally spaced knots covering the data."""
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    inner = np.linspace(lo, hi, n_basis - degree + 1)
    h = inner[1] - inner[0]
    knots = np.r_[lo - h * np.arange(degree, 0, -1), inner, hi + h * np.arange(1, degree + 1)]
    return BSpline.design_matrix(np.clip(x, lo, hi), knots, degree, extrapolate=False).toarray()

def difference_penalty(d, order=2):
    """K = D'D for the order-th difference matrix D: the P-spline penalty."""
    D = np.diff(np.eye(d), n=order, axis=0)
    return D.T @ D

def spline_term(x, n_basis=15, order=2, weight=None):
    """A P-spline term: its design, its penalty, and the map from new covariate values to
    design rows. The constraint sum_i f(x_i) = 0 is absorbed into the basis, and
    weight = u turns the term into the varying-coefficient term f(x) u."""
    lo, hi = x.min(), x.max()
    B = bspline_basis(x, n_basis, lo=lo, hi=hi)
    Q, _ = np.linalg.qr(B.sum(axis=0)[:, None], mode="complete")
    U = Q[:, 1:]                                  # a basis of {gamma : sum_i f(x_i) = 0}
    Z = B if weight is None else B * np.asarray(weight)[:, None]

    def design(xnew, wnew=None):
        Bn = bspline_basis(np.asarray(xnew, float), n_basis, lo=lo, hi=hi)
        return (Bn if wnew is None else Bn * np.asarray(wnew)[:, None]) @ U

    return Z @ U, U.T @ difference_penalty(n_basis, order) @ U, design

def spline_block(x, n_basis=15, order=2, weight=None):
    """The (design, penalty) pair of spline_term, without the prediction map."""
    return spline_term(x, n_basis, order, weight)[:2]
```

The data of @exm-add-backfit are simulated from two dependent covariates:

```{.python .run #cell-backfitting-data}
def simulate(n, rho, seed=20440):
    """Two dependent covariates on [0, 1] and an additive mean with normal errors."""
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 2)) @ np.array([[1.0, rho], [0.0, np.sqrt(1 - rho**2)]])
    z1, z2 = 0.5 * (1 + np.tanh(0.8 * u[:, 0])), 0.5 * (1 + np.tanh(0.8 * u[:, 1]))
    f1 = 16.0 * (z1 - 0.5) ** 3 - 0.6 * np.sin(4.5 * z1)
    f2 = 1.4 * np.exp(-25 * (z2 - 0.35) ** 2) - 0.5 * z2
    f1, f2 = f1 - f1.mean(), f2 - f2.mean()
    y = 3.0 + f1 + f2 + rng.normal(scale=0.4, size=n)
    return z1, z2, f1, f2, y
```

Backfitting is a dozen lines, and the direct solve another dozen:

```{.python .run #cell-backfitting-backfit}
def backfit(y, blocks, lambdas, tol=1e-12, maxit=600):
    """Backfitting: cycle over the terms, refitting each to the current partial residual."""
    A = [Z.T @ Z + lam * K for (Z, K), lam in zip(blocks, lambdas)]
    f = [np.zeros_like(y) for _ in blocks]
    b0, history = y.mean(), []
    for _ in range(maxit):
        change = 0.0
        for j, (Z, _) in enumerate(blocks):
            r = y - b0 - sum(f[k] for k in range(len(blocks)) if k != j)
            new = Z @ np.linalg.solve(A[j], Z.T @ r)      # smooth the partial residual
            change = max(change, np.abs(new - f[j]).max())
            f[j] = new
        history.append(change)
        if change < tol:
            break
    return b0, f, np.array(history)

def direct_fit(y, blocks, lambdas):
    """One solve of the penalized normal equations, with an intercept column in front."""
    Z = np.column_stack([np.ones(len(y))] + [Zj for Zj, _ in blocks])
    K = np.zeros((Z.shape[1], Z.shape[1]))
    start = 1
    for (_, Kj), lam in zip(blocks, lambdas):
        d = Kj.shape[0]
        K[start:start + d, start:start + d] = lam * Kj
        start += d
    C = np.linalg.solve(Z.T @ Z + K, Z.T)                # gamma = C y
    edf = [np.trace(Zj @ C[s:s + Zj.shape[1]]) for Zj, s in
           zip([Zj for Zj, _ in blocks], np.cumsum([1] + [Zj.shape[1] for Zj, _ in blocks]))]
    return Z @ (C @ y), np.trace(Z @ C), edf
```

The two agree, and @prp-add-df falls out of the direct solve:

```{.python .run #cell-backfitting-compare}
n, lam = 400, 8.0
z1, z2, f1_true, f2_true, y = simulate(n, rho=0.5)
blocks = [spline_block(z1), spline_block(z2)]
lambdas = [lam, lam]

b0, f, history = backfit(y, blocks, lambdas)
fitted, edf_total, edf = direct_fit(y, blocks, lambdas)
print("cycles:", len(history))
print("largest disagreement:", np.abs(fitted - (b0 + f[0] + f[1])).max())
print(f"edf {edf_total:.3f} = 1 + {edf[0]:.3f} + {edf[1]:.3f}")
```

The direct solve costs \( O(nd^2+d^3) \) and hands over \( \A^{-1} \), hence @prp-add-df
and the covariances, free. Backfitting costs \( O(\sum_jnd_j^2) \) a cycle, never forms a
\( d\times d \) matrix and accepts smoothers no penalty produces, but delivers only the fit
and takes more cycles the more the terms depend on one another. Modern software solves
directly with sparse factorizations; backfitting survives where \( d \) is enormous and
where the terms are unusual. It is also close kin to the boosting algorithms of
[Chapter 30](../ch30-regularization-boosting/index.html), which fit a block to the current
residual as backfitting does but damp the update and pick the block greedily instead of
cycling (@thm-reg-boosting).

## Concurvity

The rate in @exm-add-backfit worsens with \( \rho \) for the reason collinearity hurts a
linear model: two terms that can each explain the same feature compete for it. Here the
competition is between *function spaces*, not columns, and it has its own name.

::: {#def-add-concurvity-index}
[Concurvity]

Let \( \mathcal{V}_0=\C(\X) \) and \( \mathcal{V}_j=\C(\Z_j) \) be the model spaces of the
terms, all centred so that \( \mathcal{V}_j\subseteq\bone\perpc \) for \( j\ge1 \). The
terms exhibit **exact concurvity** if there exist \( \bv_j\in\mathcal{V}_j \), not all
zero, with \( \sum_{j\ge0}\bv_j=\bzero \). The **concurvity index** of term \( j \) is
\[
\kappa_j=\max\Bigl\{\frac{\norm{\M_{-j}\bv}^2}{\norm{\bv}^2}:
\bv\in\mathcal{V}_j,\ \bv\ne\bzero\Bigr\}\in[0,1],
\]{#eq-add-concurvity}

where \( \M_{-j} \) is the orthogonal projection onto
\( \sum_{k\ne j}\mathcal{V}_k \). Equivalently \( \kappa_j \) is the largest squared
canonical correlation between \( \mathcal{V}_j \) and \( \sum_{k\ne j}\mathcal{V}_k \), and
\( \kappa_j=1 \) iff there is exact concurvity involving term \( j \).
:::

::: {#prp-add-concurvity}
[What concurvity does]

With the notation of @def-add-concurvity-index and \( \lambda_j=0 \) for all \( j \):

::: {.enumerate options="label=(\alph*)"}
1. \( \A=\Z\T\Z \) is singular iff there is exact concurvity. The fitted values are
   unaffected — they are the projection onto \( \sum_j\mathcal{V}_j \) — but the
   decomposition into components is not unique.

2. Assume in addition that there is no exact concurvity and that each \( \Z_j \) has full
   column rank. Then there is a direction \( \bv\in\mathcal{V}_j \), \( \bv\ne\bzero \),
   with
   \[
   \Var\bigl(\bv\T\hat{\mathbf{f}}_j\bigr)=\frac{\sigma^2\norm{\bv}^2}{1-\kappa_j} ,
   \]{#eq-add-variance-inflation}

   against the value \( \sigma^2\norm{\bv}^2 \) that the same functional would have if
   \( \mathcal{V}_j \) were orthogonal to the other spaces. The inflation factor is
   \( 1/(1-\kappa_j) \), and in the single-column case
   \( \mathcal{V}_j=\spn\{\x_j\} \) it is the variance inflation factor of @def-col-vif,
   with \( \Var(\hat\beta_j)=\sigma^2/\{\norm{\x_j}^2(1-\kappa_j)\} \).

3. *(What a penalty repairs.)* With \( \lambda_j>0 \) for every \( j \) and
   \( \Null(\Z_j)\cap\Null(\bP_j)=\{\bzero\} \), the matrix \( \A \) is nonsingular
   and the fitted components are unique, whatever \( \kappa_j \) may be.
:::

:::

::: {.proof}
(a) \( \Z\T\Z \) is singular iff \( \Z\bgamma=\bzero \) for some \( \bgamma\ne\bzero \),
that is iff the vectors \( \bv_j=\Z_j\bgamma_j \) sum to zero without all being zero (if
some \( \bgamma_j\ne\bzero \) gives \( \bv_j=\bzero \), then \( \Z_j \) itself is rank
deficient, a within-term problem that the constraint of @prp-add-identifiability(c) and the
choice of basis exclude). The fitted vector is the projection onto \( \C(\Z) \), which does
not depend on the basis (@thm-proj-M-formula).

(b) Put \( \A_j=\Z_j\T\Z_j \), \( \tilde{\Z}_j=(\I-\M_{-j})\Z_j \) and
\( \tilde{\A}_j=\tilde{\Z}_j\T\tilde{\Z}_j \), both positive definite. Writing
\( \bv=\Z_j\mathbf{c} \), @eq-add-concurvity says that \( 1-\kappa_j \) is the minimum
of \( \mathbf{c}\T\tilde{\A}_j\mathbf{c}/\mathbf{c}\T\A_j\mathbf{c} \), attained at
a generalized eigenvector \( \mathbf{c} \) satisfying
\( \tilde{\A}_j\mathbf{c}=(1-\kappa_j)\A_j\mathbf{c} \) by @cor-mat-generalized-rayleigh.
Now @thm-proj-fwl with \( \X_1=\Z_j \) and \( \X_2 \)
the other blocks gives \( \hat{\bgamma}_j=\tilde{\A}_j^{-1}\tilde{\Z}_j\T\y \) and
\( \Cov(\hat{\bgamma}_j)=\sigma^2\tilde{\A}_j^{-1} \). Since
\( \bv\T\hat{\mathbf{f}}_j=\mathbf{c}\T\A_j\hat{\bgamma}_j \),
\[
\Var\bigl(\bv\T\hat{\mathbf{f}}_j\bigr)
=\sigma^2\,\mathbf{c}\T\A_j\tilde{\A}_j^{-1}\A_j\mathbf{c}
=\frac{\sigma^2}{1-\kappa_j}\,\mathbf{c}\T\A_j\mathbf{c}
=\frac{\sigma^2\norm{\bv}^2}{1-\kappa_j},
\]
because the eigenvector relation gives
\( \tilde{\A}_j^{-1}\A_j\mathbf{c}=\mathbf{c}/(1-\kappa_j) \). If
\( \mathcal{V}_j\perp\sum_{k\ne j}\mathcal{V}_k \), then \( \tilde{\A}_j=\A_j \),
\( \kappa_j=0 \) and the variance is \( \sigma^2\norm{\bv}^2 \). For a single column,
\( \tilde{\A}_j=\norm{(\I-\M_{-j})\x_j}^2=\norm{\x_j}^2(1-\kappa_j) \),
which is @exr-proj-vif.

(c) Under the stated condition, \( \bgamma\T\A\bgamma=\norm{\Z\bgamma}^2+\bgamma\T\bP\bgamma=0 \)
forces \( \Z_j\bgamma_j=\bzero \) and \( \bP_j\bgamma_j=\bzero \) for every
\( j \), hence \( \bgamma=\bzero \).
:::

Part (b) is what "concurvity is collinearity for function spaces" means precisely: the
worst direction inside the \( j \)th model space pays exactly the inflation factor
\( 1/(1-\kappa_j) \), and \( \kappa_j \) is computed from the design alone, before any
data are seen.

::: {.remark}
[A penalty restores uniqueness, not precision]

Part (c) says nothing about variances, and the usual informal claim — that they blow up as
\( \kappa_j\to1 \) — is not a theorem: with \( \lambda_j>0 \) fixed, \( \A \) stays
nonsingular and the variances of @prp-add-df(c) stay bounded. What the penalty does is take
over the identification: the closer \( \kappa_j \) is to one, the less of
\( \hat{\mathbf{f}}_j \) is settled by the data and the more by \( \bP_j \).
Exact variances cost one matrix inverse; compute them and \( \kappa_j \), and report
both.
:::

::: {#exm-add-concurvity}
[Two covariates, one relation]

Let \( z_1 \) be uniform on \( [0,1] \) and \( z_2=z_1^{1.5}+\text{noise} \), and fit
two P-spline terms with \( \lambda=5.0 \) to \( n=300 \) observations. The
design is fixed and the criterion quadratic, so @prp-add-df(c) gives the exact standard
deviation of each component without simulation. As the noise falls from \( 0.30 \) through
\( 0.10 \) to \( 0.03 \), the concurvity index of the first term climbs from
\( 0.609 \) to \( 0.933 \) to
\( 0.994 \), and the standard deviation of the first component at
\( z_1=0.3 \) climbs from \( 0.049 \) to
\( 0.087 \) to \( 0.207 \), an inflation of
\( 4.2 \); over the same range the standard deviation of the *fitted mean*
there moves only from \( 0.094 \) to \( 0.080 \). The fit
is fine; its attribution to the two covariates is not.
:::

::: {when-format="html"}
![**Figure 44.2.2.** Concurvity. (a) and (b) twenty-five replicate estimates of the first component
(truth dashed) under weak and strong dependence. (c) its exact standard deviation at
\( z_1=0.3 \) against the concurvity index @eq-add-concurvity.](concurvity.svg){#fig-add-concurvity width=100%}
:::

::: {when-format="pdf"}
![Concurvity. (a) and (b) twenty-five replicate estimates of the first component
(truth dashed) under weak and strong dependence. (c) its exact standard deviation at
\( z_1=0.3 \) against the concurvity index @eq-add-concurvity.](concurvity.pdf){width=100%}
:::

```{.python .run #cell-concurvity-setup}
def design(n, noise, seed=44012):
    """z2 is a noisy monotone function of z1: small noise means strong concurvity."""
    rng = np.random.default_rng(seed)
    z1 = rng.uniform(0, 1, n)
    z2 = np.clip(z1**1.5 + noise * rng.normal(size=n), 0.001, 0.999)
    return z1, z2

def concurvity(Z1, Z2):
    """The largest squared canonical correlation between the two column spaces."""
    Q1, _ = np.linalg.qr(Z1)
    Q2, _ = np.linalg.qr(Z2)
    return np.linalg.svd(Q1.T @ Q2, compute_uv=False)[0] ** 2
```

::: {.warning}
Concurvity is harder to spot than collinearity, because it is invisible in pairwise
scatterplots and correlations: two covariates can be uncorrelated and still have model
spaces that nearly intersect, since \( \kappa_j \) involves *all* smooth functions of each,
not just the linear ones (@exr-add-concurvity-uncorrelated). Report \( \kappa_j \)
alongside the components, and treat a value near one as a warning that the individual
curves, and any story told about them, are not identified by the data even though the
predictions are.
:::

## Exercises

### A. Check your understanding

::: {#exr-add-backfit-one}
[A1]

Run backfitting with \( q=1 \). How many cycles does it take, and why?
:::

::: {.solution}
One. With a single block the update sets \( \bgamma_1=\A_{11}^{-1}\Z_1\T\y \), which
already solves @eq-add-normal; a second cycle reproduces it. In @thm-add-backfitting(b),
\( \bL=\bzero \), so \( \mathbf{T}=\I-\A^{-1}\A=\bzero \).
:::

### B. Practice

::: {#exr-add-two-block-rate}
[B1]

For \( q=2 \) with no parametric part, show that the Gauss–Seidel iteration matrix
\( \mathbf{T} \) of @eq-add-gauss-seidel has spectral radius equal to the largest
eigenvalue of \( \A_{11}^{-1}\A_{12}\A_{22}^{-1}\A_{21} \), and that for unpenalized
one-column blocks this is the squared correlation between the two columns.
:::

::: {.solution}
With \( \bD+\bL=\begin{psmallmatrix}\A_{11}&\bzero\\\A_{21}&\A_{22}\end{psmallmatrix} \)
and \( \bL\T=\begin{psmallmatrix}\bzero&\A_{12}\\\bzero&\bzero\end{psmallmatrix} \),
\[
\mathbf{T}=-(\bD+\bL)^{-1}\bL\T
=-\begin{pmatrix}\bzero&\A_{11}^{-1}\A_{12}\\
\bzero&-\A_{22}^{-1}\A_{21}\A_{11}^{-1}\A_{12}\end{pmatrix},
\]
block triangular with eigenvalues \( 0 \) and those of
\( \A_{22}^{-1}\A_{21}\A_{11}^{-1}\A_{12} \), whose nonzero eigenvalues are those of
\( \A_{11}^{-1}\A_{12}\A_{22}^{-1}\A_{21} \). With one column per block and no penalty this
is \( (\x_1\T\x_2)^2/(\norm{\x_1}^2\norm{\x_2}^2) \), the squared correlation of centred
columns.
:::

::: {#exr-add-modified}
[B2]

*Modified backfitting.* Suppose the smoothers \( \bS_j \) are not of the form
@eq-add-blocksmoother but are symmetric with eigenvalues in \( [0,1] \) and satisfy
\( \bS_j\bone=\bone \). Show that replacing \( \bS_j \) by
\( \bS_j^{c}=(\I-n^{-1}\bone\bone\T)\bS_j \) makes every fitted component sum to zero at
each step, and that \( \A \) in @thm-add-backfitting is then replaced by a matrix with the
same structure acting on \( \bone\perpc \).
:::

::: {#exr-add-concurvity-uncorrelated}
[B3]

Let \( z_1 \) be uniform on \( [-1,1] \) and \( z_2=z_1^2 \). Show that
\( \Cov(z_1,z_2)=0 \), so the two covariates are uncorrelated, and yet the two model spaces
of quadratic functions exhibit exact concurvity. What does this say about diagnosing
concurvity by a correlation matrix?
:::

::: {.solution}
\( \E(Z_1Z_1^2)=\E(Z_1^3)=0=\E(Z_1)\E(Z_1^2) \), so the covariance is zero. But
\( g_1(z_1)=z_1^2 \) lies in the first space, \( g_2(z_2)=z_2 \) in the second, and
\( g_1(z_{i1})-g_2(z_{i2})=0 \) for every \( i \): a nontrivial additive relation, so
\( \kappa_1=\kappa_2=1 \). A correlation sees only the linear part of the relation between
covariates; concurvity sees all of it, so it must be measured on the model spaces.
:::

### C. Going deeper

::: {#exr-add-backfit-population}
[C1]

*Population backfitting.* Let \( Z_1,Z_2 \) be random with \( \E(Y^2)<\infty \), and
consider minimizing \( \E\{Y-g_1(Z_1)-g_2(Z_2)\}^2 \) over centred square-integrable
\( g_1,g_2 \). Show that the minimizers satisfy
\( g_1(z)=\E\{Y-g_2(Z_2)\mid Z_1=z\}-c \) and symmetrically, and that the sample algorithm
of this section is the empirical version with conditional expectation replaced by a
smoother.
:::

