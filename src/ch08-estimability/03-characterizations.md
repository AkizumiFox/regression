# Characterizations of estimability

@thm-est-estimable-identifiable reduced estimability to a single condition,
\( \blambda\in\C(\X\T) \). In practice that condition appears in several forms, depending on
what is at hand. We may have the rows of \( \X \), the matrix \( \X\T\X \) and one of its
generalized inverses, a collection of least squares solutions, or an SVD. This section
collects the equivalent forms in one theorem, derives the estimator of an estimable function
with its variance, and then turns to the practical question of deciding estimability in
floating-point arithmetic, where exact rank is not available.

## Eight equivalent conditions

::: {#thm-est-characterization}
[Characterizations of estimability]

Let \( \X \) be \( n\times p \), let \( \G \) denote a generalized inverse of \( \X\T\X \), and let
\( \blambda\in\Real^p \). The following are equivalent:

::: {.enumerate options="label=(\alph*)"}
1. \( \blambda\T\bbeta \) is estimable;

2. \( \blambda\in\C(\X\T) \): there is \( \boldsymbol{\uprho}\in\Real^n \) with \( \blambda\T=\boldsymbol{\uprho}\T\X \);

3. \( \blambda\perp\Null(\X) \): \( \blambda\T\bv=0 \) whenever \( \X\bv=\bzero \);

4. \( \blambda\in\C(\X\T\X) \): the system \( \X\T\X\br=\blambda \) has a solution;

5. \( \blambda\T\G\X\T\X=\blambda\T \) for some generalized inverse \( \G \);

6. \( \blambda\T\G\X\T\X=\blambda\T \) for every generalized inverse \( \G \);

7. for every \( \y \), \( \blambda\T\bb \) has the same value for all solutions \( \bb \) of the normal
           equations \( \X\T\X\bb=\X\T\y \) (@thm-proj-normal-equations);

8. \( \rank\begin{pmatrix}\X\\\blambda\T\end{pmatrix}=\rank(\X) \).
:::

:::

::: {.proof}
(a)\( \Leftrightarrow \)(b) is @thm-est-estimable-identifiable. (b)\( \Leftrightarrow \)(c) holds because
\( \C(\X\T)=\Null(\X)\perpc \) (@cor-proj-gram-colspace(d)). (b)\( \Leftrightarrow \)(d) holds because
\( \C(\X\T\X)=\C(\X\T) \) (@cor-proj-gram-colspace(b)).
(d)\( \Rightarrow \)(f): if \( \blambda=\X\T\X\br \), then for every generalized inverse
\( \blambda\T\G\X\T\X=\br\T(\X\T\X\G\X\T\X)=\br\T\X\T\X=\blambda\T \).
(f)\( \Rightarrow \)(e) holds because generalized inverses exist (@prp-proj-ginverse-exists).
(e)\( \Rightarrow \)(d): transposing \( \blambda\T\G\X\T\X=\blambda\T \) gives
\( \X\T\X(\G\T\blambda)=\blambda \), so \( \br=\G\T\blambda \) solves the system.
(b)\( \Leftrightarrow \)(g) is @thm-proj-invariant-functions.
(b)\( \Leftrightarrow \)(h): the rank of the stacked matrix is the dimension of its row space,
\( \C([\X\T,\blambda]) \). This contains \( \C(\X\T) \), and it has the same dimension iff it equals
\( \C(\X\T) \), that is, iff \( \blambda\in\C(\X\T) \).
:::

Conditions (e) and (f) are @cor-mat-invariant applied to the normal equations. Condition (g)
is the one that matters for interpretation: an estimable function is one whose estimate does not
depend on the arbitrary choice of solution. Condition (c) is the one to use when proving that
something is *not* estimable, since a single null vector \( \bv \) with \( \blambda\T\bv\ne0 \) settles
the question. Condition (h) is the most direct computational test, but in floating point it
inherits all the difficulties of deciding a rank.

## The estimator and its variance

For an estimable function, every least squares solution produces the same estimate, and that
estimate has a clean description.

::: {#prp-est-ls-estimator}
[Least squares estimation of estimable functions]

Let \( \blambda\T\bbeta \) be estimable, with \( \blambda=\X\T\boldsymbol{\uprho}=\X\T\X\br \), let \( \hbeta \) be any
least squares solution and \( \G \) any generalized inverse of \( \X\T\X \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \blambda\T\hbeta=\boldsymbol{\uprho}\T\M\y=\br\T\X\T\y=\blambda\T\G\X\T\y \);

2. \( \E(\blambda\T\hbeta)=\blambda\T\bbeta \);

3. if \( \Cov(\Y)=\sigma^2\I \), then
    \[
    \Var(\blambda\T\hbeta)=\sigma^2\blambda\T\G\blambda=\sigma^2\br\T\blambda=\sigma^2\norm{\M\boldsymbol{\uprho}}^2 ,
    \]
    the same for every \( \G \);

4. more generally, if \( \bLambda\T\bbeta \) is estimable, then
           \( \Cov(\bLambda\T\hbeta)=\sigma^2\bLambda\T\G\bLambda \) for every \( \G \).
:::

:::

::: {.proof}
(a) The first equality is @thm-proj-invariant-functions. For the second, the normal equations give
\( \br\T\X\T\y=\br\T\X\T\X\hbeta=\blambda\T\hbeta \). The last is the second with the solution
\( \hbeta=\G\X\T\y \) and \( \br=\G\T\blambda \) (see the proof of @thm-est-characterization).
(b)–(d) were proved in [Chapter 7](../ch07-optimality/index.html): unbiasedness and
\( \Var(\blambda\T\hbeta)=\sigma^2\blambda\T\G\blambda \) are @thm-opt-gauss-markov(a), and the
matrix form \( \Cov(\bLambda\T\hbeta)=\sigma^2\bLambda\T\G\bLambda \) is part (d) of the same theorem.
Neither depends on \( \G \), because \( \M=\X\G\X\T \) for every \( \G \) (@thm-proj-M-formula).
The second form in (c) follows from
\( \blambda\T\G\blambda=\br\T\X\T\X\G\X\T\X\br=\br\T\X\T\X\br=\br\T\blambda \).
The third is \( \boldsymbol{\uprho}\T\M\boldsymbol{\uprho}=\norm{\M\boldsymbol{\uprho}}^2 \).
:::

The formula \( \sigma^2\blambda\T\G\blambda \) is the rank-deficient version of the familiar
\( \sigma^2\blambda\T(\X\T\X)^{-1}\blambda \) (@thm-lm-moments). The generalized inverse is arbitrary,
but the quadratic form in an estimable \( \blambda \) is not. The condition
\( \blambda\in\C(\X\T) \) under which [Chapter 7](../ch07-optimality/index.html) proved optimality is
exactly estimability, so the Gauss–Markov theorem and the normal theory (@thm-opt-sampling) apply to
every estimable \( \blambda\T\hbeta \) unchanged.

## What a generalized inverse reports for a nonestimable function

Software that solves the normal equations with a particular generalized inverse will
print \( \blambda\T\G\X\T\y \) and \( \hat{\sigma}^2\blambda\T\G\blambda \) for any \( \blambda \) whatever. To see what
these numbers mean when \( \blambda \) is not estimable, look at the matrix
\[
\bH_{\G}=\G\X\T\X .
\]{#eq-est-H}

::: {#prp-est-H}
[The matrix \( \G\X\T\X \)]

For any generalized inverse \( \G \) of \( \X\T\X \), with \( \bH_{\G} \) as in @eq-est-H:

::: {.enumerate options="label=(\alph*)"}
1. \( \bH_{\G} \) is idempotent, \( \Null(\bH_{\G})=\Null(\X) \), \( \rank(\bH_{\G})=r \), and the rows of
           \( \bH_{\G} \) span \( \C(\X\T) \);

2. \( \blambda\T\bbeta \) is estimable iff \( \blambda\T\bH_{\G}=\blambda\T \);

3. the solution \( \hbeta_{\G}=\G\X\T\Y \) has \( \E(\hbeta_{\G})=\bH_{\G}\bbeta \). So for *every* \( \blambda \),
           \( \blambda\T\hbeta_{\G} \) is unbiased for the estimable function \( \blambda\T\bH_{\G}\bbeta \),
           which equals \( \blambda\T\bbeta \) iff \( \blambda\T\bbeta \) is estimable;

4. for the Moore–Penrose inverse \( \G=(\X\T\X)^+ \), \( \bH_{\G} \) is the orthogonal projection onto
           \( \C(\X\T) \).
:::

:::

::: {.proof}
(a) \( \bH_{\G}^2=\G(\X\T\X\G\X\T\X)=\G\X\T\X=\bH_{\G} \). By @prp-mat-ginverse-props(a),
\( \Null(\G\X\T\X)=\Null(\X\T\X) \), which is \( \Null(\X) \) by @cor-proj-gram-colspace(a).
Rank–nullity gives rank \( r \). The row space of \( \bH_{\G} \) is
\( \C(\bH_{\G}\T)=\Null(\bH_{\G})\perpc=\Null(\X)\perpc=\C(\X\T) \).
(b) is @thm-est-characterization(e),(f).
(c) \( \E(\G\X\T\Y)=\G\X\T\X\bbeta \). The coefficient vector \( \bH_{\G}\T\blambda \) lies in the row space
of \( \bH_{\G} \), which is \( \C(\X\T) \), so the function is estimable. It equals \( \blambda\T\bbeta \) for
all \( \bbeta \) iff \( \bH_{\G}\T\blambda=\blambda \), which is (b).
(d) By the Penrose conditions (@def-mat-moore-penrose), \( (\X\T\X)^+\X\T\X \) is symmetric and
idempotent. It is therefore an orthogonal projection (@thm-proj-sym-idem), and its null space is
\( \Null(\X) \) by (a). So it projects onto \( \Null(\X)\perpc=\C(\X\T) \).
:::

Part (c) answers the question. For a nonestimable \( \blambda \), the printed number is a perfectly good
unbiased estimate, but of a *different* function, \( \blambda\T\bH_{\G}\bbeta \), and which function
depends on the software's choice of \( \G \). The printed “variance” \( \sigma^2\blambda\T\G\blambda \)
is the variance of that statistic when \( \G \) is symmetric and reflexive, as the Moore–Penrose
inverse is, but in general it need not be the variance of anything (@exr-est-any-variance). Part (d) gives the cleanest picture. With the
Moore–Penrose inverse, \( \bbeta=\bH\bbeta+(\I-\bH)\bbeta \) splits the coefficient vector into
its component in the estimable space, which the data determine, and its component in
\( \Null(\X) \), which they cannot see. The minimum-norm solution estimates the first
component and sets the second to zero.

::: {#exm-est-oneway-checks}
[Checking estimability in a one-way layout]

Eighteen observations fall into four groups of sizes \( 4,6,3,5 \), with
\( \X=[\bone,\bz_1,\bz_2,\bz_3,\bz_4] \) and \( \bbeta=(\mu,\alpha_1,\dots,\alpha_4)\T \), so
\( \rank(\X)=4 \) and \( \Null(\X) \) is spanned by \( (1,-1,-1,-1,-1)\T \). By condition (c), \( \blambda\T\bbeta \)
is estimable iff \( \lambda_0=\lambda_1+\lambda_2+\lambda_3+\lambda_4 \). The listing tests five
vectors in two ways: condition (c) through the SVD of \( \X \), and condition (e) with three
generalized inverses (Moore–Penrose, the one that deletes the intercept row and column, and a
non-symmetric one built as in @prp-mat-ginverse-props(e)). The two tests always agree, and
condition (e) holds for all three inverses or for none, as (f) requires.

| \( \blambda\T\bbeta \) | estimable | \( \blambda\T\G\blambda \): Moore–Penrose | drop intercept | non-symmetric |
|---|---|---|---|---|
| \( \mu+\alpha_1 \) | yes | 0.250 | 0.250 | 0.250 |
| \( \alpha_1-\alpha_2 \) | yes | 0.417 | 0.417 | 0.417 |
| \( \alpha_1+\alpha_2-2\alpha_3 \) | yes | 1.750 | 1.750 | 1.750 |
| \( \alpha_1 \) | no | 0.188 | 0.250 | 0.630 |
| \( \alpha_1+\alpha_2+\alpha_3+\alpha_4 \) | no | 0.038 | 0.950 | 2.767 |

For the estimable rows the quadratic form is the same for all three and matches the direct
computation: \( \Var(\bar{y}_1-\bar{y}_2)=\sigma^2(\tfrac14+\tfrac16) \), for instance. For
\( \alpha_1 \) the “variances” disagree. With the drop-intercept inverse, the solution is
\( \hbeta_{\G}=(0,\bar{y}_1,\dots,\bar{y}_4)\T \), so the reported “\( \hat{\alpha}_1 \)” is \( \bar{y}_1 \), which
estimates \( \mu+\alpha_1 \). Its reported variance \( 0.250\,\sigma^2 \) is the variance
of \( \bar{y}_1 \), which is correct, but \( \bar{y}_1 \) estimates \( \mu+\alpha_1 \), not \( \alpha_1 \).
:::

```{.python .run #cell-estimability-test-layout}
import numpy as np

groups = np.repeat([0, 1, 2, 3], [4, 6, 3, 5])        # four groups, sizes 4, 6, 3, 5
X = np.column_stack([np.ones(groups.size), np.eye(4)[groups]])  # [1, z1, z2, z3, z4]
A = X.T @ X
```

```{.python .run #cell-estimability-test-tests}
def is_estimable(lam, X, tol=None):
    """lambda is estimable iff it is orthogonal to N(X); use the SVD of X."""
    U, s, Vt = np.linalg.svd(X)
    if tol is None:
        tol = max(X.shape) * np.finfo(float).eps * s[0]
    r = int(np.sum(s > tol))
    V0 = Vt[r:].T                                  # orthonormal basis of N(X)
    return np.linalg.norm(V0.T @ lam) <= 1e-8 * max(1.0, np.linalg.norm(lam))

G_mp = np.linalg.pinv(A)                            # Moore-Penrose inverse of X^T X
G_ref = np.zeros((5, 5))
G_ref[1:, 1:] = np.linalg.inv(A[1:, 1:])            # delete the intercept row and column
rng = np.random.default_rng(3)
G_odd = G_mp + (np.eye(5) - G_mp @ A) @ rng.normal(size=(5, 5))   # non-symmetric g-inverse

candidates = {"mu + alpha_1": [1, 1, 0, 0, 0],
              "alpha_1 - alpha_2": [0, 1, -1, 0, 0],
              "alpha_1 + alpha_2 - 2 alpha_3": [0, 1, 1, -2, 0],
              "alpha_1": [0, 1, 0, 0, 0],
              "alpha_1 + ... + alpha_4": [0, 1, 1, 1, 1]}
for name, lam in candidates.items():
    lam = np.array(lam, dtype=float)
    H_checks = [np.allclose(lam @ G @ A, lam) for G in (G_mp, G_ref, G_odd)]
    variances = [lam @ G @ lam for G in (G_mp, G_ref, G_odd)]
    print(f"{name:30s} SVD test: {is_estimable(lam, X)!s:5s}",
          f" lambda'GX'X = lambda': {H_checks}",
          f" lambda'G lambda: {np.round(variances, 3)}")
```

## Deciding estimability in floating point

Every condition in @thm-est-characterization is a statement of exact linear algebra. On a
computer, the rank of \( \X \) must be decided with a tolerance, and the answer to “is \( \blambda \)
estimable?” inherits that tolerance. The most reliable test works with \( \X \) itself, not with
\( \X\T\X \). Compute the singular value decomposition \( \X=\bU\bD\V\T \). Declare the singular values
below a threshold \( \tau \) to be zero. With the usual choice
\( \tau=\max(n,p)\,\epsilon\,d_1 \), where \( d_1 \) is the largest singular value and
\( \epsilon\approx2.2\times10^{-16} \) is the machine precision, the columns \( \V_0 \) of \( \V \) that belong
to the discarded singular values form an orthonormal basis of the numerical null space. Then
\( \blambda \) is declared estimable iff \( \norm{\V_0\T\blambda} \) is negligible relative to
\( \norm{\blambda} \). This is condition (c), and it is what the function `is_estimable` in the listing
does. [Chapter 10](../ch10-computation/index.html) treats the SVD and its rank decisions in detail (@thm-cmp-svd-ls).

Two things can go wrong, and both matter.

*Exact dependence can be blurred by rounding.* If a column is computed, say as a sum of other
columns, rounding may leave a singular value of order \( 10^{-15} \) rather than zero. The tolerance
exists to absorb this.

*Near dependence is not dependence.* If the columns are nearly but not exactly dependent, every
\( \blambda \) is estimable in the exact sense, but some are estimable only with enormous variance.
Estimability is then a matter of degree, not a yes-or-no property.

::: {#exm-est-near-estimable}
[Estimable in principle, not in practice]

Take twelve equally spaced values \( t_i \) on \( [-1,1] \) and the model matrix with columns
\( \bone \), \( \mathbf{t} \) and \( \mathbf{t}+\delta\bw \), where \( w_i=\cos(\pi t_i) \) and \( \delta \) is a small perturbation.
At \( \delta=0 \) the last two columns coincide. Then \( \beta_1+\beta_2 \) is estimable and \( \beta_1 \)
is not. For \( \delta>0 \) the matrix has full rank, so everything is estimable. But
[Figure 8.3.1](#fig-est-near) shows how the variance of \( \hat{\beta}_1 \) grows as \( \delta\to0 \):
in proportion to \( 1/\delta^2 \) (the fitted slope on the log–log scale is
\( -2.00 \)). At \( \delta=10^{-4} \), the smallest singular value of \( \X \) is
\( 1.8\times 10^{-4} \), and \( \Var(\hat{\beta}_1)=1.56\times 10^{7}\,\sigma^2 \), while
\( \Var(\hat{\beta}_1+\hat{\beta}_2)=0.212\,\sigma^2 \) is unaffected. The function that is
estimable at \( \delta=0 \) stays well estimated, and the one that is not becomes useless long before
it becomes nonestimable.

The listing adds a warning about computation. At \( \delta=10^{-8} \), the SVD of \( \X \) still sees
three nonzero singular values, so \( \beta_1 \) is (barely) estimable. But the matrix \( \X\T\X \) has
eigenvalues equal to the *squared* singular values. Its smallest eigenvalue is about
\( 3\times10^{-16} \), below any sensible tolerance, and `pinv(X.T @ X)` silently treats \( \X\T\X \)
as having rank \( 2 \). It then reports \( 0.0529 \) as the “variance” of \( \hat{\beta}_1 \), which is the
Moore–Penrose quadratic form of a nonestimable function. The true variance, computed from the
SVD of \( \X \) itself, is more than \( 10^{16} \) times larger. Forming \( \X\T\X \) squares the condition number, as
[Section 6.10](../ch06-projections/10-computation.html) showed, and here that turns a
nearly nonestimable function into a numerically nonestimable one.
:::

::: {when-format="html"}
![**Figure 8.3.1.** Variances of three linear functions of \( \hbeta \) in the design with
columns \( \bone,\mathbf{t},\mathbf{t}+\delta\bw \). As \( \delta\to0 \), the functions that are not estimable at
\( \delta=0 \) (\( \beta_1 \) alone and \( \beta_1-\beta_2 \)) have variances growing like
\( \delta^{-2} \), while \( \beta_1+\beta_2 \), which is estimable at \( \delta=0 \), is
unaffected.](near_estimable.svg){#fig-est-near width=62%}
:::

::: {when-format="pdf"}
![Variances of three linear functions of \( \hbeta \) in the design with
columns \( \bone,\mathbf{t},\mathbf{t}+\delta\bw \). As \( \delta\to0 \), the functions that are not estimable at
\( \delta=0 \) (\( \beta_1 \) alone and \( \beta_1-\beta_2 \)) have variances growing like
\( \delta^{-2} \), while \( \beta_1+\beta_2 \), which is estimable at \( \delta=0 \), is
unaffected.](near_estimable.pdf){width=62%}
:::

```{.python .run #cell-estimability-test-near}
t = np.linspace(-1, 1, 12)
w = np.cos(np.pi * t)                                # a fixed direction to perturb along
def near_design(delta):
    """Columns 1, t and t + delta*w: exactly rank deficient at delta = 0."""
    return np.column_stack([np.ones_like(t), t, t + delta * w])

for delta in [1.0, 1e-2, 1e-4, 1e-8, 0.0]:
    Xd = near_design(delta)
    s = np.linalg.svd(Xd, compute_uv=False)
    G = np.linalg.pinv(Xd.T @ Xd)
    v_sum = np.array([0, 1, 1]) @ G @ [0, 1, 1]     # beta_1 + beta_2
    v_one = np.array([0, 1, 0]) @ G @ [0, 1, 0]     # beta_1 alone
    print(f"delta={delta:7.0e}  smallest singular value {s[-1]:8.1e}",
          f" Var(b1+b2)/s2 = {v_sum:8.3f}  Var(b1)/s2 = {v_one:9.2e}",
          f" beta_1 estimable? {is_estimable(np.array([0., 1, 0]), Xd)}")
```

The practical rule has two parts. Decide estimability with the SVD of \( \X \), never with
\( \X\T\X \). Then look at the variance, not only the verdict. A function that passes the
test with a variance of \( 10^{7}\sigma^2 \) has, for practical purposes, not been estimated.
Collinearity, in [Chapter 26](../ch26-collinearity/index.html), is the study of this middle ground.

## Exercises

### A. Check your understanding

::: {#exr-est-oneway-list}
[A1]

In @exm-est-oneway-checks, decide which of the following are estimable, using condition (c):
\( 2\mu+\alpha_1+\alpha_2 \); \( \alpha_1+\alpha_2-\alpha_3-\alpha_4 \); \( \mu \);
\( 4\mu+\alpha_1+\alpha_2+\alpha_3+\alpha_4 \). For each estimable one give its least squares estimator in
terms of the group means and its variance.
:::

### B. Practice

::: {#exr-est-what-is-reported}
[B1]

In @exm-est-oneway-checks, compute \( \bH_{\G} \) for the drop-intercept generalized inverse, and
identify the estimable function \( \blambda\T\bH_{\G}\bbeta \) that the software is really estimating when
it reports “\( \hat{\alpha}_k \)”. Do the same for the generalized inverse that deletes the row and column of
\( \alpha_4 \).
:::

::: {.solution}
With \( \G \) holding \( \diag(1/n_1,\dots,1/n_4) \) in its lower right block,
\( \G\X\T\X \) has first row zero and row \( k+1 \) equal to \( (1,\mathbf{e}_k\T) \), because the row of \( \X\T\X \) for
\( \alpha_k \) is \( n_k(1,\mathbf{e}_k\T) \). So \( \E(\hat{\alpha}_k)=\mu+\alpha_k \): the reported “effect” estimates the
group mean. If instead \( \alpha_4 \) is deleted, the solution sets \( \hat{\alpha}_4=0 \),
\( \hat{\mu}=\bar{y}_4 \) and \( \hat{\alpha}_k=\bar{y}_k-\bar{y}_4 \), so \( \E(\hat{\mu})=\mu+\alpha_4 \) and
\( \E(\hat{\alpha}_k)=\alpha_k-\alpha_4 \) for \( k\le3 \).
:::

::: {#exr-est-cov-pd}
[B2]

Let \( \bLambda\T\bbeta \) be estimable with \( \bLambda \) of full column rank \( q \), and suppose
\( \Cov(\Y)=\sigma^2\I \) with \( \sigma^2>0 \). Show that \( \Cov(\bLambda\T\hbeta)=\sigma^2\bLambda\T\G\bLambda \) is
positive definite, whichever generalized inverse is used. *Hint:* write \( \bLambda=\X\T\mathbf{P} \) and
show \( \M\mathbf{P} \) has full column rank.
:::

### C. Going deeper

::: {#exr-est-any-variance}
[C1]

Let \( \blambda\T\bbeta \) be nonestimable. Show that as \( \G \) ranges over all generalized inverses of
\( \X\T\X \), the number \( \blambda\T\G\blambda \) takes *every* real value, including negative ones. So the
“variance” that software attaches to a nonestimable function can be anything at all.
:::

::: {.solution}
Let \( \A=\X\T\X \) and \( \G_0 \) be one generalized inverse. By
@prp-mat-ginverse-props(e), \( \G=\G_0+(\I-\G_0\A)\bU \) is a generalized inverse for every \( \bU \), and
\[
\blambda\T\G\blambda=\blambda\T\G_0\blambda+\mathbf{u}\T\bU\blambda,\qquad \mathbf{u}=(\I-\G_0\A)\T\blambda .
\]
Because \( \blambda \) is not estimable, \( \blambda\T\G_0\A\ne\blambda\T \) (by (a)\( \Leftrightarrow \)(e) of @thm-est-characterization), so
\( \mathbf{u}\ne\bzero \). Also \( \blambda\ne\bzero \), since \( \bzero \) is in \( \C(\X\T) \). Take
\( \bU=t\,\mathbf{u}\blambda\T/(\norm{\mathbf{u}}^2\norm{\blambda}^2) \). Then
\( \blambda\T\G\blambda=\blambda\T\G_0\blambda+t \), which takes every real value as \( t \) varies.
:::

::: {#exr-est-projection-split}
[C2]

Let \( \bP \) be the orthogonal projection onto \( \C(\X\T) \). Show that the minimum-norm least squares
solution \( \hbeta^+ \) satisfies \( \E(\hbeta^+)=\bP\bbeta \), and that for any \( \blambda \),
\( \blambda\T\hbeta^+ \) is the least squares estimator of the estimable function
\( (\bP\blambda)\T\bbeta \). Interpret \( \bP\blambda \) as the estimable function “closest” to
\( \blambda\T\bbeta \).
:::
