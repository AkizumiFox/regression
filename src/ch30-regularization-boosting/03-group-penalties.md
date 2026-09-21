# Group penalties

A factor with \( k \) levels enters a linear model as \( k-1 \) indicator columns
([Section 8.5](../ch08-estimability/05-factor-coding.html), @def-est-coding), and a smooth effect as several basis columns. The question is whether
*the factor* matters. The lasso on the indicator columns can keep some levels and drop others, and its answer depends on the
coding. Yuan and Lin (2006) proposed a penalty that treats a block of columns as a unit.

## The group lasso

Partition the columns of \( \X \) into groups \( g=1,\dots,G \), with \( \X_g \) the \( n\times p_g \) block of group \( g \) and
\( \bb_g \) the matching block of coefficients, so that \( \X\bb=\sum_g\X_g\bb_g \).

::: {#def-reg-group-lasso}
[Group lasso]

For weights \( w_g>0 \), the **group lasso** criterion is
\[
Q(\bb)=\tfrac12\Bigl\lVert\y-\sum_{g=1}^G\X_g\bb_g\Bigr\rVert^2+\lambda\sum_{g=1}^Gw_g\norm{\bb_g},
\]{#eq-reg-group-lasso}

with Euclidean norms \( \norm{\bb_g} \), not squared. The usual weights are \( w_g=\sqrt{p_g} \). The columns of each group are
usually centred and **orthonormalized within the group** before fitting, so that \( \X_g\T\X_g=\I_{p_g} \). Different groups are
not made orthogonal to each other.
:::

With singleton groups and unit weights this is the lasso. In general it acts like the lasso *across* groups and like ridge
*within* them: \( \norm{\bb_g} \) has a corner only at \( \bb_g=\bzero \), so a group is either absent or present with, typically,
all its coefficients nonzero.

::: {#prp-reg-group-lasso}
[Group lasso and block soft thresholding]

::: {.enumerate options="label=(\alph*)"}
1. A minimizer of @eq-reg-group-lasso exists, and all minimizers have the same fitted vector.

2. *(KKT conditions.)* \( \hbeta \) is a minimizer iff, with \( \br=\y-\X\hbeta \), for every group \( g \),
   \[
   \X_g\T\br=\lambda w_g\,\frac{\hbeta_g}{\norm{\hbeta_g}}\quad\text{if }\hbeta_g\neq\bzero,
   \qquad \norm{\X_g\T\br}\le\lambda w_g\quad\text{if }\hbeta_g=\bzero .
   \]{#eq-reg-group-kkt}

3. Let \( \lambda>0 \). Then \( \hbeta=\bzero \) iff \( \lambda\ge\lambda_{\max}=\max_g\norm{\X_g\T\y}/w_g \).

4. *(Orthonormal groups.)* If \( \X_g\T\X_h=\I_{p_g} \) for \( g=h \) and \( \bzero \) for \( g\neq h \), the unique minimizer is **block soft
   thresholding** of \( \mathbf{z}_g=\X_g\T\y \):
   \[
   \hbeta_g=\Bigl(1-\frac{\lambda w_g}{\norm{\mathbf{z}_g}}\Bigr)_+\mathbf{z}_g .
   \]{#eq-reg-block-soft}

5. *(Block update; compare @prp-reg-coordinate.)* If \( \X_g\T\X_g=\I_{p_g} \), then minimizing @eq-reg-group-lasso over \( \bb_g \) with the other
   blocks fixed gives @eq-reg-block-soft with \( \mathbf{z}_g \) replaced by \( \mathbf{s}_g=\X_g\T\br_{(g)} \), where
   \( \br_{(g)}=\y-\sum_{h\neq g}\X_h\bb_h \) is the partial residual.
:::

:::

::: {.proof}
(a) The penalty is convex and continuous. Since \( \sum_g\norm{\bb_g}\ge\norm{\bb} \), it is at least
\( (\min_gw_g)\norm{\bb} \). Apply @thm-reg-existence(a)(ii) and (b).

(b) Combine @lem-reg-optimality with @lem-reg-norm-subgradient, using the groups as the blocks.

(c) At \( \hbeta=\bzero \), the residual is \( \y \), and @eq-reg-group-kkt says that \( \norm{\X_g\T\y}\le\lambda w_g \) for every \( g \).
If this holds, \( \bzero \) is a minimizer. Every minimizer then has fitted vector \( \bzero \) and, by @thm-reg-existence(b), penalty
\( 0 \), so it is \( \bzero \). If it fails, \( \bzero \) is not a minimizer.

(d) Under the assumption, \( \norm{\y-\X\bb}^2=\norm{\y}^2-2\sum_g\mathbf{z}_g\T\bb_g+\sum_g\norm{\bb_g}^2 \). So @eq-reg-group-lasso
separates into the problems \( \min_{\bb_g}\{\tfrac12\norm{\mathbf{z}_g-\bb_g}^2+\lambda w_g\norm{\bb_g}\} \), each strictly convex.
Let \( c=\lambda w_g \) and let \( \bb_g \) be given by @eq-reg-block-soft. If \( \norm{\mathbf{z}_g}\le c \), then \( \bb_g=\bzero \) and
\( \norm{\mathbf{z}_g-\bzero}\le c \). If \( \norm{\mathbf{z}_g}>c \), then \( \bb_g=(\norm{\mathbf{z}_g}-c)\mathbf{z}_g/\norm{\mathbf{z}_g}\neq\bzero \) and
\( \mathbf{z}_g-\bb_g=c\,\mathbf{z}_g/\norm{\mathbf{z}_g}=c\,\bb_g/\norm{\bb_g} \). In both cases the one-group version of @eq-reg-group-kkt holds.

(e) As a function of \( \bb_g \), the criterion is \( \tfrac12\norm{\br_{(g)}-\X_g\bb_g}^2+\lambda w_g\norm{\bb_g} \) plus a constant. With
\( \X_g\T\X_g=\I \), this is the one-group problem of (d) with \( \mathbf{z}_g \) replaced by \( \mathbf{s}_g \).
:::

Block soft thresholding keeps the direction of the least squares vector \( \mathbf{z}_g \) and reduces its length by
\( \lambda w_g \): the lasso rule applied to a length. The weights come from (c): for pure noise and orthonormal groups,
\( \norm{\X_g\T\y}^2/\sigma^2\sim\chi^2(p_g) \) has mean \( p_g \), so without \( \sqrt{p_g} \) large groups would enter merely for being
large (@exr-reg-group-weights). Block coordinate descent with update (e) has every limit point a minimizer because the penalty is
separable across blocks (Tseng 2001), so it converges when the minimizer is unique; without orthonormality the block problem has no closed form (@exr-reg-group-nonorthonormal).

## Factors as groups

Orthonormalizing within groups also makes the group lasso independent of how a factor is coded.

::: {#prp-reg-group-invariance}
[Invariance to the coding of a group]

::: {.enumerate options="label=(\alph*)"}
1. Let \( \mathbf{U}_g \) be \( p_g\times p_g \) orthogonal matrices, and replace every \( \X_g \) by \( \X_g\mathbf{U}_g \). Then \( \hbeta \) minimizes
   the original criterion iff the vector with blocks \( \mathbf{U}_g\T\hbeta_g \) minimizes the new one. The fitted contributions
   \( \X_g\hbeta_g \) are the same.

2. Suppose each group is orthonormalized within itself, so that \( \X_g \) is an orthonormal basis of a subspace
   \( \mathcal V_g \). Then the minimizers, expressed through the fitted contributions \( \X_g\hbeta_g\in\mathcal V_g \), depend on the groups
   only through the subspaces \( \mathcal V_g \). In particular, for a factor they do not depend on the choice of reference level
   or of contrasts.
:::

:::

::: {.proof}
(a) With \( \bb_g=\mathbf{U}_g\mathbf{c}_g \), we have \( \X_g\mathbf{U}_g\mathbf{c}_g=\X_g\bb_g \) and \( \norm{\mathbf{c}_g}=\norm{\bb_g} \). So the new criterion at
\( \mathbf{c} \) equals the old criterion at \( \bb \), and \( \bb\leftrightarrow\mathbf{c} \) is a bijection.

(b) Two orthonormal bases \( \X_g \) and \( \X_g' \) of \( \mathcal V_g \) satisfy \( \X_g'=\X_g\mathbf{U}_g \) with
\( \mathbf{U}_g=\X_g\T\X_g' \) orthogonal. Apply (a). For a factor with \( k \) levels, the centred indicator columns of any reference
coding, and the centred columns of any other full-rank contrast coding, span the same \( (k-1) \)-dimensional subspace: the
vectors that are constant within levels and orthogonal to \( \bone \).
:::

With orthonormalized groups \( \norm{\bb_g}=\norm{\X_g\bb_g} \): the penalty is the total length of the groups' fitted
contributions, a coding-free quantity. Simon and Tibshirani (2012) argue for this version.

::: {#exm-reg-factors}
[Selecting factors]

The script simulates \( n=120 \) observations of four factors \( A,B,C,D \) (with \( 4 \), \( 3 \), \( 5 \) and \( 3 \) levels) and three
numeric regressors; the response depends on \( A \), weakly on \( B \), and on \( x_1 \). The group lasso path starts at
\( \lambda_{\max}=6.876 \), and groups enter in the order \( x_1 \), \( A \), \( B \), \( D \). At \( \lambda=2.0 \) only \( A \) and
\( x_1 \) are selected, with \( \norm{\X_A\hbeta_A}=2.834 \) and \( \norm{\X_{x_1}\hbeta_{x_1}}=4.861 \). Another reference level
gives the same fitted contributions to machine precision. The lasso on the raw indicator columns does not: at
\( \lambda=6.0 \) its fit moves by \( 8.6 \) per cent of its length when the reference levels change. Along its path the columns
of \( A \) enter at values of \( \lambda \) spread over a factor of \( 2.8 \), and it admits its first column of a null factor
(\( C \) or \( D \)) at a residual sum of squares of \( 112.0 \), a worse fit than the \( 110.7 \) that the group lasso reaches with both
null factors still excluded ([Figure 30.3.1](#fig-reg-group-paths)).
:::

::: {when-format="html"}
![**Figure 30.3.1.** (a) Group lasso path for @exm-reg-factors: the length of each group's fitted contribution
\( \norm{\X_g\hbeta_g} \) against \( \log_{10}\lambda \). Factors enter as wholes. (b) The lasso path for the same data, with each
indicator column penalized separately, coloured by factor. The dashed curves are the numeric
regressors.](group_paths.svg){#fig-reg-group-paths width=100%}
:::

::: {when-format="pdf"}
![(a) Group lasso path for @exm-reg-factors: the length of each group's fitted contribution
\( \norm{\X_g\hbeta_g} \) against \( \log_{10}\lambda \). Factors enter as wholes. (b) The lasso path for the same data, with each
indicator column penalized separately, coloured by factor. The dashed curves are the numeric
regressors.](group_paths.pdf){width=100%}
:::

```{.python .run #cell-group-lasso-fit}
import numpy as np

rng = np.random.default_rng(3003)
n = 120
levels = {"A": 4, "B": 3, "C": 5, "D": 3}
factors = {f: rng.integers(0, k, size=n) for f, k in levels.items()}
numeric = rng.normal(size=(n, 3))
effect_A = np.array([0.0, 1.0, -0.5, 0.8])
effect_B = np.array([0.0, 0.3, -0.3])
y = effect_A[factors["A"]] + effect_B[factors["B"]] + 0.6 * numeric[:, 0] + rng.normal(size=n)
y -= y.mean()

def indicator_columns(codes, k, reference=0):
    """Indicator columns for the levels other than the reference level."""
    return np.column_stack([(codes == l).astype(float) for l in range(k) if l != reference])

def orthonormal_group(Z):
    """Centre the columns, then an orthonormal basis of their span (Z = Q R)."""
    Z = Z - Z.mean(axis=0)
    Q, _ = np.linalg.qr(Z)
    return Q

def build(reference=0):
    blocks = [orthonormal_group(indicator_columns(factors[f], k, reference))
              for f, k in levels.items()]
    blocks += [orthonormal_group(numeric[:, [j]]) for j in range(3)]
    return blocks

def group_lasso(blocks, y, lam, b=None, tol=1e-12, max_sweeps=100000):
    """(1/2)||y - sum_g X_g b_g||^2 + lam sum_g sqrt(p_g) ||b_g||, with X_g'X_g = I."""
    b = [np.zeros(Xg.shape[1]) for Xg in blocks] if b is None else [v.copy() for v in b]
    r = y - sum(Xg @ bg for Xg, bg in zip(blocks, b))
    for _ in range(max_sweeps):
        largest = 0.0
        for g, Xg in enumerate(blocks):
            s = Xg.T @ r + b[g]                        # block least squares on the partial residual
            w = np.sqrt(Xg.shape[1])
            norm_s = np.linalg.norm(s)
            new = max(0.0, 1 - lam * w / norm_s) * s if norm_s > 0 else 0 * s
            r += Xg @ (b[g] - new)
            largest = max(largest, np.max(np.abs(new - b[g])))
            b[g] = new
        if largest < tol:
            return b
    raise RuntimeError("block coordinate descent did not converge")

blocks = build()
names = list(levels) + ["x1", "x2", "x3"]
lam = 2.0
b = group_lasso(blocks, y, lam)
for name, Xg, bg in zip(names, blocks, b):
    print(f"{name}: {Xg.shape[1]} columns, ||X_g b_g|| = {np.linalg.norm(Xg @ bg):.3f}")
```

The **sparse group lasso** adds an \( \ell_1 \) term so that a selected group can have zero coefficients (Simon, Friedman,
Hastie and Tibshirani, 2013).

## Ordered coefficients: the fused lasso

When coefficients have a natural order (levels of an ordered factor, successive times), neighbours are expected to be
similar. The **fused lasso** of Tibshirani, Saunders, Rosset, Zhu and Knight (2005) penalizes their differences:
\[
Q(\bb)=\tfrac12\norm{\y-\X\bb}^2+\lambda_1\sum_{j=1}^p\lvert b_j\rvert+\lambda_2\sum_{j=2}^p\lvert b_j-b_{j-1}\rvert .
\]{#eq-reg-fused}

It produces exact zeros and exact ties \( \hat\beta_j=\hat\beta_{j-1} \), so the coefficient sequence is piecewise constant.
The penalty is convex but not separable, so coordinate descent can stall at a tie. For an ordered factor, fusion merges adjacent levels the data cannot distinguish (Gertheiss and Tutz 2010).
@exr-reg-fused-two solves the smallest case.

## Exercises

### A. Check your understanding

::: {#exr-reg-group-weights}
[A1]

Under orthonormal groups and \( \y\sim\Normal_n(\bzero,\sigma^2\I) \), show that \( P(\hbeta_g\neq\bzero)=P\{\chi^2(p_g)>\lambda^2w_g^2/\sigma^2\} \).
With \( w_g=1 \), compare these probabilities for \( p_g=1 \) and \( p_g=10 \) when \( \lambda=2\sigma \). What do the weights
\( \sqrt{p_g} \) achieve, and what do they not achieve?
:::

::: {.solution}
By @eq-reg-block-soft, \( \hbeta_g\neq\bzero \) iff \( \norm{\mathbf{z}_g}>\lambda w_g \). Here \( \mathbf{z}_g=\X_g\T\y\sim\Normal_{p_g}(\bzero,\sigma^2\I) \), so
\( \norm{\mathbf{z}_g}^2/\sigma^2\sim\chi^2(p_g) \). With \( w_g=1 \) and \( \lambda=2\sigma \), the probabilities are \( P\{\chi^2(1)>4\}=0.046 \) and
\( P\{\chi^2(10)>4\}=0.947 \), so a large null group almost always enters. With \( w_g=\sqrt{p_g} \) the threshold is
\( 4p_g \), four times the mean of \( \chi^2(p_g) \). This removes the gross imbalance, but it does not equalize the probabilities:
they are \( 0.046 \) and \( P\{\chi^2(10)>40\}=1.7\times 10^{-5} \). Equal null entry probabilities would need thresholds at a common quantile of
\( \chi^2(p_g) \).
:::

### B. Practice

::: {#exr-reg-fused-two}
[B1]

Solve the fused lasso @eq-reg-fused with \( \X=\I_2 \), \( \lambda_1=0 \) and \( \lambda_2=\lambda \). Show that
\( \hat\beta_1=\hat\beta_2=\tfrac12(y_1+y_2) \) if \( \lvert y_1-y_2\rvert\le2\lambda \). Otherwise, each estimate moves a distance
\( \lambda \) from its observation towards the other.
:::

::: {.solution}
The penalty is \( \lambda\lvert b_1-b_2\rvert=\lambda\lvert\mathbf{a}\T\bb\rvert \) with \( \mathbf{a}=(1,-1)\T \). Its subgradients are \( \lambda s\mathbf{a} \) with
\( s\in\partial\lvert\cdot\rvert(b_1-b_2) \). By @lem-reg-optimality, \( y_1-b_1=\lambda s \) and \( y_2-b_2=-\lambda s \). If \( b_1=b_2=b \), adding
the equations gives \( b=\tfrac12(y_1+y_2) \), and then \( s=(y_1-y_2)/(2\lambda) \) must lie in \( [-1,1] \). If
\( \lvert y_1-y_2\rvert>2\lambda \), then \( b_1\neq b_2 \), \( s=\operatorname{sign}(b_1-b_2) \), and \( b_1=y_1-\lambda s \), \( b_2=y_2+\lambda s \). This is
consistent when \( s=\operatorname{sign}(y_1-y_2) \), since then \( b_1-b_2=y_1-y_2-2\lambda s \) has the sign of \( s \).
:::

::: {#exr-reg-group-vs-lasso}
[B2]

Under an orthonormal design, compare the group lasso (one group of size \( 3 \), \( w=\sqrt3 \)) with the lasso on the same three
columns, for \( \mathbf{z}=(2,0.3,-0.2)\T \) and \( \lambda=1 \). Which coefficients are zero in each? Show that, under an orthonormal design,
the group lasso never sets a single coordinate of a selected group to zero unless the corresponding \( z_j=0 \).
:::

### C. Going deeper

::: {#exr-reg-group-nonorthonormal}
[C1]

Drop the assumption \( \X_g\T\X_g=\I \). Show that a nonzero solution of the block problem
\( \min_{\bb}\tfrac12\norm{\br-\X_g\bb}^2+c\norm{\bb} \) is \( \bb=(\X_g\T\X_g+(c/t)\I)^{-1}\X_g\T\br \), where \( t=\norm{\bb} \) is the root
of a one-dimensional equation. Show that this equation has a unique root when \( \norm{\X_g\T\br}>c \), so the block update
is a ridge regression (@thm-shr-ridge) with a data-dependent parameter.
:::

