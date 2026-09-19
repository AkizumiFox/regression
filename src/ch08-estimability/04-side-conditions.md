# Reparameterization and side conditions

A rank-deficient model can be handled in three ways. One can work with it as it stands,
estimating only estimable functions, which is what the last two sections did. One can
*reparameterize*, replacing \( \X \) by a full-rank matrix with the same column space. Or one can
impose *side conditions*, linear restrictions on \( \bbeta \) that select a single coefficient vector
from each set \( \bbeta+\Null(\X) \). This section shows that the second and third are two
descriptions of the same act: choosing coordinates for the mean vector. Neither changes any
estimable quantity. It also records what the most common software does by default.

## Full-rank reparameterizations

By @cor-est-estimable-space, at most \( r \) estimable functions are linearly independent. A set of
exactly \( r \) of them can serve as a new parameter vector.

::: {#thm-est-reparameterization}
[Reparameterizing to full rank]

Let \( \rank(\X)=r \). Let \( \bU \) be an \( r\times p \) matrix of rank \( r \) whose rows lie in
\( \C(\X\T) \), so that \( \bgamma=\bU\bbeta \) is a vector of \( r \) linearly independent estimable functions,
and put
\[
\Z=\X\bU\T(\bU\bU\T)^{-1}.
\]{#eq-est-reparameterization-matrix}

Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \X=\Z\bU \), \( \Z \) has full column rank \( r \), and \( \C(\Z)=\C(\X) \);

2. the models \( \E(\Y)=\X\bbeta \), \( \bbeta\in\Real^p \), and \( \E(\Y)=\Z\bgamma \), \( \bgamma\in\Real^r \), describe
           the same set of mean vectors, with \( \bgamma=\bU\bbeta \);

3. the least squares estimate in the new model is \( \hat{\bgamma}=(\Z\T\Z)^{-1}\Z\T\y=\bU\hbeta \) for every least
           squares solution \( \hbeta \) of the old one;

4. \( \blambda\T\bbeta \) is estimable iff \( \blambda=\bU\T\mathbf{c} \) for some \( \mathbf{c}\in\Real^r \), which is then
           unique, and in that case \( \blambda\T\hbeta=\mathbf{c}\T\hat{\bgamma} \);

5. conversely, if \( \Z \) is any \( n\times r \) matrix of full column rank with \( \C(\Z)=\C(\X) \), then
           \( \X=\Z\bU \) for a unique \( \bU \), and \( \bU \) has rank \( r \) and estimable rows.
:::

:::

::: {.proof}
Since \( \C(\bU\T)\subseteq\C(\X\T) \) and both have dimension \( r \), \( \C(\bU\T)=\C(\X\T) \). The matrix
\( \mathbf{Q}=\bU\T(\bU\bU\T)^{-1}\bU \) is symmetric and idempotent with column space \( \C(\bU\T) \), so it is the
orthogonal projection onto \( \C(\X\T) \) (@thm-proj-sym-idem). Hence \( \I-\mathbf{Q} \) projects onto
\( \C(\X\T)\perpc=\Null(\X) \), and \( \X(\I-\mathbf{Q})=\bzero \), that is, \( \X=\X\mathbf{Q}=\Z\bU \).
(a) Now \( r=\rank(\X)=\rank(\Z\bU)\le\rank(\Z)\le r \), so \( \Z \) has full column rank.
\( \C(\Z)\subseteq\C(\X) \) because \( \Z=\X(\cdot) \), and \( \C(\X)\subseteq\C(\Z) \) because \( \X=\Z\bU \).
(b) \( \X\bbeta=\Z(\bU\bbeta) \), and \( \bbeta\mapsto\bU\bbeta \) maps \( \Real^p \) onto \( \Real^r \) because \( \bU \)
has full row rank. So both models describe \( \C(\X) \).
(c) \( \Z\hat{\bgamma}=\M\y=\X\hbeta=\Z(\bU\hbeta) \) by @thm-proj-ls-projection, and \( \Z \) has full column rank.
(d) \( \blambda \) is estimable iff \( \blambda\in\C(\X\T)=\C(\bU\T) \). The vector \( \mathbf{c} \) is unique because
\( \bU\T \) has full column rank, and \( \blambda\T\hbeta=\mathbf{c}\T\bU\hbeta=\mathbf{c}\T\hat{\bgamma} \).
(e) Each column of \( \X \) lies in \( \C(\Z) \), so \( \X=\Z\bU \) with \( \bU=(\Z\T\Z)^{-1}\Z\T\X \), unique because
\( \Z \) has full column rank. Then \( \bU\bbeta=(\Z\T\Z)^{-1}\Z\T(\X\bbeta) \) is a linear function of \( \X\bbeta \), so each
row of \( \bU \) is estimable. Finally \( r=\rank(\X)\le\rank(\bU)\le r \).
:::

So reparameterization is the same thing as choosing a basis for the estimable space.
@thm-proj-reparam of [Chapter 6](../ch06-projections/index.html) already showed that fitted values
and residuals do not care which basis is chosen. Part (d) adds that estimable functions do not
care either: each is a fixed linear combination of the new parameters, and its estimate is the same
combination of their estimates. Everything in the full-rank theory of
[Chapter 5](../ch05-model-and-least-squares/index.html) now applies to \( \hat{\bgamma} \), including
@thm-lm-moments and the variance estimator of @thm-lm-sigma2, with \( p \) replaced by \( r \).

For the one-way layout with \( \bbeta=(\mu,\alpha_1,\dots,\alpha_g)\T \), two choices of \( \bU \) stand out.
The rows \( (1,\mathbf{e}_k\T) \), \( k=1,\dots,g \), give \( \bgamma=(\mu+\alpha_1,\dots,\mu+\alpha_g)\T \), the
cell means, and \( \Z \) is the indicator matrix. The rows \( (1,\mathbf{e}_1\T) \) and \( (0,\mathbf{e}_k\T-\mathbf{e}_1\T) \),
\( k=2,\dots,g \), give the first group's mean and the differences from it. This is the reference coding of
[Section 8.5](05-factor-coding.html).

## Side conditions

The other classical device leaves \( \X \) alone and adds equations. A **side condition** is a set of
linear restrictions \( \bT\bb=\bzero \), with \( \bT \) a \( q\times p \) matrix, imposed on the solutions of the
normal equations. The aim is to single out one solution. Examples in the one-way layout are
\( \alpha_1=0 \), \( \sum_k\alpha_k=0 \) and \( \mu=0 \). The theorem says exactly which restrictions
achieve this.

::: {#thm-est-side-conditions}
[Side conditions]

Let \( \X \) be \( n\times p \) of rank \( r \) and \( \bT \) be \( q\times p \). The following are equivalent:

::: {.enumerate options="label=(\alph*)"}
1. for every \( \y\in\Real^n \) there is exactly one \( \bb \) with \( \X\T\X\bb=\X\T\y \) and \( \bT\bb=\bzero \);

2. for every \( \bmu\in\C(\X) \) there is exactly one \( \bb \) with \( \X\bb=\bmu \) and \( \bT\bb=\bzero \);

3. \( \Real^p=\Null(\X)\dirsum\Null(\bT) \);

4. \( \rank\begin{pmatrix}\X\\\bT\end{pmatrix}=p \) and \( \C(\bT\T)\cap\C(\X\T)=\{\bzero\} \);

5. \( \rank(\bT)=p-r \), and no nonzero linear combination of the rows of \( \bT \) is estimable.
:::

When these hold, \( \X\T\X+\bT\T\bT \) is nonsingular, the solution in (a) is
\[
\hbeta_{\bT}=(\X\T\X+\bT\T\bT)^{-1}\X\T\y ,
\]{#eq-est-side-solution}

and \( (\X\T\X+\bT\T\bT)^{-1} \) is a generalized inverse of \( \X\T\X \).
:::

::: {.proof}
(a)\( \Leftrightarrow \)(b). The solutions of the normal equations for \( \y \) are the \( \bb \) with
\( \X\bb=\M\y \) (@thm-proj-normal-equations), and \( \M\y \) runs over all of \( \C(\X) \) as \( \y \) does.

(b)\( \Leftrightarrow \)(c). The set \( \{\bb:\X\bb=\bmu\} \) is a coset \( \bb_0+\Null(\X) \), and every \( \bb_0\in\Real^p \)
arises (take \( \bmu=\X\bb_0 \)). So (b) says that every coset of \( \Null(\X) \) meets \( \Null(\bT) \) in exactly
one point. Some point \( \bb_0+\bv \) of the coset lies in \( \Null(\bT) \) iff
\( \bb_0\in\Null(\bT)+\Null(\X) \), so existence for every coset means \( \Real^p=\Null(\X)+\Null(\bT) \). Two
points of the same coset that both lie in \( \Null(\bT) \) differ by an element of
\( \Null(\X)\cap\Null(\bT) \), and a nonzero element of that intersection added to one such point gives a
second. So uniqueness means \( \Null(\X)\cap\Null(\bT)=\{\bzero\} \). Together, these say the sum is direct.

(c)\( \Leftrightarrow \)(d). By @lem-proj-null-colspace and @lem-proj-complement-intersection,
\[
\bigl(\Null(\X)+\Null(\bT)\bigr)\perpc=\C(\X\T)\cap\C(\bT\T),\qquad
\bigl(\Null(\X)\cap\Null(\bT)\bigr)\perpc=\C(\X\T)+\C(\bT\T).
\]
So \( \Null(\X)+\Null(\bT)=\Real^p \) iff \( \C(\X\T)\cap\C(\bT\T)=\{\bzero\} \), and
\( \Null(\X)\cap\Null(\bT)=\{\bzero\} \) iff \( \C(\X\T)+\C(\bT\T)=\Real^p \). The last space is the row space of
the stacked matrix, so this holds iff that matrix has rank \( p \).

(d)\( \Leftrightarrow \)(e). A vector in \( \C(\bT\T) \) is a combination of the rows of \( \bT \), and it is estimable
iff it lies in \( \C(\X\T) \). So the second halves of (d) and (e) say the same thing. When they hold,
the sum \( \C(\X\T)+\C(\bT\T) \) is direct, with dimension \( r+\rank(\bT) \). This equals \( p \) iff
\( \rank(\bT)=p-r \).

For the formula, suppose (d). If \( (\X\T\X+\bT\T\bT)\bv=\bzero \), then
\( 0=\bv\T\X\T\X\bv+\bv\T\bT\T\bT\bv=\norm{\X\bv}^2+\norm{\bT\bv}^2 \), so \( \bv \) lies in
\( \Null(\X)\cap\Null(\bT)=\{\bzero\} \). Hence the matrix is nonsingular. The solution \( \bb \) of (a) satisfies
\( (\X\T\X+\bT\T\bT)\bb=\X\T\y+\bT\T\bzero \), which gives @eq-est-side-solution. Finally, write
\( \G_{\bT}=(\X\T\X+\bT\T\bT)^{-1} \). For every \( \y \), \( \G_{\bT}\X\T\y \) solves the normal equations, so
\( \X\T\X\G_{\bT}\X\T=\X\T \). Multiplying on the right by \( \X \) gives \( \X\T\X\G_{\bT}\X\T\X=\X\T\X \).
:::

Condition (e) is the familiar rule of thumb, made exact: a side condition must supply exactly as many
restrictions as the rank deficiency, and none of them, singly or in combination, may be estimable. The
reason for the second requirement is clear from (d). An estimable restriction constrains the mean
vector itself, and so it changes the model instead of merely choosing coordinates.

::: {.remark}
[Existence and uniqueness separately]

The proof shows more than the theorem states. The condition \( \C(\bT\T)\cap\C(\X\T)=\{\bzero\} \) alone
guarantees that the restricted normal equations always have *at least* one solution. The condition
\( \rank\begin{psmallmatrix}\X\\\bT\end{psmallmatrix}=p \) alone guarantees *at most* one. A constraint that is
too weak, such as a single restriction when \( p-r=2 \), leaves many solutions. One with an estimable
part leaves none for most \( \y \) (@exr-est-bad-side).
:::

## What the constrained coefficients mean

Under a side condition every coefficient becomes identifiable, and it acquires a meaning.

::: {#prp-est-side-meaning}
[Constrained coefficients are estimable functions]

Suppose \( \bT \) satisfies the conditions of @thm-est-side-conditions, and let
\( \bP_{\bT}=(\X\T\X+\bT\T\bT)^{-1}\X\T\X \). For each \( \bbeta \), \( \bP_{\bT}\bbeta \) is the unique vector
\( \bb \) with \( \X\bb=\X\bbeta \) and \( \bT\bb=\bzero \). The matrix \( \bP_{\bT} \) is idempotent, with range
\( \Null(\bT) \) and null space \( \Null(\X) \). Every row of \( \bP_{\bT} \) lies in \( \C(\X\T) \). So each
constrained coefficient \( \mathbf{e}_j\T\bP_{\bT}\bbeta \) is an estimable function of the original parameters,
and \( \hbeta_{\bT} \) estimates it without bias.
:::

::: {.proof}
With \( \y=\X\bbeta \), the vector \( \bP_{\bT}\bbeta=\G_{\bT}\X\T(\X\bbeta) \) is the unique solution of
@thm-est-side-conditions(a) for that \( \y \). It satisfies \( \X\bb=\M\X\bbeta=\X\bbeta \) and
\( \bT\bb=\bzero \). The matrix \( \bP_{\bT} \) is \( \bH_{\G} \) of @eq-est-H for the generalized inverse
\( \G_{\bT} \), so @prp-est-H gives idempotence, the null space \( \Null(\X) \) and rows in \( \C(\X\T) \).
Its range lies in \( \Null(\bT) \). The range has dimension \( r \), and so does \( \Null(\bT) \), because
\( \rank(\bT)=p-r \). So the two coincide. Unbiasedness is @prp-est-H(c).
:::

In words, \( \bP_{\bT} \) is the projection of \( \Real^p \) onto \( \Null(\bT) \) along \( \Null(\X) \), generally an
oblique one. It sends the true \( \bbeta \) to the unique representative of \( \bbeta+\Null(\X) \) that
satisfies the side condition. That representative's coordinates are functions of the mean, so they
can be estimated. Their *meaning* is fixed by \( \bT \). In the one-way layout, the condition
\( \sum_k\alpha_k=0 \) makes \( \mu \) the unweighted average of the \( g \) cell means and \( \alpha_k \) the
deviation of cell \( k \) from it. The condition \( \sum_kn_k\alpha_k=0 \) makes \( \mu \) the size-weighted
average instead. The condition \( \alpha_1=0 \) makes \( \mu \) the first cell mean. The data do not
choose among these. The analyst chooses, by deciding which question the coefficients should answer.

One side condition deserves special mention, because it is the one many programs impose without saying
so.

::: {#prp-est-min-norm-side}
[Minimum norm as a side condition]

Let the columns of the \( p\times(p-r) \) matrix \( \V_0 \) be a basis of \( \Null(\X) \). Then \( \bT=\V_0\T \)
satisfies the conditions of @thm-est-side-conditions, and the resulting solution \( \hbeta_{\bT} \) is the
minimum-norm least squares solution \( \X^+\y \).
:::

::: {.proof}
\( \C(\bT\T)=\C(\V_0)=\Null(\X) \), which has dimension \( p-r \) and meets \( \C(\X\T)=\Null(\X)\perpc \) only in
\( \bzero \). So (d) holds. The least squares solutions are \( \hbeta_{\bT}+\Null(\X) \), and
\( \hbeta_{\bT}\in\Null(\bT)=\C(\V_0)\perpc=\Null(\X)\perpc \). For \( \bv\in\Null(\X) \), Pythagoras gives
\( \norm{\hbeta_{\bT}+\bv}^2=\norm{\hbeta_{\bT}}^2+\norm{\bv}^2 \), so \( \hbeta_{\bT} \) has the smallest norm. By
@prp-proj-min-norm this is \( \X^+\y \).
:::

In the one-way layout the null space is spanned by \( (1,-1,\dots,-1)\T \), so the minimum-norm solution
is the one satisfying \( \mu=\alpha_1+\dots+\alpha_g \). No one would choose this restriction on purpose,
and its coefficients have no natural interpretation.

::: {#exm-est-five-conditions}
[Five side conditions, one fit]

Eighteen simulated observations fall into four groups of sizes \( 5,3,6,4 \), with group means
\( 11.111 \), \( 15.204 \), \( 10.558 \) and \( 13.624 \). In the
model \( \E(y_{kj})=\mu+\alpha_k \), \( p-r=1 \), so a side condition is a single nonestimable restriction.
The listing solves @eq-est-side-solution for five of them and computes the minimum-norm solution.

| Side condition | \( \hat{\mu} \) | \( \hat{\alpha}_1 \) | \( \hat{\alpha}_2 \) | \( \hat{\alpha}_3 \) | \( \hat{\alpha}_4 \) |
|---|---|---|---|---|---|
| \( \alpha_1=0 \) | 11.111 | 0.000 | 4.093 | -0.553 | 2.513 |
| \( \alpha_4=0 \) | 13.624 | -2.513 | 1.580 | -3.066 | 0.000 |
| \( \sum_k\alpha_k=0 \) | 12.624 | -1.513 | 2.580 | -2.066 | 1.000 |
| \( \sum_kn_k\alpha_k=0 \) | 12.167 | -1.056 | 3.037 | -1.609 | 1.457 |
| \( \mu=0 \) | 0.000 | 11.111 | 15.204 | 10.558 | 13.624 |
| minimum norm (\( \mu=\sum_k\alpha_k \)) | 10.100 | 1.012 | 5.105 | 0.459 | 3.525 |

Every row solves the normal equations. In every row \( \hat{\mu}+\hat{\alpha}_k \) is the \( k \)th group mean, and
\( \hat{\alpha}_1-\hat{\alpha}_2 \) is the same. The value of \( \hat{\mu} \) is, in turn, the first group mean, the last
group mean, the average of the four group means (\( 12.624 \)), the overall mean
\( \bar{y}=12.167 \), zero, and an uninterpretable \( 10.100 \).

Compare a restriction that *is* estimable, \( \alpha_1=\alpha_2 \). It is not a side condition. Imposing it
(by solving the Lagrange system for least squares under a linear restriction) forces the first two group means
to be equal. The residual sum of squares rises from \( 14.665 \) to
\( 46.078 \). The fit has changed, which is what testing the hypothesis
\( \alpha_1=\alpha_2 \) in [Chapter 11](../ch11-general-linear-hypothesis/index.html) measures.
:::

```{.python .run #cell-side-conditions-layout}
import numpy as np

sizes = [5, 3, 6, 4]
groups = np.repeat(np.arange(4), sizes)
rng = np.random.default_rng(84)
y = np.array([12.0, 15.0, 11.0, 14.0])[groups] + rng.normal(scale=1.0, size=groups.size)
X = np.column_stack([np.ones(groups.size), np.eye(4)[groups]])   # (mu, alpha_1..alpha_4)
A, c = X.T @ X, X.T @ y
```

```{.python .run #cell-side-conditions-sideconditions}
n_k = np.array(sizes, dtype=float)
side = {"alpha_1 = 0":         [0, 1, 0, 0, 0],
        "alpha_4 = 0":         [0, 0, 0, 0, 1],
        "sum alpha_k = 0":     [0, 1, 1, 1, 1],
        "sum n_k alpha_k = 0": [0, *n_k],
        "mu = 0":              [1, 0, 0, 0, 0]}
solutions = {}
for name, t in side.items():
    T = np.array([t], dtype=float)
    G = np.linalg.inv(A + T.T @ T)          # (X'X + T'T)^{-1}: a g-inverse of X'X
    b = G @ c
    solutions[name] = b
    print(f"{name:20s} b = {np.round(b, 3)}   T b = {(T @ b)[0]:.1e}")
b_mn = np.linalg.pinv(X) @ y
print(f"{'minimum norm':20s} b = {np.round(b_mn, 3)}",
      f"  mu - sum(alpha) = {b_mn[0] - b_mn[1:].sum():.1e}")
```

```{.python .run #cell-side-conditions-estimablerestriction}
# "alpha_1 = alpha_2" is a restriction on an estimable function: it changes the model
T_bad = np.array([[0, 1, -1, 0, 0]], dtype=float)
K = np.block([[A, T_bad.T], [T_bad, np.zeros((1, 1))]])        # Lagrange (bordered) system
b_bad = np.linalg.lstsq(K, np.append(c, 0.0), rcond=None)[0][:5]
sse_free = np.sum((y - X @ solutions["alpha_1 = 0"]) ** 2)
sse_bad = np.sum((y - X @ b_bad) ** 2)
print(f"SSE without restriction {sse_free:.3f};  with alpha_1 = alpha_2 {sse_bad:.3f}")
```

Reparameterization and side conditions are two languages for the same choice. Deleting the
column of \( \alpha_1 \) from \( \X \) and fitting the full-rank remainder gives exactly the solution with
\( \alpha_1=0 \). More generally, deleting a set of columns whose removal leaves a full-rank matrix with the same
column space is the side condition that sets those coefficients to zero, and it satisfies
@thm-est-side-conditions(c) (@exr-est-drop-columns).

## What software does

Programs treat a rank-deficient \( \X \) differently.

- *Code the factor at full rank before fitting.* Formula interfaces such as R's `lm()` and statsmodels'
  formula API expand a factor with \( g \) levels into an intercept and \( g-1 \) columns, by default the reference
  coding of [Section 8.5](05-factor-coding.html). The overparameterized matrix never appears. This is a
  reparameterization in the sense of @thm-est-reparameterization.
- *Drop columns found to be dependent.* If the columns supplied are dependent anyway, R's `lm()`
  computes a QR decomposition in which a column found to be (numerically) a combination of earlier ones is
  moved to the end and left out. Its coefficient is reported as `NA`. This is the side condition setting that
  coefficient to zero, and which column is dropped depends on the order of the columns. SAS's GLM procedure
  similarly solves the normal equations with a generalized inverse that zeroes the coefficients of dependent
  columns, and it flags the resulting estimates as not unique.
- *Return the minimum-norm solution.* `numpy.linalg.lstsq` returns \( \X^+\y \) silently (it reports the
  numerical rank only as a second output), and statsmodels' `OLS`, whose `fit` uses the pseudoinverse by
  default, returns the same vector with a `SingularMatrixWarning`. This is the side condition of
  @prp-est-min-norm-side, and it does not depend on the order of the columns.
- *Pivot greedily.* A QR decomposition with column pivoting ([Chapter 10](../ch10-computation/index.html), @thm-cmp-rank-revealing)
  chooses at each step the remaining column with the largest component orthogonal to those already chosen,
  and the columns left at the end are dropped. In @exm-est-five-conditions it takes the intercept, then
  \( \alpha_3 \) (the largest group), \( \alpha_1 \) and \( \alpha_4 \), and drops \( \alpha_2 \), the smallest group, a
  choice that depends on the group sizes.

```{.python .run #cell-side-conditions-software}
import scipy.linalg
import statsmodels.api as sm
b_sm = sm.OLS(y, X).fit().params                 # statsmodels default: method="pinv"
b_np = np.linalg.lstsq(X, y, rcond=None)[0]      # LAPACK gelsd: minimum-norm solution
_, R, piv = scipy.linalg.qr(X, mode="economic", pivoting=True)
print("statsmodels:", np.round(b_sm, 3))
print("numpy lstsq:", np.round(b_np, 3))
print("pivoted QR column order:", piv, " |R_kk| =", np.round(np.abs(np.diag(R)), 3))
```

The lesson of this section is that all these outputs are correct, and none is more correct than the others.
They are different coordinate systems for one fitted mean vector. Any conclusion that changes when the
software changes is a conclusion about a nonestimable function.

::: {.warning}
Before interpreting a coefficient from a rank-deficient fit, find out which side condition produced it. A
coefficient reported as “the effect of level 3” may be a difference from level 1, a deviation from an average
of cell means, a deviation from the overall mean, or, under the minimum-norm convention, a quantity with no
interpretation at all.
:::

## Exercises

### A. Check your understanding

::: {#exr-est-valid-sides}
[A1]

In the one-way model with four groups and \( \bbeta=(\mu,\alpha_1,\dots,\alpha_4)\T \), decide which of the
following are side conditions in the sense of @thm-est-side-conditions:
(i) \( \alpha_2=0 \); (ii) \( \mu=\alpha_3 \); (iii) \( \alpha_1-\alpha_2=0 \); (iv) \( \alpha_1=0 \) together with
\( \alpha_2=0 \); (v) \( 2\mu+\alpha_1+\alpha_2=0 \).
:::

### B. Practice

::: {#exr-est-drop-columns}
[B1]

Let \( \X=[\X_1,\X_2] \), where \( \X_1 \) has full column rank and \( \C(\X_1)=\C(\X) \). Show that
the least squares fit using \( \X_1 \) alone, padded with zeros for the coefficients of \( \X_2 \), is the
solution of the side condition \( \bT=[\bzero,\I] \). Check the conditions of @thm-est-side-conditions
directly.
:::

::: {.solution}
The padded vector \( (\hbeta_1\T,\bzero\T)\T \) satisfies
\( \X(\hbeta_1\T,\bzero\T)\T=\X_1\hbeta_1=\M\y \), because \( \C(\X_1)=\C(\X) \). So it solves the normal equations
and satisfies \( \bT\bb=\bzero \). For the conditions: \( \rank(\bT)=p_2 \), and \( \rank(\X)=\rank(\X_1)=p_1 \), so
\( \rank(\bT)=p-r \). If \( \bT\T\mathbf{c}=(\bzero\T,\mathbf{c}\T)\T \) were in \( \C(\X\T) \), say
\( \X_1\T\boldsymbol{\rho}=\bzero \) and \( \X_2\T\boldsymbol{\rho}=\mathbf{c} \), then \( \boldsymbol{\rho}\perp\C(\X_1)=\C(\X)\supseteq\C(\X_2) \), so
\( \mathbf{c}=\bzero \). Condition (e) holds.
:::

::: {#exr-est-two-way-reparameterization}
[B2]

For the additive model \( \E(y_{ij})=\mu+\alpha_i+\beta_j \) with \( i=1,2 \), \( j=1,2,3 \) and one observation
per cell, take \( \bU \) with rows giving \( \mu+\alpha_1+\beta_1 \), \( \alpha_2-\alpha_1 \), \( \beta_2-\beta_1 \) and
\( \beta_3-\beta_1 \). Compute \( \Z \) from @eq-est-reparameterization-matrix, and check that it is the reference-coded model
matrix.
:::

::: {#exr-est-lagrange}
[B3]

Let \( \bT \) satisfy the conditions of @thm-est-side-conditions and have full row rank \( p-r \). Show that
the bordered matrix
\( \begin{psmallmatrix}\X\T\X&\bT\T\\\bT&\bzero\end{psmallmatrix} \) is nonsingular, that the solution of the
bordered system with right-hand side \( (\X\T\y,\bzero) \) has first block \( \hbeta_{\bT} \) and second block
\( \bzero \), and interpret the zero Lagrange multiplier. Relate this to @exr-proj-side-conditions.
:::

### C. Going deeper

::: {#exr-est-min-norm-scale}
[C1]

Replace \( \X \) by \( \X\bD \) for a nonsingular diagonal \( \bD \) (a change of units). Show that the minimum-norm
solution for \( \X\bD \), converted back to the original units, is the solution for the side condition
\( \V_0\T\bD^{-2}\bb=\bzero \), which differs in general from \( \V_0\T\bb=\bzero \). Conclude that the minimum-norm
convention depends on the units of measurement, while no estimable function does.
:::

::: {#exr-est-bad-side}
[C2]

Suppose \( \rank(\bT)=p-r \) but \( \C(\bT\T)\cap\C(\X\T)\ne\{\bzero\} \). Show that for some \( \y \) the restricted normal
equations of @thm-est-side-conditions(a) have no solution, and that whenever they have one they have
infinitely many.
:::

::: {.solution}
By the dimension count, \( \C(\X\T)+\C(\bT\T) \) has dimension
less than \( r+(p-r)=p \). So \( \Null(\X)\cap\Null(\bT)\ne\{\bzero\} \), and any solution can be moved along a nonzero vector
of that intersection. For existence, the proof of @thm-est-side-conditions shows that a solution exists for every
\( \y \) iff \( \Null(\X)+\Null(\bT)=\Real^p \), which fails because its orthogonal complement
\( \C(\X\T)\cap\C(\bT\T) \) is nonzero. So some coset \( \bb_0+\Null(\X) \) misses \( \Null(\bT) \), and the response \( \y=\X\bb_0 \)
has no solution.
:::
