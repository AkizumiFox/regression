# Least squares via projection

We now return to the linear model \( \Y=\X\bbeta+\be \), with \( \X \) an \( n\times p \)
matrix of rank \( r\le p \), and \( \M \) the orthogonal projection onto \( \C(\X) \). Nothing in this
section uses probability. The results hold for every data vector \( \y \) and so
apply whatever the distribution of the errors.

## The characterization of least squares estimates

::: {#thm-proj-ls-projection}
[Least squares is projection]

A vector \( \hbeta \) is a least squares estimate of \( \bbeta \) iff
\[
\X\hbeta=\M\y .
\]
Consequently the fitted values \( \hY=\M\y \) and residuals \( \he=(\I-\M)\y \) are the
same for every least squares estimate.
:::

::: {.proof}
Apply the projection theorem (@thm-proj-projection-theorem) with
\( \mathcal S=\C(\X) \) and \( \bu=\X\bb \). By @eq-proj-distance-split,
\[
\norm{\y-\X\bb}^2=\norm{\y-\M\y}^2+\norm{\M\y-\X\bb}^2
\qquad\text{for every }\bb\in\Real^p .
\]
The first term does not involve \( \bb \), and the second is zero iff \( \X\bb=\M\y \).
Such a \( \bb \) exists because \( \M\y\in\C(\X) \).
:::

The orthogonality half of the projection theorem turns this into linear
equations. The residual \( \y-\X\bb \) is orthogonal to \( \C(\X) \) iff it is orthogonal
to every column of \( \X \), that is, iff \( \X\T(\y-\X\bb)=\bzero \).

::: {#thm-proj-normal-equations}
[Normal equations]

The **normal equations**
\[
\X\T\X\bb=\X\T\y
\]{#eq-proj-normal-equations}

are consistent for every \( \y \), and their solutions are exactly the least squares
estimates. One solution is \( \hbeta=(\X\T\X)\ginv\X\T\y \), for any generalized inverse,
and the set of all solutions is
\[
\hbeta+\Null(\X)=\{\hbeta+\bv:\X\bv=\bzero\}.
\]
In particular, the least squares estimate is unique iff \( \rank(\X)=p \), and then
\( \hbeta=(\X\T\X)^{-1}\X\T\y \).
:::

::: {.proof}
A vector \( \bb \) solves @eq-proj-normal-equations iff \( \y-\X\bb\perp\C(\X) \).
Since \( \X\bb\in\C(\X) \) automatically, this holds iff \( \y=\X\bb+(\y-\X\bb) \) is the
orthogonal decomposition of \( \y \) (@thm-proj-direct-sum), that is, iff \( \X\bb=\M\y \). By
@thm-proj-ls-projection this is the same as \( \bb \) being a least squares
estimate. For the particular solution, \( \X(\X\T\X)\ginv\X\T\y=\M\y \) by @thm-proj-M-formula.
Finally, \( \X\bb=\M\y=\X\hbeta \) iff
\( \X(\bb-\hbeta)=\bzero \), and \( \Null(\X)=\{\bzero\} \) iff \( \rank(\X)=p \) by @eq-proj-rank-nullity.
:::

The logic of this proof runs opposite to the usual calculus derivation, and
the difference matters. Setting a gradient to zero finds stationary points. One
then has to argue separately that they are minima and that one exists. Here
existence, minimality and the full solution set all come from the
decomposition \( \Real^n=\C(\X)\dirsum\C(\X)\perpc \). The consistency of the normal
equations, which can look surprising when \( \X\T\X \) is singular, is just the
statement \( \X\T\y\in\C(\X\T)=\C(\X\T\X) \) from @cor-proj-gram-colspace.

## Fitted values and residuals

Collecting the consequences for \( \hY=\M\y \) and \( \he=(\I-\M)\y \):

::: {#prp-proj-fit-residual}
For any data vector \( \y \):

::: {.enumerate options="label=(\alph*)"}
1. \( \X\T\he=\bzero \) and \( \hY\T\he=0 \);

2. if \( \bone\in\C(\X) \) (for example, the model has an intercept), then
           \( \sum_i\hat{\varepsilon}_i=0 \) and the mean of the fitted values is \( \bar y \);

3. \( \norm{\y}^2=\norm{\hY}^2+\norm{\he}^2 \);

4. \( \norm{\he}^2=\y\T(\I-\M)\y \) and \( \norm{\hY}^2=\y\T\M\y \);

5. if \( \y=\X\bbeta+\be \) for some \( \bbeta \), then
           \( \hY=\X\bbeta+\M\be \) and \( \he=(\I-\M)\be \).
:::

:::

::: {.proof}
(a) is the orthogonality of \( \he \) to \( \C(\X) \), and \( \hY\in\C(\X) \). For (b), \( \bone\in\C(\X) \) gives
\( \bone\T\he=0 \), and then \( \bone\T\hY=\bone\T\y \). (c) is Pythagoras. (d) uses symmetry
and idempotence: \( \norm{(\I-\M)\y}^2=\y\T(\I-\M)\T(\I-\M)\y=\y\T(\I-\M)\y \). For (e),
\( \M\X=\X \) and \( (\I-\M)\X=\bzero \).
:::

Part (e) is simple but important. The residual vector does not depend on
\( \bbeta \) at all. It is the projection of the unobservable error vector onto the
\( (n-r) \)-dimensional space \( \C(\X)\perpc \). Everything about \( \sigma^2 \) that the
data can reveal must therefore come from \( (\I-\M)\be \). If
\( \E(\be)=\bzero \) and \( \Cov(\be)=\sigma^2\I \), the expectation of a quadratic form
([Chapter 2](../ch02-random-vectors/index.html)) gives
\[
\E\norm{\he}^2=\E\bigl[\be\T(\I-\M)\be\bigr]=\sigma^2\tr(\I-\M)=\sigma^2(n-r),
\]
using @prp-proj-trace-rank. So \( \norm{\he}^2/(n-r) \) is unbiased for
\( \sigma^2 \). The divisor is the dimension of the space in which the residuals live,
not the number of observations and not \( n-p \). [Chapter 7](../ch07-optimality/index.html) develops the
sampling theory, and [Chapter 9](../ch09-sums-of-squares/index.html) develops the decomposition.

## Coefficients are not unique; fitted values are

When \( \rank(\X)<p \), @thm-proj-normal-equations gives a whole affine
subspace of least squares estimates, one for each generalized inverse and more
besides. Statistical software has to pick one. Some programs drop columns
that are linearly dependent on earlier ones. Others impose “side
conditions” such as a zero coefficient for a reference level or effects that sum to
zero. Others report the minimum-norm solution. These choices look very different
in the printed coefficient table. The geometry says they are all describing the
same fitted vector.

::: {#exm-proj-ginverse-numeric}
[Four answers, one fit]

Nine observations fall in three groups of sizes \( 3,2,4 \) with group means
\( 5.700 \), \( 7.800 \) and \( 4.300 \). The model
matrix \( \X=[\bone,\bz_1,\bz_2,\bz_3] \) has \( p=4 \) columns and rank
\( 3 \), so \( \X\T\X \) is singular. We compute \( \hbeta=\G\X\T\y \) for four
generalized inverses \( \G \) of \( \X\T\X \):
the Moore–Penrose inverse; the inverse obtained by deleting the intercept row and
column (@thm-mat-ginverse-exists); the one obtained by deleting the first
group's row and column; and a deliberately strange, non-symmetric one built from
the general formula of @exr-proj-all-ginverses.

| Generalized inverse | \( \hat{\beta}_0 \) | \( \hat{\beta}_1 \) | \( \hat{\beta}_2 \) | \( \hat{\beta}_3 \) | \( \hat{\beta}_1-\hat{\beta}_2 \) | \( \max\lvert\M-\M_{\text{MP}}\rvert \) |
|---|---|---|---|---|---|---|
| Moore–Penrose | 4.450 | 1.250 | 3.350 | -0.150 | -2.100 | 0 |
| drop intercept | 0.000 | 5.700 | 7.800 | 4.300 | -2.100 | \(1.1\times 10^{-16}\) |
| drop group 1 | 5.700 | 0.000 | 2.100 | -1.400 | -2.100 | \(1.1\times 10^{-16}\) |
| random | -9.102 | 14.802 | 16.902 | 13.402 | -2.100 | \(1.8\times 10^{-15}\) |

The coefficient vectors disagree wildly, including in sign. Yet each produces
the same projection matrix up to rounding error, and hence the same fitted values
(the group means) and residuals. The difference \( \hat{\beta}_1-\hat{\beta}_2 \)
is the same in every row, and so is \( \hat{\beta}_0+\hat{\beta}_k \) for each \( k \).
Individual coefficients are not. The trace of \( \M \) is
\( 3.000 \), the rank of \( \X \).
:::

```{.python .run #cell-ginverse-invariance-layout}
import numpy as np

groups = np.array([0, 0, 0, 1, 1, 2, 2, 2, 2])     # three groups, sizes 3, 2, 4
y = np.array([5.1, 6.3, 5.7, 8.2, 7.4, 4.0, 4.9, 3.8, 4.5])
n, g = len(y), 3
X = np.column_stack([np.ones(n), np.eye(g)[groups]])  # [1, indicator columns]
A = X.T @ X
print("n =", n, " p =", X.shape[1], " rank(X) =", np.linalg.matrix_rank(X))
```

```{.python .run #cell-ginverse-invariance-ginverses}
def g_drop(j):
    """Generalized inverse obtained by deleting row and column j of A."""
    keep = [k for k in range(A.shape[0]) if k != j]
    G = np.zeros_like(A)
    G[np.ix_(keep, keep)] = np.linalg.inv(A[np.ix_(keep, keep)])
    return G

rng = np.random.default_rng(6)
G0 = np.linalg.pinv(A)                             # Moore-Penrose inverse
U, V = rng.normal(size=A.shape), rng.normal(size=A.shape)
G_random = G0 + (np.eye(4) - G0 @ A) @ U + V @ (np.eye(4) - A @ G0)

ginverses = {
    "Moore-Penrose": G0,
    "drop intercept": g_drop(0),
    "drop group 1": g_drop(1),
    "random": G_random,
}
```

```{.python .run #cell-ginverse-invariance-compare}
M_ref = X @ G0 @ X.T
for name, G in ginverses.items():
    assert np.allclose(A @ G @ A, A)               # G really is a g-inverse
    beta = G @ X.T @ y
    M = X @ G @ X.T
    print(f"{name:15s} beta = {np.round(beta, 3)}",
          f" max|M - M_ref| = {np.abs(M - M_ref).max():.1e}",
          f" group1-group2 = {beta[1] - beta[2]:.3f}")
```

Which linear functions of \( \hbeta \) come out the same for every least squares
estimate? The answer is a first look at estimability, the subject of
[Chapter 8](../ch08-estimability/index.html).

::: {#thm-proj-invariant-functions}
For \( \blambda\in\Real^p \), the value \( \blambda\T\hbeta \) is the same for every least
squares estimate \( \hbeta \) (and every \( \y \)) iff \( \blambda\in\C(\X\T) \), that is, iff
\( \blambda\T=\boldsymbol{\uprho}\T\X \) for some \( \boldsymbol{\uprho}\in\Real^n \). In that case
\( \blambda\T\hbeta=\boldsymbol{\uprho}\T\M\y \).
:::

::: {.proof}
By @thm-proj-normal-equations the least squares estimates form
\( \hbeta+\Null(\X) \). So \( \blambda\T\hbeta \) is invariant iff \( \blambda\T\bv=0 \) for all
\( \bv\in\Null(\X) \), that is, iff \( \blambda\in\Null(\X)\perpc=\C(\X\T) \)
(@cor-proj-gram-colspace). If \( \blambda\T=\boldsymbol{\uprho}\T\X \) then
\( \blambda\T\hbeta=\boldsymbol{\uprho}\T\X\hbeta=\boldsymbol{\uprho}\T\M\y \).
:::

In @exm-proj-ginverse-numeric, \( \blambda=(0,1,-1,0)\T \) is the difference
between any row of \( \X \) from group 1, \( (1,1,0,0) \), and any row from group 2,
\( (1,0,1,0) \). So \( \blambda\in\C(\X\T) \) and \( \hat{\beta}_1-\hat{\beta}_2 \) is invariant.
For \( \blambda=(0,1,0,0)\T \), note that \( \boldsymbol{\uprho}\T\X \) has first entry
\( \sum_i\rho_i \) and remaining entries equal to the sums of \( \rho_i \) within each
group. Matching \( (0,1,0,0) \) would need the group sums to be \( 1,0,0 \) and their
total to be \( 0 \), which is impossible. So \( \hat{\beta}_1 \) is not invariant.

::: {.remark}
[Side conditions choose coordinates]

Imposing \( \beta_0=0 \) or \( \beta_1=0 \) does not add information about the data.
It picks one point from the affine set \( \hbeta+\Null(\X) \). A side condition
\( \mathbf{c}\T\bbeta=0 \) picks exactly one point iff \( \mathbf{c}\notin\C(\X\T) \). A condition
inside \( \C(\X\T) \) would constrain the fitted values and change the model rather
than just its coordinates. @thm-est-side-conditions makes this precise.
:::

## Exercises

### A. Check your understanding

::: {#exr-proj-three-points-M}
[A1]

For @exm-proj-three-points, write out the \( 3\times3 \) matrix \( \M \) with exact fractional
entries, verify \( \tr\M=2 \), and recover the fitted values and residuals quoted in the text.
:::

::: {#exr-proj-residual-sum}
[A2]

Show that the residuals sum to zero whenever \( \bone\in\C(\X) \), even if no column of \( \X \)
is constant. Give an example of a model without such a column in which this happens, and
one in which the residuals do not sum to zero.
:::

### B. Practice

::: {#exr-proj-moments-fit}
[B1]

Suppose \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\I \). Show that \( \E(\hY)=\X\bbeta \),
\( \Cov(\hY)=\sigma^2\M \), \( \Cov(\he)=\sigma^2(\I-\M) \) and \( \Cov(\hY,\he)=\bzero \). Which of
these statements survive when \( \Cov(\Y)=\sigma^2\V \) for a general positive definite \( \V \)?
:::

::: {#exr-proj-side-conditions}
[B2]

Let \( \rank(\X)=r<p \) and let \( \mathbf{C} \) be a \( (p-r)\times p \) matrix. Show that the
constrained problem “minimize \( \norm{\y-\X\bb} \) subject to \( \mathbf{C}\bb=\bzero \)” has a
unique solution for every \( \y \), and that its fitted values are \( \M\y \), iff
\( \rank\begin{psmallmatrix}\X\\\mathbf{C}\end{psmallmatrix}=p \). Show that the condition implies
\( \C(\mathbf{C}\T)\cap\C(\X\T)=\{\bzero\} \), and interpret: the constraints must not restrict
any estimable function.
:::

### C. Going deeper

::: {#exr-proj-every-solution-ginverse}
[C1]

Suppose \( \X\T\y\neq\bzero \). Show that every least squares estimate \( \bb \) can be written as
\( \G\X\T\y \) for some generalized inverse \( \G \) of \( \A=\X\T\X \). *Hint:* with
\( \mathbf{c}=\X\T\y \), try \( \G=\G_0+(\bb-\G_0\mathbf{c})\mathbf{c}\T/(\mathbf{c}\T\mathbf{c}) \). What goes wrong when
\( \X\T\y=\bzero \)?
:::
