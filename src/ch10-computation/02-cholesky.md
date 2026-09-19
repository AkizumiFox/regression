# The Cholesky factorization

If the cross-product matrix is to be used at all, the normal equations should be solved by
factoring it as a triangular matrix times its transpose. One factorization of the augmented
cross-product matrix then delivers every standard regression quantity.

## The algorithm

[Chapter 1](../ch01-matrix-algebra/index.html) proved that a positive definite matrix has
a unique factorization \( \A=\bL\bL\T \) with \( \bL \) lower triangular and positive on the diagonal
(@thm-mat-pd-characterizations(d)). Numerical work usually writes it with the upper
triangular factor \( \R=\bL\T \).

::: {#thm-cmp-cholesky}
[Cholesky factorization]

Let \( \A \) be a \( p\times p \) positive definite matrix.

::: {.enumerate options="label=(\alph*)"}
1. There is exactly one upper triangular \( \R \) with positive diagonal entries such that
           \( \A=\R\T\R \).

2. Its rows are given, for \( j=1,\dots,p \) in turn, by
           \[
           \begin{aligned}
           r_{jj}&=\Bigl(a_{jj}-\sum_{k<j}r_{kj}^2\Bigr)^{1/2},\\
           r_{ji}&=\frac1{r_{jj}}\Bigl(a_{ji}-\sum_{k<j}r_{kj}r_{ki}\Bigr),\quad i>j,
           \end{aligned}
           \]{#eq-cmp-cholesky}

   and every radicand in it is positive.

3. The computation takes \( p^3/3+O(p^2) \) flops.

4. If \( \A=\X\T\X \) and \( \X=\Q\R_0 \) is a QR factorization with \( \Q\T\Q=\I_p \) and \( \R_0 \)
           upper triangular with positive diagonal, then \( \R_0=\R \).
:::

:::

::: {.proof}
(a) is @thm-mat-pd-characterizations(d) with \( \R=\bL\T \).
(b) Entry \( (j,i) \), \( i\ge j \), of \( \R\T\R \) is \( \sum_{k\le j}r_{kj}r_{ki} \), because
\( r_{kj}=0 \) for \( k>j \). Setting it equal to \( a_{ji} \) and isolating the \( k=j \) term gives the
recursion @eq-cmp-cholesky. So the entries of the unique factor of (a) satisfy @eq-cmp-cholesky, and
@eq-cmp-cholesky determines row \( j \) from rows \( 1,\dots,j-1 \). By induction on \( j \), the rows
computed by @eq-cmp-cholesky are those of the factor in (a). In particular the \( j \)th radicand
is \( r_{jj}^2>0 \).
(c) Row \( j \) needs about \( 2(j-1) \) flops for each of its \( p-j+1 \) entries, so the total is
\( \sum_j2(j-1)(p-j+1)=p^3/3+O(p^2) \).
(d) \( \X\T\X=\R_0\T\Q\T\Q\R_0=\R_0\T\R_0 \), and uniqueness in (a).
:::

By part (d), the triangle of a QR factorization of \( \X \) and the Cholesky factor of \( \X\T\X \) are
the same matrix (see also @exr-mat-qr-unique). The two routes of this chapter compute the same
object, one from \( \X \), the other from \( \X\T\X \) after the small singular values may already have
been damaged (@prp-cmp-gram-rounding).

The algorithm itself is backward stable *for the factorization*. The computed \( \hat{\R} \)
satisfies \( \hat{\R}\T\hat{\R}=\A+\Delta\A \) with \( \lvert\Delta\A\rvert\le\gamma_{p+1}\lvert\hat{\R}\rvert\T\lvert\hat{\R}\rvert \)
entrywise, and it runs to completion whenever \( \A \) is not too close to singular
(Higham 2002, chapter 10), with no pivoting. The weakness of the normal equations lies entirely in forming \( \A=\X\T\X \), not in factoring it.

```{.python .run #cell-cholesky-cholesky}
import numpy as np
import statsmodels.api as sm

def cholesky_upper(A):
    """Upper triangular R with positive diagonal and R'R = A (A positive definite)."""
    A = np.array(A, dtype=float)
    p = A.shape[0]
    R = np.zeros_like(A)
    for j in range(p):
        d = A[j, j] - R[:j, j] @ R[:j, j]
        if d <= 0:
            raise np.linalg.LinAlgError(f"not positive definite at step {j + 1}")
        R[j, j] = np.sqrt(d)
        R[j, j + 1:] = (A[j, j + 1:] - R[:j, j] @ R[:j, j + 1:]) / R[j, j]
    return R

def solve_upper(R, b):
    """Back substitution for R x = b, R upper triangular."""
    x = np.zeros(np.shape(b))
    for i in range(len(b) - 1, -1, -1):
        x[i] = (b[i] - R[i, i + 1:] @ x[i + 1:]) / R[i, i]
    return x

def solve_lower(L, b):
    """Forward substitution for L x = b, L lower triangular."""
    x = np.zeros(np.shape(b))
    for i in range(len(b)):
        x[i] = (b[i] - L[i, :i] @ x[:i]) / L[i, i]
    return x
```

## Everything from one factorization

Solving the normal equations \( \R\T\R\bb=\X\T\y \) takes two triangular solves:
\( \R\T\bz=\X\T\y \) for \( \bz \), then \( \R\bb=\bz \). Each costs \( p^2 \) flops. It is neater
to factor the cross-product matrix of the augmented matrix \( \Z=[\X,\y] \). The factor then
contains the solution of the first triangular system and the residual sum of squares.

::: {#prp-cmp-augmented}
[The augmented factorization]

Let \( \X \) be \( n\times p \) of full column rank, suppose \( \y\notin\C(\X) \), and let
\[
\Z\T\Z=\begin{pmatrix}\X\T\X&\X\T\y\\\y\T\X&\y\T\y\end{pmatrix}=\bT\T\bT,
\qquad
\bT=\begin{pmatrix}\R&\bz\\\bzero\T&d\end{pmatrix},
\]
be the Cholesky factorization (@thm-cmp-cholesky) of the cross-product matrix of \( \Z=[\X,\y] \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \R \) is the Cholesky factor of \( \X\T\X \), the least squares estimate solves
           \( \R\hbeta=\bz \), and \( d^2=\text{SSE}=\norm{\y-\X\hbeta}^2 \);

2. for each \( k\le p \), the residual sum of squares of the regression of \( \y \) on the first \( k \)
           columns of \( \X \) is \( d^2+\sum_{j>k}z_j^2 \). The drop in residual sum of squares when
           column \( k \) joins columns \( 1,\dots,k-1 \) is \( z_k^2 \);

3. \( \det(\X\T\X)=\prod_jr_{jj}^2 \), \( (\X\T\X)^{-1}=\R^{-1}\R^{-\top} \), and the leverage of
           row \( \x_{(i)}\T \) is \( h_{ii}=\norm{\R^{-\top}\x_{(i)}}^2 \).
:::

:::

::: {.proof}
\( \Z \) has full column rank because \( \y\notin\C(\X) \), so \( \Z\T\Z \) is positive definite and
\( \bT \) exists. Multiplying out \( \bT\T\bT \) and comparing blocks with \( \Z\T\Z \) gives
\[
\R\T\R=\X\T\X,\qquad\R\T\bz=\X\T\y,\qquad\bz\T\bz+d^2=\y\T\y .
\]
The first identity and uniqueness make \( \R \) the Cholesky factor of \( \X\T\X \). The second gives
\( \R\T\R\hbeta=\R\T\bz \) for the solution of the normal equations, so \( \R\hbeta=\bz \).

For (b), and for the rest of (a), let \( \Q=\Z\bT^{-1} \). Then
\( \Q\T\Q=\bT^{-\top}\Z\T\Z\bT^{-1}=\I \), so the columns \( \mathbf{q}_1,\dots,\mathbf{q}_{p+1} \) of \( \Q \) are
orthonormal, and \( \Z=\Q\bT \). Because \( \bT \) is upper triangular with nonzero diagonal, the
first \( k \) columns of \( \X \) span the same subspace as \( \mathbf{q}_1,\dots,\mathbf{q}_k \). The last column of
\( \Z=\Q\bT \) reads \( \y=\sum_{j\le p}z_j\mathbf{q}_j+d\,\mathbf{q}_{p+1} \). By
@prp-proj-orthonormal-formula the projection of \( \y \) onto the span of the first \( k \) columns
is \( \sum_{j\le k}z_j\mathbf{q}_j \), and the residual sum of squares is
\( \norm{\y}^2-\sum_{j\le k}z_j^2=d^2+\sum_{j>k}z_j^2 \). With \( k=p \) this is \( d^2 \), which
completes (a). Taking differences gives the rest of (b).

(c) The determinant of a triangular matrix is the product of its diagonal. The inverse formula
follows from \( \X\T\X=\R\T\R \). Finally
\( h_{ii}=\x_{(i)}\T(\X\T\X)^{-1}\x_{(i)}=\x_{(i)}\T\R^{-1}\R^{-\top}\x_{(i)}=\norm{\R^{-\top}\x_{(i)}}^2 \).
:::

Part (b) says that the squares \( z_1^2,\dots,z_p^2 \) are the *sequential sums of squares*
of the regressors in the order of the columns (@def-ss-sequential). A single
factorization gives the whole nested sequence of fits. This is the content of
@thm-proj-nested, now in coordinates. Part (c) costs about \( np^2 \) flops for all \( n \)
leverages, each a triangular solve of \( p^2 \) flops.

::: {#exm-cmp-state-cholesky}
[Murder rates by the augmented factorization]

For the state data of @exm-proj-fwl-crime (\( n=50 \), with an intercept and the
percentages in poverty, in single-parent households and urban, so \( p=4 \)), the
augmented factorization (@prp-cmp-augmented) gives the coefficients
\( -9.0946 \), \( 0.2538 \), \( 0.3980 \) and \( 0.0046 \), and the corner entry
\( d=10.0976 \), so \( \text{SSE}=d^2=101.96 \). The condition number of \( \X \) is only
\( 509 \), and the solution agrees with a QR solution to a relative difference of
\( 2\times 10^{-14} \). The largest leverage is \( 0.260 \), and
\( \log_{10}\det(\X\T\X)=11.14 \). The script also checks the standard errors
\( s\,\norm{\mathbf{e}_j\T\R^{-1}} \) against statsmodels.
:::

```{.python .run #cell-cholesky-data}
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
```

```{.python .run #cell-cholesky-augmented}
Z = np.column_stack([X, y])              # augmented matrix [X, y]
S = Z.T @ Z                              # (p+1) x (p+1) cross-product matrix
T = cholesky_upper(S)                    # T = [[R, z], [0, d]]
R, z, d = T[:p, :p], T[:p, p], T[p, p]

beta = solve_upper(R, z)                 # R b = z
sse = d ** 2                             # the last diagonal entry is sqrt(SSE)
K = solve_upper(R, np.eye(p))            # R^{-1}
XtX_inv = K @ K.T                        # (X'X)^{-1} = R^{-1} R^{-T}
h = np.array([np.sum(solve_lower(R.T, x) ** 2) for x in X])   # leverages
print("coefficients:", np.round(beta, 4))
print(f"SSE = {sse:.4f},  d = {d:.4f}")
```

## When the normal equations are good enough

@prp-cmp-gram-rounding and @eq-cmp-rule predict a relative error of order \( \kappa^2u \) in \( \hbeta \) from
the normal equations, where, as [Section 10.5](05-conditioning.html) shows, \( \kappa \) is the condition
number of \( \X \) after its columns are scaled to equal length. With \( \kappa\le10^4 \) about eight digits
survive, more than any coefficient is known statistically. The normal equations are therefore sound
for well-conditioned designs, when only cross-products are available, and inside iterations such as the
reweighted least squares of Chapter 34. Nearly collinear observational data are where the \( \kappa^2 \)
loss bites.

## The sweep operator

Before QR became standard, regression programs worked on the augmented cross-product matrix
with an operation that adds one variable at a time to the fit. It is still used for fitting
many subsets of variables.

::: {#def-cmp-sweep}
[Sweep]

Let \( \A \) be symmetric with \( a_{kk}\ne0 \). **Sweeping** \( \A \) on pivot \( k \) produces the
symmetric matrix \( \B \) with
\[
b_{kk}=-\frac1{a_{kk}},\qquad b_{ik}=b_{ki}=\frac{a_{ik}}{a_{kk}},\qquad
b_{ij}=a_{ij}-\frac{a_{ik}a_{kj}}{a_{kk}}\quad(i,j\ne k).
\]
:::

::: {#prp-cmp-sweep}
[What sweeping computes]

Let \( \bS=\Z\T\Z \) be the cross-product matrix of \( \Z=[\bz_1,\dots,\bz_m] \). Let \( K \) be a set of
column indices such that \( \Z_K \) has full column rank, and let \( L \) be the remaining
indices. Sweeping \( \bS \) on the pivots in \( K \), one at a time and in any order, meets only
positive pivots and produces the matrix \( \W \) with blocks
\[
\W_{KK}=-(\Z_K\T\Z_K)^{-1},\qquad
\W_{KL}=\W_{LK}\T=(\Z_K\T\Z_K)^{-1}\Z_K\T\Z_L,\qquad
\W_{LL}=\Z_L\T(\I-\M_K)\Z_L ,
\]
where \( \M_K \) is the projection onto \( \C(\Z_K) \). Column \( j\in L \) of \( \W_{KL} \) holds the
least squares coefficients of \( \bz_j \) regressed on \( \Z_K \), and \( \W_{LL} \) holds the residual
cross-products. With \( \Z=[\X,\y] \) and \( K \) the columns of \( \X \), the last column of \( \W \)
contains \( \hbeta \) and, in its corner, the residual sum of squares.
:::

::: {.proof}
We induct on the size of \( K \); the case \( K=\emptyset \) is \( \W=\bS \). Add a pivot \( m\in L \) with
\( \Z_{K\cup\{m\}} \) of full column rank, and write \( \G=(\Z_K\T\Z_K)^{-1} \), \( \B_j=\G\Z_K\T\bz_j \) and
\( \bw=(\I-\M_K)\bz_m \). The pivot is \( \bz_m\T(\I-\M_K)\bz_m=\norm{\bw}^2>0 \). For \( i,j\in L\setminus\{m\} \) the new
entries are \( \bz_i\T(\I-\M_K-\bw\bw\T/\norm{\bw}^2)\bz_j=\bz_i\T(\I-\M_{K\cup\{m\}})\bz_j \) by @lem-proj-fwl-split.
The new entries \( \bw\T\bz_j/\norm{\bw}^2 \) in row \( m \), and \( \B_j-\B_m\bw\T\bz_j/\norm{\bw}^2 \) in the rows of
\( K \), are the coefficients of \( \bz_j \) on \( [\Z_K,\bz_m] \) by @thm-proj-fwl. The partitioned inverse
(@thm-mat-partitioned-inverse), with Schur complement \( \norm{\bw}^2 \), shows that
\( (\Z_{K\cup\{m\}}\T\Z_{K\cup\{m\}})^{-1} \) has blocks \( \G+\B_m\B_m\T/\norm{\bw}^2 \), \( -\B_m/\norm{\bw}^2 \) and
\( 1/\norm{\bw}^2 \), whose negatives are what the sweep produces: \( -\G-\B_m\B_m\T/\norm{\bw}^2 \),
\( \B_m/\norm{\bw}^2 \) and \( -1/\norm{\bw}^2 \). The result depends only on the set \( K \), so the
order of the sweeps does not matter.
:::

Sweeping a pivot costs about \( (p+1)^2 \) flops, so a new variable enters a fitted model
for \( O(p^2) \) work. The price is the same as for the normal equations: the sweep works with
cross-products and inherits their \( \kappa^2 \) sensitivity. For the state data, sweeping all
four columns reproduces the coefficients and SSE of @exm-cmp-state-cholesky. Sweeping only
the intercept and the single-parent column puts the SSE of that two-column model,
\( 121.92 \), in the corner.

```{.python .run #cell-cholesky-sweep}
def sweep(A, k):
    """Sweep the symmetric matrix A on pivot k."""
    A = np.array(A, dtype=float)
    a = A[k, k]
    B = A - np.outer(A[:, k], A[k, :]) / a
    B[k, :] = A[k, :] / a
    B[:, k] = A[:, k] / a
    B[k, k] = -1 / a
    return B

W = S.copy()
for k in range(p):                       # sweep the p columns of X, one at a time
    W = sweep(W, k)
print("coefficients from the sweep:", np.round(W[:p, p], 4))
print(f"SSE from the sweep: {W[p, p]:.4f}")
```

## Forming cross-products: the centring problem

With an intercept in the model, the cross-product matrix is usually computed in centred form,
\( \sum_i(x_{ij}-\bar x_j)(x_{ik}-\bar x_k) \), which by @thm-proj-reparam changes the
parameterization but not the fit. How it is computed matters a great deal. The
"textbook" formula \( \sum_ix_i^2-n\bar x^2 \) subtracts two large, nearly equal numbers
whenever the mean is large compared with the spread. The loss is governed by the ratio
\[
\kappa_{\x}=\frac{\norm{\x}}{\norm{\x-\bar x\bone}}=\Bigl(1+\frac{n\bar x^2}{\sum_i(x_i-\bar x)^2}\Bigr)^{1/2},
\]
which is the condition number of the column, and it is large when the coefficient of
variation is small. Three algorithms are in common use:

- *textbook*: accumulate \( \sum x_i \) and \( \sum x_i^2 \) in one pass, then subtract;
- *two-pass*: compute \( \bar x \) in a first pass, then sum \( (x_i-\bar x)^2 \) in a second;
- *updating*: one pass, maintaining the running mean and the running centred sum by the
  recursions below.

::: {#prp-cmp-welford}
[Updating means and centred cross-products]

For vectors \( \x_1,\x_2,\dots \) in \( \Real^q \) let \( \mathbf{m}_k \) be the mean of the first \( k \) and
\( \mathbf{C}_k=\sum_{i\le k}(\x_i-\mathbf{m}_k)(\x_i-\mathbf{m}_k)\T \). With \( \mathbf{d}_k=\x_k-\mathbf{m}_{k-1} \),
\[
\mathbf{m}_k=\mathbf{m}_{k-1}+\frac{\mathbf{d}_k}{k},\qquad
\mathbf{C}_k=\mathbf{C}_{k-1}+\frac{k-1}{k}\,\mathbf{d}_k\mathbf{d}_k\T .
\]
:::

::: {.proof}
The first identity is \( k\mathbf{m}_k=(k-1)\mathbf{m}_{k-1}+\x_k \). For the second, write
\( \mathbf{C}_k=\sum_{i\le k}\x_i\x_i\T-k\,\mathbf{m}_k\mathbf{m}_k\T \) and put \( \mathbf{m}=\mathbf{m}_{k-1} \),
\( \mathbf{d}=\mathbf{d}_k \). Then
\[
\mathbf{C}_k-\mathbf{C}_{k-1}=\x_k\x_k\T-k\,\mathbf{m}_k\mathbf{m}_k\T+(k-1)\mathbf{m}\mathbf{m}\T .
\]
Substituting \( \x_k=\mathbf{m}+\mathbf{d} \) and \( k\mathbf{m}_k\mathbf{m}_k\T=k\mathbf{m}\mathbf{m}\T+\mathbf{m}\mathbf{d}\T+\mathbf{d}\mathbf{m}\T+\mathbf{d}\mathbf{d}\T/k \),
all the terms in \( \mathbf{m} \) cancel and \( \mathbf{d}\mathbf{d}\T(1-1/k) \) remains.
:::

The update adds a *centred* rank-one term, so it never subtracts large quantities. It is the
cross-product analogue of the row updates of [Section 10.6](06-updating.html). Chan, Golub and
LeVeque (1983) give error bounds of roughly \( n\kappa_{\x}^2u \) for the textbook formula,
\( n\kappa_{\x}u \) for updating and \( nu+n^2\kappa_{\x}^2u^2 \) for the two-pass algorithm.

::: {#exm-cmp-centring}
[Centring a column with a large mean]

A sample of \( 1000 \) standard normal values is shifted by \( c=1,10,\dots,10^8 \), and the centred
sum of squares is computed by the three algorithms and compared with its exact value (found in
rational arithmetic). At \( c=10^6 \), where \( \kappa_{\x}=1.0\times 10^{6} \), the relative errors are
\( 8.9\times 10^{-5} \) for the textbook formula, \( 8.1\times 10^{-12} \) for updating and
\( 2.7\times 10^{-17} \) for two passes. At \( c=10^8 \) the textbook formula has relative error
\( 3.1 \): not a single correct digit, and possibly a negative "sum of squares". Updating
gives \( 7.4\times 10^{-10} \) and two passes \( 2.3\times 10^{-17} \).
[Figure 10.2.1](#fig-cmp-centring) shows the three growth rates \( \kappa^2u \), \( \kappa u \)
and \( u \).
:::

::: {when-format="html"}
![**Figure 10.2.1.** Relative error of the centred sum of squares of a column of \( 1000 \) values
with standard deviation one and mean up to \( 10^8 \), against the column condition number
\( \kappa_{\x} \). The one-pass textbook formula loses accuracy like \( \kappa^2u \), the one-pass
updating formula like \( \kappa u \), and the two-pass algorithm not at all.](centring_accuracy.svg){#fig-cmp-centring width=62%}
:::

::: {when-format="pdf"}
![Relative error of the centred sum of squares of a column of \( 1000 \) values
with standard deviation one and mean up to \( 10^8 \), against the column condition number
\( \kappa_{\x} \). The one-pass textbook formula loses accuracy like \( \kappa^2u \), the one-pass
updating formula like \( \kappa u \), and the two-pass algorithm not at all.](centring_accuracy.pdf){width=62%}
:::

```{.python .run #cell-cholesky-centring}
def ss_textbook(x):
    return np.sum(x * x) - len(x) * np.mean(x) ** 2        # one pass, then subtract

def ss_two_pass(x):
    m = np.mean(x)                                         # first pass: the mean
    return np.sum((x - m) ** 2)                            # second pass: deviations

def ss_updating(x):
    m, s = 0.0, 0.0
    for k, xk in enumerate(x, start=1):                    # one pass, no cancellation
        d = xk - m
        m = m + d / k
        s = s + (k - 1) / k * d * d
    return s
```

Calendar years, population totals and price indices have large means and small spread. A program
that forms cross-products should centre them with two passes (or update), solve for the slopes, and
recover the intercept as \( \bar y-\bar{\x}\T\hbeta_{\text{slopes}} \).

## Exercises

### A. Check your understanding

::: {#exr-cmp-chol-hand}
[A1]

Compute the Cholesky factor of
\( \A=\begin{psmallmatrix}4&2&-2\\2&10&2\\-2&2&6\end{psmallmatrix} \) by @eq-cmp-cholesky, and
use it to solve \( \A\bb=(2,14,6)\T \) and to find \( \det\A \).
:::

::: {.solution}
Row 1: \( r_{11}=2 \), \( r_{12}=1 \), \( r_{13}=-1 \). Row 2: \( r_{22}=\sqrt{10-1}=3 \),
\( r_{23}=(2-(1)(-1))/3=1 \). Row 3: \( r_{33}=\sqrt{6-1-1}=2 \). So
\( \R=\begin{psmallmatrix}2&1&-1\\0&3&1\\0&0&2\end{psmallmatrix} \) and
\( \det\A=(2\cdot3\cdot2)^2=144 \). Forward substitution in \( \R\T\bz=(2,14,6)\T \) gives
\( z_1=1 \), \( z_2=(14-1)/3=13/3 \), \( z_3=(6+1-13/3)/2=4/3 \). Back substitution in \( \R\bb=\bz \)
gives \( b_3=2/3 \), \( b_2=(13/3-2/3)/3=11/9 \) and \( b_1=(1-11/9+2/3)/2=2/9 \). Check:
\( \A\bb=(8/9+22/9-12/9,\ 4/9+110/9+12/9,\ -4/9+22/9+36/9)=(2,14,6) \).
:::

::: {#exr-cmp-chol-simple}
[A2]

For simple regression, \( \X=[\bone,\x] \), find the Cholesky factor (@thm-cmp-cholesky) of \( \X\T\X \) and of the
augmented matrix, and show that it reproduces the slope \( S_{xy}/S_{xx} \) of @thm-lm-simple-ls.
Where does the centring happen?
:::

::: {.solution}
\( \X\T\X=\begin{psmallmatrix}n&n\bar x\\n\bar x&\sum x_i^2\end{psmallmatrix} \) gives
\( r_{11}=\sqrt n \), \( r_{12}=\sqrt n\,\bar x \) and
\( r_{22}=(\sum x_i^2-n\bar x^2)^{1/2}=S_{xx}^{1/2} \). The augmented column gives
\( z_1=n\bar y/\sqrt n=\sqrt n\,\bar y \) and \( z_2=(\sum x_iy_i-n\bar x\bar y)/S_{xx}^{1/2}=S_{xy}/S_{xx}^{1/2} \).
Back substitution gives the slope \( b_2=z_2/r_{22}=S_{xy}/S_{xx} \) and the intercept
\( b_1=(z_1-r_{12}b_2)/r_{11}=\bar y-\bar xb_2 \). The centring happens in the subtraction
\( a_{22}-r_{12}^2 \), which is the textbook formula of @exm-cmp-centring. So the Cholesky
factorization of an uncentred cross-product matrix inherits its cancellation.
:::

### B. Practice

::: {#exr-cmp-se-from-R}
[B1]

Show that the standard error of \( \hat\beta_j \) is \( s\,\norm{\mathbf{e}_j\T\R^{-1}} \), where
\( s^2=\text{SSE}/(n-p) \) (@thm-lm-sigma2). Show that all \( p \) standard errors cost \( p^3/3+O(p^2) \)
flops, and all \( n \) leverages \( np^2+O(np) \) flops. Why is computing \( \M \) itself a bad idea?
:::

::: {.solution}
By @thm-lm-moments the estimated variance of \( \hat\beta_j \) is \( s^2[(\X\T\X)^{-1}]_{jj} \), and
\( [(\X\T\X)^{-1}]_{jj}=\mathbf{e}_j\T\R^{-1}\R^{-\top}\mathbf{e}_j=\norm{\mathbf{e}_j\T\R^{-1}}^2 \). The inverse of a
triangular matrix is triangular, and column \( j \) of \( \R^{-1} \) is found by back substitution on
its leading \( j\times j \) block, at about \( j^2 \) flops. The total is \( \sum_jj^2\approx p^3/3 \). Each
leverage is a triangular solve with \( \R\T \) (\( p^2 \) flops) and a squared norm. The projection
\( \M \) has \( n^2 \) entries, which for \( n=10^6 \) is too many to store, and only its diagonal is ever
needed.
:::

::: {#exr-cmp-pooling}
[B2]

Two groups of observations have sizes \( n_a,n_b \), means \( \mathbf{m}_a,\mathbf{m}_b \) and centred
cross-product matrices \( \mathbf{C}_a,\mathbf{C}_b \). Show that the pooled data have centred
cross-product matrix
\[
\mathbf{C}=\mathbf{C}_a+\mathbf{C}_b+\frac{n_an_b}{n_a+n_b}(\mathbf{m}_a-\mathbf{m}_b)(\mathbf{m}_a-\mathbf{m}_b)\T .
\]
Show that @prp-cmp-welford is the special case \( n_b=1 \), and explain how the formula allows
cross-products to be computed in parallel on separate blocks of data.
:::

### C. Going deeper

::: {#exr-cmp-pivoted-cholesky}
[C1]

Let \( \A \) be nonnegative definite of rank \( r \). Show that the Cholesky recursion with
*symmetric pivoting*, which at step \( j \) moves the largest remaining diagonal entry of the
current Schur complement to position \( j \), produces a permutation \( \boldsymbol{\Pi} \) and
\( \boldsymbol{\Pi}\T\A\boldsymbol{\Pi}=\R\T\R \) with \( \R=\begin{psmallmatrix}\R_{11}&\R_{12}\\\bzero&\bzero\end{psmallmatrix} \)
and \( \R_{11} \) an \( r\times r \) nonsingular upper triangle. (Use @prp-mat-pd-properties(a).)
How is this related to the pivoted QR of [Section 10.7](07-rank-revealing.html)?
:::
