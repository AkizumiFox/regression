# Rank-revealing factorizations

A model matrix with an intercept and indicators for every level of a factor has dependent, or
**aliased**, columns. Software must notice them and decide which to set aside, but in floating point it
can only ask whether a column is *nearly* a combination of others, relative to a tolerance. The SVD
answers that ([Section 10.4](04-svd.html)) but does not say which columns to drop, and it costs more. The
standard tool is QR with column pivoting.

## QR with column pivoting

At step \( k \) of Householder QR, the trailing columns (rows \( k \) to \( n \), columns \( k \) to \( p \) of the
partly reduced matrix) are what remains of the original columns after the parts explained by the
first \( k-1 \) chosen columns have been removed. **Column pivoting** (Businger and Golub 1965)
moves the trailing column of largest norm, the one least explained by the columns already chosen,
into position \( k \) before reflecting. The result is
\[
\X\boldsymbol{\Pi}=\Q\begin{pmatrix}\R\\\bzero\end{pmatrix},\qquad
\R=\begin{pmatrix}\R_{11}&\R_{12}\\\bzero&\R_{22}\end{pmatrix},
\]
with \( \boldsymbol{\Pi} \) a permutation matrix, \( \R_{11} \) the leading \( k\times k \) block and \( \R_{22} \) the trailing
\( (p-k)\times(p-k) \) block, for any \( k \). We need one fact about submatrices.

::: {#lem-cmp-column-subset}
[Singular values of a submatrix]

If \( \X_1 \) is formed from some of the columns of \( \X \), then \( \sigma_i(\X_1)\le\sigma_i(\X) \) for every \( i \) up
to the number of columns of \( \X_1 \). The same holds for a submatrix formed from some of the rows.
:::

::: {.proof}
Write \( \X_1=\X\bS \), where \( \bS \) holds the corresponding columns of the identity, so \( \norm{\bS}_2=1 \). Let
\( \X_{i-1} \) be the truncated SVD with \( i-1 \) terms, or \( \X \) itself if \( \rank\X<i \), so that
\( \norm{\X-\X_{i-1}}_2=\sigma_i(\X) \). The matrix \( \X_{i-1}\bS \) has rank at most \( i-1 \), so by the
Eckart–Young inequality, in the form used in the proof of @prp-cmp-weyl,
\[
\sigma_i(\X_1)\le\norm{\X\bS-\X_{i-1}\bS}_2\le\norm{\X-\X_{i-1}}_2\norm{\bS}_2=\sigma_i(\X).
\]
For rows, apply this to \( \X\T \), which has the same singular values.
:::

::: {#thm-cmp-rank-revealing}
[Column-pivoted QR]

Let \( \X \) be \( n\times p \) with \( n\ge p \), and let \( \X\boldsymbol{\Pi}=\Q\begin{psmallmatrix}\R\\\bzero\end{psmallmatrix} \) be
computed by Householder QR (@thm-cmp-householder) with column pivoting.

::: {.enumerate options="label=(\alph*)"}
1. For every \( k\le j \), \( r_{kk}^2\ge\sum_{i=k}^jr_{ij}^2 \). In particular
           \( \lvert r_{11}\rvert\ge\lvert r_{22}\rvert\ge\dots\ge\lvert r_{pp}\rvert \).

2. For every \( k<p \),
           \[
           \sigma_{k+1}(\X)\le\norm{\R_{22}}_2\le\sqrt{p-k}\,\lvert r_{k+1,k+1}\rvert
           \quad\text{and}\quad
           \sigma_{\min}(\R_{11})\le\sigma_k(\X).
           \]

3. In exact arithmetic, if \( \rank\X=r \), then \( r_{kk}\ne0 \) for \( k\le r \) and \( \R_{22}=\bzero \) for
           \( k=r \). The first \( r \) columns of \( \X\boldsymbol{\Pi} \) are a basis of \( \C(\X) \), and
           \( \X\boldsymbol{\Pi}=\Q\begin{psmallmatrix}\R_{11}&\R_{12}\\\bzero&\bzero\end{psmallmatrix} \) with \( \R_{11} \)
           nonsingular.

4. The cost is that of the unpivoted factorization, \( 2np^2-2p^3/3 \) flops, plus \( O(np) \) for keeping track of
           the column norms.
:::

:::

::: {.proof}
(a) Reflections preserve the length of each column segment they act on. When column \( k \) is chosen,
its segment in rows \( k,\dots,n \) has the largest norm among the trailing columns, and the reflection
maps it to \( (r_{kk},0,\dots,0)\T \), so \( \lvert r_{kk}\rvert \) equals that norm. Every column that ends in a
position \( j>k \) was among the trailing columns at step \( k \), so its segment in rows \( k,\dots,n \) had
norm at most \( \lvert r_{kk}\rvert \). The reflections at steps \( k,\dots,p \) act on rows \( k,\dots,n \) only, so
they preserve this norm, and at the end the segment is \( (r_{kj},\dots,r_{jj},0,\dots,0)\T \). Hence
\( \sum_{i=k}^jr_{ij}^2\le r_{kk}^2 \). With \( j=k+1 \) this gives \( r_{k+1,k+1}^2\le r_{kk}^2 \).

(b) \( \X\boldsymbol{\Pi} \) has the same singular values as \( \X \), because \( \boldsymbol{\Pi} \) is orthogonal. The matrix
\( \Q\begin{psmallmatrix}\R_{11}&\R_{12}\\\bzero&\bzero\\\bzero&\bzero\end{psmallmatrix} \) has rank at most \( k \) and differs from
\( \X\boldsymbol{\Pi} \) by \( \Q\begin{psmallmatrix}\bzero&\bzero\\\bzero&\R_{22}\\\bzero&\bzero\end{psmallmatrix} \), whose norm is
\( \norm{\R_{22}}_2 \). The Eckart–Young inequality (@prp-mat-svd-norms(c)) gives
\( \sigma_{k+1}(\X)\le\norm{\R_{22}}_2 \). Next,
\( \norm{\R_{22}}_2\le\norm{\R_{22}}_F\le\sqrt{p-k}\max_j\norm{\text{column }j\text{ of }\R_{22}} \), and by (a) with
\( k+1 \) in place of \( k \), each column of \( \R_{22} \) has norm at most \( \lvert r_{k+1,k+1}\rvert \). Finally, the first
\( k \) columns of \( \X\boldsymbol{\Pi} \) equal \( \Q_1\begin{psmallmatrix}\R_{11}\\\bzero\end{psmallmatrix} \), so their singular values are
those of \( \R_{11} \). Then @lem-cmp-column-subset gives \( \sigma_k(\R_{11})\le\sigma_k(\X) \).

(c) Because \( \R_{11} \) is triangular, \( \rank\R=k+\rank\R_{22} \) whenever \( r_{11},\dots,r_{kk} \) are nonzero (subtract
multiples of the first \( k \) columns from the others to clear \( \R_{12} \)). Also \( \rank\R=\rank\X=r \). So while
\( k<r \), \( \R_{22}\ne\bzero \), some trailing column is nonzero, and the next pivot \( r_{k+1,k+1} \) is nonzero. At
\( k=r \), \( \rank\R_{22}=0 \). The first \( r \) columns of \( \X\boldsymbol{\Pi} \) are \( \Q_1 \) times the independent columns
\( \begin{psmallmatrix}\R_{11}\\\bzero\end{psmallmatrix} \), so they are independent, and there are \( r=\dim\C(\X) \) of them.

(d) The norms of the trailing column segments are updated at each step by subtracting \( r_{kj}^2 \), at
\( O(p) \) flops per step. When cancellation makes the update inaccurate, the norm is recomputed.
:::

Part (b) is the rank-revealing property, and it is one-sided. A small pivot \( \lvert r_{k+1,k+1}\rvert \)
*certifies* a small singular value \( \sigma_{k+1} \): pivoted QR never reports a near dependence that is
not there. The converse can fail. A matrix can have a tiny singular value although no pivot is
small.

::: {#exm-cmp-kahan}
[Kahan's matrix]

Kahan (1966) gave the upper triangular matrix with diagonal \( 1,s,s^2,\dots,s^{m-1} \) and entries
\( -c\,s^{i-1} \) above the diagonal in row \( i \), where \( c^2+s^2=1 \). Its columns are already in pivoted
order (a tiny perturbation of the column scaling breaks ties in favour of that order), so pivoted
QR leaves it unchanged. With \( m=100 \) and \( c=0.2 \), the smallest diagonal entry is
\( 0.133 \), yet the smallest singular value is \( 3.7\times 10^{-9} \). No
diagonal entry hints at this ([Figure 10.7.1](#fig-cmp-pivoted)(a)). Such matrices are rare in
practice. Chan (1987) gave the first rank-revealing QR algorithm; factorizations whose bounds in both
directions are polynomial in \( p \) exist (Hong and Pan 1992) and can be computed efficiently (Gu and
Eisenstat 1996). Panel (b) shows the typical case: a \( 50\times10 \) matrix of
rank six plus noise of relative size \( 10^{-9} \). The pivoted diagonal drops by eight orders of magnitude
exactly where the singular values do.
:::

::: {when-format="html"}
![**Figure 10.7.1.** Diagonal of the pivoted triangle, \( \lvert r_{kk}\rvert \), and singular values \( \sigma_k \).
(a) Kahan's matrix: the diagonal decreases slowly while the last singular value collapses. (b) A matrix
of numerical rank six: both drop at \( k=7 \), well below the tolerance.](pivoted_qr.svg){#fig-cmp-pivoted width=100%}
:::

::: {when-format="pdf"}
![Diagonal of the pivoted triangle, \( \lvert r_{kk}\rvert \), and singular values \( \sigma_k \).
(a) Kahan's matrix: the diagonal decreases slowly while the last singular value collapses. (b) A matrix
of numerical rank six: both drop at \( k=7 \), well below the tolerance.](pivoted_qr.pdf){width=100%}
:::

```{.python .run #cell-pivoted-qr-pivoted}
import numpy as np
import statsmodels.api as sm

def pivoted_qr(X):
    """Householder QR with column pivoting: X[:, perm] = Q R, |r_11| >= |r_22| >= ..."""
    A = np.array(X, dtype=float)
    n, p = A.shape
    perm = np.arange(p)
    for k in range(min(n, p)):
        norms = np.sum(A[k:, k:] ** 2, axis=0)          # squared norms of the trailing columns
        j = k + int(np.argmax(norms))
        A[:, [k, j]], perm[[k, j]] = A[:, [j, k]], perm[[j, k]]
        x = A[k:, k]
        alpha = -np.copysign(np.linalg.norm(x), x[0])
        if alpha == 0.0:
            break                                        # the trailing block is exactly zero
        v = x.copy()
        v[0] -= alpha
        v /= np.linalg.norm(v)
        A[k:, k:] -= 2.0 * np.outer(v, v @ A[k:, k:])
    return np.triu(A[:p, :]), perm

def numerical_rank(R, tol):
    """Number of diagonal entries of the pivoted triangle above tol * |r_11|."""
    d = np.abs(np.diag(R))
    return int(np.sum(d > tol * d[0]))
```

```{.python .run #cell-pivoted-qr-kahan}
m, c = 100, 0.2
sn = np.sqrt(1 - c ** 2)
K = np.diag(sn ** np.arange(m)) @ (np.eye(m) - c * np.triu(np.ones((m, m)), 1))
K = K @ np.diag(1 - 1e-10 * np.arange(m))           # tiny tilt so that pivoting keeps the order
RK, permK = pivoted_qr(K)
sK = np.linalg.svd(K, compute_uv=False)
print("pivot order unchanged:", np.array_equal(permK, np.arange(m)))
print(f"smallest |r_kk| = {np.abs(RK[-1, -1]):.3f},  smallest singular value = {sK[-1]:.2e}")
```

## Basic solutions and aliased columns

Suppose pivoted QR finds numerical rank \( r<p \), so that
\( \X\boldsymbol{\Pi}\approx\Q\begin{psmallmatrix}\R_{11}&\R_{12}\\\bzero&\bzero\end{psmallmatrix} \). Write
\( \Q\T\y=(\mathbf{c}_1\T,\mathbf{c}_2\T)\T \) with \( \mathbf{c}_1\in\Real^r \), and split the permuted coefficients as
\( (\bb_1,\bb_2) \). As in @prp-cmp-qr-quantities, the least squares estimates are the solutions of
\[
\R_{11}\bb_1+\R_{12}\bb_2=\mathbf{c}_1 ,
\]
with \( \bb_2 \) free. Setting \( \bb_2=\bzero \) gives the **basic solution**, which simply drops the last
\( p-r \) pivoted columns. It has at most \( r \) nonzero coefficients. It differs from the minimum-norm solution
of @thm-cmp-svd-ls. That solution can also be computed at QR cost, by a second orthogonal reduction
from the right that turns \( [\R_{11},\R_{12}] \) into \( [\bT,\bzero] \) (a *complete orthogonal decomposition*;
Golub and Van Loan 2013). Every least squares solution gives the same fitted values. By
@thm-proj-invariant-functions, and in the language of [Chapter 8](../ch08-estimability/index.html) (@thm-est-characterization), every
estimable function \( \blambda\T\bbeta \) has the same estimate under every convention. The coefficients of aliased columns are not estimable.

Software differs in the choice. R's `lm` takes the columns in order, moves to the end any column whose
remaining norm falls below a relative tolerance (by default \( 10^{-7} \)), and reports its coefficient as
`NA`. statsmodels' `OLS` uses the pseudoinverse and so reports the minimum-norm solution. Programs
built on the sweep operator (@prp-cmp-sweep), such as SAS's GLM procedure, skip numerically zero pivots
and set those coefficients to zero, another basic solution.

::: {#exm-cmp-grunfeld}
[Firm effects in investment data]

Grunfeld's investment data (\( n=220 \): 11 US firms over 20 years, public domain in
statsmodels) are modelled with an intercept, the firm's market value and capital stock, and indicators
for all 11 firms: \( p=14 \) columns of rank 13. Pivoted QR with a tolerance of \( 10^{-10} \) finds rank 13
and leaves out the indicator for General Motors. Three conventions give three different
intercepts: \( -50.67 \) for the minimum-norm solution (as statsmodels reports),
\( -20.58 \) for the in-order rule that marks the last indicator as aliased (as R would),
and \( -70.30 \) for the basic solution from full pivoting. They agree to every printed digit on
everything estimable: the coefficient \( 0.1101 \) of market value, \( 0.3100 \) of
capital, the difference \( -172.20 \) between the General Motors and US Steel effects, and the
residual sum of squares \( 523719 \).
:::

```{.python .run #cell-pivoted-qr-aliased}
gr = sm.datasets.grunfeld.load_pandas().data
firms = list(dict.fromkeys(gr["firm"]))
D = np.column_stack([(gr["firm"] == f).to_numpy(float) for f in firms])   # 11 indicators
y = gr["invest"].to_numpy()
X = np.column_stack([np.ones(len(y)), gr[["value", "capital"]], D])      # 14 columns, rank 13
n, p = X.shape
R, perm = pivoted_qr(X)
r = numerical_rank(R, 1e-10)
print("p =", p, " numerical rank =", r, " column left out:", perm[r:])
```

```{.python .run #cell-pivoted-qr-conventions}
def basic_solution(X, y, tol=1e-10):
    """Drop the columns that pivoted QR puts last; fit the rest (zeros for dropped columns)."""
    R, perm = pivoted_qr(X)
    r = numerical_rank(R, tol)
    b = np.zeros(X.shape[1])
    b[perm[:r]] = np.linalg.lstsq(X[:, perm[:r]], y, rcond=None)[0]
    return b

def in_order_solution(X, y, tol=1e-7):
    """Keep a column only if it is not (nearly) a combination of the columns kept before it."""
    kept = []
    for j in range(X.shape[1]):
        Z = X[:, kept + [j]]
        Q, R = np.linalg.qr(Z)
        if abs(R[-1, -1]) > tol * np.linalg.norm(X[:, j]):
            kept.append(j)
    b = np.full(X.shape[1], np.nan)                  # aliased coefficients reported as missing
    b[kept] = np.linalg.lstsq(X[:, kept], y, rcond=None)[0]
    return b

solutions = {
    "minimum norm (pinv)": np.linalg.pinv(X) @ y,
    "pivoted QR, basic": basic_solution(X, y),
    "in order, aliased = NA": in_order_solution(X, y),
}
for name, b in solutions.items():
    fitted = X @ np.nan_to_num(b)
    print(f"{name:24s} value {b[1]:.5f}  capital {b[2]:.5f}  const {b[0]:9.3f}"
          f"  GM - US Steel {b[3] - b[4]:9.3f}  SSE {np.sum((y - fitted) ** 2):.1f}")
```

## Choosing a tolerance

Exact aliasing, as in the example, is easy: the offending pivot is at the level of rounding error, far
below any reasonable tolerance.
The hard case is a column that is nearly, but not exactly, a combination of others. The numerical rank then
depends on the tolerance, and so do the coefficients, dramatically.

::: {#exm-cmp-tolerance}
[A nearly aliased column]

Add to the Grunfeld regressors (intercept, value, capital) a fourth column equal to value plus capital
plus noise of relative size \( 10^{-9} \). Its smallest singular value, relative to the largest, is
\( 4.0\times 10^{-10} \). The default tolerance of `numpy.linalg.lstsq`,
\( \max(n,p)\,\epsilon_M \) relative to \( \sigma_1 \), treats the matrix as having full rank. It returns coefficients as
large as \( 8.1\times 10^{5} \), which spend one more degree of freedom to reduce the residual sum of
squares from \( 1768678.4 \) to \( 1768252.5 \), by \( 0.02 \) per cent. With the
tolerance set to \( 10^{-8} \), the fourth column is treated as aliased, and the largest coefficient is
\( 38.41 \).
:::

```{.python .run #cell-pivoted-qr-tolerance}
rng = np.random.default_rng(108)
z = X[:, 1] + X[:, 2] + 1e-9 * np.linalg.norm(X[:, 1]) * rng.normal(size=n) / np.sqrt(n)
Xz = np.column_stack([X[:, :3], z])              # constant, value, capital, value + capital + tiny
s = np.linalg.svd(Xz, compute_uv=False)
print("relative singular values:", np.array2string(s / s[0], precision=2))
for tol in [None, 1e-8]:
    b, _, rank, _ = np.linalg.lstsq(Xz, y, rcond=tol)
    print(f"rcond = {tol}: rank {rank}, largest |b_j| = {np.abs(b).max():.3g},"
          f" SSE = {np.sum((y - Xz @ b) ** 2):.1f}")
```

Neither answer is wrong as arithmetic. But a column that differs from a combination of the others
only at the level of \( 10^{-9} \) of its size carries no information a real measurement could support.
Three rules follow. Scale the columns to comparable lengths before deciding the rank (@prp-cmp-van-der-sluis); R's
relative test does this implicitly. Tie the tolerance to the accuracy of the data, not of the machine:
regressors recorded to six significant figures cannot resolve near dependences below about \( 10^{-6} \) of
the column lengths (@prp-cmp-weyl, applied to the recording error). And look for a gap, as in
[Figure 10.7.1](#fig-cmp-pivoted)(b). Without one the rank decision is fragile, the data are showing
collinearity ([Chapter 26](../ch26-collinearity/index.html)), and shrinkage is often a better response than a sharp cut-off ([Chapter 27](../ch27-shrinkage/index.html)).

::: {.idea}
Pivoted QR brings the most independent columns to the front. A small trailing pivot proves a near
dependence, and setting those coefficients aside changes nothing estimable. The tolerance that decides
"small" is a statement about the accuracy of the data.
:::

## Exercises

### A. Check your understanding

::: {#exr-cmp-pivot-hand}
[A1]

Apply column-pivoted QR (@thm-cmp-rank-revealing) by hand to the matrix with columns \( \x_1=(1,1,0,0)\T \), \( \x_2=(0,0,1,1)\T \) and
\( \x_3=(1,1,1,1)\T \) (use Gram–Schmidt on the pivoted columns, which gives the same \( \lvert r_{ij}\rvert \)).
Find the pivot order, \( \R \), and the rank. Give the basic solution for \( \y=(1,2,3,6)\T \).
:::

::: {.solution}
The norms are \( \sqrt2,\sqrt2,2 \), so \( \x_3 \) is chosen first, with \( \lvert r_{11}\rvert=2 \) and
\( \mathbf{q}_1=\x_3/2 \). Then \( \mathbf{q}_1\T\x_1=\mathbf{q}_1\T\x_2=1 \), and the remainders are
\( \x_1-\mathbf{q}_1=(\tfrac12,\tfrac12,-\tfrac12,-\tfrac12)\T \) and \( \x_2-\mathbf{q}_1 \), its negative. Both have norm \( 1 \),
so the tie goes to \( \x_1 \), with \( \lvert r_{22}\rvert=1 \). The remainder of \( \x_2 \) is then zero, so
\( r_{33}=0 \) and the rank is \( 2 \). The pivot order is \( (3,1,2) \), and up to signs
\( \R=\begin{psmallmatrix}2&1&1\\0&1&-1\\0&0&0\end{psmallmatrix} \). The basic solution regresses \( \y \) on
\( \x_3 \) and \( \x_1 \) with \( b_2=0 \). The fit is the pair of group means, \( 1.5 \) for the first two observations
and \( 4.5 \) for the last two. So \( b_3=4.5 \) and \( b_1=1.5-4.5=-3 \).
:::

### B. Practice

::: {#exr-cmp-basic-vs-minnorm}
[B1]

In the setting of this section, show that the basic solution and the minimum-norm solution have the same
fitted values. Show that the minimum-norm solution is orthogonal to \( \Null(\X) \), while the basic solution in
general is not. For the one-way layout of @exm-proj-ginverse-numeric, which generalized inverse in that
example gives the basic solution that drops the intercept?
:::

::: {.solution}
Both satisfy \( \R_{11}\bb_1+\R_{12}\bb_2=\mathbf{c}_1 \), so both are least squares estimates (@prp-cmp-qr-quantities),
and by @thm-proj-ls-projection both have fitted values \( \M\y \). The minimum-norm solution lies in
\( \C(\X\T)=\Null(\X)\perpc \) (@thm-cmp-svd-ls(b)). The basic solution has zeros in the positions of the dropped
columns. A null vector with a nonzero entry in one of those positions is generally not orthogonal to it.
For the one-way layout, dropping the intercept column and fitting the three indicators is the "drop
intercept" generalized inverse of @exm-proj-ginverse-numeric, which gives the group means as coefficients
and a zero intercept.
:::

::: {#exr-cmp-r11-bounds}
[B2]

Show that \( \lvert r_{11}\rvert=\max_j\norm{\x_j} \), and deduce \( \lvert r_{11}\rvert\le\sigma_1(\X)\le\sqrt p\,\lvert r_{11}\rvert \).
Combine with @thm-cmp-rank-revealing(b) to show that \( \kappa(\X)\ge\lvert r_{11}/r_{pp}\rvert \), so that a large ratio
of pivots proves ill-conditioning.
:::

### C. Going deeper

::: {#exr-cmp-kahan-null}
[C1]

For Kahan's matrix \( \mathbf{K}_m(c) \), find a vector \( \x \) with \( \norm{\x}=1 \) and \( \norm{\mathbf{K}_m\x} \) exponentially small in \( m \)
when \( c \) is fixed. (Solve \( \mathbf{K}_m\x=\mathbf{e}_m\,s^{m-1} \) by back substitution and examine the growth of the entries.) Why
does column pivoting not detect this?
:::
