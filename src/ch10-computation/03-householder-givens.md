# Householder and Givens QR

Orthogonal matrices preserve lengths (@prp-mat-orthogonal(a)). For any orthogonal \( \Q \),
\[
\norm{\y-\X\bb}=\norm{\Q\T\y-\Q\T\X\bb}\qquad\text{for every }\bb .
\]
So a least squares problem can be replaced by an equivalent one in which \( \Q\T\X \) is
upper triangular, without ever forming \( \X\T\X \). Such a \( \Q \) is built from reflections, which zero
a whole column below the diagonal at once, or from rotations, which zero one entry at a time.

## Householder reflections

A reflection \( \bH=\I-2\bv\bv\T/(\bv\T\bv) \) is symmetric and orthogonal, and it maps \( \bv \) to \( -\bv \)
and fixes every vector orthogonal to \( \bv \) (@exr-mat-householder). The whole numerical
content lies in choosing \( \bv \) so that a given vector is mapped onto a coordinate axis, and in
choosing its sign so that no cancellation occurs.

::: {#lem-cmp-reflector}
[Householder vector]

Let \( \x\in\Real^m \), \( \x\ne\bzero \), let \( \alpha=-\operatorname{sign}(x_1)\norm{\x} \) (with
\( \operatorname{sign}(0)=1 \)), and let \( \bv=\x-\alpha\bm e_1 \). Then \( \bv\ne\bzero \), and
\( \bH=\I-2\bv\bv\T/(\bv\T\bv) \) satisfies \( \bH\x=\alpha\bm e_1 \). The first entry
\( v_1=x_1-\alpha \) is a sum of two numbers of the same sign, so
\( \lvert v_1\rvert=\lvert x_1\rvert+\norm{\x} \). Applying \( \bH \) to a vector \( \bw \) as
\( \bw-2\bv(\bv\T\bw)/(\bv\T\bv) \) costs about \( 4m \) flops, and \( \bH \) is never formed.
:::

::: {.proof}
\( \alpha \) and \( x_1 \) have opposite signs (or \( x_1=0 \)), so \( v_1=x_1-\alpha \) adds two numbers of
the same sign and \( \lvert v_1\rvert=\lvert x_1\rvert+\norm{\x}>0 \). Next,
\( \bv\T\x=\norm{\x}^2-\alpha x_1 \) and \( \bv\T\bv=\norm{\x}^2-2\alpha x_1+\alpha^2=2(\norm{\x}^2-\alpha x_1) \),
using \( \alpha^2=\norm{\x}^2 \). Hence \( 2\bv\T\x/\bv\T\bv=1 \) and \( \bH\x=\x-\bv=\alpha\bm e_1 \).
:::

The opposite sign choice, \( \alpha=+\operatorname{sign}(x_1)\norm{\x} \), also gives a reflection
onto the axis. But then \( v_1=x_1-\alpha \) subtracts two nearly equal numbers whenever \( \x \) is
nearly parallel to \( \bm e_1 \). This is the cancellation of
[Section 10.1](01-why-not-invert.html), and it can make \( \bv \) inaccurate in every digit.

::: {#thm-cmp-householder}
[Householder QR]

Let \( \X \) be \( n\times p \) with \( n\ge p \).

::: {.enumerate options="label=(\alph*)"}
1. There are reflections \( \bH_1,\dots,\bH_p \), with \( \bH_j \) acting only on coordinates
           \( j,\dots,n \) (or equal to \( \I \)), such that
           \[
           \bH_p\cdots\bH_1\X=\begin{pmatrix}\R\\\bzero\end{pmatrix}
           \]
           with \( \R \) upper triangular of size \( p\times p \). Thus \( \X=\Q\begin{psmallmatrix}\R\\\bzero\end{psmallmatrix} \)
           with \( \Q=\bH_1\cdots\bH_p \) orthogonal, and \( \X=\Q_1\R \) where \( \Q_1 \) holds the first \( p \)
           columns of \( \Q \). If \( \rank\X=p \), then \( \R \) is nonsingular, and changing the signs of the
           rows of \( \R \) with negative diagonal entries gives the Cholesky factor of \( \X\T\X \).

2. The reduction takes \( 2np^2-2p^3/3 \) flops to leading order. Applying \( \Q\T \) (or \( \Q \)) to a
           vector takes \( 4np-2p^2 \), and forming \( \Q_1 \) explicitly takes another \( 2np^2-2p^3/3 \).

3. (Backward stability.) In floating-point arithmetic, the computed \( \hat{\R} \) is the exact
           triangular factor of \( \X+\Delta\X \), where each column satisfies
           \( \norm{\Delta\x_j}\le c_{np}u\norm{\x_j} \). The least squares solution
           \( \hat{\bb} \) computed from \( \hat{\R} \) and the computed \( \Q\T\y \) is the exact least squares
           solution for data \( (\X+\Delta\X,\ \y+\Delta\y) \) with perturbations of the same kind,
           \( \norm{\Delta\y}\le c_{np}u\norm{\y} \). Here \( c_{np} \) is a modest multiple of \( np \).
:::

:::

::: {.proof}
(a) Induction on the column. Suppose \( \bH_{j-1}\cdots\bH_1\X \) has zeros below the diagonal in its
first \( j-1 \) columns. Let \( \x \) be the part of column \( j \) in rows \( j,\dots,n \). If
\( \x=\bzero \), take \( \bH_j=\I \). Otherwise take the reflection of @lem-cmp-reflector for \( \x \),
acting on rows \( j,\dots,n \). It zeroes column \( j \) below the diagonal. It leaves rows
\( 1,\dots,j-1 \) alone. It maps the first \( j-1 \) columns, which are zero in rows \( j,\dots,n \),
to themselves. After \( p \) steps the matrix is \( \begin{psmallmatrix}\R\\\bzero\end{psmallmatrix} \), and \( \Q \) is a
product of orthogonal matrices. Orthogonal factors do not change rank, so
\( \rank\R=\rank\X \). If \( \bD \) is the diagonal matrix of signs of the diagonal of \( \R \), then
\( \X=(\Q_1\bD)(\bD\R) \) is a QR factorization with positive diagonal, and \( \bD\R \) is the Cholesky
factor by @thm-cmp-cholesky(d).

(b) Step \( j \) builds \( \bv \) (about \( 2(n-j) \) flops) and applies the reflection to the
\( p-j \) remaining columns of length \( n-j+1 \), at about \( 4(n-j)(p-j) \) flops. Summing,
\[
\sum_{j=1}^p4(n-j)(p-j)\approx4\int_0^p(n-t)(p-t)\,dt=2np^2-\tfrac23p^3 .
\]
Applying the \( p \) reflections to a vector costs \( \sum_j4(n-j+1)\approx4np-2p^2 \). Forming \( \Q_1 \)
means applying the reflections in reverse order to the first \( p \) columns of \( \I_n \).
Reflection \( j \) then acts on columns \( j,\dots,p \) only, which gives the same integral again.

(c) This is a theorem of rounding-error analysis, not of linear algebra, and we state it
without proof. Its essence is that a reflection applied in floating point is an exact
reflection applied to a slightly perturbed matrix, that orthogonal matrices do not amplify
perturbations, and that the \( p \) steps add their perturbations. See Higham (2002,
chapters 19 and 20) and Björck (1996, chapter 2).
:::

Part (c) is why QR is the default. Its backward error is small relative to each column separately,
which [Section 10.5](05-conditioning.html) uses to explain why column scaling does not change the
accuracy of QR. Software never forms \( \Q \); it stores the vectors \( \bv_j \) in the zeroed part of \( \X \).

```{.python .run #cell-householder-givens-householder}
import numpy as np
import statsmodels.api as sm

def house(x):
    """Unit vector v with (I - 2 v v') x = alpha e_1, alpha = -sign(x_1) ||x||."""
    alpha = -np.copysign(np.linalg.norm(x), x[0])
    v = x.astype(float).copy()
    v[0] -= alpha                        # x_1 - alpha adds two numbers of the same sign
    return v / np.linalg.norm(v), alpha

def householder_qr(X):
    """Return the reflection vectors and the p x p triangle R of X = Q [R; 0]."""
    A = np.array(X, dtype=float)
    n, p = A.shape
    vs = []
    for j in range(p):
        v, alpha = house(A[j:, j])
        A[j:, j:] -= 2.0 * np.outer(v, v @ A[j:, j:])    # apply I - 2vv' to the trailing block
        vs.append(v)
    return vs, np.triu(A[:p, :])

def apply_Qt(vs, y):
    """Q'y, applying the reflections H_1, ..., H_p in turn (Q itself is never formed)."""
    z = np.array(y, dtype=float)
    for j, v in enumerate(vs):
        z[j:] -= 2.0 * v * (v @ z[j:])
    return z

def apply_Q(vs, z):
    """Q z, applying the reflections in reverse order."""
    w = np.array(z, dtype=float)
    for j in range(len(vs) - 1, -1, -1):
        v = vs[j]
        w[j:] -= 2.0 * v * (v @ w[j:])
    return w
```

## Regression quantities from the factorization

[Section 6.10](../ch06-projections/10-computation.html) derived the regression quantities from the thin
factorization \( \X=\Q_1\R \) (@eq-proj-qr-ls). The full factorization, with \( \Q \) square, adds one
thing: the part of \( \Q\T\y \) that the model cannot reach, which *is* the residual in rotated
coordinates.

::: {#prp-cmp-qr-quantities}
[Least squares from QR]

Let \( \X=\Q\begin{psmallmatrix}\R\\\bzero\end{psmallmatrix} \) with \( \Q \) orthogonal and \( \rank\X=p \),
and partition \( \bm c=\Q\T\y \) as \( (\bm c_1\T,\bm c_2\T)\T \) with \( \bm c_1\in\Real^p \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta \) is the solution of the triangular system \( \R\bb=\bm c_1 \), and
           \( \text{SSE}=\norm{\bm c_2}^2 \);

2. \( \hY=\Q\begin{psmallmatrix}\bm c_1\\\bzero\end{psmallmatrix} \) and
           \( \he=\Q\begin{psmallmatrix}\bzero\\\bm c_2\end{psmallmatrix} \);

3. \( (\X\T\X)^{-1}=\R^{-1}\R^{-\top} \) and \( h_{ii}=\norm{\bm e_i\T\Q_1}^2 \);

4. \( c_k^2 \) is the drop in residual sum of squares when column \( k \) is added to columns
           \( 1,\dots,k-1 \), the sequential sum of squares of @def-ss-sequential.
:::

:::

::: {.proof}
(a) With \( \Q_1 \) the first \( p \) columns of \( \Q \), \( \X=\Q_1\R \) is a thin QR factorization
and \( \bm c_1=\Q_1\T\y \), so \( \R\hbeta=\bm c_1 \) and \( \text{SSE}=\norm{\y}^2-\norm{\bm c_1}^2 \)
are @eq-proj-qr-ls. What the full \( \Q \) adds is \( \norm{\bm c}=\norm{\y} \), hence
\( \text{SSE}=\norm{\bm c_2}^2 \).
(b) \( \X\hbeta=\Q\begin{psmallmatrix}\R\hbeta\\\bzero\end{psmallmatrix}=\Q\begin{psmallmatrix}\bm c_1\\\bzero\end{psmallmatrix} \),
and \( \y=\Q\bm c \).
(c) \( \X\T\X=\R\T\Q_1\T\Q_1\R=\R\T\R \), and the leverage formula is @eq-proj-qr-ls.
(d) The first \( k \) columns of \( \X \) equal the first \( k \) columns of \( \Q \) times the leading
\( k\times k \) block of \( \R \), which is nonsingular. So they span the same subspace as
\( \bm q_1,\dots,\bm q_k \), and the squared length of the projection of \( \y \) onto it is
\( \sum_{j\le k}(\bm q_j\T\y)^2=\sum_{j\le k}c_j^2 \) (@prp-proj-orthonormal-formula).
:::

Part (d) is the QR form of @prp-cmp-augmented(b): by uniqueness of the triangle, \( \bm c_1 \) equals
the vector \( \bz \) of that proposition up to signs. The vector \( \bm c_2 \) replaces the \( n \) correlated
residuals by \( n-p \) uncorrelated ones with the same sum of squares (@exr-cmp-uncorrelated-residuals).

::: {#exm-cmp-state-qr}
[Sequential sums of squares from one factorization]

For the state data of @exm-cmp-state-cholesky, Householder QR (@thm-cmp-householder) gives
\( \R \) with diagonal \( -7.071 \), \( 21.511 \), \( 21.609 \) and
\( 113.413 \). The first entry is \( -\sqrt{50} \): the sign convention of
@lem-cmp-reflector makes it negative. The squares of the first four entries of \( \Q\T\y \) are
the sequential sums of squares: \( 1018.81 \) for the intercept (which is
\( n\bar y^2 \)), \( 102.49 \) for poverty, \( 80.04 \) for single-parent households
given poverty, and \( 0.27 \) for urbanization given the other two. The remaining
\( 46 \) entries have squared length \( \text{SSE}=101.96 \). The
urbanization variable adds almost nothing once the other two are in the model, which is what
its small coefficient in @exm-proj-fwl-crime suggested. The script checks every number against
nested fits, and checks the leverages against the squared row lengths of \( \Q_1 \).
:::

```{.python .run #cell-householder-givens-statecrime}
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape

vs, R = householder_qr(X)
c = apply_Qt(vs, y)                      # c = Q'y = (c1, c2)
c1, c2 = c[:p], c[p:]
beta = np.linalg.solve(R, c1)            # triangular system R b = c1
sse = c2 @ c2                            # SSE = ||c2||^2
resid = apply_Q(vs, np.concatenate([np.zeros(p), c2]))
seq_ss = c1 ** 2                         # sequential sums of squares, in column order
print("R diagonal:", np.round(np.diag(R), 4))
print("coefficients:", np.round(beta, 4), f"  SSE = {sse:.4f}")
print("sequential SS (intercept, poverty, single, urban):", np.round(seq_ss, 3))
```

## Givens rotations

A reflection zeroes a whole column at once. Sometimes it is better to zero one entry at a time,
touching only two rows. This is the case when \( \X \) is sparse, when rows arrive one by one, or
when an existing factorization must be modified.

::: {#thm-cmp-givens}
[Givens rotations]

::: {.enumerate options="label=(\alph*)"}
1. For \( (a,b)\ne(0,0) \) let \( r=(a^2+b^2)^{1/2} \), \( c=a/r \) and \( s=b/r \). The rotation
           \( \bm G=\begin{psmallmatrix}c&s\\-s&c\end{psmallmatrix} \) is orthogonal and
           \( \bm G(a,b)\T=(r,0)\T \). Embedded in rows \( i \) and \( k \) of the identity, it changes only
           rows \( i \) and \( k \) of a matrix it multiplies, at \( 6 \) flops per column.

2. Zeroing the subdiagonal of \( \X \) column by column, each entry against the row above it,
           produces \( \bm G\X=\begin{psmallmatrix}\R\\\bzero\end{psmallmatrix} \) with \( \bm G \) orthogonal. It takes
           \( np-p(p+1)/2 \) rotations and \( 3np^2-p^3 \) flops to leading order.

3. (Row-wise QR.) Let \( \R \) be upper triangular and \( \x\in\Real^p \). There are \( p \) rotations, the
           \( k \)th acting on row \( k \) of \( \R \) and on the new row, that reduce
           \( \begin{psmallmatrix}\R\\\x\T\end{psmallmatrix} \) to \( \begin{psmallmatrix}\R_+\\\bzero\T\end{psmallmatrix} \)
           with \( \R_+ \) upper triangular and \( \R_+\T\R_+=\R\T\R+\x\x\T \). The cost is \( 3p^2 \) flops. If
           \( \R \) and \( \x \) are augmented by the response, as in @prp-cmp-augmented, the new row
           ends as \( (0,\dots,0,t) \), and the residual sum of squares increases by \( t^2 \).

4. A rotation of rows \( i \) and \( k \) creates a nonzero in column \( j \) only if row \( i \) or row \( k \)
           already had one there.
:::

:::

::: {.proof}
(a) \( c^2+s^2=1 \) gives \( \bm G\T\bm G=\I \). Also \( ca+sb=(a^2+b^2)/r=r \) and \( -sa+cb=0 \).
Each column of the product needs four multiplications and two additions.
(b) When entry \( (i,j) \) is zeroed against row \( i-1 \), working upwards from row \( n \), both rows are
already zero in columns \( 1,\dots,j-1 \), so no earlier zero is destroyed. Column \( j \) needs \( n-j \)
rotations, each acting on \( p-j+1 \) columns, so the count is
\( \sum_j6(n-j)(p-j+1)\approx6\int_0^p(n-t)(p-t)\,dt=3np^2-p^3 \).
(c) Before rotation \( k \), the new row is zero in columns \( 1,\dots,k-1 \), and so is row \( k \) of
\( \R \). The rotation chosen by (a) to zero the \( k \)th entry of the new row keeps both rows zero in
those columns. After \( p \) rotations the new row is zero, and row \( k \) of \( \R \) has been replaced by a
row that is zero before column \( k \). With \( \bm G \) the product of the rotations,
\[
\R_+\T\R_+=\begin{pmatrix}\R\\\x\T\end{pmatrix}\T\bm G\T\bm G\begin{pmatrix}\R\\\x\T\end{pmatrix}=\R\T\R+\x\x\T .
\]
Rotation \( k \) acts on \( p-k+1 \) columns, which gives \( \sum_k6(p-k+1)\approx3p^2 \). In the augmented
case, apply \( p \) rotations to the augmented \( (p+1)\times(p+1) \) triangle \( \bT \) and the row
\( (\x\T,y) \). The new row becomes \( (\bzero\T,t) \), and \( \bT_+\T\bT_+=\bT\T\bT+(\x\T,y)\T(\x\T,y)-(\bzero\T,t)\T(\bzero\T,t) \).
The corner entry of \( \bT_+ \) is still \( d \), so the augmented triangle with the new row included is
obtained by one more rotation, of \( d \) against \( t \). This gives a corner \( (d^2+t^2)^{1/2} \), and
the claim follows from @prp-cmp-augmented(a).
(d) Each new row is a combination of the two old ones.
:::

Part (c) turns QR into a streaming algorithm (Gentleman 1973): feed the rows of \( [\X,\y] \) one at a
time into \( \R=\bzero \). Memory is \( O(p^2) \), one pass suffices, and the result is as accurate as
Householder QR, so the one practical advantage of accumulating cross-products is kept without the
\( \kappa^2 \) loss. Part (d) keeps \( \R \) sparse for sparse designs, such as factors with many levels.

```{.python .run #cell-householder-givens-givens}
def givens(a, b):
    """c, s, r with [[c, s], [-s, c]] @ [a, b] = [r, 0]."""
    if b == 0.0:
        return 1.0, 0.0, a
    r = np.hypot(a, b)                   # sqrt(a^2 + b^2) without overflow
    return a / r, b / r, r
```

```{.python .run #cell-householder-givens-rowwise}
def add_row(R, z, x, y):
    """Rotate one observation (x, y) into the triangle R and the vector z = (Q'y)[:p].

    Returns the new R, z and the entry that leaves: its square is the increase in SSE.
    """
    R, z, x, y = R.copy(), z.copy(), np.array(x, dtype=float), float(y)
    for k in range(len(x)):
        c, s, r = givens(R[k, k], x[k])
        Rk, xk = R[k, k:].copy(), x[k:].copy()
        R[k, k:], x[k:] = c * Rk + s * xk, -s * Rk + c * xk
        z[k], y = c * z[k] + s * y, -s * z[k] + c * y
    return R, z, y

def rowwise_ls(X, y):
    """Least squares with one pass over the rows and O(p^2) memory."""
    p = X.shape[1]
    R, z, sse = np.zeros((p, p)), np.zeros(p), 0.0
    for xi, yi in zip(X, y):
        R, z, leftover = add_row(R, z, xi, yi)
        sse += leftover ** 2
    return R, z, sse

R_row, z_row, sse_row = rowwise_ls(X, y)
print("coefficients (row-wise Givens):", np.round(np.linalg.solve(R_row, z_row), 4))
print(f"SSE (row-wise Givens): {sse_row:.4f}")
```

The rotation in `givens` uses `hypot`, which computes \( (a^2+b^2)^{1/2} \) without overflow or
underflow. Square-root-free "fast" variants (Gentleman 1973) matter less on modern hardware.

## Gram–Schmidt revisited

Gram–Schmidt (@thm-mat-qr) also computes a QR factorization, at \( 2np^2 \) flops, and forms
\( \Q_1 \) explicitly, which gives leverages directly. The classical form computes \( r_{kj}=\bm q_k\T\x_j \) from the *original* column; the modified form
subtracts each new direction from all remaining columns at once. The two agree in exact arithmetic
and differ greatly in floating point.

::: {#exm-cmp-orthogonality}
[Loss of orthogonality]

For \( 100\times12 \) matrices with prescribed condition numbers from \( 10 \) to \( 10^{15} \),
[Figure 10.3.1](#fig-cmp-orthogonality) plots \( \norm{\hat{\Q}\T\hat{\Q}-\I}_2 \) for the computed
\( \hat{\Q} \). At \( \kappa=10^8 \) classical Gram–Schmidt gives
\( 3.5\times 10^{-2} \), and its columns are no longer an orthonormal basis in any useful
sense. Modified Gram–Schmidt gives \( 4.2\times 10^{-9} \), in line with the bound
\( O(\kappa u) \) of Björck (1967). The \( \Q \) accumulated from Householder reflections gives
\( 6.3\times 10^{-16} \), orthogonal to working precision whatever the condition number.
:::

::: {when-format="html"}
![**Figure 10.3.1.** Departure from orthogonality of the computed \( \Q \) against the condition
number of \( \X \), for classical and modified Gram–Schmidt and for Householder QR. The dotted
line is \( \kappa u \).](orthogonality_loss.svg){#fig-cmp-orthogonality width=62%}
:::

::: {when-format="pdf"}
![Departure from orthogonality of the computed \( \Q \) against the condition
number of \( \X \), for classical and modified Gram–Schmidt and for Householder QR. The dotted
line is \( \kappa u \).](orthogonality_loss.pdf){width=62%}
:::

```{.python .run #cell-householder-givens-gram-schmidt}
def cgs(X):
    """Classical Gram-Schmidt: inner products with the original column."""
    n, p = X.shape
    Q = np.zeros((n, p))
    for j in range(p):
        w = X[:, j] - Q[:, :j] @ (Q[:, :j].T @ X[:, j])
        Q[:, j] = w / np.linalg.norm(w)
    return Q

def mgs(X):
    """Modified Gram-Schmidt: subtract each direction from the working vector at once."""
    Q = np.array(X, dtype=float)
    n, p = Q.shape
    for j in range(p):
        Q[:, j] /= np.linalg.norm(Q[:, j])
        Q[:, j + 1:] -= np.outer(Q[:, j], Q[:, j] @ Q[:, j + 1:])
    return Q
```

Loss of orthogonality does not make modified Gram–Schmidt useless. It is numerically equivalent to
Householder QR applied to \( \X \) with \( p \) rows of zeros on top (Björck and Paige 1992), so applied to
\( [\X,\y] \), with \( \Q_1\T\y \) taken from the last column, it gives a backward stable least squares
solution. Leverages computed from the rows of \( \hat{\Q}_1 \) are not reliable when \( \kappa \) is large.

::: {.idea}
Orthogonal transformations reduce \( \X \) to a triangle while preserving every residual length.
Householder reflections do it fastest; Givens rotations do it one entry at a time, which suits
streaming, updating and sparse data.
:::

## Exercises

### A. Check your understanding

::: {#exr-cmp-house-hand}
[A1]

Find the Householder vector and \( \alpha \) of @lem-cmp-reflector for \( \x=(3,4)\T \), and verify
\( \bH\x=\alpha\bm e_1 \). For \( \x=(1,10^{-9})\T \), compute \( v_1 \) with both sign choices in
double precision, and explain what goes wrong with the unstable one.
:::

::: {.solution}
\( \norm{\x}=5 \), \( \alpha=-5 \), \( \bv=(8,4)\T \), \( \bv\T\bv=80 \) and \( \bv\T\x=40 \), so
\( \bH\x=\x-2\cdot\tfrac{40}{80}\bv=(3,4)-(8,4)=(-5,0) \). For \( \x=(1,10^{-9})\T \),
\( \norm{\x}=(1+10^{-18})^{1/2} \) rounds to \( 1 \). The stable choice gives \( v_1=1+1=2 \). The
unstable one gives \( v_1=1-1=0 \) in floating point, while the true value is
\( 1-\norm{\x}\approx-5\times10^{-19} \). With \( \bv=(0,10^{-9}) \) the "reflection" swaps the
sign of the second coordinate and maps \( \x \) to \( (1,-10^{-9}) \), not onto the axis. An
algebraically equivalent formula, \( v_1=-(x_2^2+\dots+x_m^2)/(x_1+\norm{\x}) \), rescues that
choice.
:::

### B. Practice

::: {#exr-cmp-uncorrelated-residuals}
[B1]

In the setting of @prp-cmp-qr-quantities, let \( \Q=[\Q_1,\Q_2] \). Show that under the linear model
with \( \Cov(\be)=\sigma^2\I \), the vector \( \bm c_2=\Q_2\T\y \) has mean \( \bzero \) and covariance
\( \sigma^2\I_{n-p} \), and that \( \he=\Q_2\bm c_2 \). Under normality, deduce from \( \bm c_2 \) alone that
\( \text{SSE}/\sigma^2\sim\chi^2(n-p) \) independently of \( \hbeta \) (compare @thm-opt-sampling). Is
\( \bm c_2 \) unique?
:::

::: {.solution}
\( \Q_2\T\X=\bzero \), because \( \Q\T\X=\begin{psmallmatrix}\R\\\bzero\end{psmallmatrix} \). So
\( \bm c_2=\Q_2\T\X\bbeta+\Q_2\T\be=\Q_2\T\be \), with mean \( \bzero \) and covariance
\( \sigma^2\Q_2\T\Q_2=\sigma^2\I \). By @prp-cmp-qr-quantities(b),
\( \he=\Q\begin{psmallmatrix}\bzero\\\bm c_2\end{psmallmatrix}=\Q_2\bm c_2 \). Under normality,
\( \Q\T\be\sim\Normal_n(\bzero,\sigma^2\I) \) (@thm-mvn-linear), so \( \bm c_1-\R\bbeta=\Q_1\T\be \) and
\( \bm c_2=\Q_2\T\be \) are independent normal vectors with identity covariance times \( \sigma^2 \).
Then \( \text{SSE}/\sigma^2=\norm{\bm c_2/\sigma}^2\sim\chi^2(n-p) \), and
\( \hbeta=\R^{-1}\bm c_1 \) is a function of \( \bm c_1 \) alone. \( \bm c_2 \) is not unique: \( \Q_2 \) can be
replaced by \( \Q_2\bm O \) for any orthogonal \( \bm O \), which changes \( \bm c_2 \) but not its length or its
distribution.
:::

::: {#exr-cmp-wls-condition}
[B2]

Weighted least squares with a diagonal weight matrix \( \W \) of positive weights is ordinary least
squares for \( \W^{1/2}\X \) and \( \W^{1/2}\y \). Show that
\( \kappa(\W^{1/2}\X)\le\kappa(\W^{1/2})\,\kappa(\X) \). Give an example with two very different
weights in which \( \kappa(\W^{1/2}\X) \) is large although \( \kappa(\X) \) is small, and one in which
it is small although \( \kappa(\W^{1/2}) \) is large.
:::

::: {.solution}
\( \sigma_{\max}(\W^{1/2}\X)\le\sigma_{\max}(\W^{1/2})\sigma_{\max}(\X) \) and
\( \norm{\W^{1/2}\X\bv}\ge\sigma_{\min}(\W^{1/2})\norm{\X\bv}\ge\sigma_{\min}(\W^{1/2})\sigma_{\min}(\X) \)
for unit \( \bv \), by @prp-mat-svd-norms. Divide. For the first example take
\( \X=\begin{psmallmatrix}1&1\\1&-1\end{psmallmatrix} \) (with \( \kappa=1 \)) and \( \W=\diag(1,w^2) \) with
\( w \) large. Then \( \W^{1/2}\X \) has rows \( (1,1) \) and \( (w,-w) \), whose singular values are
\( \sqrt2 \) and \( \sqrt2w \), so \( \kappa=w \). For the second take \( \X=\begin{psmallmatrix}1&0\\0&1/w\end{psmallmatrix} \)
with the same \( \W \). Then \( \W^{1/2}\X=\I \).
:::

### C. Going deeper

::: {#exr-cmp-mgs-householder}
[C1]

Show that, in exact arithmetic, Householder QR (@thm-cmp-householder) applied to
\( \begin{psmallmatrix}\bzero_{p\times p}\\\X\end{psmallmatrix} \) produces the same \( \R \) as modified
Gram–Schmidt applied to \( \X \) up to the signs of its rows, and that with the sign convention of
@lem-cmp-reflector its reflection vectors are
\( \bv_j=\begin{psmallmatrix}\bm e_j\\\bm q_j\end{psmallmatrix}/\sqrt2 \). This is the
observation of Björck and Paige (1992), who take \( \alpha=+\norm{\x} \) and so obtain
\( \begin{psmallmatrix}-\bm e_j\\\bm q_j\end{psmallmatrix}/\sqrt2 \).
:::
