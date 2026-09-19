# Updating and downdating

A fitted regression often has to be changed slightly: new observations arrive, a suspicious case
is set aside, a window slides along a time series, a variable is added or removed. Refitting from scratch costs \( O(np^2) \) each time. Modifying the factorization costs
\( O(p^2) \) for a change of one row, and \( O(np) \) for a change of one column. This section develops
those modifications and uses them to compute every leave-one-out quantity from a single fit.

## Adding and deleting observations

Adding a row is part (c) of @thm-cmp-givens: \( p \) rotations fold the new row into the triangle.
Removing a row is harder, because orthogonal transformations can only combine rows that are
present. The trick is to run the addition backwards.

::: {#thm-cmp-updating}
[Updating and downdating a triangular factor]

Let \( \bT \) be an \( m\times m \) nonsingular upper triangular matrix with \( \bT\T\bT=\A \), and let
\( \bz\in\Real^m \).

::: {.enumerate options="label=(\alph*)"}
1. (*Adding a row.*) There are \( m \) rotations, the \( k \)th acting on row \( k \) of \( \bT \) and on
           an extra row, that turn \( \begin{psmallmatrix}\bT\\\bz\T\end{psmallmatrix} \) into
           \( \begin{psmallmatrix}\bT_+\\\bzero\T\end{psmallmatrix} \), with \( \bT_+ \) upper triangular and
           \( \bT_+\T\bT_+=\A+\bz\bz\T \).

2. (*Deleting a row.*) Let \( \mathbf{a} \) solve \( \bT\T\mathbf{a}=\bz \). The matrix \( \A-\bz\bz\T \) is positive
           definite iff \( \norm{\mathbf{a}}<1 \). In that case put \( \alpha=(1-\norm{\mathbf{a}}^2)^{1/2} \). Choose
           rotations \( \mathbf{G}_m,\dots,\mathbf{G}_1 \), with \( \mathbf{G}_k \) acting on coordinates \( k \) and \( m+1 \), such
           that \( \mathbf{G}_1\cdots\mathbf{G}_m(\mathbf{a}\T,\alpha)\T=(\bzero\T,1)\T \). Then
           \[
           \mathbf{G}_1\cdots\mathbf{G}_m\begin{pmatrix}\bT\\\bzero\T\end{pmatrix}=\begin{pmatrix}\bT_-\\\bz\T\end{pmatrix}
           \]
           with \( \bT_- \) upper triangular and \( \bT_-\T\bT_-=\A-\bz\bz\T \).

3. Each operation costs \( O(m^2) \) flops: \( 3m^2 \) for (a), and \( 3m^2 \) plus a triangular solve
           of \( m^2 \) for (b).
:::

:::

::: {.proof}
(a) is @thm-cmp-givens(c).

(b) \( \A-\bz\bz\T=\bT\T(\I-\mathbf{a}\mathbf{a}\T)\bT \). The matrix \( \I-\mathbf{a}\mathbf{a}\T \) has eigenvalue
\( 1-\norm{\mathbf{a}}^2 \) on \( \mathbf{a} \) and \( 1 \) on \( \mathbf{a}\perpc \). So it is positive definite iff \( \norm{\mathbf{a}}<1 \), and then
so is \( \A-\bz\bz\T \) (@prp-mat-pd-properties(c)). The rotations exist: \( \mathbf{G}_m \) combines \( a_m \) with
\( \alpha \) and leaves \( 0 \) in position \( m \) and \( (a_m^2+\alpha^2)^{1/2} \) in position \( m+1 \). Then \( \mathbf{G}_{m-1} \)
does the same with \( a_{m-1} \), and so on (@thm-cmp-givens(a)). The final vector has norm
\( (\norm{\mathbf{a}}^2+\alpha^2)^{1/2}=1 \) and is zero except in the last position. That entry is positive,
so it is \( 1 \).

Apply the same rotations, in the same order, to \( \begin{psmallmatrix}\bT\\\bzero\T\end{psmallmatrix} \).
Before \( \mathbf{G}_k \) is applied, the extra row is zero in columns \( 1,\dots,k \), because it has only been
mixed with rows \( k+1,\dots,m \) of \( \bT \), which are zero there. Row \( k \) of \( \bT \) is zero in columns
\( 1,\dots,k-1 \). So \( \mathbf{G}_k \) produces a new row \( k \) that is zero before column \( k \), and the result is
\( \begin{psmallmatrix}\bT_-\\\bw\T\end{psmallmatrix} \) with \( \bT_- \) upper triangular. Write \( \mathbf{G}=\mathbf{G}_1\cdots\mathbf{G}_m \),
which is orthogonal. Then
\[
\bz=\bT\T\mathbf{a}=\begin{pmatrix}\bT\\\bzero\T\end{pmatrix}\T\begin{pmatrix}\mathbf{a}\\\alpha\end{pmatrix}
=\Bigl(\mathbf{G}\begin{pmatrix}\bT\\\bzero\T\end{pmatrix}\Bigr)\T\mathbf{G}\begin{pmatrix}\mathbf{a}\\\alpha\end{pmatrix}
=\begin{pmatrix}\bT_-\\\bw\T\end{pmatrix}\T\begin{pmatrix}\bzero\\1\end{pmatrix}=\bw ,
\]
and in the same way \( \bT\T\bT=\bT_-\T\bT_-+\bw\bw\T \). Hence \( \bT_-\T\bT_-=\A-\bz\bz\T \).

(c) The rotations act on rows of length at most \( m \), and there are \( m \) of them.
:::

If \( \A=\X\T\X \) and \( \bz=\x_{(i)} \), then \( \norm{\mathbf{a}}^2=\bz\T\A^{-1}\bz \) is the leverage \( h_{ii} \)
(@def-proj-leverage): a case can be deleted exactly when the others still determine every coefficient.
Downdating the augmented triangle of \( [\X,\y] \) (@prp-cmp-augmented) gives the fit without case \( i \), with
the square root of the new SSE in the corner. Downdating is less stable than updating: \( 1-\norm{\mathbf{a}}^2 \)
cancels when \( h_{ii}\approx1 \), and the relative error can reach about \( u/(1-h_{ii}) \) (Stewart 1979).

```{.python .run #cell-updating-rotations}
import numpy as np
import statsmodels.api as sm

def givens(a, b):
    """c, s, r with [[c, s], [-s, c]] @ [a, b] = [r, 0]."""
    if b == 0.0:
        return 1.0, 0.0, a
    r = np.hypot(a, b)
    return a / r, b / r, r

def update(T, z):
    """Triangle of [A; z'] from the triangle T of A: rotate the new row into T."""
    T, z = T.copy(), np.array(z, dtype=float)
    for k in range(len(z)):
        c, s, r = givens(T[k, k], z[k])
        Tk, zk = T[k, k:].copy(), z[k:].copy()
        T[k, k:], z[k:] = c * Tk + s * zk, -s * Tk + c * zk
    return T

def downdate(T, z):
    """Triangle of A with the row z' removed, from the triangle T of A (T'T = A'A)."""
    m = len(z)
    a = np.linalg.solve(T.T, z)                   # T'a = z; ||a||^2 is the leverage of z
    alpha2 = 1.0 - a @ a
    if alpha2 <= 0:
        raise np.linalg.LinAlgError("the row cannot be removed (leverage >= 1)")
    alpha = np.sqrt(alpha2)
    T, w = T.copy(), np.zeros(m)                  # w: the extra row, initially zero
    for k in range(m - 1, -1, -1):                # rotate (a_k, alpha) into (0, new alpha)
        c, s, r = givens(alpha, a[k])
        alpha = r
        Tk = T[k, k:].copy()
        T[k, k:], w[k:] = c * Tk - s * w[k:], s * Tk + c * w[k:]
    return T                                      # on exit w equals z
```

The alternative, **recursive least squares** (Plackett 1950), updates \( \mathbf{K}=(\X\T\X)^{-1} \) by the
Sherman–Morrison formula (@thm-mat-woodbury, @exm-mat-deletion): adding \( (\x,y) \) gives
\[
\mathbf{K}_+=\mathbf{K}-\frac{\mathbf{K}\x\x\T\mathbf{K}}{1+\x\T\mathbf{K}\x},\qquad
\hbeta_+=\hbeta+\frac{\mathbf{K}\x\,(y-\x\T\hbeta)}{1+\x\T\mathbf{K}\x}.
\]
It also costs \( O(p^2) \), but it carries an explicit inverse of the cross-product matrix, with its \( \kappa^2 \)
sensitivity, and nothing corrects its accumulated errors (@exm-cmp-moving-window).

## Leave-one-out quantities from one fit

Many diagnostics ask what would happen if observation \( i \) were left out. Computing them by \( n \)
separate fits costs \( O(n^2p^2) \), but all of them follow from one fit.

::: {#prp-cmp-loo}
[Deletion formulas]

Let \( \X \) have full column rank, let \( h_{ii}<1 \), and let \( \hbeta_{(i)} \) and \( \text{SSE}_{(i)} \) denote the
least squares estimate and residual sum of squares computed without observation \( i \). Then

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta-\hbeta_{(i)}=(\X\T\X)^{-1}\x_{(i)}\,\hat{\varepsilon}_i/(1-h_{ii}) \);

2. \( y_i-\x_{(i)}\T\hbeta_{(i)}=\hat{\varepsilon}_i/(1-h_{ii}) \);

3. \( \text{SSE}_{(i)}=\text{SSE}-\hat{\varepsilon}_i^2/(1-h_{ii}) \).
:::

:::

::: {.proof}
(a) is @exr-mat-leave-one-out, which follows from the deletion formula of @exm-mat-deletion.
(b) Multiply (a) by \( \x_{(i)}\T \): \( \x_{(i)}\T\hbeta_{(i)}=\x_{(i)}\T\hbeta-h_{ii}\hat{\varepsilon}_i/(1-h_{ii}) \), so
\( y_i-\x_{(i)}\T\hbeta_{(i)}=\hat{\varepsilon}_i+h_{ii}\hat{\varepsilon}_i/(1-h_{ii})=\hat{\varepsilon}_i/(1-h_{ii}) \).
(c) Put \( g=\hat{\varepsilon}_i/(1-h_{ii}) \) and let \( h_{ji}=\x_{(j)}\T(\X\T\X)^{-1}\x_{(i)} \) be the entries of
\( \M \). By (a), \( y_j-\x_{(j)}\T\hbeta_{(i)}=\hat{\varepsilon}_j+h_{ji}\,g \) for every \( j \), including \( j=i \), where
it gives \( g \) by (b). Summing squares over all \( j \),
\[
\sum_j(\hat{\varepsilon}_j+h_{ji}g)^2=\text{SSE}+2g\,(\M\he)_i+g^2(\M^2)_{ii}=\text{SSE}+g^2h_{ii},
\]
because \( \M\he=\bzero \) and \( \M^2=\M \). Removing the \( j=i \) term, which is \( g^2 \), leaves
\( \text{SSE}_{(i)}=\text{SSE}-g^2(1-h_{ii})=\text{SSE}-\hat{\varepsilon}_i^2/(1-h_{ii}) \).
:::

With a QR factorization the leverages cost \( O(np) \) (@prp-cmp-qr-quantities(c)), and the rest \( O(n) \).
This gives the PRESS residuals \( \hat{\varepsilon}_i/(1-h_{ii}) \), whose sum of squares is leave-one-out
cross-validation (Allen 1974; Chapter 29), and the deleted variances \( s_{(i)}^2=\text{SSE}_{(i)}/(n-p-1) \)
behind the externally studentized residuals \( t_i=\hat{\varepsilon}_i/(s_{(i)}\sqrt{1-h_{ii}}) \) of Chapter 20.

::: {#exm-cmp-dc}
[Leaving out the District of Columbia]

Refit the state murder-rate regression with all \( 51 \) jurisdictions, including the
District of Columbia, which [Chapter 6](../ch06-projections/index.html) set aside. The residual sum of squares is
\( 142.53 \) and the PRESS statistic is \( 230.52 \). The largest
PRESS residual belongs to the District of Columbia. Its ordinary residual is
\( 4.475 \) and its leverage \( 0.506 \), so its prediction error when it is
left out is \( \hat{\varepsilon}_i/(1-h_{ii})=9.066 \). Its externally
studentized residual is \( 4.28 \). By @prp-cmp-loo(c), deleting it reduces the residual
sum of squares to \( 101.96 \), the value found for the 50 states
in @exm-cmp-state-cholesky. The script confirms all \( 51 \) deletions three ways: by the
formulas, by downdating the augmented triangle, and by refitting.
:::

```{.python .run #cell-updating-loo}
data = sm.datasets.statecrime.load_pandas().data           # all 50 states and DC
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape

Q, R = np.linalg.qr(X)
beta = np.linalg.solve(R, Q.T @ y)
e = y - X @ beta
sse = e @ e
h = np.sum(Q ** 2, axis=1)                          # leverages
press_resid = e / (1 - h)                           # y_i minus its prediction without case i
sse_loo = sse - e ** 2 / (1 - h)                    # SSE after deleting case i
s_loo = np.sqrt(sse_loo / (n - 1 - p))
t_ext = e / (s_loo * np.sqrt(1 - h))                # externally studentized residuals
i = int(np.argmax(np.abs(press_resid)))
print(f"largest PRESS residual: {data.index[i]}, e = {e[i]:.3f}, h = {h[i]:.3f},"
      f" e/(1-h) = {press_resid[i]:.3f}, t = {t_ext[i]:.2f}")
print(f"PRESS = {np.sum(press_resid ** 2):.2f}, SSE = {sse:.2f}")
```

## Adding and deleting variables

Changing the columns of \( \X \) is a Gram–Schmidt step, or its reverse.

::: {#prp-cmp-add-column}
[Adding a column]

Let \( \X=\Q_1\R \) be a thin QR factorization of a full-column-rank \( \X \), and let \( \bz\notin\C(\X) \).
Put \( \br=\Q_1\T\bz \), \( \bw=\bz-\Q_1\br \), \( \rho=\norm{\bw} \) and \( \mathbf{q}=\bw/\rho \). Then
\[
[\X,\bz]=[\Q_1,\mathbf{q}]\begin{pmatrix}\R&\br\\\bzero\T&\rho\end{pmatrix}
\]
is a thin QR factorization. The residual sum of squares falls by \( (\mathbf{q}\T\y)^2 \) when \( \bz \) is added,
and the coefficient of \( \bz \) in the enlarged model is \( \mathbf{q}\T\y/\rho=\bw\T\y/\bw\T\bw \).
:::

::: {.proof}
\( \bw=(\I-\M)\bz\ne\bzero \) because \( \bz\notin\C(\X) \), and \( \mathbf{q} \) is a unit vector orthogonal to \( \C(\Q_1) \).
The product of the factors is \( [\Q_1\R,\ \Q_1\br+\rho\mathbf{q}]=[\X,\bz] \). By @prp-cmp-qr-quantities(d), applied
to the enlarged factorization, the new column's sequential sum of squares is the square of its
coordinate \( \mathbf{q}\T\y \). The last row of the triangular system of @prp-cmp-qr-quantities(a) reads
\( \rho b_z=\mathbf{q}\T\y \). This is also the Frisch–Waugh–Lovell coefficient \( \bw\T\y/\bw\T\bw \) of @thm-proj-fwl.
:::

If \( \bz \) is nearly in \( \C(\X) \), the computed \( \bw \) loses orthogonality to \( \Q_1 \); repeating the step once
restores it (Björck 1996, chapter 2). Deleting column \( j \) of \( \R \) leaves an upper Hessenberg block that
\( p-j \) rotations restore to triangular form, at \( O(p^2) \) cost (@exr-cmp-delete-column). So stepwise
procedures move between models at \( O(p^2) \) or \( O(np) \) cost per step.

```{.python .run #cell-updating-add-column}
def add_column(Q1, R, x):
    """Thin QR of [X, x] from the thin QR of X (one Gram-Schmidt step, repeated once)."""
    r = Q1.T @ x
    w = x - Q1 @ r
    r2 = Q1.T @ w                                 # reorthogonalize: guards against cancellation
    w, r = w - Q1 @ r2, r + r2
    rho = np.linalg.norm(w)
    Rn = np.block([[R, r[:, None]], [np.zeros((1, R.shape[1])), np.array([[rho]])]])
    return np.column_stack([Q1, w / rho]), Rn
```

## A moving window

::: {#exm-cmp-moving-window}
[Ten-year regressions of consumption]

Quarterly US real consumption is regressed on an intercept, real disposable income, real
investment, population and unemployment (the public-domain macroeconomic
series of @exm-mat-macro-svd), in a window of 40 quarters that slides forward one quarter at a time.
Each step adds one row and deletes one. Over the \( 163 \) steps, whose windows have
condition numbers up to \( 2.9\times 10^{6} \), the fit is maintained in two ways. The
first applies two Sherman–Morrison updates to \( (\X\T\X)^{-1} \). The second updates and then downdates
the augmented triangle. After each step both are compared with a fresh QR fit of the window
([Figure 10.6.1](#fig-cmp-moving-window)). Inverse updating reaches a relative error of
\( 6.2\times 10^{-7} \) and ends at \( 8.4\times 10^{-8} \). Updating the triangle stays
below \( 1.1\times 10^{-10} \) and ends at \( 3.6\times 10^{-11} \). Both are cheap. Only one keeps QR's accuracy.
:::

::: {when-format="html"}
![**Figure 10.6.1.** A moving-window regression maintained by rank-one updates of
\( (\X\T\X)^{-1} \) and by Givens updating and downdating of the augmented triangle. The vertical axis is
the relative difference in the coefficients from a fresh QR fit of each window.](moving_window.svg){#fig-cmp-moving-window width=62%}
:::

::: {when-format="pdf"}
![A moving-window regression maintained by rank-one updates of
\( (\X\T\X)^{-1} \) and by Givens updating and downdating of the augmented triangle. The vertical axis is
the relative difference in the coefficients from a fresh QR fit of each window.](moving_window.pdf){width=62%}
:::

```{.python .run #cell-updating-window-setup}
macro = sm.datasets.macrodata.load_pandas().data
yw = macro["realcons"].to_numpy()
Xw = np.column_stack([np.ones(len(yw)), macro[["realdpi", "realinv", "pop", "unemp"]]])
w = 40                                            # ten years of quarterly data
print("kappa of the full design:", f"{np.linalg.cond(Xw):.2e}")
```

```{.python .run #cell-updating-window}
K = np.linalg.inv(Xw[:w].T @ Xw[:w])              # (a) rank-one updates of the inverse
g = Xw[:w].T @ yw[:w]
Zw = np.column_stack([Xw, yw])
Tw = np.linalg.qr(Zw[:w], mode="r")               # (b) triangle of [X, y]
p1 = Xw.shape[1]
err_sm, err_giv = [], []
for t in range(w, len(yw)):
    new, old = Xw[t], Xw[t - w]
    Kx = K @ new                                  # add the new quarter (Sherman-Morrison)
    K = K - np.outer(Kx, Kx) / (1 + new @ Kx)
    Kx = K @ old                                  # delete the oldest quarter
    K = K + np.outer(Kx, Kx) / (1 - old @ Kx)
    g = g + new * yw[t] - old * yw[t - w]
    b_sm = K @ g
    Tw = downdate(update(Tw, Zw[t]), Zw[t - w])
    b_giv = np.linalg.solve(Tw[:p1, :p1], Tw[:p1, p1])
    b_ref = np.linalg.lstsq(Xw[t - w + 1:t + 1], yw[t - w + 1:t + 1], rcond=None)[0]   # (c) refit
    err_sm.append(np.linalg.norm(b_sm - b_ref) / np.linalg.norm(b_ref))
    err_giv.append(np.linalg.norm(b_giv - b_ref) / np.linalg.norm(b_ref))
print(f"final relative error: inverse updating {err_sm[-1]:.1e}, Givens {err_giv[-1]:.1e}")
```

::: {.idea}
A change of one row or column costs \( O(p^2) \) or \( O(np) \) on the triangular factor, not \( O(np^2) \).
All \( n \) leave-one-out fits cost one, and the leverage says whether, and how accurately, a case can
be removed.
:::

## Exercises

### A. Check your understanding

::: {#exr-cmp-press-balanced}
[A1]

Show that \( \text{PRESS}=\sum_i\hat{\varepsilon}_i^2/(1-h_{ii})^2\ge\text{SSE} \). Show that if all leverages equal
\( p/n \), then \( \text{PRESS}=\text{SSE}\,n^2/(n-p)^2 \). Evaluate this for the state data with \( n=51 \) and
\( p=4 \), and compare with the actual ratio in @exm-cmp-dc.
:::

### B. Practice

::: {#exr-cmp-augmented-leverage}
[B1]

Show that the leverage of row \( i \) of the augmented matrix \( \Z=[\X,\y] \) is
\( h_{ii}+\hat{\varepsilon}_i^2/\text{SSE} \). Deduce the condition under which the augmented triangle can be
downdated, and show that it fails only if \( h_{ii}=1 \) or deleting case \( i \) leaves a model that fits the
remaining data exactly.
:::

::: {.solution}
\( \C(\Z)=\C(\X)\dirsum\spn(\he) \), with the two summands orthogonal, so by @thm-proj-sum the projection onto
\( \C(\Z) \) is \( \M+\he\he\T/\norm{\he}^2 \). Its \( i \)th diagonal entry is \( h_{ii}+\hat{\varepsilon}_i^2/\text{SSE} \). Downdating
requires this to be below \( 1 \). By @prp-cmp-loo(c),
\( 1-h_{ii}-\hat{\varepsilon}_i^2/\text{SSE}=(1-h_{ii})\,\text{SSE}_{(i)}/\text{SSE} \), which vanishes iff \( h_{ii}=1 \) or
\( \text{SSE}_{(i)}=0 \).
:::

::: {#exr-cmp-delete-column}
[B2]

Let \( \X=\Q_1\R \) and let \( \X_{-j} \) be \( \X \) without column \( j \). Show that \( \X_{-j}=\Q_1\mathbf{H} \), where \( \mathbf{H} \) is \( \R \)
without column \( j \), which is upper Hessenberg in its last \( p-j \) columns. Show that \( p-j \) rotations of
adjacent rows reduce \( \mathbf{H} \) to \( \begin{psmallmatrix}\R_{-j}\\\bzero\T\end{psmallmatrix} \), and that \( \R_{-j} \) is the
triangle of \( \X_{-j} \). How much does the residual sum of squares rise?
:::

::: {.solution}
Columns of \( \X \) are \( \Q_1 \) times columns of \( \R \), so deleting a column of both gives \( \X_{-j}=\Q_1\mathbf{H} \).
Column \( k\ge j \) of \( \mathbf{H} \) is column \( k+1 \) of \( \R \), which has nonzeros down to row \( k+1 \): one below the
diagonal. Rotate rows \( j,j+1 \) to zero entry \( (j+1,j) \). This mixes two rows that are zero before column
\( j \), so nothing earlier fills in. Continue with rows \( (j+1,j+2) \) and so on. After \( p-j \) rotations,
\( \mathbf{G}\mathbf{H}=\begin{psmallmatrix}\R_{-j}\\\bzero\T\end{psmallmatrix} \) with \( \R_{-j} \) upper triangular, and
\( \X_{-j}=(\Q_1\mathbf{G}\T)\begin{psmallmatrix}\R_{-j}\\\bzero\T\end{psmallmatrix} \), where \( \Q_1\mathbf{G}\T \) has orthonormal
columns. Apply the same rotations to \( \mathbf{c}_1=\Q_1\T\y \). The last entry of the rotated vector is the
coordinate of \( \y \) along the direction that has dropped out of the model, and its square is the rise in the
residual sum of squares. This is the partial sum of squares of column \( j \) (@def-ss-partial).
:::

### C. Going deeper

::: {#exr-cmp-recursive-residuals}
[C1]

Feed the rows of \( [\X,\y] \) one at a time into the augmented triangle, as in @thm-cmp-givens(c). For
\( k>p \), let \( t_k \) be the entry left in the new row when case \( k \) is added. Show that
\[
t_k=\pm\frac{y_k-\x_{(k)}\T\hbeta_{k-1}}{(1+\x_{(k)}\T(\X_{k-1}\T\X_{k-1})^{-1}\x_{(k)})^{1/2}},
\]
where \( \hbeta_{k-1} \) and \( \X_{k-1} \) use the first \( k-1 \) cases. Show that under the normal linear model the
\( t_k \) are independent \( \Normal(0,\sigma^2) \). These are the recursive residuals of Brown, Durbin and Evans (1975),
used to test the constancy of a regression over time.
:::
