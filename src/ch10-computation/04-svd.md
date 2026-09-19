# The singular value decomposition

QR is the workhorse when \( \X \) has full column rank. When it does not, or nearly does not,
the singular value decomposition tells the whole story: every least squares solution, the one of
minimum norm, and, in the smallest singular values, the combinations of the columns that the data
cannot resolve. @thm-mat-svd established the decomposition, and @prp-proj-min-norm previewed its
use for least squares.

## Least squares through the SVD

Throughout, \( \X \) is \( n\times p \) of rank \( r\ge1 \), with singular value decomposition
\( \X=\sum_{i\le r}\sigma_i\bu_i\bv_i\T \) (@thm-mat-svd). We complete \( \bu_1,\dots,\bu_r \) and
\( \bv_1,\dots,\bv_r \) to orthonormal bases \( \bu_1,\dots,\bu_n \) of \( \Real^n \) and \( \bv_1,\dots,\bv_p \) of
\( \Real^p \), and write \( c_i=\bu_i\T\y \) for the coordinates of the data in the left singular basis.

::: {#thm-cmp-svd-ls}
[Least squares via the SVD]

::: {.enumerate options="label=(\alph*)"}
1. \( \bb \) is a least squares estimate iff \( \bv_i\T\bb=c_i/\sigma_i \) for \( i=1,\dots,r \). The
           coordinates \( \bv_i\T\bb \) for \( i>r \) are arbitrary, and the set of least squares
           estimates is \( \X^+\y+\Null(\X) \) with \( \Null(\X)=\spn(\bv_{r+1},\dots,\bv_p) \).

2. \( \X^+\y=\sum_{i\le r}(c_i/\sigma_i)\,\bv_i \) is the unique least squares estimate of smallest
           Euclidean norm.

3. \( \hY=\sum_{i\le r}c_i\bu_i \) and \( \text{SSE}=\sum_{i>r}c_i^2=\norm{\y}^2-\sum_{i\le r}c_i^2 \).

4. If \( r=p \) and \( \Cov(\Y)=\sigma^2\I \), the estimate \( \hbeta=\X^+\Y \) has
           \( \Cov(\hbeta)=\sigma^2\sum_{i\le p}\sigma_i^{-2}\bv_i\bv_i\T \). In particular
           \( \Var(\bv_i\T\hbeta)=\sigma^2/\sigma_i^2 \), and the total variance is
           \( \tr\Cov(\hbeta)=\sigma^2\sum_i\sigma_i^{-2} \).
:::

:::

::: {.proof}
Let \( \bU=[\bu_1,\dots,\bu_n] \), an orthogonal matrix. For any \( \bb \), \( \X\bb=\sum_{i\le r}\sigma_i(\bv_i\T\bb)\bu_i \),
so the coordinates of \( \X\bb \) in the basis \( \bu_1,\dots,\bu_n \) are \( \sigma_i\bv_i\T\bb \) for \( i\le r \) and
zero for \( i>r \). By orthogonal invariance,
\[
\norm{\y-\X\bb}^2=\norm{\bU\T\y-\bU\T\X\bb}^2=\sum_{i\le r}\bigl(c_i-\sigma_i\bv_i\T\bb\bigr)^2+\sum_{i>r}c_i^2 .
\]{#eq-cmp-svd-rss}

The second sum does not involve \( \bb \), and the first vanishes exactly when
\( \bv_i\T\bb=c_i/\sigma_i \) for \( i\le r \). This proves (a) and (c), since \( \bb=\sum_i(\bv_i\T\bb)\bv_i \)
and \( \Null(\X) \) is spanned by the \( \bv_i \) with \( i>r \) (@thm-mat-svd). Part (b) is @prp-proj-min-norm. In
coordinates, \( \norm{\bb}^2=\sum_{i\le p}(\bv_i\T\bb)^2 \), so the minimum-norm estimate is the one
whose free coordinates are zero, \( \sum_{i\le r}(c_i/\sigma_i)\bv_i \). For (d),
\( \hbeta=\sum_i\sigma_i^{-1}\bv_i\,\bu_i\T\Y \), and \( \Cov(\bu_i\T\Y,\bu_j\T\Y)=\sigma^2\bu_i\T\bu_j \) is
\( \sigma^2 \) or \( 0 \). By @thm-rv-linear,
\( \Cov(\hbeta)=\sigma^2\sum_i\sigma_i^{-2}\bv_i\bv_i\T \). Its trace is \( \sigma^2\sum_i\sigma_i^{-2} \)
because \( \norm{\bv_i}=1 \).
:::

Part (d) is the statistical face of ill-conditioning. The combination \( \bv_p\T\bbeta \) of the
coefficients, along the direction in which the columns of \( \X \) are nearly dependent, is estimated
with standard deviation \( \sigma/\sigma_p \). When \( \sigma_p \) is tiny, the data carry almost no
information about that combination, and no algorithm can change this.

```{.python .run #cell-svd-rank-svd-ls}
import numpy as np

def svd_ls(X, y, tol=None):
    """Minimum-norm least squares solution, treating singular values <= tol as zero."""
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    if tol is None:
        tol = max(X.shape) * np.finfo(float).eps * s[0]   # the usual default
    r = int(np.sum(s > tol))                               # numerical rank
    c = U[:, :r].T @ y                                     # coordinates of the fit
    b = Vt[:r].T @ (c / s[:r])
    sse = y @ y - c @ c
    return b, r, sse, s
```

## Computing the SVD

The SVD is not computed from the eigenvalues of \( \X\T\X \), which has already lost the small singular
values ([Section 10.1](01-why-not-invert.html)). The standard method (Golub and Kahan 1965; Golub and
Reinsch 1970) first applies Householder reflections alternately from the left and the right to reduce
\( \X \) to upper bidiagonal form, at about twice the cost of QR, and then drives the superdiagonal to zero
by an implicitly shifted QR-type iteration that never forms a cross-product. When \( n\gg p \), it pays to
compute \( \X=\Q_1\R \) first and then an SVD \( \R=\bU_R\bSigma\V\T \) of the triangle, since
\( \X=(\Q_1\bU_R)\bSigma\V\T \) (@exr-cmp-rsvd). Least squares needs only \( \bU\T\y \), never \( \bU \) itself, which
gives the count \( 2np^2+11p^3 \) of [Section 10.1](01-why-not-invert.html) (Golub and Van Loan 2013).

The computed singular values are those of a nearby matrix \( \X+\bm E \) with
\( \norm{\bm E}_2\le c\,u\norm{\X}_2 \). The next result turns this into accuracy.

::: {#prp-cmp-weyl}
[Singular values are well conditioned]

For any matrices \( \A \) and \( \bm E \) of the same size and every \( k \),
\[
\lvert\sigma_k(\A+\bm E)-\sigma_k(\A)\rvert\le\norm{\bm E}_2 .
\]
:::

::: {.proof}
For any matrix \( \bm C \) and any \( \B \) of rank at most \( k-1 \), \( \sigma_k(\bm C)\le\norm{\bm C-\B}_2 \).
When \( k-1<\rank\bm C \) this is the Eckart–Young inequality (@prp-mat-svd-norms(c)), and
otherwise \( \sigma_k(\bm C)=0 \). Let \( \A_{k-1} \) be the truncated SVD of \( \A \) with \( k-1 \) terms, or \( \A \)
itself if \( \rank\A<k \), so that \( \norm{\A-\A_{k-1}}_2=\sigma_k(\A) \). Then
\[
\sigma_k(\A+\bm E)\le\norm{\A+\bm E-\A_{k-1}}_2\le\sigma_k(\A)+\norm{\bm E}_2 .
\]
Exchanging the roles of \( \A \) and \( \A+\bm E \) gives the other inequality.
:::

So a backward stable SVD computes every singular value with *absolute* error of order \( u\sigma_1 \),
and a singular value of that size or below is indistinguishable from zero. An exactly rank-deficient
design, such as a one-way layout with an intercept and all its indicators, therefore shows a smallest
singular value near \( 10^{-16}\sigma_1 \) rather than zero: \( 5\times 10^{-17} \) times the largest in the
script's example.

## Numerical rank

In floating point, and with data measured to finite precision, the meaningful question is how
many singular values are *clearly* nonzero.

::: {#def-cmp-numerical-rank}
[Numerical rank]

For a tolerance \( \tau>0 \), the **numerical rank** of \( \X \) is
\( r_\tau=\#\{i:\sigma_i>\tau\} \). The **truncated SVD** with \( k \) terms is
\( \X_k=\sum_{i\le k}\sigma_i\bu_i\bv_i\T \).
:::

By the Eckart–Young inequality, \( r_\tau \) is the least rank of any matrix within spectral distance \( \tau \)
of \( \X \) (@exr-cmp-numerical-rank-char), so \( \tau \) should reflect how well \( \X \) is known. The default in
numpy and LAPACK-based software, \( \tau=\max(n,p)\,\epsilon_M\,\sigma_1 \), accounts for rounding only; data
recorded to six significant figures cannot resolve singular values below about \( 10^{-6}\sigma_1 \). The
choice matters, because the minimum-norm solution is not a continuous function of \( \X \).

::: {#exm-cmp-pinv-discontinuity}
[The pseudoinverse jumps]

Let \( \x_1=(1,2,3,4)\T \), let \( \bm d=(1,-1,-1,1)\T \) be orthogonal to it, and let
\( \X_t=[\x_1,\ \x_1+t\bm d] \), with \( \y=(1,3,2,5)\T \). For every \( t\ne0 \) the column space is the
same plane \( \spn(\x_1,\bm d) \), so the fitted values do not depend on \( t \). The coefficients do.
At \( t=10^{-6} \) they are \( -249999 \) and \( 250000 \), because the
coefficients must reproduce a fixed multiple of \( \bm d \) using the small difference
\( t\bm d \) of the two columns. At \( t=0 \) the rank drops to one. The minimum-norm solution
jumps to \( 0.550 \) for both columns, and the fitted values jump too, because the
\( \bm d \)-direction has disappeared from the model. Treating singular values below a tolerance as zero makes the two answers agree. Whether it should is a statistical
question: is the difference between the columns signal, or recording error?
:::

```{.python .run #cell-svd-rank-discontinuity}
x1 = np.array([1.0, 2.0, 3.0, 4.0])
d = np.array([1.0, -1.0, -1.0, 1.0])            # orthogonal to x1
yd = np.array([1.0, 3.0, 2.0, 5.0])
for t in [1e-2, 1e-6, 1e-10, 0.0]:
    Xt = np.column_stack([x1, x1 + t * d])        # second column tends to the first
    b = np.linalg.pinv(Xt, rcond=1e-15) @ yd
    print(f"t = {t:7.0e}   b = {np.round(b, 3)}   fit = {np.round(Xt @ b, 3)}")
```

## Truncated SVD as an estimator

The truncated solution \( \bb_k=\sum_{i\le k}(c_i/\sigma_i)\bv_i \) is more than a numerical
safeguard. It is a statistical estimator in its own right, and it trades bias for variance.

::: {#prp-cmp-tsvd}
[Bias and variance of truncated SVD]

Suppose \( \rank\X=p \), \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\I \), and let
\( \bb_k=\sum_{i\le k}(\bu_i\T\Y/\sigma_i)\bv_i \) for \( k\le p \). Then
\[
\E(\bb_k)=\sum_{i\le k}(\bv_i\T\bbeta)\bv_i,\qquad
\Cov(\bb_k)=\sigma^2\sum_{i\le k}\frac{\bv_i\bv_i\T}{\sigma_i^2},
\]
\[
\E\norm{\bb_k-\bbeta}^2=\sum_{i>k}(\bv_i\T\bbeta)^2+\sigma^2\sum_{i\le k}\frac1{\sigma_i^2}.
\]{#eq-cmp-tsvd-mse}

Hence \( \bb_k \) has smaller mean squared error than the least squares estimate \( \bb_p \) iff
\( \sum_{i>k}(\bv_i\T\bbeta)^2<\sigma^2\sum_{i>k}\sigma_i^{-2} \).
:::

::: {.proof}
\( \E(\bu_i\T\Y)=\bu_i\T\X\bbeta=\sigma_i\bv_i\T\bbeta \), which gives the mean. The covariance is part (d)
of @thm-cmp-svd-ls with the sum stopped at \( k \). For the mean squared error, write
\( \E\norm{\bb_k-\bbeta}^2=\norm{\E\bb_k-\bbeta}^2+\tr\Cov(\bb_k) \). The bias \( \E\bb_k-\bbeta \) is
\( -\sum_{i>k}(\bv_i\T\bbeta)\bv_i \), with squared length \( \sum_{i>k}(\bv_i\T\bbeta)^2 \). Comparing
@eq-cmp-tsvd-mse for \( k \) and \( p \) gives the last claim.
:::

Each dropped component costs its squared coefficient \( (\bv_i\T\bbeta)^2 \) in bias and saves
\( \sigma^2/\sigma_i^2 \) in variance. It is worth dropping exactly when the signal along \( \bv_i \) is
weaker than the noise in its estimate, that is, when \( \sigma_i\lvert\bv_i\T\bbeta\rvert<\sigma \). The
estimator \( \bb_k \) is principal components regression in its uncentred form. Ridge regression replaces
the hard cut-off by weights \( \sigma_i^2/(\sigma_i^2+\lambda) \). Both are developed in Chapter 27,
and choosing \( k \) from data is a model-selection problem (Chapter 29).

::: {#exm-cmp-tsvd}
[Truncating a polynomial design]

The powers \( 1,t,\dots,t^7 \) at \( 60 \) equally spaced points in \( [0,1] \) form a design with
\( \kappa(\X)=1.1\times 10^{5} \) and smallest singular value \( 9.0\times 10^{-5} \).
Take \( \sigma=0.1 \) and a coefficient vector whose components along \( \bv_1,\dots,\bv_8 \) decay
(\( 2,-1.5,1,-0.6,0.3,-0.1,0.05,-0.02 \)), as they do for smooth regression functions. By
@eq-cmp-tsvd-mse, least squares has mean squared error \( 1231055 \). The best truncation keeps
\( 4 \) components and has mean squared error \( 0.21 \), made up of squared bias
\( 0.10 \) and variance \( 0.11 \). A simulation with \( 4000 \) replicates
agrees with the formula at every \( k \) ([Figure 10.4.1](#fig-cmp-tsvd)). The default numerical
tolerance treats the design as having full rank, which it does to working precision. The statistical
"effective rank" is set by the noise level instead.
:::

::: {when-format="html"}
![**Figure 10.4.1.** Truncated SVD for a degree-7 polynomial design. (a) Singular values, with the
default numerical-rank tolerance and the noise level \( \sigma \). (b) Squared bias, variance and mean
squared error of the truncated estimate \( \bb_k \) (@eq-cmp-tsvd-mse), with simulated values.
Least squares is \( k=8 \).](truncated_svd.svg){#fig-cmp-tsvd width=100%}
:::

::: {when-format="pdf"}
![Truncated SVD for a degree-7 polynomial design. (a) Singular values, with the
default numerical-rank tolerance and the noise level \( \sigma \). (b) Squared bias, variance and mean
squared error of the truncated estimate \( \bb_k \) (@eq-cmp-tsvd-mse), with simulated values.
Least squares is \( k=8 \).](truncated_svd.pdf){width=100%}
:::

```{.python .run #cell-svd-rank-tsvd-setup}
n, p, sigma = 60, 8, 0.1
t = np.linspace(0, 1, n)
X = np.vander(t, p, increasing=True)            # 1, t, ..., t^7: badly conditioned
U, s, Vt = np.linalg.svd(X, full_matrices=False)
a = np.array([2.0, -1.5, 1.0, -0.6, 0.3, -0.1, 0.05, -0.02])   # decaying components v_i' beta
beta = Vt.T @ a
print("condition number:", s[0] / s[-1])
print("singular values:", np.array2string(s, precision=2))
```

```{.python .run #cell-svd-rank-tsvd-mse}
bias2 = np.array([np.sum((Vt[k:] @ beta) ** 2) for k in range(1, p + 1)])   # ||(I - V_k V_k') beta||^2
var = np.array([sigma ** 2 * np.sum(1 / s[:k] ** 2) for k in range(1, p + 1)])
mse = bias2 + var
print("k    bias^2      variance     MSE")
for k in range(p):
    print(f"{k + 1}  {bias2[k]:10.3e}  {var[k]:10.3e}  {mse[k]:10.3e}")
```

::: {.idea}
The SVD gives every least squares solution, the precision \( \sigma/\sigma_i \) of each combination
\( \bv_i\T\bbeta \), and, relative to a tolerance, how many combinations are determined at all.
:::

## Exercises

### A. Check your understanding

::: {#exr-cmp-svd-hand}
[A1]

Let \( \X=\begin{psmallmatrix}1&1\\1&1\\0&0\end{psmallmatrix} \) and \( \y=(1,3,5)\T \). Find the SVD of \( \X \),
the set of least squares estimates, the minimum-norm estimate and the SSE.
:::

::: {.solution}
\( \X=2\,\bu_1\bv_1\T \) with \( \bu_1=(1,1,0)\T/\sqrt2 \), \( \bv_1=(1,1)\T/\sqrt2 \), \( \sigma_1=2 \) and
\( r=1 \). Here \( c_1=\bu_1\T\y=4/\sqrt2 \), so the least squares estimates are the \( \bb \) with
\( \bv_1\T\bb=c_1/\sigma_1=\sqrt2 \), that is, \( b_1+b_2=2 \). The minimum-norm estimate is
\( \sqrt2\,\bv_1=(1,1)\T \). The fitted vector is \( c_1\bu_1=(2,2,0)\T \), and
\( \text{SSE}=\norm{\y}^2-c_1^2=35-8=27 \), which is \( (1-2)^2+(3-2)^2+5^2 \).
:::

### B. Practice

::: {#exr-cmp-numerical-rank-char}
[B1]

Show that if \( \tau \) is not a singular value of \( \X \), then
\( r_\tau=\min\{\rank\B:\norm{\X-\B}_2\le\tau\} \), and that the minimum is attained by the truncated
SVD \( \X_{r_\tau} \).
:::

::: {.solution}
\( \norm{\X-\X_k}_2=\sigma_{k+1} \). For \( k=r_\tau \), \( \sigma_{k+1}\le\tau \) (with \( \sigma_{k+1}=0 \) if
\( k=\rank\X \)), so \( \X_{r_\tau} \) is within \( \tau \) and has rank \( r_\tau \). If \( \rank\B=k<r_\tau \), then by
@prp-mat-svd-norms(c) \( \norm{\X-\B}_2\ge\sigma_{k+1}\ge\sigma_{r_\tau}>\tau \). So no matrix of smaller
rank is within \( \tau \).
:::

::: {#exr-cmp-rsvd}
[B2]

Let \( \X=\Q_1\R \) be a thin QR and \( \R=\bU_R\bSigma\V\T \) an SVD of the triangle. Show that
\( \X=(\Q_1\bU_R)\bSigma\V\T \) is an SVD of \( \X \), and that the minimum-norm least squares solution can be
computed from \( \R \) and \( \Q_1\T\y \) alone. Why is this attractive when \( n \) is much larger than \( p \)?
:::

::: {#exr-cmp-tsvd-constrained}
[B3]

Show that \( \bb_k \) minimizes \( \norm{\y-\X\bb} \) over \( \bb\in\spn(\bv_1,\dots,\bv_k) \), and that it is the
minimum-norm least squares estimate for the matrix \( \X_k \). Interpret \( \bb_k \) as least squares after
deleting the directions in which \( \X \) is nearly degenerate.
:::

### C. Going deeper

::: {#exr-cmp-ridge-svd}
[C1]

Show that the ridge estimate \( (\X\T\X+\lambda\I)^{-1}\X\T\y \) equals
\( \sum_i\frac{\sigma_i}{\sigma_i^2+\lambda}c_i\bv_i \). Derive its bias and variance in the form of
@prp-cmp-tsvd, and show that for every \( \bbeta \) there is a \( \lambda>0 \) (depending on \( \bbeta \)
and \( \sigma^2 \)) for which its mean squared error is below that of least squares.
:::
