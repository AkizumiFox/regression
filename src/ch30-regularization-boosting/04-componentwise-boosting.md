# Componentwise boosting

Boosting is an algorithm, not a criterion: it improves the fit by many small steps, each changing the single coefficient
that helps most, and the number of steps plays the part of a penalty. Begun as a way to combine classifiers (Freund and
Schapire, 1997) and recast as gradient descent on a loss (Friedman, Hastie and Tibshirani, 2000; Friedman, 2001), with squared
error and single regressors it is the \( L_2 \) boosting of Bühlmann and Yu (2003).

## The algorithm

Let \( \y \) and the columns of \( \X \) be centred, and let no column be zero. For a vector \( \br \), the least squares fit of \( \br \) on
the single column \( \x_j \) has coefficient \( \x_j\T\br/\norm{\x_j}^2 \). It reduces the sum of squares by
\( (\x_j\T\br)^2/\norm{\x_j}^2 \), by the one-column case of @thm-proj-ls-projection.

::: {#def-reg-l2boost}
[Componentwise \( L_2 \) boosting]

Fix a **step size** \( \nu\in(0,2) \). Start from \( \bb^{(0)}=\bzero \) and \( \br^{(0)}=\y \). For \( m=0,1,2,\dots \):

1. choose \( j_m=\operatorname{arg\,max}_j(\x_j\T\br^{(m)})^2/\norm{\x_j}^2 \), the column whose simple regression fits the current residual best;
2. set \( b^{(m+1)}_{j_m}=b^{(m)}_{j_m}+\nu\,\x_{j_m}\T\br^{(m)}/\norm{\x_{j_m}}^2 \), leave the other coordinates unchanged, and put
   \( \br^{(m+1)}=\y-\X\bb^{(m+1)} \).

The estimate after \( m \) steps is \( \bb^{(m)} \), and the **stopping iteration** \( m \) is the tuning parameter.
:::

With unit-length columns, step 1 picks the column most correlated with the residual. The residual is the negative
gradient of \( \tfrac12\norm{\y-\mathbf{f}}^2 \) in the fitted vector \( \mathbf{f} \), so each step fits the negative gradient by
the best single column and moves a fraction \( \nu \) of the way. The theory below allows \( \nu<2 \); with \( \nu>1 \) a
step overshoots the simple regression fit, and in practice \( \nu\le1 \), often \( \nu=0.1 \). With \( \nu=1 \) this is *matching pursuit* (Mallat and
Zhang, 1993), which is not forward selection, since earlier coefficients are never refitted.

## Convergence to least squares

Boosting converges to the least squares fit at a geometric rate, slow when the design is badly conditioned.

::: {#thm-reg-boosting}
[Convergence of componentwise \( L_2 \) boosting]

Let \( \nu\in(0,2) \), let \( Q(\bb)=\tfrac12\norm{\y-\X\bb}^2 \) with minimum \( Q^* \), let \( s>0 \) be the smallest nonzero singular value of
\( \X \), and put
\[
\kappa=\frac{\nu(2-\nu)\,s^2}{p\max_j\norm{\x_j}^2}.
\]
Then \( 0<\kappa\le1 \), and:

::: {.enumerate options="label=(\alph*)"}
1. each step lowers the residual sum of squares by exactly
   \[
   Q(\bb^{(m)})-Q(\bb^{(m+1)})=\nu\Bigl(1-\frac\nu2\Bigr)\max_j\frac{(\x_j\T\br^{(m)})^2}{\norm{\x_j}^2};
   \]{#eq-reg-boost-drop}

2. \( Q(\bb^{(m)})-Q^*\le(1-\kappa)^m(Q(\bzero)-Q^*) \). Hence the fitted values converge geometrically to the least
   squares fit: \( \norm{\X\bb^{(m)}-\M\y}^2=2\{Q(\bb^{(m)})-Q^*\}\to0 \);

3. the coefficients \( \bb^{(m)} \) converge to a least squares estimate. If \( \X \) has full column rank, the limit is \( \hbeta \).
:::

:::

::: {.proof}
(a) Write \( \br=\br^{(m)} \), \( j=j_m \) and \( c=\x_j\T\br/\norm{\x_j}^2 \). Then \( \br^{(m+1)}=\br-\nu c\x_j \), and
\[
\begin{aligned}
\norm{\br-\nu c\x_j}^2&=\norm{\br}^2-2\nu c\,\x_j\T\br+\nu^2c^2\norm{\x_j}^2\\
&=\norm{\br}^2-(2\nu-\nu^2)\frac{(\x_j\T\br)^2}{\norm{\x_j}^2}.
\end{aligned}
\]
Halving gives @eq-reg-boost-drop, since \( j \) maximizes the ratio.

(b) Let \( \bu_m=\M\br^{(m)}=\M\y-\X\bb^{(m)} \), which lies in \( \C(\X) \). The residual splits orthogonally as
\( \br^{(m)}=(\I-\M)\y+\bu_m \), so \( Q(\bb^{(m)})-Q^*=\tfrac12\norm{\bu_m}^2 \) and \( \X\T\br^{(m)}=\X\T\bu_m \). By the singular value
decomposition (@thm-mat-svd), \( \norm{\X\T\bu}\ge s\norm{\bu} \) for every \( \bu\in\C(\X) \). A maximum is at least an average, so
\[
\begin{aligned}
\max_j\frac{(\x_j\T\br^{(m)})^2}{\norm{\x_j}^2}&\ge\frac{\sum_j(\x_j\T\bu_m)^2}{p\max_j\norm{\x_j}^2}\\
&=\frac{\norm{\X\T\bu_m}^2}{p\max_j\norm{\x_j}^2}\ge\frac{s^2\norm{\bu_m}^2}{p\max_j\norm{\x_j}^2}.
\end{aligned}
\]
With (a), \( Q(\bb^{(m)})-Q(\bb^{(m+1)})\ge\kappa\{Q(\bb^{(m)})-Q^*\} \), that is,
\( Q(\bb^{(m+1)})-Q^*\le(1-\kappa)\{Q(\bb^{(m)})-Q^*\} \). Iterating gives (b). Finally \( \kappa>0 \), and \( \kappa\le1 \) because \( \nu(2-\nu)\le1 \) and
\( s^2\le\tr(\X\T\X)=\sum_j\norm{\x_j}^2\le p\max_j\norm{\x_j}^2 \).

(c) The \( m \)th step changes one coordinate, by
\( \nu\lvert\x_{j_m}\T\bu_m\rvert/\norm{\x_{j_m}}^2\le\nu\norm{\x_{j_m}}\norm{\bu_m}/\norm{\x_{j_m}}^2\le\nu\norm{\bu_m}/\min_j\norm{\x_j} \), using Cauchy–Schwarz.
By (b), \( \norm{\bu_m}\le(1-\kappa)^{m/2}\norm{\bu_0} \), which is summable. So \( \sum_m\norm{\bb^{(m+1)}-\bb^{(m)}}<\infty \), the sequence is Cauchy,
and it converges to some \( \bb^\infty \). By continuity, \( \X\bb^\infty=\lim\X\bb^{(m)}=\M\y \), so \( \bb^\infty \) is a least squares estimate (@thm-proj-ls-projection). With full column rank there is only one.
:::

The bound is crude, but it shows what controls the speed: the step size and the conditioning of \( \X \).

::: {.remark}
[When \( p>n \)]

If \( \rank\X=n \), as is typical when \( p>n \), boosting run long enough *interpolates* the data, so stopping early is a
necessity. Its limit depends on the order in which columns were chosen and is in general not the least squares estimate of
smallest norm (@prp-proj-min-norm), which full gradient descent from zero does reach (@exr-reg-landweber; this is the ridgeless
estimator of @thm-hd-min-norm).
:::

## Early stopping as shrinkage

For a base learner linear in the data, the effect of stopping early is exact. Let an \( n\times n \) matrix \( \bS \) map a
residual vector to its fitted values; boosting updates \( \mathbf{f}^{(m+1)}=\mathbf{f}^{(m)}+\nu\bS(\y-\mathbf{f}^{(m)}) \) from
\( \mathbf{f}^{(0)}=\bzero \).

::: {#prp-reg-boosting-operator}
[The boosting operator]

::: {.enumerate options="label=(\alph*)"}
1. \( \mathbf{f}^{(m)}=\mathbf{B}_m\y \) with \( \mathbf{B}_m=\I-(\I-\nu\bS)^m \).

2. Let \( \bS \) be symmetric with eigenvalues \( d_k\ge0 \) and orthonormal eigenvectors \( \bu_k \), and let \( \nu d_k\le1 \) for every \( k \). Then
   \[
   \mathbf{B}_m=\sum_k\bigl\{1-(1-\nu d_k)^m\bigr\}\bu_k\bu_k\T .
   \]
   Each factor \( 1-(1-\nu d_k)^m \) lies in \( [0,1] \), is nondecreasing in \( m \), and tends to \( 1 \) when \( d_k>0 \). The trace
   \( \tr\mathbf{B}_m \) increases with \( m \) towards \( \rank\bS \).

3. *(Least squares learner.)* If \( \bS=\M \), then \( \mathbf{B}_m=\{1-(1-\nu)^m\}\M \). Early stopping shrinks the least squares fit by a
   single factor.

4. *(Gradient descent.)* The iteration \( \bb^{(m+1)}=\bb^{(m)}+\nu\X\T(\y-\X\bb^{(m)}) \) from \( \bb^{(0)}=\bzero \) is the case
   \( \bS=\X\X\T \). If \( \X \) has singular values \( \sigma_k \) and \( \nu\sigma_k^2\le1 \), its fitted values shrink the \( k \)th singular component
   of the least squares fit by \( 1-(1-\nu\sigma_k^2)^m \). Ridge regression shrinks the same component by \( \sigma_k^2/(\sigma_k^2+\lambda) \).
:::

:::

::: {.proof}
(a) \( \y-\mathbf{f}^{(m+1)}=(\I-\nu\bS)(\y-\mathbf{f}^{(m)}) \), so \( \y-\mathbf{f}^{(m)}=(\I-\nu\bS)^m\y \). (b) \( (\I-\nu\bS)^m=\sum_k(1-\nu d_k)^m\bu_k\bu_k\T \)
by the spectral theorem (@thm-mat-spectral), and \( 0\le1-\nu d_k\le1 \). (c) \( \M \) has eigenvalues \( 0 \) and \( 1 \), and it is idempotent.
(d) \( \X\bb^{(m+1)}=\X\bb^{(m)}+\nu\X\X\T(\y-\X\bb^{(m)}) \). The eigenvalues of \( \X\X\T \) are the \( \sigma_k^2 \) and zeros, and the
eigenvectors for the \( \sigma_k^2 \) are the left singular vectors (@thm-mat-svd). For ridge, see @thm-shr-ridge.
:::

For small \( \sigma_k^2 \), gradient descent keeps the fraction
\( 1-(1-\nu\sigma_k^2)^m\approx m\nu\sigma_k^2 \) of a component and ridge keeps \( \approx\sigma_k^2/\lambda \); for large \( \sigma_k^2 \) both keep nearly
everything. So \( m \) steps behave like ridge with \( \lambda\approx1/(m\nu) \). Part (c) is the proportional
shrinkage of @prp-opt-shrinkage(c).

**Degrees of freedom.** For a linear smoother the complexity is \( \tr\mathbf{B}_m \), which equals the covariance penalty
\( \sigma^{-2}\sum_i\Cov(\hat y_i,y_i) \) (@thm-sel-optimism). Componentwise boosting is not linear, since the chosen columns
depend on \( \y \), but for a fixed sequence \( j_0,\dots,j_{m-1} \) its fitted values are \( \mathbf{B}_m\y \) with
\[
\mathbf{B}_m=\I-(\I-\nu\bH_{j_{m-1}})\cdots(\I-\nu\bH_{j_0}),\qquad \bH_j=\frac{\x_j\x_j\T}{\norm{\x_j}^2},
\]{#eq-reg-boost-operator}

by @prp-reg-boosting-operator(a) applied step by step. Bühlmann (2006) used \( \tr\mathbf{B}_m \) along the observed sequence as
the degrees of freedom in a corrected Akaike criterion for choosing \( m \). That treats the selection as fixed in advance, and
the next example shows how much it can understate.

::: {#exm-reg-boosting}
[Boosting a sparse regression]

The script draws \( n=60 \) observations of \( p=30 \) normal regressors with correlations \( 0.5^{\lvert j-k\rvert} \), five nonzero
coefficients and \( \sigma=1.5 \), and boosts on the standardized regressors. The expected squared error on a new case (the
*risk*) is computed exactly from the design law; its floor is \( \sigma^2=2.25 \).

*Convergence.* For \( \nu=0.1 \), \( \kappa=3.18\times 10^{-4} \), and \( \norm{\bb^{(m)}-\hbeta}/\norm{\hbeta} \) is \( 0.497 \), \( 0.161 \)
and \( 2.7\times 10^{-8} \) after \( 100 \), \( 1000 \) and \( 50000 \) steps; the script checks @eq-reg-boost-drop and the bound (b)
at every step.

*Early stopping.* Least squares has risk \( 4.046 \). With \( \nu=0.1 \) the risk is smallest, \( 2.981 \), after \( 108 \) steps (with
\( 14 \) nonzero coefficients), and within 5 per cent of that from \( m=72 \) to \( m=305 \); five-fold cross-validation chooses
\( m=229 \), with risk \( 3.047 \). With \( \nu=1 \) the best risk, \( 2.893 \), occurs at \( m=4 \) and only there is the risk within 5 per
cent of it: a needle that a data-based stopping rule is unlikely to hit ([Figure 30.4.1](#fig-reg-boosting)).

*Degrees of freedom.* Over \( 2000 \) responses with the same design and mean, the covariance degrees of freedom after
\( m=20 \) steps are \( 3.1 \) against an average trace of \( 1.5 \), and after \( m=120 \) steps \( 12.6 \) against \( 6.2 \).
The trace does not see the adaptation of the chosen columns to the noise.
:::

::: {when-format="html"}
![**Figure 30.4.1.** Componentwise boosting for @exm-reg-boosting. (a) Coefficient paths with \( \nu=0.1 \) (true
nonzero coefficients in blue); the dashed line is the cross-validated stopping iteration. (b) Exact risk along the path
for \( \nu=0.1 \) and \( \nu=1 \), with the least squares risk dashed. Both horizontal axes are linear up to \( 10 \) and logarithmic
beyond. (c) Two degrees of freedom: the average trace of the boosting operator @eq-reg-boost-operator, and the covariance degrees of
freedom estimated from \( 2000 \) simulated responses.](boosting.svg){#fig-reg-boosting width=100%}
:::

::: {when-format="pdf"}
![Componentwise boosting for @exm-reg-boosting. (a) Coefficient paths with \( \nu=0.1 \) (true
nonzero coefficients in blue); the dashed line is the cross-validated stopping iteration. (b) Exact risk along the path
for \( \nu=0.1 \) and \( \nu=1 \), with the least squares risk dashed. Both horizontal axes are linear up to \( 10 \) and logarithmic
beyond. (c) Two degrees of freedom: the average trace of the boosting operator @eq-reg-boost-operator, and the covariance degrees of
freedom estimated from \( 2000 \) simulated responses.](boosting.pdf){width=100%}
:::

```{.python .run #cell-boosting-test}
import numpy as np

def l2_boost(X, y, nu, steps):
    """Componentwise L2 boosting from b = 0; returns the coefficient path (steps + 1) x p."""
    n, p = X.shape
    norms2 = np.sum(X ** 2, axis=0)
    b = np.zeros(p)
    r = y.astype(float).copy()
    path = [b.copy()]
    for _ in range(steps):
        corr = X.T @ r
        j = np.argmax(corr ** 2 / norms2)               # largest drop in the residual sum of squares
        step = nu * corr[j] / norms2[j]
        b[j] += step
        r -= step * X[:, j]
        path.append(b.copy())
    return np.array(path)

rng = np.random.default_rng(3004)
n, p, sigma = 60, 30, 1.5
C = 0.5 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
L = np.linalg.cholesky(C)
beta = np.zeros(p)
beta[[0, 4, 9, 14, 19]] = [2.0, -1.5, 1.0, 1.0, -0.5]

def draw(m):
    Z = rng.normal(size=(m, p)) @ L.T
    return Z, Z @ beta + sigma * rng.normal(size=m)

X, y = draw(n)
mx, sx, my = X.mean(axis=0), X.std(axis=0), y.mean()
Xs, ys = (X - mx) / sx, y - my                          # centre and scale by the training data

def risk(path):
    """Expected squared error on a new case, for each row b of path (exact, from the design law)."""
    d = beta[None, :] - path / sx                       # coefficient error on the original scale
    offset = path @ (mx / sx) - my                      # error in the fitted intercept
    return np.sum((d @ C) * d, axis=1) + offset ** 2 + sigma ** 2

b_ls = np.linalg.lstsq(Xs, ys, rcond=None)[0]
M = 1500
test_err = {}
for step in [0.1, 1.0]:
    test_err[step] = risk(l2_boost(Xs, ys, step, M))
err_ls = risk(b_ls[None, :])[0]
m_best = int(np.argmin(test_err[0.1]))
print(f"least squares: risk {err_ls:.3f}")
print(f"nu = 0.1: best m = {m_best}, risk {test_err[0.1][m_best]:.3f}")
print(f"nu = 1.0: best m = {np.argmin(test_err[1.0])}, risk {test_err[1.0].min():.3f}")
```

The next cell repeats the comparison with \( 300 \) responses and \( m=20 \); the script's run of it gives \( 3.2 \) against \( 1.5 \).

```{.python .run #cell-boosting-df-small}
mu = Xs @ (beta * sx)                                    # the true mean, in the centred scale
def operator_trace(X, path, nu):
    """Trace of B_m = I - prod_k (I - nu H_{j_k}), with the selected columns held fixed."""
    n = X.shape[0]
    Rm = np.eye(n)
    for k in range(len(path) - 1):
        j = int(np.flatnonzero(path[k + 1] != path[k])[0])
        x = X[:, j]
        Rm -= nu * np.outer(x, x @ Rm) / (x @ x)
    return n - np.trace(Rm)

reps, m_df = 300, 20
rng_df = np.random.default_rng(30041)
Y = mu + sigma * rng_df.normal(size=(reps, n))            # repeated responses, same design and mean
paths = [l2_boost(Xs, Y[r], 0.1, m_df) for r in range(reps)]
F = np.array([Xs @ P[-1] for P in paths])                 # fitted values after m_df steps
cov_df = np.sum(np.mean((Y - Y.mean(axis=0)) * (F - F.mean(axis=0)), axis=0)) / sigma ** 2
trace_df = np.mean([operator_trace(Xs, P, 0.1) for P in paths])
print(f"m = {m_df}: mean trace {trace_df:.1f}, covariance df {cov_df:.1f}")
```

The gap is the optimism of selection (@prp-sel-selection-bias). For the lasso, by contrast, the covariance degrees of freedom
equal the expected number of nonzero coefficients (Zou, Hastie and Tibshirani, 2007; we do not reproduce the proof): shrinkage
offsets the cost of selection (compare @exr-reg-boost-df-orthonormal). In practice \( m \) is chosen by cross-validation
([Chapter 29](../ch29-model-selection/index.html)), which needs no degrees of freedom.

## Exercises

### A. Check your understanding

::: {#exr-reg-boost-orthonormal-step}
[A1]

Let \( \X\T\X=\I \), \( \mathbf{z}=\X\T\y=(3,-1,0.5)\T \) and \( \nu=0.1 \). Carry out the first three boosting steps by hand. Show that under an
orthonormal design the current residual correlations are \( \mathbf{z}-\bb^{(m)} \), so boosting always moves the coordinate whose
remaining correlation is largest in absolute value, and removes one tenth of that correlation.
:::

::: {.solution}
\( \X\T\br^{(m)}=\X\T\y-\X\T\X\bb^{(m)}=\mathbf{z}-\bb^{(m)} \), with unit column lengths. Step 1: coordinate 1, \( b_1=0.3 \), residual
correlation \( 2.7 \). Step 2: coordinate 1 again, \( b_1=0.3+0.27=0.57 \), correlation \( 2.43 \). Step 3: coordinate 1, \( b_1=0.57+0.243=0.813 \).
Coordinate 1 keeps being chosen until its remaining correlation \( 3(0.9)^m \) falls below \( 1 \), which first happens at
\( m=11 \). This is the approach to soft thresholding that [Section 30.5](05-boosting-as-regularization.html) makes
precise.
:::

### B. Practice

::: {#exr-reg-boost-drop-sum}
[B1]

Use @eq-reg-boost-drop to show that \( \sum_{m\ge0}\max_j(\x_j\T\br^{(m)})^2/\norm{\x_j}^2\le\norm{\y}^2/\{\nu(2-\nu)\} \). Deduce, without the
singular value argument, that \( \X\T\br^{(m)}\to\bzero \) for every \( \nu\in(0,2) \).
:::

::: {.solution}
Summing @eq-reg-boost-drop over \( m=0,\dots,M-1 \) gives
\( \nu(1-\nu/2)\sum_{m<M}\max_j(\cdots)=Q(\bzero)-Q(\bb^{(M)})\le\tfrac12\norm{\y}^2 \). Multiply by \( 2/\{\nu(2-\nu)\} \) and let \( M\to\infty \). The terms of a
convergent series tend to zero. So \( \max_j\lvert\x_j\T\br^{(m)}\rvert/\norm{\x_j}\to0 \), and hence every coordinate of \( \X\T\br^{(m)} \) does.
:::

::: {#exr-reg-landweber}
[B2]

For the gradient descent iteration of @prp-reg-boosting-operator(d) with \( 0<\nu<2/\sigma_1^2 \), show that every \( \bb^{(m)} \) lies in
\( \C(\X\T) \), and that \( \bb^{(m)}\to\X^+\y \), the minimum-norm least squares estimate (@prp-proj-min-norm). Why does this argument fail for
componentwise boosting?
:::

::: {.solution}
Each increment \( \nu\X\T(\y-\X\bb^{(m)}) \) lies in \( \C(\X\T) \), and \( \bb^{(0)}=\bzero \). In the singular value basis,
\( \bb^{(m)}=\sum_k\{1-(1-\nu\sigma_k^2)^m\}\sigma_k^{-1}(\bu_k\T\y)\bv_k \), because the iteration acts on each component separately. With
\( \lvert1-\nu\sigma_k^2\rvert<1 \) this converges to \( \sum_k\sigma_k^{-1}(\bu_k\T\y)\bv_k=\X^+\y \). Componentwise boosting adds multiples of
coordinate vectors \( \mathbf{e}_j \), which need not lie in \( \C(\X\T) \) when \( \rank\X<p \). So its limit can have a component in \( \Null(\X) \).
:::

### C. Going deeper

::: {#exr-reg-boost-df-orthonormal}
[C1]

Under an orthonormal design with \( \y\sim\Normal_n(\X\bbeta,\sigma^2\I) \), the lasso fit with parameter \( \lambda \) has covariance
degrees of freedom \( \E\#\{j:\lvert z_j\rvert>\lambda\} \) (by Stein's lemma, @lem-shr-stein, applied coordinatewise). Show that
hard thresholding at the same \( \lambda \) has covariance degrees of freedom
\( \E\#\{j:\lvert z_j\rvert>\lambda\}+\sum_j(\lambda/\sigma)\{\phi((\lambda-\beta_j)/\sigma)+\phi((\lambda+\beta_j)/\sigma)\} \), where \( \phi \) is the
standard normal density. The extra term is the price of selecting without shrinking.
:::

