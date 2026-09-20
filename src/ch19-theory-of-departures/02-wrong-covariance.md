# An incorrect covariance matrix

Suppose the mean is right but the errors are correlated or have unequal variances. Write
\[
\E(\Y)=\X\bbeta,\qquad \Cov(\Y)=\sigma^2\V ,
\]{#eq-dep-cov-model}

with \( \V \) positive definite and \( \X \) of full rank \( p<n \). The analyst fits ordinary least squares and
reports \( s^2(\X\T\X)^{-1} \) as if \( \V=\I \). Two things go wrong: \( \hbeta \) is no longer the best linear unbiased
estimator (@cor-opt-aitken), and \( s^2(\X\T\X)^{-1} \) no longer estimates its covariance. The second failure is
usually much more serious.

## What happens to the estimator and its standard error

The bounds below need a fact about traces of projections, a special case of an inequality of
Ky Fan.

::: {#lem-dep-ky-fan}
[Trace of a compressed matrix]

Let \( \A \) be symmetric with eigenvalues \( \lambda_1\ge\dots\ge\lambda_n \), and let \( \bP \) be a symmetric idempotent
matrix of rank \( k \). Then
\[
\sum_{i=n-k+1}^{n}\lambda_i\ \le\ \tr(\bP\A)\ \le\ \sum_{i=1}^{k}\lambda_i .
\]
:::

::: {.proof}
By @thm-mat-spectral, \( \A=\sum_i\lambda_i\mathbf{q}_i\mathbf{q}_i\T \) with orthonormal \( \mathbf{q}_i \). So
\( \tr(\bP\A)=\sum_i\lambda_iw_i \) with \( w_i=\mathbf{q}_i\T\bP\mathbf{q}_i=\norm{\bP\mathbf{q}_i}^2 \). Each \( w_i \) lies in
\( [0,1] \) because the eigenvalues of \( \bP \) are \( 0 \) or \( 1 \) (@prp-proj-trace-rank), so
\( 0\le\mathbf{q}_i\T\bP\mathbf{q}_i\le\mathbf{q}_i\T\mathbf{q}_i=1 \), and
\( \sum_iw_i=\tr\bigl(\bP\sum_i\mathbf{q}_i\mathbf{q}_i\T\bigr)=\tr\bP=k \). Then
\[
\begin{aligned}
\sum_i\lambda_iw_i-\sum_{i\le k}\lambda_i&=\sum_{i\le k}\lambda_i(w_i-1)+\sum_{i>k}\lambda_iw_i\\
&\le\lambda_k\sum_{i\le k}(w_i-1)+\lambda_k\sum_{i>k}w_i=\lambda_k\Bigl(\sum_iw_i-k\Bigr)=0,
\end{aligned}
\]
since \( w_i-1\le0 \) and \( \lambda_i\ge\lambda_k \) for \( i\le k \), while \( w_i\ge0 \) and \( \lambda_i\le\lambda_k \) for
\( i>k \). The lower bound follows by applying this to \( -\A \).
:::

::: {#thm-dep-covariance}
[Least squares under a wrong covariance]

Under @eq-dep-cov-model, let \( \mathbf{a}\in\Real^p \) be nonzero and let \( v=s^2\mathbf{a}\T(\X\T\X)^{-1}\mathbf{a} \) be the usual
estimate of \( \Var(\mathbf{a}\T\hbeta) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\hbeta)=\bbeta \) and \( \Cov(\hbeta)=\sigma^2(\X\T\X)^{-1}\X\T\V\X(\X\T\X)^{-1} \).

2. \( \E(s^2)=\sigma^2\tr\bigl((\I-\M)\V\bigr)/(n-p) \).

3. With \( \bu=\X(\X\T\X)^{-1}\mathbf{a}/\norm{\X(\X\T\X)^{-1}\mathbf{a}} \), a unit vector in \( \C(\X) \),
   \[
   \frac{\E(v)}{\Var(\mathbf{a}\T\hbeta)}=\frac{\tr\bigl((\I-\M)\V\bigr)/(n-p)}{\bu\T\V\bu} .
   \]{#eq-dep-se-ratio}

4. Let \( \lambda_1\ge\dots\ge\lambda_n>0 \) be the eigenvalues of \( \V \), let \( \bar{\lambda}_- \) be the mean of the
   \( n-p \) smallest and \( \bar{\lambda}_+ \) the mean of the \( n-p \) largest. Then
   \[
   \frac{\bar{\lambda}_-}{\lambda_1}\le\frac{\E(v)}{\Var(\mathbf{a}\T\hbeta)}\le\frac{\bar{\lambda}_+}{\lambda_n} .
   \]
:::

:::

::: {.proof}
(a) is @prp-lm-misspecified(b) with \( \bSigma=\sigma^2\V \). (b) is @thm-rv-quadform-mean with the mean term
vanishing, since \( (\I-\M)\X\bbeta=\bzero \). (c) Put \( \mathbf{c}=\X(\X\T\X)^{-1}\mathbf{a} \). By (a),
\( \Var(\mathbf{a}\T\hbeta)=\sigma^2\mathbf{c}\T\V\mathbf{c} \), and
\( \mathbf{c}\T\mathbf{c}=\mathbf{a}\T(\X\T\X)^{-1}\mathbf{a} \). By (b),
\( \E(v)=\sigma^2\mathbf{c}\T\mathbf{c}\,\tr\bigl((\I-\M)\V\bigr)/(n-p) \). Divide, and write \( \bu=\mathbf{c}/\norm{\mathbf{c}} \).
(d) \( \I-\M \) is symmetric idempotent of rank \( n-p \), so by @lem-dep-ky-fan the numerator of
@eq-dep-se-ratio lies in \( [\bar{\lambda}_-,\bar{\lambda}_+] \). By @thm-mat-extremal-rayleigh the denominator lies in
\( [\lambda_n,\lambda_1] \).
:::

Part (c) shows what the usual standard error gets wrong. The numerator of @eq-dep-se-ratio averages \( \V \) over
the residual space, which is what \( s^2 \) measures. The denominator evaluates \( \V \) in the one direction \( \bu \)
that \( \mathbf{a}\T\hbeta \) uses. When the errors vary more along \( \bu \) than on average, the usual standard error is too
small. The bounds in (d), due to Swindel (1968), hold for every design and linear function and can be attained (@exr-dep-swindel-attained). They are typically wide.

## When ordinary least squares loses little

Generalized least squares has variance \( \sigma^2\mathbf{a}\T(\X\T\V^{-1}\X)^{-1}\mathbf{a} \) for \( \mathbf{a}\T\bbeta \) (@cor-opt-aitken). How much larger can the variance of \( \mathbf{a}\T\hbeta \) be? The answer again depends only
on the extreme eigenvalues of \( \V \).

::: {#thm-dep-efficiency}
[Efficiency of least squares]

Under @eq-dep-cov-model, with \( \lambda_1 \) and \( \lambda_n \) the largest and smallest eigenvalues of \( \V \),
\[
\frac{\mathbf{a}\T(\X\T\V^{-1}\X)^{-1}\mathbf{a}}{\mathbf{a}\T(\X\T\X)^{-1}\X\T\V\X(\X\T\X)^{-1}\mathbf{a}}
\ \ge\ \frac{4\lambda_1\lambda_n}{(\lambda_1+\lambda_n)^2}
\]
for every \( \mathbf{a}\ne\bzero \).
The bound is attained for \( p=1 \) and \( \X=\x=(\mathbf{q}_1+\mathbf{q}_n)/\sqrt2 \), where \( \mathbf{q}_1,\mathbf{q}_n \) are unit eigenvectors for
\( \lambda_1,\lambda_n \).
:::

::: {.proof}
Write \( \X=\Q\R \) with \( \Q \) having orthonormal columns and \( \R \) nonsingular (@thm-mat-qr), and put
\( \bb=\R^{-\top}\mathbf{a} \), \( \A=\Q\T\V\Q \) and \( \B=\Q\T\V^{-1}\Q \), both positive definite. Then
\( (\X\T\X)^{-1}\X\T\V\X(\X\T\X)^{-1}=\R^{-1}\A\R^{-\top} \) and \( (\X\T\V^{-1}\X)^{-1}=\R^{-1}\B^{-1}\R^{-\top} \), so the
ratio is \( \bb\T\B^{-1}\bb/\bb\T\A\bb \). It suffices to show \( \A\le\kappa\B^{-1} \) in the nonnegative definite
order, with \( \kappa=(\lambda_1+\lambda_n)^2/(4\lambda_1\lambda_n) \).

Every eigenvalue \( \lambda \) of \( \V \) satisfies \( (\lambda-\lambda_1)(\lambda-\lambda_n)\le0 \), that is,
\( \lambda+\lambda_1\lambda_n/\lambda\le\lambda_1+\lambda_n \). By the spectral decomposition this says
\( \V+\lambda_1\lambda_n\V^{-1}\le(\lambda_1+\lambda_n)\I \). Multiplying by \( \Q\T \) on the left and \( \Q \) on the right
preserves the order, so
\[
\A\le(\lambda_1+\lambda_n)\I-\lambda_1\lambda_n\B .
\]
Now let \( \mu>0 \) be any eigenvalue of \( \B \). By the inequality between arithmetic and geometric means,
\( \kappa/\mu+\lambda_1\lambda_n\mu\ge2\sqrt{\kappa\lambda_1\lambda_n}=\lambda_1+\lambda_n \). Applied along the
eigenvectors of \( \B \), this gives \( (\lambda_1+\lambda_n)\I-\lambda_1\lambda_n\B\le\kappa\B^{-1} \), and hence \( \A\le\kappa\B^{-1} \).

For the attainment, with \( \x=(\mathbf{q}_1+\mathbf{q}_n)/\sqrt2 \) we have \( \x\T\x=1 \),
\( \x\T\V\x=(\lambda_1+\lambda_n)/2 \) and \( \x\T\V^{-1}\x=(1/\lambda_1+1/\lambda_n)/2 \). The ratio is
\( 1/\bigl((\x\T\V\x)(\x\T\V^{-1}\x)\bigr)=4\lambda_1\lambda_n/(\lambda_1+\lambda_n)^2 \).
:::

The scalar inequality behind the proof,
\( (\x\T\V\x)(\x\T\V^{-1}\x)\le\kappa(\x\T\x)^2 \), is Kantorovich's inequality, and the proof above is a matrix form of it. Bloomfield and Watson (1975)
and Knott (1975) found the corresponding sharp bound when efficiency is measured by the determinant
of the whole covariance matrix.

The bound is a worst case over all designs. For particular designs ordinary least squares can be fully efficient.
Kruskal's theorem (@thm-proj-kruskal) says that \( \hbeta \) *equals* the generalized least squares estimate for every
\( \Y \) iff \( \C(\V\X)\subseteq\C(\X) \), and then the efficiency is one. Even then the standard errors may be wrong.

::: {#exm-dep-equicorrelated}
[Equicorrelated errors]

Let \( \V=(1-\rho)\I+\rho\bone\bone\T \), with \( -1/(n-1)<\rho<1 \) so that \( \V \) is positive definite, and let the model
contain an intercept, \( \X=[\bone,\X_1] \). It arises when all observations share one random
shock. Since \( \V\X=(1-\rho)\X+\rho\bone(\bone\T\X) \) has columns in \( \C(\X) \),
Kruskal's condition holds and \( \hbeta \) is the best linear unbiased estimator. Its covariance is not
\( \sigma^2(\X\T\X)^{-1} \), however. Because \( \X\mathbf{e}_1=\bone \) for the first standard basis vector \( \mathbf{e}_1 \), we have
\( (\X\T\X)^{-1}\X\T\bone=\mathbf{e}_1 \), and @thm-dep-covariance(a) gives
\[
\Cov(\hbeta)=\sigma^2(1-\rho)(\X\T\X)^{-1}+\sigma^2\rho\,\mathbf{e}_1\mathbf{e}_1\T .
\]
By @exr-rv-rss-correlated, \( \E(s^2)=\sigma^2(1-\rho) \). So \( s^2(\X\T\X)^{-1} \) is exactly unbiased for the
covariance of the slopes, whose information comes from comparisons in which the common shock cancels.
For the intercept it misses the term \( \sigma^2\rho \). When \( \rho>0 \) this term does not shrink as \( n \) grows,
so the intercept is not even consistent, and its reported standard error goes to zero anyway (@exr-dep-equicorrelated-intercept). The same holds for any fitted mean \( \x_0\T\hbeta \), since
\( \x_0\T\mathbf{e}_1=1 \) whenever \( \x_0 \) has a leading \( 1 \).
:::

## Autocorrelated errors

::: {#exm-dep-ar1}
[A trend with autocorrelated errors]

Take \( n=20 \) equally spaced times and errors with \( v_{st}=\rho^{\lvert s-t\rvert} \), the covariance of a stationary
first-order autoregression, with \( \rho=0.6 \). Fit a straight line in two designs, a linear trend
\( x_t=t \) and an alternating regressor \( x_t=(-1)^t \). For the trend,
the expected usual variance of the slope is \( 0.256 \) times its true variance, so the
root of the expected usual variance is \( 0.506 \) of the true standard error. Yet ordinary least squares is
nearly efficient: generalized least squares has \( 0.905 \) of its variance. For the
alternating regressor the error goes the other way. The ratio is \( 3.271 \), so the standard
error is too large by a factor \( 1.809 \), and the efficiency is
\( 0.950 \). Swindel's bounds for this \( \V \) are \( 0.193 \) and
\( 4.31 \), so both designs lie well towards the extremes. The Kantorovich bound,
\( 0.2351 \) from \( \lambda_1=3.757 \) and
\( \lambda_n=0.2514 \), is far below either efficiency.

In \( 100000 \) simulated series with normal errors, the nominal \( 95\% \) \( t \) interval for the slope
covered the true value in \( 0.693 \) of them for the trend and \( 0.9976 \) for the
alternating regressor. [Figure 19.2.1](#fig-dep-covariance) shows the whole range of \( \rho \).
:::

```{.python .run #cell-wrong-covariance-ratio}
import numpy as np
n = 20
t = np.arange(1, n + 1)
designs = {"trend": t.astype(float), "alternating": (-1.0) ** t}

def ar1(rho):
    return rho ** np.abs(t[:, None] - t[None, :])

def slope_summary(x, V):
    """True variance, mean of the usual estimate, GLS variance of the slope (sigma = 1)."""
    X = np.column_stack([np.ones(n), x])
    G = np.linalg.inv(X.T @ X)
    A = G @ X.T                                     # beta_hat = A y
    true_var = (A @ V @ A.T)[1, 1]                  # the sandwich
    M = X @ A
    Es2 = np.trace((np.eye(n) - M) @ V) / (n - 2)   # E(s^2)
    usual = Es2 * G[1, 1]                           # E of s^2 [(X'X)^{-1}]_{11}
    gls = np.linalg.inv(X.T @ np.linalg.solve(V, X))[1, 1]
    return true_var, usual, gls

rho = 0.6
for name, x in designs.items():
    true_var, usual, gls = slope_summary(x, ar1(rho))
    print(f"{name:12s} usual/true = {usual / true_var:.3f}   GLS/OLS variance = {gls / true_var:.3f}")
```

```{.python .run #cell-wrong-covariance-coverage}
from scipy import stats
rng = np.random.default_rng(1902)
reps = 100_000
L = np.linalg.cholesky(ar1(rho))
E = rng.normal(size=(reps, n)) @ L.T              # rows are AR(1) error vectors
tq = stats.t.ppf(0.975, n - 2)
for name, x in designs.items():
    X = np.column_stack([np.ones(n), x])
    G = np.linalg.inv(X.T @ X)
    B = E @ (G @ X.T).T                             # beta_hat - beta, one row per data set
    s2 = np.sum((E - B @ X.T) ** 2, axis=1) / (n - 2)
    cover = np.mean(np.abs(B[:, 1]) <= tq * np.sqrt(s2 * G[1, 1]))
    print(f"{name:12s} coverage of the nominal 95% interval: {cover:.4f}")
```

::: {when-format="html"}
![**Figure 19.2.1.** Straight-line regression with \( n=20 \) and first-order autoregressive errors.
(a) The square root of the ratio of the expected usual variance of the slope to its true variance, for a trend
and an alternating regressor, with Swindel's bounds (shaded). (b) The efficiency of ordinary least squares
relative to generalized least squares, with the Kantorovich bound of @thm-dep-efficiency.](wrong_covariance.svg){#fig-dep-covariance width=100%}
:::

::: {when-format="pdf"}
![Straight-line regression with \( n=20 \) and first-order autoregressive errors.
(a) The square root of the ratio of the expected usual variance of the slope to its true variance, for a trend
and an alternating regressor, with Swindel's bounds (shaded). (b) The efficiency of ordinary least squares
relative to generalized least squares, with the Kantorovich bound of @thm-dep-efficiency.](wrong_covariance.pdf){width=100%}
:::

@eq-dep-se-ratio explains the pattern. A smooth regressor puts \( \bu \) along slowly varying directions, where
positively correlated errors have their largest variance; an alternating one puts it where positive correlation
cancels. Most time-series regressors are smooth, which is why positive autocorrelation usually makes ordinary
standard errors too small.

::: {.idea}
A wrong covariance matrix is mainly a problem for *standard errors*, not *estimates*: \( s^2(\X\T\X)^{-1} \) can be wrong by a large factor in either direction, and
nothing in the fitted model reveals this. [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) shows how to estimate the correct covariance @eq-lm-sandwich without knowing
\( \V \), and how to detect serial correlation. [Chapter 31](../ch31-general-gauss-markov/index.html) treats generalized least squares when \( \V \) is known up to a few
parameters (@def-ggm-feasible).
:::

## Exercises

### A. Check your understanding

::: {#exr-dep-diagonal-s2}
[A1]

Let \( \V=\diag(v_1,\dots,v_n) \). Show that \( \E(s^2)=\sigma^2\sum_i(1-h_{ii})v_i/(n-p) \), a weighted average of the
variances \( \sigma^2v_i \).
:::

::: {#exr-dep-swindel-trivial}
[A2]

Show that the bounds of @thm-dep-covariance(d) do not change when \( \V \) is replaced by \( c\V \) for \( c>0 \), and that
both equal \( 1 \) iff \( \V \) is a multiple of \( \I \).
:::

### B. Practice

::: {#exr-dep-hetero-origin}
[B1]

For the regression through the origin \( \E(Y_i)=\beta x_i \) with \( \Var(Y_i)=\sigma^2v_i \), show that
\( \Var(\hat{\beta})=\sigma^2\sum_ih_iv_i/\sum_ix_i^2 \) and \( \E(v)=\sigma^2\sum_i(1-h_i)v_i/\bigl((n-1)\sum_ix_i^2\bigr) \),
where \( h_i=x_i^2/\sum_jx_j^2 \) are the leverages. Conclude that the usual standard error is too small on average iff the
leverage-weighted average of the \( v_i \) exceeds the \( (1-h_i) \)-weighted average.
:::

::: {.solution}
\( \hat{\beta}=\sum_ix_iY_i/S \) with \( S=\sum_jx_j^2 \), so \( \Var(\hat{\beta})=\sigma^2\sum_ix_i^2v_i/S^2=\sigma^2\sum_ih_iv_i/S \).
@exr-dep-diagonal-s2 with \( p=1 \) gives \( \E(s^2) \), and \( v=s^2/S \). The weights \( h_i \) sum to \( 1 \), and so do
\( (1-h_i)/(n-1) \). The usual variance is too small on average iff \( \sum_ih_iv_i>\sum_i(1-h_i)v_i/(n-1) \). In
particular, when the error variance grows with \( \lvert x_i\rvert \), the points that dominate \( \hat{\beta} \) are the noisiest,
and the usual standard error understates the truth.
:::

::: {#exr-dep-equicorrelated-intercept}
[B2]

In @exm-dep-equicorrelated with \( \rho>0 \), show that \( \Var(\hat{\beta}_0)\ge\sigma^2\rho \) for every \( n \), while
\( \E(v)\to0 \) for the intercept when \( \bigl[(\X\T\X)^{-1}\bigr]_{11}\to0 \). Explain the result by writing
\( \varepsilon_i=\sqrt{\rho}\,\sigma Z+\sqrt{1-\rho}\,\sigma U_i \) with independent standardized \( Z,U_1,\dots,U_n \).
:::

::: {.solution}
From the covariance in the example, \( \Var(\hat{\beta}_0)=\sigma^2(1-\rho)\bigl[(\X\T\X)^{-1}\bigr]_{11}+\sigma^2\rho\ge\sigma^2\rho \),
while \( \E(v)=\sigma^2(1-\rho)\bigl[(\X\T\X)^{-1}\bigr]_{11} \). The representation has covariance
\( \sigma^2\rho\bone\bone\T+\sigma^2(1-\rho)\I=\sigma^2\V \). The shared term \( \sqrt{\rho}\sigma Z\bone \) lies in \( \C(\X) \) and is absorbed
entirely by the intercept, since \( (\X\T\X)^{-1}\X\T\bone=\mathbf{e}_1 \). It never averages out, and it leaves no trace in the
residuals, so \( s^2 \) cannot see it.
:::

### C. Going deeper

::: {#exr-dep-swindel-attained}
[C1]

Show that both bounds of @thm-dep-covariance(d) are attained. *Hint:* take \( \C(\X) \) spanned by eigenvectors of \( \V \),
and choose \( \mathbf{a} \) so that \( \bu \) is one of them.
:::

::: {.solution}
Let \( \mathbf{q}_1,\dots,\mathbf{q}_n \) be orthonormal eigenvectors for \( \lambda_1\ge\dots\ge\lambda_n \). If
\( \X=[\mathbf{q}_1,\dots,\mathbf{q}_p] \), then \( \I-\M=\sum_{i>p}\mathbf{q}_i\mathbf{q}_i\T \) and
\( \tr\bigl((\I-\M)\V\bigr)=\sum_{i>p}\lambda_i=(n-p)\bar{\lambda}_- \). With \( \mathbf{a}=\mathbf{e}_1 \), \( \bu=\mathbf{q}_1 \) and
\( \bu\T\V\bu=\lambda_1 \), so the ratio @eq-dep-se-ratio equals \( \bar{\lambda}_-/\lambda_1 \). Symmetrically,
\( \X=[\mathbf{q}_{n-p+1},\dots,\mathbf{q}_n] \) and \( \bu=\mathbf{q}_n \) give \( \bar{\lambda}_+/\lambda_n \).
:::
