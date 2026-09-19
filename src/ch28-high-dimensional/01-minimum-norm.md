# Minimum-norm interpolation

Gene expression studies measure tens of thousands of genes on a few hundred patients, and a
model with many splines, interactions or random features easily has more columns than rows. This
section asks what least squares does when \( p>n \). It does not fail, but it no longer picks out one
answer, and we study the answer that the SVD and the common algorithms return.

## Least squares when p exceeds n

Suppose \( \rank(\X)=n<p \), the usual case (a matrix with continuously distributed entries has full
rank with probability one). Three consequences follow from Parts I and II.

First, \( \C(\X)=\Real^n \), so \( \M=\I \): every least squares fit reproduces the data, the residuals are zero,
\( R^2=1 \), and the estimator \( s^2 \) of @thm-lm-sigma2 has \( n-r=0 \) in its denominator. The data contain no
internal information about \( \sigma^2 \).

Second, \( \Null(\X) \) has dimension \( p-n\ge1 \), and the least squares estimates form the affine set
\( \X^+\y+\Null(\X) \) (@thm-cmp-svd-ls). Coefficient vectors that differ by an element of \( \Null(\X) \) give the same
mean, so the data cannot tell them apart.

Third, \( \blambda\T\bbeta \) is estimable iff \( \blambda\in\C(\X\T) \) (@thm-est-estimable-identifiable), and a coordinate
vector \( \mathbf{e}_j \) lies in the \( n \)-dimensional row space only if it is orthogonal to \( \Null(\X) \), which fails for
generic designs. For a generic \( \X \) with \( p>n \), *no single coefficient is estimable*. Any method that
answers questions about individual coefficients is adding information, explicitly or through its algorithm.

::: {#exm-hd-two-interpolators}
[Two interpolators]

Take \( n=4 \) observations on \( p=6 \) regressors, with the design of the listing below and
\( \y=(1,2,0,-1)\T \). The minimum-norm solution \( \hbeta^+=\X^+\y \) (@prp-mat-min-norm, @prp-proj-min-norm) has length
\( 1.322 \). Adding twice a unit vector of \( \Null(\X) \) gives a second solution of
length \( 2.398 \). Both reproduce \( \y \) exactly. They disagree about the coefficients: the
first, for instance, is \( -0.002 \) in one and \( -0.059 \) in the other.
They also disagree about predictions at a point that was not observed: at \( \x_0=\bone \)
they predict \( 2.131 \) and \( 2.588 \). The data give no reason to prefer
either. The minimum-norm solution is singled out by a choice, not by the data.
:::

```{.python .run #cell-min-norm-small}
import numpy as np

rng = np.random.default_rng(2801)
n, p = 4, 6
X = np.round(rng.normal(size=(n, p)), 1)            # a 4 x 6 design, rank 4
y = np.array([1.0, 2.0, 0.0, -1.0])

b_plus = X.T @ np.linalg.solve(X @ X.T, y)          # X^T (X X^T)^{-1} y
print("min-norm solution ", np.round(b_plus, 3))
print("residual norm     ", np.linalg.norm(y - X @ b_plus))

# another interpolator: add a null-space vector of X
_, _, Vt = np.linalg.svd(X)
b_other = b_plus + 2.0 * Vt[-1]                     # Vt[-1] spans part of N(X)
print("other solution    ", np.round(b_other, 3))
print("residual norm     ", np.linalg.norm(y - X @ b_other))
print("norms", np.linalg.norm(b_plus), np.linalg.norm(b_other))

x0 = np.ones(p)                                      # a new point
print("predictions at x0:", x0 @ b_plus, x0 @ b_other)
```

## The minimum-norm solution

In this chapter the singular values of \( \X \) are written \( d_1\ge\dots\ge d_r>0 \), so that
\( \sigma \) stays free for the error standard deviation. By @thm-mat-svd,
\[
\X=\sum_{i=1}^r d_i\,\bu_i\bv_i\T,\qquad \X^+=\sum_{i=1}^r d_i^{-1}\,\bv_i\bu_i\T ,
\]
with orthonormal \( \bu_1,\dots,\bu_r\in\Real^n \) and \( \bv_1,\dots,\bv_r\in\Real^p \). The matrix
\[
\bP=\X^+\X=\sum_{i=1}^r\bv_i\bv_i\T
\]{#eq-hd-row-projection}

is the orthogonal projection onto the row space \( \C(\X\T)=\spn(\bv_1,\dots,\bv_r) \). The next
theorem collects what is known about \( \hbeta^+=\X^+\y \) when \( p>n \). Part (a) restates earlier
results in the form needed here. Parts (b) and (c) show that \( \hbeta^+ \) is what two common
procedures produce when they are pushed to the limit. Part (d) gives its bias and variance.

::: {#thm-hd-min-norm}
[The minimum-norm interpolator]

Let \( \X \) be \( n\times p \) with rank \( r \) and singular values \( d_1\ge\dots\ge d_r>0 \), and let
\( \hbeta^+=\X^+\y \).

::: {.enumerate options="label=(\alph*)"}
1. *(Characterization.)* The least squares estimates are the vectors \( \hbeta^++\bv \) with
   \( \bv\in\Null(\X) \). Among them \( \hbeta^+ \) is the only one in \( \C(\X\T) \), and it has the
   smallest Euclidean norm. If \( r=n \), then \( \hbeta^+=\X\T(\X\X\T)^{-1}\y \) and every least squares
   estimate interpolates: \( \X\hbeta^+=\y \).

2. *(Ridgeless limit.)* For \( \lambda>0 \) let \( \hbeta_\lambda=(\X\T\X+\lambda\I)^{-1}\X\T\y \) be the
   ridge estimate. Then
   \( \norm{\hbeta_\lambda-\hbeta^+}\le(\lambda/d_r^2)\,\norm{\hbeta^+} \), so \( \hbeta_\lambda\to\hbeta^+ \) as
   \( \lambda\downarrow0 \).

3. *(Gradient descent.)* Let \( \bb_0\in\Real^p \), \( 0<\eta<2/d_1^2 \) and
   \( \bb_{t+1}=\bb_t+\eta\,\X\T(\y-\X\bb_t) \). Then \( \bb_t\to\hbeta^++(\I-\bP)\bb_0 \). In particular,
   gradient descent started at \( \bb_0=\bzero \) converges to \( \hbeta^+ \).

4. *(Bias and variance.)* If \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\I \), then
   \( \E(\X^+\Y)=\bP\bbeta \) and \( \Cov(\X^+\Y)=\sigma^2\sum_{i\le r}d_i^{-2}\bv_i\bv_i\T \). For every
   \( \mathbf{a}\in\Real^p \),
   \[
\E\bigl(\mathbf{a}\T\X^+\Y-\mathbf{a}\T\bbeta\bigr)^2
=\bigl(\mathbf{a}\T(\I-\bP)\bbeta\bigr)^2+\sigma^2\sum_{i=1}^r\frac{(\mathbf{a}\T\bv_i)^2}{d_i^2},
\]{#eq-hd-min-norm-mse}

   and \( \E\norm{\X^+\Y-\bbeta}^2=\norm{(\I-\bP)\bbeta}^2+\sigma^2\sum_{i\le r}d_i^{-2} \).
:::

:::

::: {.proof}
(a) The first two statements are @thm-cmp-svd-ls, and @prp-mat-min-norm applied to the normal
equations of @thm-proj-normal-equations; \( \hbeta^+=\sum_i(\bu_i\T\y/d_i)\bv_i \) lies in \( \spn(\bv_1,\dots,\bv_r)=\C(\X\T) \). If
\( r=n \), then \( \X\X\T \) is a nonsingular \( n\times n \) matrix, and \( \X\T(\X\X\T)^{-1} \) satisfies the
four conditions of @def-mat-moore-penrose, so it equals \( \X^+ \). Moreover \( \C(\X)=\Real^n \), so the
projection onto \( \C(\X) \) is \( \I \) and every least squares fit equals \( \y \).

(b) Write \( c_i=\bu_i\T\y \). By @thm-shr-ridge(a),
\( \hbeta_\lambda=\sum_{i\le r}d_ic_i\bv_i/(d_i^2+\lambda) \), and @thm-shr-ridge(d) gives the limit
\( \hbeta^+=\sum_{i\le r}(c_i/d_i)\bv_i \); only the rate is new. Subtracting,
\( \hbeta^+-\hbeta_\lambda=\sum_{i\le r}\{\lambda/(d_i^2+\lambda)\}(c_i/d_i)\bv_i \), so by orthonormality
\( \norm{\hbeta^+-\hbeta_\lambda}^2\le(\lambda/d_r^2)^2\sum_i(c_i/d_i)^2=(\lambda/d_r^2)^2\norm{\hbeta^+}^2 \).

(c) Complete \( \bv_1,\dots,\bv_r \) to an orthonormal basis of \( \Real^p \) with vectors
\( \bv_{r+1},\dots,\bv_p \) spanning \( \Null(\X) \). Since \( \X\T(\y-\X\bb)=\sum_{i\le r}d_i(c_i-d_i\bv_i\T\bb)\bv_i \),
the coordinates of \( \bb_t \) evolve separately:
\[
\bv_i\T\bb_{t+1}=(1-\eta d_i^2)\,\bv_i\T\bb_t+\eta d_ic_i\quad(i\le r),\qquad
\bv_i\T\bb_{t+1}=\bv_i\T\bb_t\quad(i>r).
\]
For \( i\le r \), subtracting the fixed point \( c_i/d_i \) gives
\( \bv_i\T\bb_t-c_i/d_i=(1-\eta d_i^2)^t(\bv_i\T\bb_0-c_i/d_i)\to0 \), because \( 0<\eta d_i^2<2 \). The
coordinates with \( i>r \) never move. So \( \bb_t\to\sum_{i\le r}(c_i/d_i)\bv_i+\sum_{i>r}(\bv_i\T\bb_0)\bv_i
=\hbeta^++(\I-\bP)\bb_0 \).

(d) \( \E(\X^+\Y)=\X^+\X\bbeta=\bP\bbeta \), and
\( \Cov(\X^+\Y)=\sigma^2\X^+(\X^+)\T=\sigma^2\sum_id_i^{-2}\bv_i\bv_i\T \) (@thm-rv-linear). The mean
squared error of \( \mathbf{a}\T\X^+\Y \) is its squared bias \( \bigl(\mathbf{a}\T(\bP-\I)\bbeta\bigr)^2 \) plus its variance,
which is @eq-hd-min-norm-mse. Summing @eq-hd-min-norm-mse over \( \mathbf{a}=\mathbf{e}_1,\dots,\mathbf{e}_p \) gives
the last formula.
:::

Parts (b) and (c) explain why \( \hbeta^+ \) matters although nobody chooses it on purpose. Gradient
iterates started at zero never leave \( \C(\X\T) \), because every gradient \( \X\T(\y-\X\bb) \) lies in it.
Algorithms for large models are typically started near zero and run until the training error
vanishes, and in the linear case they return the minimum-norm interpolator. This is called the
*implicit regularization* of the algorithm.

```{.python .run #cell-min-norm-ridgeless}
for lam in [1.0, 0.1, 0.01, 0.001]:
    b_ridge = np.linalg.solve(X.T @ X + lam * np.eye(p), X.T @ y)
    print(f"lambda = {lam:6.3f}   distance to min-norm {np.linalg.norm(b_ridge - b_plus):.2e}")

# gradient descent on ||y - X b||^2 / 2 started at zero
step = 1.0 / np.linalg.norm(X, 2) ** 2
b = np.zeros(p)
for t in range(5000):
    b = b + step * X.T @ (y - X @ b)
print("gradient descent: distance to min-norm", np.linalg.norm(b - b_plus))
```

For the design of @exm-hd-two-interpolators the ridge estimate at \( \lambda=0.01 \) is within
\( 0.0132 \) of \( \hbeta^+ \), and the distance shrinks in proportion to \( \lambda \), as part (b)
predicts.

In part (d), the bias \( (\I-\bP)\bbeta \) is the component of \( \bbeta \) in the null space, which the data
cannot see and the minimum-norm solution sets to zero. The variance is the familiar \( \sigma^2/d_i^2 \) of
@thm-cmp-svd-ls, summed over the directions the data do see. So \( \mathbf{a}\T\hbeta^+ \) is unbiased for all
\( \bbeta \) exactly when \( \mathbf{a}\T\bbeta \) is estimable (@exr-hd-unbiased-functions). For any other \( \mathbf{a} \), its bias
depends on the invisible part of \( \bbeta \) and can be anything.

## The risk with random regressors

To say whether interpolation predicts well, the design must be random, since prediction is at
new regressor values. The cleanest case has independent standard normal regressors, where
both terms of part (d) can be averaged exactly.

::: {#prp-hd-isotropic-risk}
[Risk under isotropic regressors]

Let the rows of \( \X \) be independent \( \Normal_p(\bzero,\I_p) \) vectors, let
\( \Y=\X\bbeta+\be \) with \( \be\sim\Normal_n(\bzero,\sigma^2\I) \) independent of \( \X \), and let
\( \x_0\sim\Normal_p(\bzero,\I_p) \) be independent of both. The prediction risk
\( \mathcal{R}=\E(\x_0\T\X^+\Y-\x_0\T\bbeta)^2 \) is
\[
\mathcal{R}=
\begin{cases}
\sigma^2\,\dfrac{p}{n-p-1}, & p\le n-2,\\[2ex]
\norm{\bbeta}^2\Bigl(1-\dfrac np\Bigr)+\sigma^2\,\dfrac{n}{p-n-1}, & p\ge n+2 .
\end{cases}
\]{#eq-hd-isotropic-risk}

:::

::: {.proof}
Given \( \X \) and \( \Y \), the vector \( \hbeta^+=\X^+\Y \) is fixed and
\( \E\{(\x_0\T(\hbeta^+-\bbeta))^2\mid\X,\Y\}=\norm{\hbeta^+-\bbeta}^2 \), because \( \Cov(\x_0)=\I \). So
\( \mathcal{R}=\E\norm{\hbeta^+-\bbeta}^2 \), and by @thm-hd-min-norm(d) applied conditionally on \( \X \),
\[
\mathcal{R}=\E\norm{(\I-\bP)\bbeta}^2+\sigma^2\,\E\sum_{i=1}^r d_i^{-2}.
\]

*Case \( p\le n-2 \).* \( \X\T\X=\sum_i\x_{(i)}\x_{(i)}\T \) has the Wishart law \( W_p(n,\I) \) (@def-cor-wishart),
so by @lem-cor-inverse-wishart it is nonsingular with probability one and
\( \E(\X\T\X)^{-1}=\I/(n-p-1) \). Then \( r=p \), \( \bP=\I \), the bias term vanishes, and
\( \sum_id_i^{-2}=\tr(\X\T\X)^{-1} \) has mean \( p/(n-p-1) \).

*Case \( p\ge n+2 \).* The \( p \) columns of \( \X \) are independent \( \Normal_n(\bzero,\I_n) \) vectors, so
\( \X\X\T\sim W_n(p,\I_n) \). By @lem-cor-inverse-wishart again, \( \X\X\T \) is nonsingular with probability
one, so \( r=n \). The numbers \( d_i^2 \) are also the eigenvalues of \( \X\X\T \), so
\( \sum_id_i^{-2}=\tr(\X\X\T)^{-1} \), with mean \( n/(p-n-1) \). For the bias, let \( \Q \) be any
orthogonal \( p\times p \) matrix. The rows of \( \X\Q \) are independent \( \Normal_p(\bzero,\Q\T\Q)=\Normal_p(\bzero,\I) \)
vectors (@thm-mvn-linear), so \( \X\Q \) has the same law as \( \X \). One checks from
@def-mat-moore-penrose that \( (\X\Q)^+=\Q\T\X^+ \), so the row-space projection of \( \X\Q \) is \( \Q\T\bP\Q \).
Hence \( \bar{\bP}=\E\bP \) satisfies \( \bar{\bP}=\Q\T\bar{\bP}\Q \) for every orthogonal \( \Q \). Then
\( \bu\T\bar{\bP}\bu \) takes the same value \( c \) at every unit vector \( \bu \), since any unit vector is \( \Q\mathbf{e}_1 \)
for some \( \Q \). A symmetric matrix whose quadratic form is constant on the unit sphere is \( c\I \), by @thm-mat-extremal-rayleigh.
Taking traces, \( pc=\E\tr\bP=\E r=n \). Therefore
\( \E\norm{(\I-\bP)\bbeta}^2=\bbeta\T(\I-\E\bP)\bbeta=\norm{\bbeta}^2(1-n/p) \).
:::

Below \( p=n \) the fit is unbiased and its variance grows like \( p/(n-p-1) \), which explodes as
\( p \) approaches \( n \). Above \( p=n \) the variance *falls* again, like \( n/(p-n-1) \). The
fit has to interpolate \( n \) noisy values, and with many more directions available the
interpolator can spread the noise over them thinly. The price is the bias: the fit recovers only
the fraction \( n/p \) of \( \norm{\bbeta}^2 \). The formulas do not cover \( |p-n|\le1 \).

## Double descent

In applications \( p \) is usually a choice, such as the number of basis functions. To compare choices of
\( p \) on one problem, let the response depend on a pool of \( D \) features and let the fit use only the
first \( p \). The features left out act as extra noise.

::: {#cor-hd-double-descent}
[Fitting the first p of D features]

Let \( \x\sim\Normal_D(\bzero,\I_D) \), \( Y=\x\T\bbeta+\varepsilon \) with \( \varepsilon\sim\Normal(0,\sigma^2) \)
independent of \( \x \), and let \( n \) independent copies be observed. For \( p\le D \), split
\( \bbeta=(\bbeta_{(1)}\T,\bbeta_{(2)}\T)\T \) into its first \( p \) and last \( D-p \) entries, write
\( \omega^2=\norm{\bbeta_{(2)}}^2 \), fit the minimum-norm least squares estimate \( \hbeta^+ \) (@prp-proj-min-norm) to the
first \( p \) features, and predict with \( (\hbeta^{+\T},\bzero\T)\T \). With \( \sigma_p^2=\sigma^2+\omega^2 \), the excess
prediction risk at an independent new \( \x_0 \) is \( \omega^2 \) plus the right-hand side of
@eq-hd-isotropic-risk with \( \bbeta_{(1)} \) in place of \( \bbeta \) and \( \sigma_p^2 \) in place of \( \sigma^2 \).
:::

::: {.proof}
Split \( \X=[\X_1,\X_2] \) in the same way. Then \( \Y=\X_1\bbeta_{(1)}+(\X_2\bbeta_{(2)}+\be) \). The entries of
\( \X_2 \) are independent of \( \X_1 \), so the bracket is a \( \Normal_n(\bzero,\sigma_p^2\I) \) vector
independent of \( \X_1 \). @prp-hd-isotropic-risk applies with \( \X_1 \), \( \bbeta_{(1)} \) and \( \sigma_p^2 \),
and gives \( \E\norm{\hbeta^+-\bbeta_{(1)}}^2 \). At the new point, the prediction error
\( \x_{0,1}\T(\hbeta^+-\bbeta_{(1)})-\x_{0,2}\T\bbeta_{(2)} \) is a sum of two uncorrelated terms with mean zero,
so its mean square is \( \E\norm{\hbeta^+-\bbeta_{(1)}}^2+\omega^2 \).
:::

::: {#exm-hd-double-descent}
[Two shapes of signal]

Take \( n=50 \) observations, a pool of \( D=200 \) features, \( \sigma^2=0.25 \) and
\( \norm{\bbeta}^2=1 \). [Figure 28.1.1](#fig-hd-double-descent) plots the exact risk of
@cor-hd-double-descent against \( p \) for two coefficient vectors, with simulated values as a
check. It also shows the risk of ridge regression (@thm-shr-ridge) with the best penalty for each \( p \).

In panel (a) all \( 200 \) coefficients are equal, so each feature left out adds to the noise. No model
with \( p<n \) beats predicting zero, which has risk \( 1.000 \). The risk climbs to \( 12.31 \) at \( p=45 \), is
\( 12.94 \) at \( p=55 \), and then falls. The overall minimum is at \( p=200 \):
\( 0.834 \). Here the best prediction comes from interpolating with every feature.

In panel (b) the coefficients decay geometrically, and the classical picture holds: a small model with
\( p=9 \) has risk \( 0.078 \), and nothing to the right of the peak comes close. At \( p=D \) the risk is
\( 0.834 \) again. When all features are used, @eq-hd-isotropic-risk depends on \( \bbeta \)
only through \( \norm{\bbeta} \).
:::

::: {when-format="html"}
![**Figure 28.1.1.** Excess prediction risk of the minimum-norm least squares fit on the first
\( p \) of \( D=200 \) isotropic features, with \( n=50 \) (vertical line) and \( \sigma^2=0.25 \):
the exact formula of @cor-hd-double-descent (line) and simulation (dots), and ridge regression
with the best penalty for each \( p \). The dotted line is the risk of predicting zero. (a) Equal
coefficients. (b) Geometrically decaying coefficients.](double_descent.svg){#fig-hd-double-descent width=100%}
:::

::: {when-format="pdf"}
![Excess prediction risk of the minimum-norm least squares fit on the first
\( p \) of \( D=200 \) isotropic features, with \( n=50 \) (vertical line) and \( \sigma^2=0.25 \):
the exact formula of @cor-hd-double-descent (line) and simulation (dots), and ridge regression
with the best penalty for each \( p \). The dotted line is the risk of predicting zero. (a) Equal
coefficients. (b) Geometrically decaying coefficients.](double_descent.pdf){width=100%}
:::

A risk curve that falls, peaks at \( p=n \) and falls again is called **double descent**. Three
cautions keep it in proportion.

- *The peak belongs to the estimator, not to the problem.* It comes from the variance term
  \( \sigma^2\sum_id_i^{-2} \), large when the smallest nonzero singular value is near zero, as for \( p\approx n \).
  Ridge with the best penalty has no peak, and its best values, \( 0.811 \) and \( 0.073 \), are
  slightly below the best minimum-norm values, as @thm-hd-min-norm(b) requires.
- *The second descent need not win.* That depends on how the signal is spread, as the panels show.
- *Isotropic features are a special case.* Hastie, Montanari, Rosset and Tibshirani (2022) treat general
  covariances as \( n,p\to\infty \) with \( p/n \) fixed, and Bartlett, Long, Lugosi and Tsigler (2020) find when
  interpolation is *benign*, with risk still tending to that of the best predictor.

```{.python .run #cell-min-norm-risk-small}
import numpy as np

def exact_risk(p, n, beta, sigma2):
    """Excess prediction risk of min-norm least squares on the first p of D isotropic features."""
    kept, left = np.sum(beta[:p] ** 2), np.sum(beta[p:] ** 2)
    noise = sigma2 + left                           # omitted features act as extra noise
    if p == 0:
        return left
    if p <= n - 2:
        return left + noise * p / (n - p - 1)
    if p >= n + 2:
        return left + kept * (1 - n / p) + noise * n / (p - n - 1)
    return np.inf                                   # |p - n| <= 1: not covered by the formula

def simulated_risk(p, n, beta, sigma2, reps, rng, lams=()):
    """Average over designs of the conditional risk given X; also ridge for each penalty."""
    bS, left = beta[:p], np.sum(beta[p:] ** 2)
    noise = sigma2 + left
    out, ridge = [], np.zeros(len(lams))
    for _ in range(reps):
        X = rng.normal(size=(n, p))
        _, d, Vt = np.linalg.svd(X, full_matrices=False)
        c = Vt @ bS                                 # coordinates of beta_S in the row space
        bias2 = bS @ bS - c @ c                     # ||(I - P) beta_S||^2
        out.append(left + bias2 + noise * np.sum(1 / d ** 2))
        if len(lams):
            shrink = d ** 2 / (d ** 2 + np.asarray(lams)[:, None])     # one row per penalty
            ridge += left + bias2 + np.sum(((1 - shrink) * c) ** 2, axis=1) \
                + noise * np.sum(shrink ** 2 / d ** 2, axis=1)
    return np.mean(out), np.std(out) / np.sqrt(reps), ridge / reps

n, D, sigma2 = 20, 100, 0.25
beta = np.full(D, 1 / np.sqrt(D))                  # signal spread evenly, ||beta||^2 = 1
rng = np.random.default_rng(2802)
for p in [5, 10, 15, 25, 40, 100]:
    sim, se, _ = simulated_risk(p, n, beta, sigma2, reps=200, rng=rng)
    print(f"p = {p:3d}: formula {exact_risk(p, n, beta, sigma2):.3f}   simulation {sim:.3f}")
```

For the rest of the chapter we take a different route. Instead of letting the algorithm decide
which of the many solutions to report, we *assume* that the truth is simple: only a few
coefficients are nonzero. That assumption is strong, but it is explicit, and it can be exploited
through the lasso.

## Exercises

### A. Check your understanding

::: {#exr-hd-interpolation-facts}
[A1]

Let \( \X \) be \( n\times p \) with an intercept column and \( \rank(\X)=n<p \). Show that every least
squares fit has \( R^2=1 \), and explain why neither the \( F \) test of @thm-glh-f-test nor the
residual-based diagnostics of [Chapter 20](../ch20-residuals-leverage-influence/index.html) can be used.
What are the leverages \( h_{ii} \)?
:::

::: {#exr-hd-unbiased-functions}
[A2]

In the setting of @thm-hd-min-norm(d), show that \( \mathbf{a}\T\X^+\Y \) is unbiased for \( \mathbf{a}\T\bbeta \) for
every \( \bbeta \) iff \( \mathbf{a}\in\C(\X\T) \). Show also that for such \( \mathbf{a} \) its variance is
\( \sigma^2\mathbf{a}\T(\X\T\X)^+\mathbf{a} \), and that it equals the least squares estimator \( \mathbf{a}\T\hbeta \) of
@thm-proj-invariant-functions for any least squares estimate \( \hbeta \).
:::

::: {.solution}
The bias is \( \mathbf{a}\T(\bP-\I)\bbeta \), which vanishes for all \( \bbeta \) iff \( (\I-\bP)\mathbf{a}=\bzero \), that
is, iff \( \mathbf{a}\in\C(\bP)=\C(\X\T) \). The variance is \( \sigma^2\sum_i(\mathbf{a}\T\bv_i)^2/d_i^2 \), and
\( \sum_id_i^{-2}\bv_i\bv_i\T \) is the Moore–Penrose inverse of \( \X\T\X=\sum_id_i^2\bv_i\bv_i\T \) (@thm-mat-svd).
Finally, any least squares estimate is \( \hbeta^++\bv \) with \( \bv\in\Null(\X) \), and
\( \mathbf{a}\T\bv=0 \) because \( \mathbf{a}\in\C(\X\T)=\Null(\X)\perpc \), so \( \mathbf{a}\T\hbeta=\mathbf{a}\T\hbeta^+ \).
:::

### B. Practice

::: {#exr-hd-gd-start}
[B1]

Show directly from the update in @thm-hd-min-norm(c) that \( (\I-\bP)\bb_t=(\I-\bP)\bb_0 \) for all \( t \).
Deduce that gradient descent started at \( \bb_0 \) returns the least squares estimate closest to \( \bb_0 \).
:::

::: {#exr-hd-ridge-beats-ridgeless}
[B2]

Fix \( \X \) of rank \( r \) and let \( m(\lambda)=\E\norm{\hbeta_\lambda-\bbeta}^2 \) for the ridge estimate, with
\( m(0) \) the risk of \( \hbeta^+ \) in @thm-hd-min-norm(d). Show that
\[
m(\lambda)=\norm{(\I-\bP)\bbeta}^2+\sum_{i=1}^r\Bigl[\Bigl(\frac{\lambda}{d_i^2+\lambda}\Bigr)^2(\bv_i\T\bbeta)^2
+\sigma^2\frac{d_i^2}{(d_i^2+\lambda)^2}\Bigr],
\]
and that \( m'(0)<0 \) whenever \( \sigma>0 \). So some ridge estimate beats the minimum-norm
interpolator for every \( \bbeta \). Compare @thm-shr-ridge-dominates.
:::

::: {.solution}
The linear representation in @thm-shr-ridge(a) holds for any rank \( r \). Applied to
\( \E(\Y)=\X\bbeta \), whose coordinates are \( \bu_i\T\X\bbeta=d_i\bv_i\T\bbeta \), it gives
\( \E\hbeta_\lambda=\sum_{i\le r}\{d_i^2/(d_i^2+\lambda)\}(\bv_i\T\bbeta)\bv_i \), and the covariance is
\( \sigma^2\sum_{i\le r}d_i^2(d_i^2+\lambda)^{-2}\bv_i\bv_i\T \), as in @thm-shr-ridge(c). What \( r<p \) adds is the
bias component \( -(\I-\bP)\bbeta \) in the null space, orthogonal to the rest of the bias
\( -\sum_i\{\lambda/(d_i^2+\lambda)\}(\bv_i\T\bbeta)\bv_i \); it contributes \( \norm{(\I-\bP)\bbeta}^2 \). Adding the
squared bias and the trace of the covariance gives \( m(\lambda) \).
Differentiating term by term, the first bracket has derivative \( 2\lambda d_i^2(\bv_i\T\bbeta)^2/(d_i^2+\lambda)^3 \), which is
zero at \( \lambda=0 \). The second has derivative \( -2\sigma^2d_i^2/(d_i^2+\lambda)^3 \), which equals \( -2\sigma^2/d_i^4<0 \) at
\( \lambda=0 \). Hence \( m'(0)=-2\sigma^2\sum_id_i^{-4}<0 \), and \( m(\lambda)<m(0) \) for small \( \lambda>0 \).
:::

### C. Going deeper

::: {#exr-hd-one-observation}
[C1]

Check @eq-hd-isotropic-risk for \( n=1 \) directly. With one observation \( (\x,y) \), show that
\( \hbeta^+=\x y/\norm{\x}^2 \), and use \( \E(1/\chi^2(p))=1/(p-2) \) for \( p\ge3 \) to show that
\( \mathcal{R}=\norm{\bbeta}^2(1-1/p)+\sigma^2/(p-2) \).
:::

::: {.solution}
Here \( \X=\x\T \) and \( \X\X\T=\norm{\x}^2 \), so \( \hbeta^+=\x y/\norm{\x}^2 \) and
\( \hbeta^+-\bbeta=-(\I-\x\x\T/\norm{\x}^2)\bbeta+\x\varepsilon/\norm{\x}^2 \), a sum of orthogonal vectors. So
\[
\norm{\hbeta^+-\bbeta}^2=\norm{\bbeta}^2-\frac{(\x\T\bbeta)^2}{\norm{\x}^2}+\frac{\varepsilon^2}{\norm{\x}^2}.
\]
The matrix
\( \bP=\x\x\T/\norm{\x}^2 \) is the row-space projection, and the rotation argument in the proof of
@prp-hd-isotropic-risk gives \( \E\bP=\I/p \), so \( \E(\x\T\bbeta)^2/\norm{\x}^2=\bbeta\T(\E\bP)\bbeta=\norm{\bbeta}^2/p \).
And \( \norm{\x}^2\sim\chi^2(p) \) is independent of \( \varepsilon \), so
\( \E\varepsilon^2/\norm{\x}^2=\sigma^2/(p-2) \). This agrees with @eq-hd-isotropic-risk at \( n=1 \).
:::

::: {#exr-hd-proportional}
[C2]

In @prp-hd-isotropic-risk let \( n,p\to\infty \) with \( n/p\to\gamma\in(0,1) \) and \( \norm{\bbeta}^2 \) fixed.
Show that \( \mathcal{R}\to\norm{\bbeta}^2(1-\gamma)+\sigma^2\gamma/(1-\gamma) \), and that this is below the
risk \( \norm{\bbeta}^2 \) of predicting zero iff \( \norm{\bbeta}^2/\sigma^2>1/(1-\gamma) \). What happens to the
limit as \( \gamma\to0 \)?
:::
