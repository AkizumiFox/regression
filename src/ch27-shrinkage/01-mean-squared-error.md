# The bias–variance trade-off

[Chapter 26](../ch26-collinearity/index.html) showed where collinearity does its damage (@thm-col-variance). The least squares estimate remains unbiased, but along a direction
\( \bv \) in which the regressors barely vary, the variance of \( \bv\T\hbeta \) is
\( \sigma^2 \) divided by a tiny eigenvalue; for a single coefficient the same damage appears as the factor
\( 1/(1-R_j^2) \) of @eq-proj-vif-preview, the variance inflation factor of @def-col-vif. No assumption of Parts II and III is violated. The Gauss–Markov theorem (@thm-opt-gauss-markov)
concerns the class of *linear unbiased* estimators, within which least
squares cannot be improved, so improving on it means leaving the class. This section sets up the
criterion that makes biased and unbiased estimators comparable.

## The mean squared error matrix

::: {#def-shr-mse}
[Mean squared error]

Let \( \tilde{\boldsymbol{\uptheta}} \) be an estimator of \( \boldsymbol{\uptheta}\in\Real^q \)
with finite second moments. Its **mean squared error matrix** is
\[
\text{MSE}(\tilde{\boldsymbol{\uptheta}})
=\E\bigl[(\tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta})(\tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta})\T\bigr],
\]
its **bias** is \( \boldsymbol{\updelta}=\E(\tilde{\boldsymbol{\uptheta}})-\boldsymbol{\uptheta} \), and its
**total mean squared error** is
\( \tr\text{MSE}(\tilde{\boldsymbol{\uptheta}})=\E\norm{\tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}}^2 \).
For estimators \( \tilde{\boldsymbol{\uptheta}}_1,\tilde{\boldsymbol{\uptheta}}_2 \) we say that
\( \tilde{\boldsymbol{\uptheta}}_1 \) is **at least as good in the matrix sense** if
\( \text{MSE}(\tilde{\boldsymbol{\uptheta}}_2)-\text{MSE}(\tilde{\boldsymbol{\uptheta}}_1) \) is
nonnegative definite.
:::

All three quantities depend on the unknown parameter, so one estimator may be better for some
parameter values and worse for others. For unbiased estimators the mean squared error matrix is
the covariance matrix, and the matrix sense is the ordering in which Gauss–Markov declares least
squares the winner (@thm-opt-gauss-markov(d)).

::: {#thm-shr-bias-variance}
[Bias–variance decomposition for linear estimators]

::: {.enumerate options="label=(\alph*)"}
1. For any estimator with finite second moments,
   \( \text{MSE}(\tilde{\boldsymbol{\uptheta}})=\Cov(\tilde{\boldsymbol{\uptheta}})+\boldsymbol{\updelta}\boldsymbol{\updelta}\T \). For
   every \( \mathbf{a}\in\Real^q \),
   \( \E(\mathbf{a}\T\tilde{\boldsymbol{\uptheta}}-\mathbf{a}\T\boldsymbol{\uptheta})^2=\mathbf{a}\T\text{MSE}(\tilde{\boldsymbol{\uptheta}})\mathbf{a} \),
   and for every nonnegative definite \( \W \),
   \( \E\bigl[(\tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta})\T\W(\tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta})\bigr]=\tr\bigl(\W\,\text{MSE}(\tilde{\boldsymbol{\uptheta}})\bigr) \).

2. Suppose \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\I \), and let \( \tilde{\bbeta}=\A\Y \) be a
   linear estimator of \( \bbeta \). Then \( \boldsymbol{\updelta}=(\A\X-\I)\bbeta \) and
   \( \text{MSE}(\tilde{\bbeta})=\sigma^2\A\A\T+(\A\X-\I)\bbeta\bbeta\T(\A\X-\I)\T \).

3. In the setting of (b), let \( \X \) have full column rank, let \( \hbeta \) be least squares, and
   put \( \boldsymbol{\Delta}_{\A}=\sigma^2(\X\T\X)^{-1}-\sigma^2\A\A\T \), the covariance that \( \tilde{\bbeta} \)
   saves. The following are equivalent:

   (i) \( \tilde{\bbeta} \) is at least as good as \( \hbeta \) in the matrix sense;

   (ii) \( \E(\mathbf{a}\T\tilde{\bbeta}-\mathbf{a}\T\bbeta)^2\le\Var(\mathbf{a}\T\hbeta) \) for every \( \mathbf{a} \);

   (iii) \( \E\bigl[(\tilde{\bbeta}-\bbeta)\T\W(\tilde{\bbeta}-\bbeta)\bigr]\le\E\bigl[(\hbeta-\bbeta)\T\W(\hbeta-\bbeta)\bigr] \) for every nonnegative definite \( \W \);

   (iv) \( \boldsymbol{\Delta}_{\A} \) is nonnegative definite, \( \boldsymbol{\updelta}\in\C(\boldsymbol{\Delta}_{\A}) \) and \( \boldsymbol{\updelta}\T\boldsymbol{\Delta}_{\A}^+\boldsymbol{\updelta}\le1 \).

   When \( \boldsymbol{\Delta}_{\A} \) is positive definite, (iv) reads \( \boldsymbol{\updelta}\T\boldsymbol{\Delta}_{\A}^{-1}\boldsymbol{\updelta}\le1 \).
:::
:::

::: {.proof}
(a) Write \( \tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}=(\tilde{\boldsymbol{\uptheta}}-\E\tilde{\boldsymbol{\uptheta}})+\boldsymbol{\updelta} \).
The cross terms have expectation zero, which gives the decomposition. The scalar identities
follow from \( \E(\mathbf{a}\T\mathbf{w})^2=\mathbf{a}\T\E(\mathbf{w}\mathbf{w}\T)\mathbf{a} \) and
\( \mathbf{w}\T\W\mathbf{w}=\tr(\W\mathbf{w}\mathbf{w}\T) \) with
\( \mathbf{w}=\tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta} \).
(b) \( \E(\A\Y)=\A\X\bbeta \) and \( \Cov(\A\Y)=\sigma^2\A\A\T \) (@thm-rv-linear); apply (a).
(c) By (a), \( \text{MSE}(\hbeta)-\text{MSE}(\tilde{\bbeta})=\boldsymbol{\Delta}_{\A}-\boldsymbol{\updelta}\boldsymbol{\updelta}\T \), and (i) says this
matrix is nonnegative definite. By (a), (ii) says
\( \mathbf{a}\T(\boldsymbol{\Delta}_{\A}-\boldsymbol{\updelta}\boldsymbol{\updelta}\T)\mathbf{a}\ge0 \) for all \( \mathbf{a} \), which is (i).
Next, (i) implies (iii) because \( \tr(\W\mathbf{N})\ge0 \) when \( \W \) and \( \mathbf{N} \) are both
nonnegative definite: with \( \W=\sum_k\omega_k\mathbf{q}_k\mathbf{q}_k\T \) (@thm-mat-spectral),
\( \tr(\W\mathbf{N})=\sum_k\omega_k\mathbf{q}_k\T\mathbf{N}\mathbf{q}_k\ge0 \). And (iii) with
\( \W=\mathbf{a}\mathbf{a}\T \) is (ii).
For (i)\( \Leftrightarrow \)(iv), put \( \mathbf{N}=\boldsymbol{\Delta}_{\A}-\boldsymbol{\updelta}\boldsymbol{\updelta}\T \). If \( \mathbf{N} \) is
nonnegative definite, so is \( \boldsymbol{\Delta}_{\A}=\mathbf{N}+\boldsymbol{\updelta}\boldsymbol{\updelta}\T \). If \( \boldsymbol{\updelta}\notin\C(\boldsymbol{\Delta}_{\A}) \), there is
\( \mathbf{w}\in\Null(\boldsymbol{\Delta}_{\A}) \) with \( \mathbf{w}\T\boldsymbol{\updelta}\ne0 \) (@prp-mat-complement), and then
\( \mathbf{w}\T\mathbf{N}\mathbf{w}=-(\mathbf{w}\T\boldsymbol{\updelta})^2<0 \). So (i) forces \( \boldsymbol{\updelta}=\boldsymbol{\Delta}_{\A}\mathbf{c} \) for
some \( \mathbf{c} \), and \( \mathbf{c}\T\mathbf{N}\mathbf{c}=\mathbf{c}\T\boldsymbol{\Delta}_{\A}\mathbf{c}\,(1-\mathbf{c}\T\boldsymbol{\Delta}_{\A}\mathbf{c})\ge0 \)
gives \( \boldsymbol{\updelta}\T\boldsymbol{\Delta}_{\A}^+\boldsymbol{\updelta}=\mathbf{c}\T\boldsymbol{\Delta}_{\A}\mathbf{c}\le1 \) (when \( \mathbf{c}\T\boldsymbol{\Delta}_{\A}\mathbf{c}=0 \) there
is nothing to prove). Here \( \boldsymbol{\updelta}\T\boldsymbol{\Delta}_{\A}^+\boldsymbol{\updelta}=\mathbf{c}\T\boldsymbol{\Delta}_{\A}\boldsymbol{\Delta}_{\A}^+\boldsymbol{\Delta}_{\A}\mathbf{c}=\mathbf{c}\T\boldsymbol{\Delta}_{\A}\mathbf{c} \)
by the first Penrose condition. Conversely, under (iv) the Cauchy–Schwarz inequality for the
semi-inner product \( (\mathbf{w},\mathbf{z})\mapsto\mathbf{w}\T\boldsymbol{\Delta}_{\A}\mathbf{z} \) gives, for every \( \mathbf{w} \),
\[
(\mathbf{w}\T\boldsymbol{\updelta})^2=(\mathbf{w}\T\boldsymbol{\Delta}_{\A}\mathbf{c})^2\le(\mathbf{w}\T\boldsymbol{\Delta}_{\A}\mathbf{w})(\mathbf{c}\T\boldsymbol{\Delta}_{\A}\mathbf{c})\le\mathbf{w}\T\boldsymbol{\Delta}_{\A}\mathbf{w},
\]
so \( \mathbf{w}\T\mathbf{N}\mathbf{w}\ge0 \). When \( \boldsymbol{\Delta}_{\A} \) is invertible, \( \boldsymbol{\Delta}_{\A}^+=\boldsymbol{\Delta}_{\A}^{-1} \);
this case is also @lem-dep-rank-one.
:::

Part (c) makes the trade-off precise: a biased linear estimator is better in every direction and
under every quadratic loss exactly when its squared bias, measured in the metric of the covariance
it saves, is at most one. Since \( \boldsymbol{\updelta}=(\A\X-\I)\bbeta \), condition (iv)
holds when \( \bbeta \) is close to what the estimator shrinks towards, and fails far from it. No linear estimator with \( \A\X\ne\I \) wins
everywhere (@exr-shr-no-linear-dominance). @thm-dep-mse was the first instance: the short
regression is a biased linear estimator, and its winning region \( \gamma\le1 \) is condition (iv)
for that choice of \( \A \).

## What unbiasedness buys, and what it costs

Unbiasedness concerns an estimator's centre, not its accuracy, and it does not survive
transformation: \( \norm{\hbeta}^2 \) overestimates \( \norm{\bbeta}^2 \) by \( \sigma^2\tr(\X\T\X)^{-1} \) on
average (@exr-shr-bias-length), which under collinearity is enormous. Least squares estimates are
systematically *too long*, and shrinking them is a correction, not a distortion. When many
parameters are estimated at once, small biases spread over many coordinates can buy a large cut in
total variance; [Section 27.5](05-james-stein.html) turns this into a theorem. The cost is in
inference: the bias depends on the unknown \( \bbeta \), so an interval "estimate \( \pm \) multiple of
standard error" centred on a shrinkage estimate has no guaranteed coverage. [Chapter 28](../ch28-high-dimensional/index.html) and
[Chapter 29](../ch29-model-selection/index.html) return to inference after shrinkage and selection.

## Shrinkage in canonical coordinates

The collinearity problem has a natural coordinate system. Let \( \X \) have full column rank and
let \( \X=\bU\bD\V\T \) be its thin singular value decomposition (@thm-mat-svd), with
\( \bU=[\bu_1,\dots,\bu_p] \) and \( \V=[\bv_1,\dots,\bv_p] \) having orthonormal columns and
\( \bD=\diag(d_1,\dots,d_p) \), \( d_1\ge\dots\ge d_p>0 \). We write \( d_j \) instead of the
\( \sigma_j \) of [Chapter 1](../ch01-matrix-algebra/index.html) because \( \sigma \) is the error standard deviation. Put
\[
\boldsymbol{\upalpha}=\V\T\bbeta,\qquad
\hat{\boldsymbol{\upalpha}}=\V\T\hbeta=\bD^{-1}\bU\T\Y,\qquad
c_j=\bu_j\T\Y=d_j\hat{\alpha}_j .
\]{#eq-shr-canonical}

The \( \alpha_j \) are the **canonical coefficients**: the coordinates of \( \bbeta \) along the
eigenvectors \( \bv_j \) of \( \X\T\X \), whose eigenvalues are \( d_j^2 \). By @thm-cmp-svd-ls(d),
the \( \hat{\alpha}_j \) are uncorrelated with \( \E(\hat{\alpha}_j)=\alpha_j \) and
\( \Var(\hat{\alpha}_j)=\sigma^2/d_j^2 \). A small \( d_j \) is a direction \( \bv_j \) along which the
regressors hardly move, and along which the data therefore say little (@prp-dep-collinear-directions).

Many shrinkage estimators act on the canonical coordinates one at a time. For numbers
\( f_1,\dots,f_p \), called **shrinkage factors** or filter factors, define
\[
\tilde{\bbeta}_{\mathbf{f}}=\sum_{j=1}^pf_j\hat{\alpha}_j\bv_j=\V\mathbf{F}\hat{\boldsymbol{\upalpha}},
\qquad \mathbf{F}=\diag(f_1,\dots,f_p).
\]
Least squares has every \( f_j=1 \). Principal component regression
([Section 27.3](03-principal-components.html)) has \( f_j\in\{0,1\} \), and ridge regression
([Section 27.2](02-ridge.html)) has \( f_j=d_j^2/(d_j^2+\lambda) \). The signal-to-noise ratio of
direction \( j \) is
\[
\tau_j^2=\frac{\alpha_j^2}{\Var(\hat{\alpha}_j)}=\frac{d_j^2\alpha_j^2}{\sigma^2},
\]
the squared mean of the standardized coordinate \( \hat{\alpha}_j/\sqrt{\Var(\hat{\alpha}_j)} \). Under
normal errors it is the noncentrality of \( c_j^2/\sigma^2\sim\chi^2(1,\tau_j^2) \).

::: {#prp-shr-canonical}
[Shrinking along the principal directions]

Let \( \E(\Y)=\X\bbeta \), \( \Cov(\Y)=\sigma^2\I \) and \( \rank(\X)=p \). For fixed factors \( f_j \):

::: {.enumerate options="label=(\alph*)"}
1. \( \E\norm{\tilde{\bbeta}_{\mathbf{f}}-\bbeta}^2=\sum_j\bigl\{f_j^2\sigma^2/d_j^2+(1-f_j)^2\alpha_j^2\bigr\} \);

2. \( \E\norm{\X\tilde{\bbeta}_{\mathbf{f}}-\X\bbeta}^2=\sum_j\bigl\{f_j^2\sigma^2+(1-f_j)^2d_j^2\alpha_j^2\bigr\} \);

3. both are minimized by the same factors \( f_j^*=\tau_j^2/(1+\tau_j^2) \), and the minima are
   \( \sum_j(\sigma^2/d_j^2)\,\tau_j^2/(1+\tau_j^2) \) and \( \sigma^2\sum_j\tau_j^2/(1+\tau_j^2) \);

4. replacing \( f_j=1 \) by some \( f_j\in[0,1) \), with the other factors fixed, lowers both
   quantities iff \( \tau_j^2<(1+f_j)/(1-f_j) \). In particular dropping direction \( j \)
   (\( f_j=0 \)) helps iff \( \tau_j^2<1 \).
:::
:::

::: {.proof}
\( \tilde{\bbeta}_{\mathbf{f}}-\bbeta=\V(\mathbf{F}\hat{\boldsymbol{\upalpha}}-\boldsymbol{\upalpha}) \) and
\( \V \) is orthogonal, so \( \norm{\tilde{\bbeta}_{\mathbf{f}}-\bbeta}^2=\sum_j(f_j\hat{\alpha}_j-\alpha_j)^2 \).
Each term has expectation \( \Var(f_j\hat{\alpha}_j)+(f_j\alpha_j-\alpha_j)^2=f_j^2\sigma^2/d_j^2+(1-f_j)^2\alpha_j^2 \),
which is (a). For (b), \( \X(\tilde{\bbeta}_{\mathbf{f}}-\bbeta)=\bU\bD(\mathbf{F}\hat{\boldsymbol{\upalpha}}-\boldsymbol{\upalpha}) \)
and \( \bU \) has orthonormal columns, so each term of (a) is multiplied by \( d_j^2 \). (c) Each
term of (a) is a convex quadratic in \( f_j \) with derivative
\( 2f_j\sigma^2/d_j^2-2(1-f_j)\alpha_j^2 \), which vanishes at
\( f_j=\alpha_j^2/(\alpha_j^2+\sigma^2/d_j^2)=\tau_j^2/(1+\tau_j^2) \). Terms of (b) are those of (a)
times \( d_j^2 \), so they have the same minimizer; substituting gives the minima. (d) The
\( j \)th term of (a) falls below its value \( \sigma^2/d_j^2 \) at \( f_j=1 \) iff
\( (1-f_j)^2\alpha_j^2<(1-f_j^2)\sigma^2/d_j^2 \). Dividing by \( (1-f_j)\sigma^2/d_j^2>0 \) gives
\( (1-f_j)\tau_j^2<1+f_j \). The same holds for (b).
:::

The ideal shrinkage depends on \( \tau_j^2 \), which combines the spread of the regressors in direction
\( j \) *and* the size of the coefficient there, so a direction with tiny \( d_j \) but large \( \alpha_j \)
should be kept. The ideal factors are the same for estimating \( \bbeta \) and \( \X\bbeta \). The oracle factors depend
on \( \bbeta \) and \( \sigma^2 \); each method of this chapter chooses factors from the data instead, ridge by
one tuning constant, principal components by a cut-off, and James–Stein by an unbiased estimate
of the risk. Part (d) with \( f_j=0 \) is the criterion of @prp-cmp-tsvd.

::: {#exm-shr-two-regressors}
[Two nearly collinear regressors]

Take \( n=50 \) observations of two centred regressors, scaled so that \( \X\T\X=n\mathbf{R} \) where
\( \mathbf{R} \) is a correlation matrix with off-diagonal entry \( 0.95 \). Let \( \sigma=1 \) and
\( \bbeta=(1,0.5)\T \). The eigenvalues of \( \X\T\X \) are \( d_1^2=n(1+0.95)=97.5 \) and
\( d_2^2=n(1-0.95)=2.5 \), with eigenvectors \( (1,1)\T/\sqrt2 \) and \( (1,-1)\T/\sqrt2 \). The
canonical coefficients have sizes \( \lvert\alpha_1\rvert=1.0607 \) and \( \lvert\alpha_2\rvert=0.3536 \), and the signal-to-noise
ratios are \( \tau_1^2=109.7 \) and \( \tau_2^2=0.3125 \). The first direction, the sum of the
two regressors, is estimated precisely. The second, their difference, is not.

| Estimator | Factors \( (f_1,f_2) \) | Total MSE |
|---|---|---|
| least squares | \( (1,1) \) | 0.4103 |
| drop the second direction | \( (1,0) \) | 0.1353 |
| halve the second direction | \( (1,0.5) \) | 0.1415 |
| oracle \( f_j^* \) | \( (0.991,0.238) \) | 0.1054 |

Dropping the difference direction cuts the total mean squared error by two thirds, because
\( \tau_2^2<1 \), and wins in the matrix sense too. Here \( \boldsymbol{\Delta}_{\A}=\sigma^2\bv_2\bv_2\T/d_2^2 \)
is singular, the bias is \( \boldsymbol{\updelta}=-\alpha_2\bv_2\in\C(\boldsymbol{\Delta}_{\A}) \), and
\( \boldsymbol{\updelta}\T\boldsymbol{\Delta}_{\A}^+\boldsymbol{\updelta}=\alpha_2^2d_2^2/\sigma^2=\tau_2^2=0.3125\le1 \). So by
@thm-shr-bias-variance(c), *every* linear combination \( \mathbf{a}\T\bbeta \) is estimated at least
as accurately. Had the truth been
\( \bbeta=(1,-0.5)\T \), with the same length but pointing along the difference, then
\( \tau_2^2=2.5\times1.125=2.81 \) and dropping would have increased the total error.
:::

```{.python .run #cell-mse-example-example}
import numpy as np
n, rho, sigma = 50, 0.95, 1.0
XtX = n * np.array([[1.0, rho], [rho, 1.0]])
beta = np.array([1.0, 0.5])

evals, V = np.linalg.eigh(XtX)
order = np.argsort(evals)[::-1]
d2, V = evals[order], V[:, order]          # d_1^2 >= d_2^2, eigenvectors v_1, v_2
alpha = V.T @ beta                         # canonical coefficients
tau2 = d2 * alpha ** 2 / sigma ** 2        # squared signal-to-noise ratio of each direction

def mse(f):
    """Total MSE of the estimator that multiplies alpha_hat_j by f_j."""
    f = np.asarray(f, dtype=float)
    return np.sum(f ** 2 * sigma ** 2 / d2 + (1 - f) ** 2 * alpha ** 2)

f_oracle = tau2 / (1 + tau2)
for label, f in [("least squares", [1, 1]), ("drop v_2", [1, 0]),
                 ("halve v_2", [1, 0.5]), ("oracle", f_oracle)]:
    print(f"{label:14s} f = {np.round(f, 3)}   MSE = {mse(f):.4f}")
print("tau^2 =", np.round(tau2, 3))
```

## Exercises

### A. Check your understanding

::: {#exr-shr-scalar}
[A1]

Let \( \bar{Y} \) be the mean of \( n \) independent observations with mean \( \mu \) and variance
\( \sigma^2 \). For a constant \( c \), find the mean squared error of \( c\bar{Y} \) as an estimator of
\( \mu \). Show that it is minimized at \( c^*=\mu^2/(\mu^2+\sigma^2/n) \) and that
\( c\bar{Y} \) beats \( \bar{Y} \) iff \( (1-c)\mu^2<(1+c)\sigma^2/n \). Why is \( c^* \) not an estimator?
:::

::: {#exr-shr-two-check}
[A2]

In @exm-shr-two-regressors, verify the least squares total mean squared error
\( \sigma^2(1/d_1^2+1/d_2^2)=0.4103 \) by hand, and check the entry for "drop the second
direction" using @prp-shr-canonical(a).
:::

### B. Practice

::: {#exr-shr-no-linear-dominance}
[B1]

Let \( \X \) have full column rank and let \( \tilde{\bbeta}=\A\Y \) with \( \A\X\ne\I \). Show that
there is a \( \bbeta \) for which \( \E\norm{\tilde{\bbeta}-\bbeta}^2>\E\norm{\hbeta-\bbeta}^2 \). So no
biased linear estimator dominates least squares.
:::

::: {.solution}
By @thm-shr-bias-variance(b), \( \E\norm{\tilde{\bbeta}-\bbeta}^2=\sigma^2\tr(\A\A\T)+\norm{(\A\X-\I)\bbeta}^2 \).
The first term does not depend on \( \bbeta \). Since \( \A\X-\I\ne\mathbf{0} \) there is a unit vector
\( \bv \) with \( (\A\X-\I)\bv\ne\bzero \). Along \( \bbeta=t\bv \) the second term is
\( t^2\norm{(\A\X-\I)\bv}^2 \), which exceeds the constant \( \sigma^2\tr(\X\T\X)^{-1} \) for large \( t \).
:::

::: {#exr-shr-smoother-mse}
[B2]

A **linear smoother** estimates the mean vector \( \bmu=\E(\Y) \) by \( \bS\Y \) for a fixed
\( n\times n \) matrix \( \bS \). Show that
\( \E\norm{\bS\Y-\bmu}^2=\sigma^2\tr(\bS\bS\T)+\norm{(\I-\bS)\bmu}^2 \) when \( \Cov(\Y)=\sigma^2\I \). For
the projection \( \bS=\M \) and \( \bmu\in\C(\X) \), recover \( \sigma^2r \) with \( r=\rank(\X) \).
:::

::: {.solution}
\( \bS\Y-\bmu=\bS(\Y-\bmu)-(\I-\bS)\bmu \). The first term has mean zero and covariance
\( \sigma^2\bS\bS\T \), so its expected squared length is \( \sigma^2\tr(\bS\bS\T) \) (@thm-rv-quadform-mean). The second term is constant, and the cross term has mean zero. For
\( \bS=\M \), \( \tr(\M\M\T)=\tr\M=r \) and \( (\I-\M)\bmu=\bzero \).
:::

::: {#exr-shr-bias-length}
[B3]

Show that \( \E\norm{\hbeta}^2=\norm{\bbeta}^2+\sigma^2\sum_jd_j^{-2} \). For which \( c \) does
\( c\hbeta \) have smaller total mean squared error than \( \hbeta \)? Express the answer through
\( \norm{\bbeta}^2 \) and \( \sigma^2\tr(\X\T\X)^{-1} \).
:::

### C. Going deeper

::: {#exr-shr-scaling}
[C1]

Least squares is equivariant under rescaling: replacing \( \X \) by \( \X\mathbf{G}^{-1} \), with
\( \mathbf{G} \) diagonal and positive, turns the estimate into \( \mathbf{G}\hbeta \). By rescaling one column
in @exm-shr-two-regressors by \( 10 \), show that \( \tilde{\bbeta}_{\mathbf{f}} \) with \( f_1=1 \), \( f_2=0 \) is not
equivariant. Why are shrinkage estimators usually applied to standardized regressors?
:::
