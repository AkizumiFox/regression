# Ridge regression

The oldest remedy for collinearity is also the simplest: when \( \X\T\X \) has eigenvalues close to
zero, add a positive constant to its diagonal before inverting. Hoerl and Kennard (1970) proposed
this in statistics; numerical analysts knew it as Tikhonov regularization. The constant bounds the
variance at the price of bias towards the origin. @exr-opt-ridge and @exr-cmp-ridge-svd met the estimator
briefly.

## Centring, scaling and the intercept

A penalty should not touch the intercept, or shifting the response by a constant would change the
fitted slopes. Leaving the intercept unpenalized is equivalent to centring, for every penalty.

::: {#prp-shr-centring}
[An unpenalized intercept]

Let \( \mathbf{X}_0 \) be an \( n\times p \) matrix of regressors with column means \( \bar{\mathbf{x}}_0 \), let
\( \X=\mathbf{X}_0-\bone\bar{\mathbf{x}}_0\T \) and \( \tilde{\y}=\y-\bar y\bone \) be the centred regressors and response,
and let \( P:\Real^p\to[0,\infty] \) be any function. Then \( (\beta_0,\bb) \) minimizes
\[
\norm{\y-\beta_0\bone-\mathbf{X}_0\bb}^2+P(\bb)
\]
iff \( \bb \) minimizes \( \norm{\tilde{\y}-\X\bb}^2+P(\bb) \) and \( \beta_0=\bar y-\bar{\mathbf{x}}_0\T\bb \).
:::

::: {.proof}
For fixed \( \bb \), the first term is the residual sum of squares of regressing
\( \y-\mathbf{X}_0\bb \) on \( \bone \), minimized uniquely by the mean \( \beta_0=\bar y-\bar{\mathbf{x}}_0\T\bb \), with
minimum \( \norm{(\I-n^{-1}\bone\bone\T)(\y-\mathbf{X}_0\bb)}^2=\norm{\tilde{\y}-\X\bb}^2 \). So the joint
minimum is the minimum over \( \bb \) of \( \norm{\tilde{\y}-\X\bb}^2+P(\bb) \), attained exactly at the
pairs described.
:::

From now on the regressors are centred and the response is centred, so \( \X \) has no intercept
column and \( \bone\T\X=\bzero\T \). A penalty also treats all coefficients as commensurable, so
the regressors should be on a common scale (@exr-shr-scaling). In this chapter's examples every
column has unit standard deviation, so \( \X\T\X=(n-1)\mathbf{R} \) with \( \mathbf{R} \) the correlation
matrix, and a coefficient is the change in mean response per standard deviation of its regressor.

## The estimator

::: {#def-shr-ridge}
[Ridge regression]

For \( \lambda\ge0 \), a **ridge estimate** \( \hbeta^{\mathrm{R}}_\lambda \) minimizes
\[
\norm{\y-\X\bb}^2+\lambda\norm{\bb}^2 .
\]{#eq-shr-ridge-criterion}

The constant \( \lambda \) is the **ridge parameter** (or tuning parameter). The **ridge trace** is
the curve \( \lambda\mapsto\hbeta^{\mathrm{R}}_\lambda \).
:::

For \( \lambda=0 \) this is least squares. For \( \lambda>0 \) it is uniquely defined whatever the rank
of \( \X \), and its behaviour is transparent in the coordinates of @eq-shr-canonical. Let
\( \X=\sum_{j\le r}d_j\bu_j\bv_j\T \) be the SVD of \( \X \), of rank \( r \), and \( c_j=\bu_j\T\y \).

::: {#thm-shr-ridge}
[Ridge regression]

Let \( \lambda>0 \).

::: {.enumerate options="label=(\alph*)"}
1. *(Closed form.)* The criterion @eq-shr-ridge-criterion has the unique minimizer
   \[
   \hbeta^{\mathrm{R}}_\lambda=(\X\T\X+\lambda\I)^{-1}\X\T\y
   =\sum_{j\le r}\frac{d_j}{d_j^2+\lambda}\,c_j\,\bv_j .
   \]
   If \( r=p \), then \( \hbeta^{\mathrm{R}}_\lambda=\sum_jf_j(\lambda)\hat{\alpha}_j\bv_j \) with **shrinkage factors**
   \( f_j(\lambda)=d_j^2/(d_j^2+\lambda)\in(0,1) \), decreasing in \( \lambda \) and smallest for the
   smallest \( d_j \).

2. *(Fitted values.)* \( \X\hbeta^{\mathrm{R}}_\lambda=\bH_\lambda\y \) with
   \( \bH_\lambda=\X(\X\T\X+\lambda\I)^{-1}\X\T=\sum_{j\le r}f_j(\lambda)\bu_j\bu_j\T \). The **effective
   degrees of freedom** \( \mathrm{df}(\lambda)=\tr\bH_\lambda=\sum_{j\le r}f_j(\lambda) \) decrease strictly
   from \( r \) (as \( \lambda\to0 \)) to \( 0 \) (as \( \lambda\to\infty \)).

3. *(Moments.)* If \( \E(\Y)=\X\bbeta \), \( \Cov(\Y)=\sigma^2\I \) and \( r=p \), then
   \[
   \E\hbeta^{\mathrm{R}}_\lambda-\bbeta=-\lambda(\X\T\X+\lambda\I)^{-1}\bbeta,
   \]
   \[
   \Cov\bigl(\hbeta^{\mathrm{R}}_\lambda\bigr)=\sigma^2\sum_j\frac{d_j^2}{(d_j^2+\lambda)^2}\bv_j\bv_j\T ,
   \]
   and, with \( \boldsymbol{\upalpha}=\V\T\bbeta \), the total mean squared error is
   \[
   \E\norm{\hbeta^{\mathrm{R}}_\lambda-\bbeta}^2=\sigma^2\sum_j\frac{d_j^2}{(d_j^2+\lambda)^2}
   +\lambda^2\sum_j\frac{\alpha_j^2}{(d_j^2+\lambda)^2}.
   \]{#eq-shr-ridge-mse}

4. *(The path.)* Unless \( \X\T\y=\bzero \), \( \norm{\hbeta^{\mathrm{R}}_\lambda} \) decreases strictly in
   \( \lambda \), with limit \( 0 \) as \( \lambda\to\infty \); and \( \hbeta^{\mathrm{R}}_\lambda\to\X^+\y \), the
   least squares estimate of smallest norm, as \( \lambda\to0 \).

5. *(Bayesian reading.)* If \( \Y\mid\bbeta\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \sigma^2 \) known and
   \( \bbeta\sim\Normal_p(\bzero,(\sigma^2/\lambda)\I) \), then
   \( \bbeta\mid\Y\sim\Normal_p\bigl(\hbeta^{\mathrm{R}}_\lambda,\sigma^2(\X\T\X+\lambda\I)^{-1}\bigr) \). With
   \( \sigma^2 \) unknown and the normal–inverse-gamma prior of @thm-opt-bayes-conjugate with
   \( \mathbf{m}_0=\bzero \) and \( \V_0=\lambda^{-1}\I \), the posterior mean of \( \bbeta \) is again
   \( \hbeta^{\mathrm{R}}_\lambda \).
:::
:::

::: {.proof}
(a) Uniqueness and the first formula are @exr-mat-ridge: \( \X\T\X+\lambda\I \) is positive definite, and
@prp-mat-quadratic-min applies. Complete \( \bv_1,\dots,\bv_r \) to an orthonormal basis
\( \bv_1,\dots,\bv_p \) of \( \Real^p \). Then \( \X\T\X+\lambda\I=\sum_{j\le r}(d_j^2+\lambda)\bv_j\bv_j\T+\sum_{j>r}\lambda\bv_j\bv_j\T \),
whose inverse has the same eigenvectors and reciprocal eigenvalues (@thm-mat-spectral). Since
\( \X\T\y=\sum_{j\le r}d_jc_j\bv_j \), the second formula follows. When \( r=p \),
\( d_jc_j/(d_j^2+\lambda)=f_j(\lambda)\,c_j/d_j=f_j(\lambda)\hat{\alpha}_j \).
(b) \( \X\bv_j=d_j\bu_j \), so \( \X\hbeta^{\mathrm{R}}_\lambda=\sum_{j\le r}f_j(\lambda)c_j\bu_j=\sum_jf_j(\lambda)\bu_j\bu_j\T\y \).
The trace of \( \bu_j\bu_j\T \) is one. Each \( f_j \) decreases strictly from \( 1 \) to \( 0 \).
(c) These are the formulas of @exr-opt-ridge in the present notation (its \( d_j \) and \( \bgamma \)
are our \( d_j^2 \) and \( \boldsymbol{\upalpha} \)); they are also @prp-shr-canonical(a) with the factors \( f_j(\lambda) \).
(d) \( \norm{\hbeta^{\mathrm{R}}_\lambda}^2=\sum_{j\le r}d_j^2c_j^2/(d_j^2+\lambda)^2 \). Every term is
nonincreasing in \( \lambda \) and tends to zero, and a term with \( c_j\ne0 \) decreases strictly. Some
\( c_j\ne0 \) iff \( \X\T\y\ne\bzero \). As \( \lambda\to0 \), \( d_j/(d_j^2+\lambda)\to1/d_j \), and
\( \sum_{j\le r}(c_j/d_j)\bv_j=\X^+\y \) by @thm-cmp-svd-ls(b).
(e) With \( \sigma^2 \) known, the posterior density is proportional to
\( \exp\{-(\norm{\y-\X\bbeta}^2+\lambda\norm{\bbeta}^2)/(2\sigma^2)\} \). Completing the square as in the
proof of @thm-opt-bayes-conjugate, with \( \mathbf{P}_0=\lambda\I \) and \( \mathbf{m}_0=\bzero \),
\[
\norm{\y-\X\bbeta}^2+\lambda\norm{\bbeta}^2=(\bbeta-\hbeta^{\mathrm{R}}_\lambda)\T(\X\T\X+\lambda\I)(\bbeta-\hbeta^{\mathrm{R}}_\lambda)+C,
\]
with \( C \) free of \( \bbeta \). This is the stated normal density. The last statement is @prp-opt-shrinkage(d).
:::

Part (a) says what ridge does. It keeps the directions \( \bv_j \) and shrinks the coordinate along each
by \( f_j(\lambda) \): directions with \( d_j^2\gg\lambda \) are hardly touched, and those with
\( d_j^2\ll\lambda \) are shrunk almost to zero, whatever the size of the coefficients. Comparing with the oracle
factors \( f_j^*=d_j^2/(d_j^2+\sigma^2/\alpha_j^2) \) of @prp-shr-canonical(c), ridge is the oracle when all
canonical coefficients have the same size, \( \alpha_j^2=\sigma^2/\lambda \). Part (e) says the same in Bayesian
terms: the prior \( \Normal(\bzero,(\sigma^2/\lambda)\I) \) treats every direction alike.

The effective degrees of freedom in (b) count each dimension by the fraction \( f_j(\lambda) \) in which
it is used; [Chapter 29](../ch29-model-selection/index.html) uses them in model-selection criteria, and @exr-shr-df-covariance gives them
a second meaning.

::: {.remark}
[Penalized and constrained forms]

For any penalty \( P\ge0 \), a minimizer \( \hat{\bb} \) of \( \norm{\y-\X\bb}^2+\lambda P(\bb) \), with
\( \lambda>0 \), also minimizes \( \norm{\y-\X\bb}^2 \) subject to \( P(\bb)\le t \), where \( t=P(\hat{\bb}) \).
Indeed, if \( P(\bb)\le t \), then
\[
\begin{aligned}
\norm{\y-\X\bb}^2&\ge\norm{\y-\X\bb}^2+\lambda\{P(\bb)-t\}\\
&\ge\norm{\y-\X\hat{\bb}}^2+\lambda\{P(\hat{\bb})-t\}=\norm{\y-\X\hat{\bb}}^2 .
\end{aligned}
\]
So ridge is least squares restricted to a ball, \( \norm{\bb}^2\le t \), whose radius is fixed by
\( \lambda \) (@exr-shr-ridge-constrained). The lasso of [Section 27.4](04-lasso.html) replaces the ball by
a cross-polytope.
:::

## Ridge beats least squares

@exr-opt-ridge showed that the total mean squared error of ridge falls as \( \lambda \) leaves \( 0 \).
Hoerl and Kennard made the range explicit, and Theobald (1974) added the matrix sense.

::: {#thm-shr-ridge-dominates}
[Hoerl–Kennard]

Let \( \X \) have full column rank, \( \E(\Y)=\X\bbeta \), \( \Cov(\Y)=\sigma^2\I \) with \( \sigma^2>0 \), and let
\( \alpha_{\max}^2=\max_j\alpha_j^2 \) with \( \boldsymbol{\upalpha}=\V\T\bbeta \). Let \( m(\lambda) \) be the total mean squared
error @eq-shr-ridge-mse, so that \( m(0)=\sigma^2\sum_jd_j^{-2} \) is that of least squares.

::: {.enumerate options="label=(\alph*)"}
1. \( m(\lambda)<m(0) \) for every \( 0<\lambda<2\sigma^2/\alpha_{\max}^2 \), and in particular for every
   \( 0<\lambda<2\sigma^2/\norm{\bbeta}^2 \). If \( \bbeta=\bzero \), for every \( \lambda>0 \).

2. \( m \) is strictly decreasing on \( [0,\sigma^2/\alpha_{\max}^2] \).

3. \( \hbeta^{\mathrm{R}}_\lambda \) is at least as good as \( \hbeta \) in the matrix sense iff
   \[
   \frac{\lambda}{\sigma^2}\sum_j\frac{d_j^2\,\alpha_j^2}{2d_j^2+\lambda}\le1 ,
   \]
   which holds whenever \( 0<\lambda\le2\sigma^2/\norm{\bbeta}^2 \).
:::
:::

::: {.proof}
(a) Compare term by term. The \( j \)th term of @eq-shr-ridge-mse is below its value \( \sigma^2/d_j^2 \)
at \( \lambda=0 \) iff
\( d_j^2(\sigma^2d_j^2+\lambda^2\alpha_j^2)<\sigma^2(d_j^2+\lambda)^2 \), that is, iff
\( \lambda^2d_j^2\alpha_j^2<2\lambda\sigma^2d_j^2+\lambda^2\sigma^2 \), or, dividing by \( \lambda>0 \),
\[
\lambda\,(d_j^2\alpha_j^2-\sigma^2)<2\sigma^2d_j^2 .
\]
If \( d_j^2\alpha_j^2\le\sigma^2 \), this holds for every \( \lambda>0 \). Otherwise it holds for
\( \lambda<2\sigma^2d_j^2/(d_j^2\alpha_j^2-\sigma^2) \), and the right side exceeds
\( 2\sigma^2d_j^2/(d_j^2\alpha_j^2)=2\sigma^2/\alpha_j^2\ge2\sigma^2/\alpha_{\max}^2 \). So for
\( 0<\lambda<2\sigma^2/\alpha_{\max}^2 \) every term decreases, and so does the sum. Finally
\( \alpha_{\max}^2\le\norm{\boldsymbol{\upalpha}}^2=\norm{\bbeta}^2 \).
(b) By @exr-opt-ridge, in the present notation,
\( m'(\lambda)=2\sum_jd_j^2(\lambda\alpha_j^2-\sigma^2)/(d_j^2+\lambda)^3 \). For \( 0\le\lambda<\sigma^2/\alpha_{\max}^2 \) every numerator is negative, so \( m'<0 \) there, and
\( m \) is strictly decreasing on the closed interval.
(c) Write \( \bS=\X\T\X \) and \( \A=(\bS+\lambda\I)^{-1}\X\T \). In the notation of
@thm-shr-bias-variance(c), all matrices involved have the eigenvectors \( \bv_j \), with eigenvalues
\[
\boldsymbol{\Delta}_{\A}:\ \sigma^2\Bigl(\frac1{d_j^2}-\frac{d_j^2}{(d_j^2+\lambda)^2}\Bigr)=\sigma^2\frac{\lambda(2d_j^2+\lambda)}{d_j^2(d_j^2+\lambda)^2}>0,
\]
and the bias is \( \boldsymbol{\updelta}=-\sum_j\lambda\alpha_j(d_j^2+\lambda)^{-1}\bv_j \). So \( \boldsymbol{\Delta}_{\A} \) is positive
definite and
\( \boldsymbol{\updelta}\T\boldsymbol{\Delta}_{\A}^{-1}\boldsymbol{\updelta}=\sum_j\lambda^2\alpha_j^2(d_j^2+\lambda)^{-2}\cdot d_j^2(d_j^2+\lambda)^2/\{\sigma^2\lambda(2d_j^2+\lambda)\} \),
which is the stated sum. Each factor \( d_j^2/(2d_j^2+\lambda) \) is below \( 1/2 \), so the sum is at most
\( \lambda\norm{\boldsymbol{\upalpha}}^2/(2\sigma^2)=\lambda\norm{\bbeta}^2/(2\sigma^2) \). Apply @thm-shr-bias-variance(c).
:::

The theorem is an existence result. The good range of \( \lambda \) depends on the unknown \( \bbeta \) and
\( \sigma^2 \), and it says nothing about a \( \lambda \) chosen from the data. Replacing \( \bbeta \) by \( \hbeta \) in the
bound is misleading, because \( \norm{\hbeta}^2 \) overestimates \( \norm{\bbeta}^2 \) by \( \sigma^2\sum_jd_j^{-2} \) on
average (@exr-shr-bias-length), which is large exactly when ridge is needed. What the theorem does
establish is that least squares is never the best member of the ridge family.

## Choosing the ridge parameter

Hoerl and Kennard suggested choosing the smallest \( \lambda \) at which the ridge trace "stabilizes", a
subjective judgement. Hoerl, Kennard and Baldwin (1975) proposed the plug-in
\( \lambda=p\hat{\sigma}^2/\norm{\hbeta}^2 \), the optimum of @exr-shr-orthonormal-ridge with estimates
substituted; since \( \norm{\hbeta}^2 \) is too large on average, it tends to shrink too little under
collinearity. A more defensible choice minimizes an estimate of prediction error. Because
ridge is a linear smoother, \( \hat{\y}=\bH_\lambda\y \), the leave-one-out prediction errors need only one
fit: the residual of case \( i \) divided by \( 1-h_{ii}(\lambda) \), with \( h_{ii}(\lambda) \) the \( i \)th diagonal entry
of \( \bH_\lambda \) (@exr-shr-ridge-loo; add \( 1/n \) for an unpenalized intercept), the ridge version of @thm-res-deletion. **Generalized cross-validation** replaces each \( h_{ii}(\lambda) \) by the average
\( \tr\bH_\lambda/n \):
\[
\mathrm{GCV}(\lambda)=\frac{n^{-1}\norm{\y-\bH_\lambda\y}^2}{\bigl(1-\tr\bH_\lambda/n\bigr)^2}
\]
(with an unpenalized intercept, \( \tr\bH_\lambda \) becomes \( 1+\mathrm{df}(\lambda) \)). [Chapter 29](../ch29-model-selection/index.html) develops both criteria for general linear smoothers (@thm-sel-loocv and @def-sel-gcv). A Bayesian alternative estimates
the prior variance \( \sigma^2/\lambda \) from the marginal likelihood, an empirical Bayes idea that returns in
[Section 27.5](05-james-stein.html).

::: {#exm-shr-longley-ridge}
[Ridge regression for Longley's data]

Longley's regression (@exm-dep-longley) relates total employment, in thousands, to six
economic series over \( n=16 \) years. With the regressors standardized, the eigenvalues
\( d_j^2 \) of \( \X\T\X \) run from \( 69.05 \) down to \( 0.0057 \). The least squares coefficients, in
thousands of persons per standard deviation, are \( 163 \) for the GNP deflator, \( -3560 \) for GNP,
\( -1888 \) for unemployment, \( -719 \) for the armed forces, \( -355 \) for population and \( 8709 \) for
the year. The GNP and YEAR coefficients are large and of opposite sign: the two series rose
together, and the data can hardly separate them.

[Figure 27.2.1](#fig-shr-ridge)(a) shows the ridge trace against the effective degrees of
freedom. As \( \lambda \) grows, the GNP and YEAR coefficients approach each other, and the GNP
coefficient changes sign at \( \lambda=0.0110 \), where \( \mathrm{df}=5.1 \). Leave-one-out
cross-validation chooses \( \lambda=0.0026 \), with \( \mathrm{df}=5.61 \), and GCV chooses
\( \lambda=0.0030 \), with \( \mathrm{df}=5.56 \). Both prefer little shrinkage. At the leave-one-out
choice the coefficients are \( -6 \), \( -1809 \), \( -1645 \), \( -669 \), \( -825 \) and \( 7415 \). The
largest changes are in the coefficients of GNP, YEAR and population, the three series that rose
together, and the deflator coefficient is essentially zero. The minimized root mean squared
leave-one-out error is \( 402 \), against \( 425 \) for least squares, but this minimum is optimistic, since
the same errors chose \( \lambda \). Choosing \( \lambda \) afresh within each leave-one-out fold (*nested*
cross-validation) gives an honest \( 441 \): on these sixteen years tuned ridge does not predict better
than least squares.

Panel (b) asks what the theorem guarantees in a world where the truth is the least squares fit,
\( \bbeta=\hbeta \), with \( \sigma=s=304.9 \). There, least squares has total mean squared error
\( 19.33 \) million, and the best ridge estimate, at \( \lambda=0.0019 \) with \( \mathrm{df}=5.69 \), has
\( 15.01 \) million, a ratio of \( 0.777 \). @thm-shr-ridge-dominates(a) guarantees an improvement
only for \( \lambda<0.00383 \), and part (c) gives dominance in the matrix sense for
\( \lambda\le0.00200 \); in fact the improvement continues up to \( \lambda=0.0058 \). The ranges are
narrow because this "truth" inherits the least squares estimate's large components along the
weak directions.
:::

::: {when-format="html"}
![**Figure 27.2.1.** Ridge regression for Longley's data. (a) The ridge trace: standardized
coefficients against the effective degrees of freedom \( \mathrm{df}(\lambda) \); the vertical line
marks the leave-one-out choice. (b) Variance, squared bias and total mean squared error of ridge,
computed exactly from @eq-shr-ridge-mse when the true coefficients are the least squares
estimates and \( \sigma=s \); the dashed line is least squares.](ridge_longley.svg){#fig-shr-ridge width=100%}
:::

::: {when-format="pdf"}
![Ridge regression for Longley's data. (a) The ridge trace: standardized
coefficients against the effective degrees of freedom \( \mathrm{df}(\lambda) \); the vertical line
marks the leave-one-out choice. (b) Variance, squared bias and total mean squared error of ridge,
computed exactly from @eq-shr-ridge-mse when the true coefficients are the least squares
estimates and \( \sigma=s \); the dashed line is least squares.](ridge_longley.pdf){width=100%}
:::

```{.python .run #cell-ridge-longley-setup}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.longley.load_pandas().data
names = ["GNPDEFL", "GNP", "UNEMP", "ARMED", "POP", "YEAR"]
Z = data[names].to_numpy()
y = data["TOTEMP"].to_numpy()
n, p = Z.shape

X = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)     # centred, unit standard deviation
yc = y - y.mean()                                     # centred response (the intercept is ybar)

U, d, Vt = np.linalg.svd(X, full_matrices=False)      # X = U diag(d) V^T
V = Vt.T
c = U.T @ yc                                          # coordinates of y along u_1, ..., u_p

def ridge(lam):
    """Ridge coefficients for the centred problem, through the SVD."""
    return V @ (d / (d ** 2 + lam) * c)

def edf(lam):
    """Effective degrees of freedom of the slopes: sum of the shrinkage factors."""
    return np.sum(d ** 2 / (d ** 2 + lam))

beta_ls = V @ (c / d)
print("singular values squared:", np.round(d ** 2, 4))
print("least squares:", np.round(beta_ls, 1))
for lam in [0.01, 0.1, 1.0]:
    print(f"lambda = {lam:5.2f}  edf = {edf(lam):.2f}  ridge =", np.round(ridge(lam), 1))
```

The second cell computes the leverages of the ridge smoother, the leave-one-out and GCV
criteria, and the two choices of \( \lambda \).

```{.python .run #cell-ridge-longley-cv}
def hat_diag(lam):
    """Leverages of the ridge smoother with an unpenalized intercept."""
    f = d ** 2 / (d ** 2 + lam)
    return 1 / n + np.sum(U ** 2 * f, axis=1)

def loo_and_gcv(lam):
    resid = yc - X @ ridge(lam)
    h = hat_diag(lam)
    loo = np.mean((resid / (1 - h)) ** 2)                  # leave-one-out shortcut
    gcv = np.mean(resid ** 2) / (1 - (1 + edf(lam)) / n) ** 2
    return loo, gcv

grid = np.logspace(-4, 2, 601)
scores = np.array([loo_and_gcv(l) for l in grid])
lam_loo = grid[np.argmin(scores[:, 0])]
lam_gcv = grid[np.argmin(scores[:, 1])]
print(f"LOO chooses lambda = {lam_loo:.4f} (edf {edf(lam_loo):.2f})")
print(f"GCV chooses lambda = {lam_gcv:.4f} (edf {edf(lam_gcv):.2f})")
print("ridge at the LOO choice:", np.round(ridge(lam_loo), 1))
```

The script also checks the closed form against a direct solve, augmented least squares (@exr-shr-augmented)
and the posterior mean of @thm-opt-bayes-conjugate, and the leave-one-out
shortcut against sixteen refits.

## Exercises

### A. Check your understanding

::: {#exr-shr-augmented}
[A1]

Show that \( \hbeta^{\mathrm{R}}_\lambda \) is the least squares estimate for the augmented data
\( \tilde{\y}=(\y\T,\bzero\T)\T \) and \( \tilde{\X}=(\X\T,\sqrt{\lambda}\,\I)\T \). Show that
\( \kappa(\tilde{\X})^2=(d_1^2+\lambda)/(d_p^2+\lambda) \), so ridge can be computed stably by QR
([Chapter 10](../ch10-computation/index.html)), and that the augmented problem has "residual sum of squares"
\( \norm{\y-\X\hbeta^{\mathrm{R}}_\lambda}^2+\lambda\norm{\hbeta^{\mathrm{R}}_\lambda}^2 \).
:::

::: {.solution}
\( \norm{\tilde{\y}-\tilde{\X}\bb}^2=\norm{\y-\X\bb}^2+\norm{\bzero-\sqrt{\lambda}\bb}^2 \), which is
@eq-shr-ridge-criterion, so the minimizers coincide and the minimum is the stated sum.
\( \tilde{\X}\T\tilde{\X}=\X\T\X+\lambda\I \) has eigenvalues \( d_j^2+\lambda \), and the singular values of
\( \tilde{\X} \) are their square roots (@thm-mat-svd). The condition number is the ratio of the
largest to the smallest (@def-mat-condition-number).
:::

::: {#exr-shr-orthonormal-ridge}
[A2]

Suppose \( \X\T\X=\I \). Show that \( \hbeta^{\mathrm{R}}_\lambda=\hbeta/(1+\lambda) \),
\( \mathrm{df}(\lambda)=p/(1+\lambda) \), and that the value of \( \lambda \) minimizing the total mean squared
error is \( p\sigma^2/\norm{\bbeta}^2 \).
:::

### B. Practice

::: {#exr-shr-ridge-loo}
[B1]

Let \( \hbeta^{\mathrm{R}}_{\lambda,(i)} \) be the ridge estimate computed without case \( i \) (no intercept),
and \( h_{ii} \) the \( i \)th diagonal entry of \( \bH_\lambda \). Show that
\[
y_i-\x_{(i)}\T\hbeta^{\mathrm{R}}_{\lambda,(i)}=\frac{y_i-\x_{(i)}\T\hbeta^{\mathrm{R}}_\lambda}{1-h_{ii}} .
\]
:::

::: {.solution}
Put \( \mathbf{K}=\X\T\X+\lambda\I \) and \( \x=\x_{(i)} \), so \( h_{ii}=\x\T\mathbf{K}^{-1}\x<1 \). Without case \( i \)
the matrix is \( \mathbf{K}-\x\x\T \) and the right side is \( \X\T\y-\x y_i \). By @eq-mat-sherman-morrison,
\( (\mathbf{K}-\x\x\T)^{-1}=\mathbf{K}^{-1}+\mathbf{K}^{-1}\x\x\T\mathbf{K}^{-1}/(1-h_{ii}) \). Hence
\[
\x\T\hbeta^{\mathrm{R}}_{\lambda,(i)}=\Bigl(\x\T\mathbf{K}^{-1}+\frac{h_{ii}\x\T\mathbf{K}^{-1}}{1-h_{ii}}\Bigr)(\X\T\y-\x y_i)
=\frac{\x\T\hbeta^{\mathrm{R}}_\lambda-h_{ii}y_i}{1-h_{ii}},
\]
and \( y_i \) minus this is \( (y_i-\x\T\hbeta^{\mathrm{R}}_\lambda)/(1-h_{ii}) \). The argument is that of
@thm-res-deletion with \( \X\T\X \) replaced by \( \mathbf{K} \); only the rank-one update is used.
:::

::: {#exr-shr-generalized-ridge}
[B2]

**Generalized ridge regression** uses a separate constant for each principal direction:
\( \hbeta_{\mathbf{k}}=(\X\T\X+\V\mathbf{K}\V\T)^{-1}\X\T\y \) with \( \mathbf{K}=\diag(k_1,\dots,k_p) \), \( k_j\ge0 \).
Show that its shrinkage factors are \( d_j^2/(d_j^2+k_j) \), and that the choice
\( k_j=\sigma^2/\alpha_j^2 \) gives the oracle factors of @prp-shr-canonical(c).
:::

::: {#exr-shr-ridge-constrained}
[B3]

Show that for every \( t \) with \( 0<t<\norm{\hbeta}^2 \) (full column rank) there is exactly one
\( \lambda>0 \) with \( \norm{\hbeta^{\mathrm{R}}_\lambda}^2=t \), and that \( \hbeta^{\mathrm{R}}_\lambda \) is then the
unique minimizer of \( \norm{\y-\X\bb}^2 \) subject to \( \norm{\bb}^2\le t \).
:::

### C. Going deeper

::: {#exr-shr-df-covariance}
[C1]

For any estimator \( \hat{\bmu}=\bS\Y \) of \( \bmu=\E(\Y) \) with \( \Cov(\Y)=\sigma^2\I \), show that
\( \sigma^{-2}\sum_i\Cov(\hat{\mu}_i,Y_i)=\tr\bS \). Deduce that \( \mathrm{df}(\lambda) \) measures how
strongly the ridge fit follows its own data, and that
\( \E\norm{\Y-\bH_\lambda\Y}^2=\E\norm{\bH_\lambda\Y-\bmu}^2+n\sigma^2-2\sigma^2\mathrm{df}(\lambda) \).
:::

::: {#exr-shr-ridge-wide}
[C2]

Let \( p>n \). Show that \( \hbeta^{\mathrm{R}}_\lambda=\X\T(\X\X\T+\lambda\I_n)^{-1}\y \), which needs only an
\( n\times n \) inverse, and that \( \hbeta^{\mathrm{R}}_\lambda\in\C(\X\T) \). What is the limit as \( \lambda\to0 \)?
[Chapter 28](../ch28-high-dimensional/index.html) studies this limit.
:::

::: {.solution}
The identity \( \X\T(\X\X\T+\lambda\I_n)=(\X\T\X+\lambda\I_p)\X\T \) gives
\( (\X\T\X+\lambda\I_p)^{-1}\X\T=\X\T(\X\X\T+\lambda\I_n)^{-1} \) after multiplying by the two inverses. The
right side is in \( \C(\X\T) \). By @thm-shr-ridge(d) the limit is \( \X^+\y \), the least squares
solution of smallest norm, which interpolates the data when \( \rank(\X)=n \).
:::
