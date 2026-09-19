# Extra sums of squares

Most questions asked of a regression compare two models. Each question names a reduced model inside
the full one, and the natural measure of what the full model adds is how much it lowers the
residual sum of squares. This section develops that measure, its computation when regressors are
added to a fitted model, and a closed form when the reduced model is defined by linear constraints.

## Notation

::: {#def-ss-extra}
[Extra sum of squares]

Let \( \X_1 \) and \( \X_2 \) be matrices with \( n \) rows, let \( \X=[\X_1,\X_2] \), and let \( \M_1 \)
and \( \M \) be the projections onto \( \C(\X_1) \) and \( \C(\X) \). The **extra sum of squares** of
\( \X_2 \) given \( \X_1 \) is
\[
\text{SS}(\X_2\mid\X_1)=\y\T(\M-\M_1)\y=\text{SSE}(\X_1)-\text{SSE}(\X),
\]
with \( \rank(\X)-\rank(\X_1) \) degrees of freedom. Here \( \text{SSE}(\cdot) \) is the residual sum
of squares of the model with the given model matrix.
:::

The two expressions agree by @eq-proj-extra-ss. One also writes \( \text{SS}(\bbeta_2\mid\bbeta_1) \),
or \( \text{SS}(B\mid\mu,A) \) in factor models; some books write \( R \), for “reduction”. Four
properties follow from @thm-proj-nested and @thm-proj-reparam.

1. *It depends only on the two column spaces.* Rescaling, recoding, centring, replacing \( \X_2 \) by
   \( \X_2+\X_1\B \), and rank deficiency change nothing.

2. *It is nonnegative*, and zero iff \( \M\y=\M_1\y \), that is, iff adding \( \X_2 \) leaves the
   fitted values unchanged.

3. *It obeys a chain rule*:
   \( \text{SS}(\X_2,\X_3\mid\X_1)=\text{SS}(\X_2\mid\X_1)+\text{SS}(\X_3\mid\X_1,\X_2) \),
   because the projections telescope.

4. *It is the sum of squares of a residualized regression.* The next lemma makes this precise.

::: {#lem-ss-residualized}
[Extra sums of squares by residualizing]

For any \( \X_1 \) and \( \X_2 \), let \( \tilde{\X}_2=(\I-\M_1)\X_2 \). Then
\( \C([\X_1,\X_2])=\C(\X_1)\dirsum\C(\tilde{\X}_2) \) with the two summands orthogonal, so
\( \M=\M_1+\bP_{\C(\tilde{\X}_2)} \) and
\[
\text{SS}(\X_2\mid\X_1)=\norm{\bP_{\C(\tilde{\X}_2)}\y}^2
=\y\T\tilde{\X}_2(\tilde{\X}_2\T\tilde{\X}_2)\ginv\tilde{\X}_2\T\y .
\]
:::

::: {.proof}
The argument for part (a) of @lem-proj-fwl-split uses no rank condition, so it applies here:
the columns of \( \tilde{\X}_2 \) lie in \( \C(\X) \) and are orthogonal to \( \C(\X_1) \), and every
\( \X_1\bb_1+\X_2\bb_2 \) equals \( \tilde{\X}_2\bb_2+(\X_1\bb_1+\M_1\X_2\bb_2) \). By @thm-proj-sum,
\( \M=\M_1+\bP_{\C(\tilde{\X}_2)} \). The formula is @thm-proj-M-formula applied to
\( \tilde{\X}_2 \).
:::

## Adding regressors to a fitted model

Suppose \( \E(\Y)=\X\bbeta \) has been fitted and we add the columns of an \( n\times t \) matrix
\( \mathbf{Z} \), giving \( \E(\Y)=\X\bbeta+\mathbf{Z}\bgamma \) with model matrix \( \W=[\X,\mathbf{Z}] \). Only a
\( t\times t \) matrix has to be inverted.

::: {#thm-ss-adding}
[Adding regressors]

Let \( \W=[\X,\mathbf{Z}] \) have full column rank \( p+t \), let \( \mathbf{C}=(\X\T\X)^{-1} \),
\( \hbeta=\mathbf{C}\X\T\y \), and \( \R=\I-\X\mathbf{C}\X\T \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \mathbf{Z}\T\R\mathbf{Z} \) is positive definite;

2. the least squares coefficients of the enlarged model are
   \[
\hat{\bgamma}=(\mathbf{Z}\T\R\mathbf{Z})^{-1}\mathbf{Z}\T\R\y,\qquad
\hbeta_W=\hbeta-\bL\hat{\bgamma},\qquad \bL=\mathbf{C}\X\T\mathbf{Z} ;
\]

3. the residual sum of squares falls by
   \[
\text{SS}(\mathbf{Z}\mid\X)=\y\T\R\mathbf{Z}(\mathbf{Z}\T\R\mathbf{Z})^{-1}\mathbf{Z}\T\R\y=\hat{\bgamma}\T\mathbf{Z}\T\R\y
=\hat{\bgamma}\T(\mathbf{Z}\T\R\mathbf{Z})\hat{\bgamma} ;
\]

4. if \( \Cov(\Y)=\sigma^2\I \), then \( \Cov(\hat{\bgamma})=\sigma^2(\mathbf{Z}\T\R\mathbf{Z})^{-1} \),
   \( \Cov(\hbeta_W,\hat{\bgamma})=-\bL\Cov(\hat{\bgamma}) \) and
   \[
\Cov(\hbeta_W)=\sigma^2\bigl[\mathbf{C}+\bL(\mathbf{Z}\T\R\mathbf{Z})^{-1}\bL\T\bigr].
\]
:::

:::

::: {.proof}
(a) \( \R \) is the symmetric idempotent residual projection of the model \( \X \), so
\( \mathbf{a}\T\mathbf{Z}\T\R\mathbf{Z}\mathbf{a}=\norm{\R\mathbf{Z}\mathbf{a}}^2\ge0 \). If \( \R\mathbf{Z}\mathbf{a}=\bzero \), then
\( \mathbf{Z}\mathbf{a}\in\C(\X) \), say \( \mathbf{Z}\mathbf{a}=\X\bb \), so \( \W(-\bb\T,\mathbf{a}\T)\T=\bzero \) and full
column rank forces \( \mathbf{a}=\bzero \).

(b) Apply the Frisch–Waugh–Lovell theorem (@thm-proj-fwl) with \( \X_1=\mathbf{Z} \) and
\( \X_2=\X \). Then \( \M_2=\I-\R \) and \( \tilde{\X}_1=\R\mathbf{Z} \), so part (a) of that theorem gives
\( \hat{\bgamma}=(\mathbf{Z}\T\R\R\mathbf{Z})^{-1}\mathbf{Z}\T\R\y \), which is the stated formula because \( \R \) is
symmetric and idempotent. Part (c) of the theorem gives
\( \hbeta_W=\mathbf{C}\X\T(\y-\mathbf{Z}\hat{\bgamma})=\hbeta-\bL\hat{\bgamma} \).

(c) By @lem-ss-residualized, the extra sum of squares is the squared length of the projection of
\( \y \) onto \( \C(\R\mathbf{Z}) \), which is
\( \y\T\R\mathbf{Z}(\mathbf{Z}\T\R\mathbf{Z})^{-1}\mathbf{Z}\T\R\y \) by @thm-proj-M-formula. Substituting (b) gives the
other two forms.

(d) Both \( \hat{\bgamma} \) and \( \hbeta \) are linear in \( \y \), and by @thm-rv-linear
\[
\Cov(\hat{\bgamma})=\sigma^2(\mathbf{Z}\T\R\mathbf{Z})^{-1}\mathbf{Z}\T\R\R\mathbf{Z}(\mathbf{Z}\T\R\mathbf{Z})^{-1}=\sigma^2(\mathbf{Z}\T\R\mathbf{Z})^{-1},
\qquad
\Cov(\hbeta,\hat{\bgamma})=\sigma^2\mathbf{C}\X\T\R\mathbf{Z}(\mathbf{Z}\T\R\mathbf{Z})^{-1}=\bzero ,
\]
since \( \X\T\R=\bzero \). With \( \Cov(\hbeta)=\sigma^2\mathbf{C} \) (@thm-lm-moments), the formula
\( \hbeta_W=\hbeta-\bL\hat{\bgamma} \) gives
\( \Cov(\hbeta_W,\hat{\bgamma})=-\bL\Cov(\hat{\bgamma}) \) and
\( \Cov(\hbeta_W)=\Cov(\hbeta)+\bL\Cov(\hat{\bgamma})\bL\T \).
:::

Since \( \bL(\mathbf{Z}\T\R\mathbf{Z})^{-1}\bL\T \) is nonnegative definite, part (d) shows that *adding
regressors never decreases the variances of the coefficients already in the model*, when both
models are correct (compare @exr-mat-adding-regressors, and see @exr-ss-no-variance-change). [Chapter 26](../ch26-collinearity/index.html) measures this price of adjustment as
variance inflation (@def-col-vif). In computation one updates a QR factorization instead of using these formulas
([Chapter 10](../ch10-computation/index.html)).

::: {#exm-ss-adding-crime}
[Adding poverty and urbanization]

In the state data, fit the murder rate on single parenthood first:
\( \hat{y}=-8.2146+0.5147\,x_{\text{single}} \), with residual sum of squares
\( 121.92 \). Now add \( \mathbf{Z}=[\x_{\text{poverty}},\x_{\text{urban}}] \). The listing computes
\( \hat{\bgamma}=(0.2538,\,0.0046) \) from the old residual projection alone, and
corrects the old coefficients to \( -9.0946 \) and \( 0.3980 \). These agree with a fresh fit of
the enlarged model. The residual sum of squares falls by \( 19.96 \), on \( 2 \) degrees of
freedom, to \( 101.96 \). Single parenthood's coefficient drops by almost a quarter once poverty
is allowed to take its share, the pattern seen in @exm-proj-sequential-order.
:::

```{.python .run #cell-extra-ss-adding}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
n = len(y)
one = np.ones(n)
cols = {"poverty": data["poverty"].to_numpy(), "single": data["single"].to_numpy(),
        "urban": data["urban"].to_numpy()}

X = np.column_stack([one, cols["single"]])                 # the model already fitted
Z = np.column_stack([cols["poverty"], cols["urban"]])       # the regressors to add
XtX_inv = np.linalg.inv(X.T @ X)
beta_X = XtX_inv @ X.T @ y
R = np.eye(n) - X @ XtX_inv @ X.T                           # residual projection of the old model

gamma = np.linalg.solve(Z.T @ R @ Z, Z.T @ R @ y)          # new coefficients
beta_new = beta_X - XtX_inv @ X.T @ Z @ gamma               # old coefficients, corrected
extra_ss = gamma @ Z.T @ R @ y                              # SS(Z | X) = SSE_X - SSE_W
print("gamma =", gamma.round(4), " corrected beta =", beta_new.round(4))
print(f"SSE_X = {y @ R @ y:.2f},  extra SS = {extra_ss:.2f}")
```

## One column at a time

With \( t=1 \), @thm-ss-adding links the extra sum of squares of a single regressor to its
coefficient and standard error.

::: {#cor-ss-single-column}
[The extra sum of squares of one column]

Let \( \X \) have full column rank \( p \), write \( c_{jj} \) for the \( j \)th diagonal entry of
\( (\X\T\X)^{-1} \), and let \( \X_{(j)} \) be \( \X \) without its \( j \)th column. Then
\[
\text{SS}(\x_j\mid\X_{(j)})=\frac{\hat{\beta}_j^2}{c_{jj}}=t_j^2\,s^2,
\qquad t_j=\frac{\hat{\beta}_j}{s\sqrt{c_{jj}}},\quad s^2=\frac{\text{SSE}}{n-p}.
\]
Consequently, the ratio of this one-degree-of-freedom sum of squares to the residual mean
square is \( t_j^2 \).
:::

::: {.proof}
Permute the columns so that \( \x_j \) comes last; this carries \( c_{jj} \) to the last diagonal
entry of the inverse. Write \( \X=[\X_{(j)},\x_j] \) and \( \R=\I-\X_{(j)}(\X_{(j)}\T\X_{(j)})^{-1}\X_{(j)}\T \)
for the residual projection of \( \X_{(j)} \). By @thm-mat-partitioned-inverse the last diagonal entry of
\( (\X\T\X)^{-1} \) is the inverse of the Schur complement
\( \x_j\T\x_j-\x_j\T\X_{(j)}(\X_{(j)}\T\X_{(j)})^{-1}\X_{(j)}\T\x_j=\x_j\T\R\x_j \), so
\( \x_j\T\R\x_j=1/c_{jj} \). Now apply @thm-ss-adding with \( \X_{(j)} \) as the old model and
\( \mathbf{Z}=\x_j \): then \( \hat{\gamma}=\hat{\beta}_j \), and part (c) gives
\( \text{SS}(\x_j\mid\X_{(j)})=\hat{\beta}_j^2\,\x_j\T\R\x_j=\hat{\beta}_j^2/c_{jj} \). The
definition of \( t_j \) gives the last form.
:::

For poverty in the three-regressor model, \( \hat{\beta}=0.2538 \) and
\( c_{jj}=0.003935 \), so the extra sum of squares is \( 0.2538^2/0.003935=16.37 \).
This is exactly the last entry for poverty in Order B of @exm-proj-sequential-order. The
\( t \) statistic is \( 2.718 \), and \( t^2s^2=7.386\times2.217=16.37 \).

```{.python .run #cell-extra-ss-single}
Xf = np.column_stack([one, cols["poverty"], cols["single"], cols["urban"]])
fit = sm.OLS(y, Xf).fit()
C = np.linalg.inv(Xf.T @ Xf)
j = 1                                                       # poverty
partial_ss = fit.params[j] ** 2 / C[j, j]
print(f"beta = {fit.params[j]:.4f}, c_jj = {C[j, j]:.6f}, beta^2/c_jj = {partial_ss:.2f}")
print(f"t = {fit.tvalues[j]:.3f},  t^2 * s^2 = {fit.tvalues[j] ** 2 * fit.scale:.2f}")
```

## Sums of squares for linear hypotheses

Reduced models often come as constraints: two coefficients are equal, a set of coefficients
vanishes, a contrast among group means is zero. [Section 6.5](../ch06-projections/05-nested.html)
identified the test space of such a constraint (@thm-proj-constraint-space). The next theorem
gives its sum of squares in terms of the least squares estimate, with no reduced model fitted.

::: {#thm-ss-hypothesis}
[Sum of squares of a linear hypothesis]

Let \( \bLambda \) be a \( p\times q \) matrix of rank \( q \) such that each entry of \( \bLambda\T\bbeta \)
is estimable, so that \( \bLambda=\X\T\bT \) for some \( n\times q \) matrix \( \bT=[\boldsymbol{\uprho}_1,\dots,\boldsymbol{\uprho}_q] \)
(@def-est-estimable, @thm-proj-invariant-functions). Let \( \G \) be any generalized inverse of
\( \X\T\X \), let \( \hbeta \) be any least squares estimate, and let \( \mathbf{d}\in\Real^q \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \bLambda\T\hbeta \) and \( \bLambda\T\G\bLambda \) do not depend on the choices of \( \hbeta \)
   and \( \G \), and \( \bLambda\T\G\bLambda \) is positive definite;

2. the subspace \( \mathcal S_0=\{\X\bb:\bLambda\T\bb=\bzero\} \) of \( \C(\X) \) has orthogonal
   complement in \( \C(\X) \) equal to \( \C(\M\bT) \), of dimension \( q \), and the projection onto
   \( \C(\M\bT) \) is \( \M\bT(\bLambda\T\G\bLambda)^{-1}\bT\T\M \);

3. the minimum of \( \norm{\y-\X\bb}^2 \) over \( \{\bb:\bLambda\T\bb=\mathbf{d}\} \) exceeds
   \( \text{SSE} \) by
   \[
\text{SS}_H=(\bLambda\T\hbeta-\mathbf{d})\T(\bLambda\T\G\bLambda)^{-1}(\bLambda\T\hbeta-\mathbf{d}).
\]{#eq-ss-hypothesis}

:::

:::

::: {.proof}
(a) \( \bLambda\T\hbeta=\bT\T\X\hbeta=\bT\T\M\y \) and
\( \bLambda\T\G\bLambda=\bT\T\X\G\X\T\bT=\bT\T\M\bT \) by @thm-proj-M-formula, and neither
expression involves the choices. Since \( \X\T\M=\X\T \), we have
\( \X\T(\M\bT)=\bLambda \), which has rank \( q \), so the \( q \) columns of \( \M\bT \) are linearly
independent. Hence \( \bT\T\M\bT=(\M\bT)\T(\M\bT) \) is positive definite.

(b) \( \bLambda\T\bb=\bT\T\X\bb \), so \( \mathcal S_0=\C(\X)\cap\Null(\bT\T) \). @thm-proj-constraint-space
identifies its orthogonal complement in \( \C(\X) \) as \( \C(\M\bT) \). The columns of \( \M\bT \) are
independent, so @thm-proj-M-formula gives the projection
\( \M\bT(\bT\T\M\bT)^{-1}\bT\T\M \), and \( \bT\T\M\bT=\bLambda\T\G\bLambda \) by (a).

(c) First let \( \mathbf{d}=\bzero \). Minimizing over \( \bLambda\T\bb=\bzero \) is minimizing
\( \norm{\y-\bmu}^2 \) over \( \bmu\in\mathcal S_0 \), so by @eq-proj-extra-ss the excess over
\( \text{SSE} \) is the squared length of the projection of \( \y \) onto
\( \C(\X)\ominus\mathcal S_0=\C(\M\bT) \), where \( \ominus \) denotes the orthogonal complement
within \( \C(\X) \). By (b) this is
\( \y\T\M\bT(\bLambda\T\G\bLambda)^{-1}\bT\T\M\y \), and \( \bT\T\M\y=\bLambda\T\hbeta \). For
general \( \mathbf{d} \), the rows of \( \bLambda\T \) are independent, so there is a \( \bb_0 \) with
\( \bLambda\T\bb_0=\mathbf{d} \). Substituting \( \bb=\bb_0+\bb_* \) turns the problem into minimizing
\( \norm{\y_*-\X\bb_*}^2 \) subject to \( \bLambda\T\bb_*=\bzero \), with \( \y_*=\y-\X\bb_0 \). The
data \( \y_* \) have the same residual sum of squares as \( \y \), and \( \hbeta-\bb_0 \) is a least
squares estimate for them, since \( \X(\hbeta-\bb_0)=\M\y_* \). Applying the case
\( \mathbf{d}=\bzero \) to \( \y_* \) gives @eq-ss-hypothesis, because
\( \bLambda\T(\hbeta-\bb_0)=\bLambda\T\hbeta-\mathbf{d} \).
:::

In @eq-ss-hypothesis, \( \bLambda\T\hbeta-\mathbf{d} \) measures how far the estimate is from the
hypothesis, and \( \sigma^2\bLambda\T\G\bLambda \) is the covariance of \( \bLambda\T\hbeta \) ([Chapter 7](../ch07-optimality/index.html)), so
\( \text{SS}_H/\sigma^2 \) is the squared Mahalanobis length of the discrepancy
(@prp-mvn-mahalanobis). No reduced model has to be fitted.

In the full-rank case the constrained minimizer itself has a closed form, which recovers the
classical restricted least squares estimator.

::: {#thm-ss-restricted}
[Restricted least squares]

Let \( \X \) have full column rank, \( \mathbf{C}=(\X\T\X)^{-1} \), and let \( \bLambda \) (\( p\times q \)) have
rank \( q \). For any \( \mathbf{d}\in\Real^q \), the unique minimizer of \( \norm{\y-\X\bb}^2 \) subject to
\( \bLambda\T\bb=\mathbf{d} \) is
\[
\hbeta_H=\hbeta-\mathbf{C}\bLambda(\bLambda\T\mathbf{C}\bLambda)^{-1}(\bLambda\T\hbeta-\mathbf{d}),
\]
and \( \norm{\y-\X\hbeta_H}^2=\text{SSE}+\norm{\X(\hbeta-\hbeta_H)}^2 \), where the last term equals
\( \text{SS}_H \) of @eq-ss-hypothesis with \( \G=\mathbf{C} \).
:::

::: {.proof}
\( \bLambda\T\mathbf{C}\bLambda \) is positive definite because \( \mathbf{C} \) is and \( \bLambda \) has full
column rank, and \( \bLambda\T\hbeta_H=\bLambda\T\hbeta-(\bLambda\T\hbeta-\mathbf{d})=\mathbf{d} \). For any
\( \bb \), @eq-proj-distance-split with \( \bu=\X\bb \) gives
\( \norm{\y-\X\bb}^2=\text{SSE}+\norm{\X(\hbeta-\bb)}^2 \). Every feasible \( \bb \) has the form
\( \hbeta_H+\bv \) with \( \bLambda\T\bv=\bzero \), and then
\[
\norm{\X(\hbeta-\bb)}^2=\norm{\X(\hbeta-\hbeta_H)}^2-2(\hbeta-\hbeta_H)\T\X\T\X\bv+\norm{\X\bv}^2 .
\]
The cross term vanishes:
\( (\hbeta-\hbeta_H)\T\X\T\X\bv=(\bLambda\T\hbeta-\mathbf{d})\T(\bLambda\T\mathbf{C}\bLambda)^{-1}\bLambda\T\mathbf{C}\X\T\X\bv
=(\bLambda\T\hbeta-\mathbf{d})\T(\bLambda\T\mathbf{C}\bLambda)^{-1}\bLambda\T\bv=0 \). So the objective is
\( \text{SSE}+\norm{\X(\hbeta-\hbeta_H)}^2+\norm{\X\bv}^2 \), minimized exactly when
\( \X\bv=\bzero \), that is, \( \bv=\bzero \). Finally
\( \norm{\X(\hbeta-\hbeta_H)}^2=(\bLambda\T\hbeta-\mathbf{d})\T(\bLambda\T\mathbf{C}\bLambda)^{-1}\bLambda\T\mathbf{C}\X\T\X\mathbf{C}\bLambda(\bLambda\T\mathbf{C}\bLambda)^{-1}(\bLambda\T\hbeta-\mathbf{d}) \),
which simplifies to \( \text{SS}_H \).
:::

The same estimator comes out of Lagrange multipliers (@prp-mat-lagrange with \( \A=\X\T\X \)).
The proof above shows instead that the constrained fit is the feasible point nearest to
\( \hbeta \) in the geometry of \( \X\T\X \).

::: {#exm-ss-hypothesis-crime}
[Equal effects of poverty and single parenthood]

Both poverty and single parenthood are percentages of the population. Is one percentage point
of either associated with the same change in the murder rate? The hypothesis is
\( \blambda\T\bbeta=0 \) with \( \blambda=(0,1,-1,0)\T \). Formula @eq-ss-hypothesis gives
\( \text{SS}_H=1.807 \) on one degree of freedom, identical to the drop in residual sum of
squares when the two regressors are replaced by their sum. Its ratio to the residual mean
square is \( 0.815 \), so the data are entirely compatible with equal effects. For a hypothesis
with a nonzero right-hand side, \( \beta_{\text{single}}=0.5 \), @thm-ss-restricted gives the
constrained coefficients \( 0.1803 \) for poverty and \( -0.0044 \) for urbanization,
and a sum of squares of \( 3.357 \), both confirmed by a direct fit of the reduced model.
:::

```{.python .run #cell-extra-ss-hypothesis}
def hypothesis_ss(X, y, Lam, d=None):
    """(L'b - d)' [L' (X'X)^- L]^{-1} (L'b - d) for estimable L'beta = d."""
    G = np.linalg.pinv(X.T @ X)
    b = G @ X.T @ y
    u = Lam.T @ b - (0 if d is None else d)
    return u @ np.linalg.solve(Lam.T @ G @ Lam, u)

def sse(Z):
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    e = y - Z @ b
    return e @ e

Lam_equal = np.array([[0.0, 1.0, -1.0, 0.0]]).T             # beta_poverty = beta_single
ss_equal = hypothesis_ss(Xf, y, Lam_equal)
X_equal = np.column_stack([one, cols["poverty"] + cols["single"], cols["urban"]])
print(f"H: equal slopes   formula {ss_equal:.3f}   SSE_0 - SSE {sse(X_equal) - sse(Xf):.3f}")
```

## Constraints that are not estimable

Without estimability, the reduced model \( \mathcal S_0=\{\X\bb:\bLambda\T\bb=\bzero\} \) is still defined,
but it may be larger than the number of constraints suggests.

::: {#prp-ss-effective-df}
[How much a constraint restricts]

For any \( p\times q \) matrix \( \bLambda \), the subspace \( \mathcal S_0=\{\X\bb:\bLambda\T\bb=\bzero\} \)
has dimension \( r-\dim\bigl(\C(\bLambda)\cap\C(\X\T)\bigr) \). So the extra sum of squares for the
constraint has \( \dim\bigl(\C(\bLambda)\cap\C(\X\T)\bigr) \) degrees of freedom. It is identically
zero iff no nonzero linear combination of the constrained functions is estimable.
:::

::: {.proof}
\( \mathcal S_0 \) is the image of \( \Null(\bLambda\T) \) under \( \bb\mapsto\X\bb \), and the kernel of
this map restricted to \( \Null(\bLambda\T) \) is \( \Null(\bLambda\T)\cap\Null(\X) \), the null space of
the stacked matrix \( [\bLambda,\X\T]\T \). By rank–nullity (@thm-mat-rank-nullity),
\[
\dim\mathcal S_0=\dim\Null(\bLambda\T)-\dim\bigl(\Null(\bLambda\T)\cap\Null(\X)\bigr)
=(p-\rank\bLambda)-\bigl(p-\rank[\bLambda,\X\T]\bigr).
\]
The dimension formula for a sum of subspaces ([Section 6.2](../ch06-projections/02-subspaces.html))
gives \( \rank[\bLambda,\X\T]=\dim\bigl(\C(\bLambda)+\C(\X\T)\bigr)=\rank\bLambda+r-\dim\bigl(\C(\bLambda)\cap\C(\X\T)\bigr) \).
Substituting, \( \dim\mathcal S_0=r-\dim\bigl(\C(\bLambda)\cap\C(\X\T)\bigr) \). The degrees of freedom
of the extra sum of squares are \( r-\dim\mathcal S_0 \). A combination \( \bLambda\mathbf{a} \) is estimable
iff it lies in \( \C(\X\T) \) (@thm-proj-invariant-functions), which gives the last statement.
:::

In the one-way layout \( \E(Y_{ij})=\mu+\alpha_i \), the constraint \( \alpha_1=\alpha_2 \) has one degree of
freedom. The constraint \( \alpha_1=0 \) has none: it does not restrict the mean vector, and only
picks one of the least squares solutions (@thm-est-side-conditions). Software asked to “test” a
nonestimable hypothesis either refuses or silently tests its estimable part. [Chapter 17](../ch17-unbalanced-data/index.html) meets this
with empty cells (@thm-ub-empty-cells).

## Exercises

### A. Check your understanding

::: {#exr-ss-simple-hypothesis}
[A1]

In simple regression, compute \( \text{SS}_H \) for \( \beta_1=0 \) from @eq-ss-hypothesis and show
that it equals \( \hat{\beta}_1^2S_{xx} \), in agreement with @exr-ss-simple-anova.
:::

### B. Practice

::: {#exr-ss-not-monotone}
[B1]

Prove the chain rule \( \text{SS}(\X_2,\X_3\mid\X_1)=\text{SS}(\X_2\mid\X_1)+\text{SS}(\X_3\mid\X_1,\X_2) \),
and deduce \( \text{SS}(\X_3\mid\X_1,\X_2)\le\text{SS}(\X_2,\X_3\mid\X_1) \). Show by an example from
@exm-ss-adding-crime or its listing that \( \text{SS}(\X_3\mid\X_1,\X_2)\le\text{SS}(\X_3\mid\X_1) \)
can fail.
:::

::: {.solution}
Let \( \M_1 \), \( \M_{12} \) and \( \M_{123} \) be the projections onto the nested column spaces. Then
\( \M_{123}-\M_1=(\M_{12}-\M_1)+(\M_{123}-\M_{12}) \), and the chain rule follows by taking
\( \y\T(\cdot)\y \). Both terms are nonnegative, which gives the inequality. For the failure, the
listing of [Section 9.3](03-sequential-partial.html) reports
\( \text{SS}(\x_{\text{urban}}\mid\bone)=6.30 \) but
\( \text{SS}(\x_{\text{urban}}\mid\bone,\x_{\text{poverty}})=29.26 \). Adjusting for poverty
more than quadruples the sum of squares for urbanization.
:::

::: {#exr-ss-restricted-moments}
[B2]

In @thm-ss-restricted with \( \Cov(\Y)=\sigma^2\I \), show that
\[
\Cov(\hbeta_H)=\sigma^2\bigl[\mathbf{C}-\mathbf{C}\bLambda(\bLambda\T\mathbf{C}\bLambda)^{-1}\bLambda\T\mathbf{C}\bigr],
\]
that \( \Cov(\hbeta)-\Cov(\hbeta_H) \) is nonnegative definite, and that \( \E(\hbeta_H)=\bbeta \) when
the constraint is true. What is \( \E(\hbeta_H) \) when it is false?
:::

::: {.solution}
Write \( \mathbf{K}=\mathbf{C}\bLambda(\bLambda\T\mathbf{C}\bLambda)^{-1} \), so \( \hbeta_H=(\I-\mathbf{K}\bLambda\T)\hbeta+\mathbf{K}\mathbf{d} \).
By @thm-rv-linear and \( \Cov(\hbeta)=\sigma^2\mathbf{C} \),
\( \Cov(\hbeta_H)=\sigma^2(\I-\mathbf{K}\bLambda\T)\mathbf{C}(\I-\bLambda\mathbf{K}\T) \). Expanding and using
\( \bLambda\T\mathbf{C}\bLambda(\bLambda\T\mathbf{C}\bLambda)^{-1}=\I \), the two cross terms and the last term
all equal \( \mathbf{K}\bLambda\T\mathbf{C}=\mathbf{C}\bLambda(\bLambda\T\mathbf{C}\bLambda)^{-1}\bLambda\T\mathbf{C} \), with
signs \( -,-,+ \), which gives the formula. The difference
\( \sigma^2\mathbf{C}\bLambda(\bLambda\T\mathbf{C}\bLambda)^{-1}\bLambda\T\mathbf{C} \) is of the form
\( \mathbf{B}\mathbf{A}\mathbf{B}\T \) with \( \mathbf{A} \) positive definite, hence nonnegative definite. Finally
\( \E(\hbeta_H)=\bbeta-\mathbf{K}(\bLambda\T\bbeta-\mathbf{d}) \), which is \( \bbeta \) iff
\( \bLambda\T\bbeta=\mathbf{d} \). A false constraint buys the smaller variance at the price of a bias.
:::

::: {#exr-ss-oneway-contrast}
[B3]

In the one-way layout with model matrix \( \X=[\bone,\bz_1,\dots,\bz_g] \) and group sizes
\( n_1,\dots,n_g \), use the generalized inverse \( \G=\diag(0,1/n_1,\dots,1/n_g) \) in
@eq-ss-hypothesis to show that the sum of squares for \( \alpha_1=\alpha_2 \) is
\( (\bar{y}_1-\bar{y}_2)^2/(1/n_1+1/n_2) \). Verify with @prp-ss-effective-df that the constraint
\( \alpha_1=0 \) has zero degrees of freedom.
:::

::: {.solution}
\( \G \) is a generalized inverse because \( \X\T\X \) has the block form with first row
\( (n,n_1,\dots,n_g) \) and lower right block \( \diag(n_1,\dots,n_g) \), and a direct multiplication
gives \( \X\T\X\G\X\T\X=\X\T\X \). With \( \blambda=(0,1,-1,0,\dots,0)\T \),
\( \G\X\T\y=(0,\bar{y}_1,\dots,\bar{y}_g)\T \), so \( \blambda\T\hbeta=\bar{y}_1-\bar{y}_2 \) and
\( \blambda\T\G\blambda=1/n_1+1/n_2 \). For \( \alpha_1=0 \), \( \blambda=(0,1,0,\dots,0)\T \) is not in
\( \C(\X\T) \) (@exm-proj-ginverse-numeric), and \( \C(\blambda) \) is a line, so
\( \C(\blambda)\cap\C(\X\T)=\{\bzero\} \).
:::

### C. Going deeper

::: {#exr-ss-extra-max}
[C1]

Using @exr-ss-max-single and @lem-ss-residualized, show that
\[
\text{SS}(\X_2\mid\X_1)=\max_{\mathbf{a}}\frac{(\mathbf{a}\T\tilde{\X}_2\T\y)^2}{\mathbf{a}\T\tilde{\X}_2\T\tilde{\X}_2\mathbf{a}},
\]
the maximum being over \( \mathbf{a} \) with \( \tilde{\X}_2\mathbf{a}\neq\bzero \). Deduce that the extra sum of
squares of a block is at least that of any single column of the block, given \( \X_1 \), and that
the sum of squares of an estimable hypothesis \( \bLambda\T\bbeta=\bzero \) equals
\( \max_{\mathbf{a}}(\mathbf{a}\T\bLambda\T\hbeta)^2/(\mathbf{a}\T\bLambda\T\G\bLambda\mathbf{a}) \).
:::

::: {#exr-ss-no-variance-change}
[C2]

@exr-mat-adding-regressors gives the matrix inequality behind the remark after @thm-ss-adding.
Using part (d) of that theorem instead, show that
\( \Var(\hat{\beta}_{W,j})=\Var(\hat{\beta}_j) \) iff the \( j \)th row of \( \bL \) is zero, and
that this holds for every \( j \) iff \( \X\T\mathbf{Z}=\bzero \). Show that when
\( \X\T\mathbf{Z}=\bzero \) the estimates \( \hbeta \) do not change at all.
:::

::: {#exr-ss-two-way-constraints}
[C3]

In the additive two-way model \( \E(Y_{ijk})=\mu+\alpha_i+\beta_j \) with all cells occupied, use
@prp-ss-effective-df to find the degrees of freedom of the extra sums of squares for the
constraints (i) \( \sum_i\alpha_i=0 \) and \( \sum_j\beta_j=0 \); (ii) \( \alpha_1=\alpha_2 \);
(iii) \( \alpha_1=0 \) and \( \alpha_1=\alpha_2 \).
:::
