# Least squares and the normal equations

The straight-line calculation of [Section 5.2](02-simple-regression.html) extends to any
model matrix. The algebra is shorter in matrix form than it was for two coefficients,
because the derivatives of linear and quadratic forms were worked out once and for
all in [Chapter 1](../ch01-matrix-algebra/index.html). In this section \( \y \) is a fixed
data vector. Nothing here uses probability. The statistical properties come in
[Section 5.5](05-properties.html).

## The least squares estimate

::: {#def-lm-least-squares}
[Least squares estimate]

For an \( n\times p \) matrix \( \X \) and \( \y\in\Real^n \), the **least squares criterion**
is
\[
S(\bb)=\norm{\y-\X\bb}^2=\sum_{i=1}^n\bigl(y_i-\x_{(i)}\T\bb\bigr)^2,\qquad\bb\in\Real^p .
\]
A **least squares estimate** of \( \bbeta \) is any \( \hbeta \) that minimizes \( S \). The
corresponding **fitted values** are \( \hY=\X\hbeta \), the **residuals** are
\( \he=\y-\X\hbeta \), and the **residual sum of squares** is
\( \text{SSE}=S(\hbeta)=\norm{\he}^2 \).
:::

Following the book's notation, we write \( \hY \) and \( \he \) in bold for the vectors,
and \( \hat{y}_i \) and \( \hat{\varepsilon}_i \) for their entries, whether the data are regarded as
fixed numbers or as random.

::: {#thm-lm-ls-full-rank}
[Least squares for a full-rank model]

Let \( \X \) be an \( n\times p \) matrix of rank \( p \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \X\T\X \) is positive definite, and in particular nonsingular;

2. the **normal equations**
   \[
   \X\T\X\bb=\X\T\y
   \]{#eq-lm-normal-equations}

   have exactly one solution,
   \[
   \hbeta=(\X\T\X)^{-1}\X\T\y ;
   \]{#eq-lm-beta-hat}

3. for every \( \bb\in\Real^p \),
   \[
   \norm{\y-\X\bb}^2=\norm{\y-\X\hbeta}^2+(\bb-\hbeta)\T\X\T\X(\bb-\hbeta),
   \]{#eq-lm-ls-decomposition}

   so \( \hbeta \) is the unique least squares estimate.
:::

If instead \( \rank(\X)<p \), least squares estimates exist but are never unique.
:::

::: {.proof}
*(a)* \( \X\T\X \) is symmetric, and \( \bb\T\X\T\X\bb=\norm{\X\bb}^2\ge0 \), with equality iff
\( \X\bb=\bzero \). Since \( \X \) has rank \( p \), rank–nullity (@thm-mat-rank-nullity) gives
\( \Null(\X)=\{\bzero\} \), so the quadratic form is positive for \( \bb\ne\bzero \). A positive
definite matrix is nonsingular (@prp-mat-pd-properties(e)).

*(b)* This is the computation of @exm-mat-ls-gradient: by @thm-mat-quadform-derivative the
gradient of \( S \) is \( -2\X\T(\y-\X\bb) \) and its Hessian is \( 2\X\T\X \).
Every minimizer of the differentiable function \( S \) is a stationary point, and the
stationary points are the solutions of @eq-lm-normal-equations. By (a) there is exactly one, @eq-lm-beta-hat.

*(c)* Put \( \bm d=\bb-\hbeta \). Then \( \y-\X\bb=\he-\X\bm d \), and
\[
\norm{\y-\X\bb}^2=\norm{\he}^2-2\bm d\T\X\T\he+\bm d\T\X\T\X\bm d .
\]
The middle term vanishes because \( \X\T\he=\X\T\y-\X\T\X\hbeta=\bzero \) by the normal
equations. The last term is positive for \( \bm d\ne\bzero \) by (a). So \( S(\bb)>S(\hbeta) \)
whenever \( \bb\ne\hbeta \).

For the last statement, suppose \( \rank(\X)<p \). The vector \( \X\T\y \) lies in
\( \C(\X\T)=\C(\X\T\X) \) (@prp-mat-rank-product(c)), so by @prp-mat-quadratic-min, applied
with \( \A=\X\T\X \) and \( S(\bb)=\y\T\y+(\bb\T\X\T\X\bb-2\bb\T\X\T\y) \), the minimizers of \( S \)
are exactly the solutions of the normal equations, and they exist. But there is a
\( \bv\ne\bzero \) with \( \X\bv=\bzero \), and then \( S(\bb+t\bv)=S(\bb) \) for every \( t \), so no
minimizer is unique.
:::

The calculus proof finds the stationary point and then needs a separate argument,
@eq-lm-ls-decomposition, to show it is the global minimum. The Hessian \( 2\X\T\X \)
tells the same story: it is positive definite, so \( S \) is a strictly convex quadratic
bowl with a single lowest point. The name "normal equations" has nothing to do
with the normal distribution. It refers to the orthogonality they express:
\( \X\T(\y-\X\hbeta)=\bzero \) says the residual vector is *normal*, that is, perpendicular,
to every column of \( \X \). [Chapter 6](../ch06-projections/index.html) makes that
geometry the starting point, derives the normal equations again without calculus
(@thm-proj-normal-equations), and treats the rank-deficient case in full.

::: {.remark}
[The rank-deficient case]

When \( \rank(\X)<p \), the proof shows that the least squares estimates form an affine
set \( \hbeta+\Null(\X) \) of dimension \( p-\rank(\X) \). Different programs report different
members of this set, according to which columns they drop or which constraints they
impose. The fitted values \( \X\hbeta \) are the same for all of them, since
\( \X(\hbeta+\bv)=\X\hbeta \) for \( \bv\in\Null(\X) \). Which linear combinations of the coefficients
share this invariance is the question of estimability, answered in
[Chapter 8](../ch08-estimability/index.html) (@def-est-estimable and @thm-est-characterization).
:::

Least squares under known linear restrictions \( \bm H\bb=\bm h \), such as shares that add to
one, has a closed form found by Lagrange multipliers in @exr-mat-constrained-ls; the increase
in the residual sum of squares it computes drives the test of the restrictions in Chapter 11.

::: {#exm-lm-simple-matrix}
[The straight line again]

For \( \X=[\bone,\x] \),
\[
\X\T\X=\begin{pmatrix}n&n\bar{x}\\n\bar{x}&\sum_ix_i^2\end{pmatrix},\qquad
\X\T\y=\begin{pmatrix}n\bar{y}\\\sum_ix_iy_i\end{pmatrix}.
\]
The determinant is \( n\sum_ix_i^2-n^2\bar{x}^2=nS_{xx} \), positive iff the \( x_i \) are not all
equal, which is exactly when \( \rank(\X)=2 \). Inverting,
\[
(\X\T\X)^{-1}=\frac{1}{nS_{xx}}\begin{pmatrix}\sum_ix_i^2&-n\bar{x}\\-n\bar{x}&n\end{pmatrix}
=\begin{pmatrix}\frac1n+\frac{\bar{x}^2}{S_{xx}}&-\frac{\bar{x}}{S_{xx}}\\[2pt]-\frac{\bar{x}}{S_{xx}}&\frac{1}{S_{xx}}\end{pmatrix},
\]{#eq-lm-simple-inverse}

using \( \sum_ix_i^2=S_{xx}+n\bar{x}^2 \). The second entry of \( (\X\T\X)^{-1}\X\T\y \) is
\( (\sum_ix_iy_i-n\bar{x}\bar{y})/S_{xx}=S_{xy}/S_{xx} \), and the first is
\( \bar{y}-\bar{x}S_{xy}/S_{xx} \). This recovers @thm-lm-simple-ls, and
@eq-lm-ls-decomposition reduces to @eq-lm-simple-decomposition.
:::

## Fitted values, residuals and the hat matrix

Substituting @eq-lm-beta-hat, the fitted values are
\[
\hY=\X\hbeta=\bH\y,\qquad \bH=\X(\X\T\X)^{-1}\X\T ,
\]{#eq-lm-hat-matrix}

and the residuals are \( \he=(\I-\bH)\y \). The \( n\times n \) matrix \( \bH \) is the
**hat matrix**, so called because it puts the hat on \( \y \). It depends only on the
regressors, never on the response.

::: {#prp-lm-fit-algebra}
[Algebra of the fit]

Let \( \X \) have rank \( p \), and let \( \bH \), \( \hY \) and \( \he \) be as above.

::: {.enumerate options="label=(\alph*)"}
1. \( \bH \) is symmetric and idempotent, \( \bH\X=\X \) and \( (\I-\bH)\X=\bzero \), and
   \( \tr\bH=p \), \( \tr(\I-\bH)=n-p \).

2. \( \X\T\he=\bzero \) and \( \hY\T\he=0 \).

3. \( \norm{\y}^2=\norm{\hY}^2+\norm{\he}^2 \), and
   \( \text{SSE}=\y\T(\I-\bH)\y=\y\T\y-\hbeta\T\X\T\y \).

4. If \( \bone=\X\bm c \) for some \( \bm c \) (for example, if \( \X \) has an intercept column), then
   \( \sum_i\hat{\varepsilon}_i=0 \), the mean of the fitted values is \( \bar{y} \), and
   \[
   \sum_i(y_i-\bar{y})^2=\sum_i(\hat{y}_i-\bar{y})^2+\sum_i\hat{\varepsilon}_i^2,
   \]
   that is, \( \text{SST}=\text{SSR}+\text{SSE} \).

5. If \( \y=\X\bbeta+\be \) for some \( \bbeta \) and \( \be \), then
   \( \hbeta=\bbeta+(\X\T\X)^{-1}\X\T\be \) and \( \he=(\I-\bH)\be \).
:::

:::

::: {.proof}
*(a)* \( \bH\T=\X\bigl((\X\T\X)^{-1}\bigr)\T\X\T=\bH \), since the inverse of a symmetric matrix is
symmetric. \( \bH\X=\X(\X\T\X)^{-1}\X\T\X=\X \), and then
\( \bH^2=\bH\X(\X\T\X)^{-1}\X\T=\X(\X\T\X)^{-1}\X\T=\bH \). By the cyclic property
(@thm-mat-trace-cyclic), \( \tr\bH=\tr\bigl((\X\T\X)^{-1}\X\T\X\bigr)=\tr\I_p=p \).
*(b)* \( \X\T\he=\X\T(\I-\bH)\y=\bigl((\I-\bH)\X\bigr)\T\y=\bzero \), and
\( \hY\T\he=\hbeta\T\X\T\he=0 \). *(c)* \( \y=\hY+\he \) with \( \hY\T\he=0 \) gives the first
identity. For the second, \( (\I-\bH) \) is symmetric and idempotent, so
\( \norm{(\I-\bH)\y}^2=\y\T(\I-\bH)\y \); and
\( \y\T\bH\y=\y\T\X\hbeta=\hbeta\T\X\T\y \). *(d)* \( \bone\T\he=\bm c\T\X\T\he=0 \) by (b). Then
\( \sum_i\hat{y}_i=\sum_iy_i \). Write \( y_i-\bar{y}=(\hat{y}_i-\bar{y})+\hat{\varepsilon}_i \). The cross
term \( \sum_i(\hat{y}_i-\bar{y})\hat{\varepsilon}_i=\hY\T\he-\bar{y}\,\bone\T\he=0 \). *(e)* Substitute
\( \y=\X\bbeta+\be \) into @eq-lm-beta-hat and into \( (\I-\bH)\y \), using (a).
:::

Part (e) will carry most of the statistical theory. The estimation error
\( \hbeta-\bbeta \) is a fixed linear function of the unobserved errors, and the residual
vector does not depend on \( \bbeta \) at all. Part (b) says that the residuals satisfy
\( p \) independent linear constraints. That is the origin of the divisor \( n-p \) in
[Section 5.6](06-error-variance.html).

For a model with an intercept, part (d) motivates the **coefficient of determination**
\[
R^2=\frac{\text{SSR}}{\text{SST}}=1-\frac{\text{SSE}}{\text{SST}},
\]
the proportion of the variation of the response about its mean that the fitted values
reproduce. It generalizes @prp-lm-simple-fit(d), and lies in \( [0,1] \). It can only
increase when a column is added to \( \X \) (@exr-lm-add-column), which is why [Chapter 9](../ch09-sums-of-squares/index.html)
introduces an adjusted version (@prp-ss-adjusted-r2). [Section 6.7](../ch06-projections/07-reparameterization.html)
interprets \( R \) as the cosine of an angle.

## Centring the regressors

For a model with an intercept, write \( \X=[\bone,\X_1] \), where \( \X_1 \) is \( n\times k \) with
column means \( \bar{\x}=\X_1\T\bone/n \). The **centred regressors** are the columns of
\[
\X_c=\X_1-\bone\bar{\x}\T=\Bigl(\I-\tfrac1n\bone\bone\T\Bigr)\X_1 ,
\]
each of which sums to zero. The next result is the matrix version of the formulas
\( \hat{\beta}_1=S_{xy}/S_{xx} \) and \( \hat{\beta}_0=\bar{y}-\hat{\beta}_1\bar{x} \).

::: {#prp-lm-centred}
[Centred form]

Let \( \X=[\bone,\X_1] \) have rank \( p=k+1 \), and partition
\( \hbeta=(\hat{\beta}_0,\hbeta_1\T)\T \). Then \( \X_c \) has rank \( k \), and
\[
\hbeta_1=(\X_c\T\X_c)^{-1}\X_c\T\y=(\X_c\T\X_c)^{-1}\X_c\T(\y-\bar{y}\bone),
\qquad \hat{\beta}_0=\bar{y}-\bar{\x}\T\hbeta_1 .
\]
:::

::: {.proof}
For \( \bb=(b_0,\bb_1\T)\T \),
\[
\X\bb=b_0\bone+\X_1\bb_1=(b_0+\bar{\x}\T\bb_1)\bone+\X_c\bb_1=[\bone,\X_c]\,\bT\bb,\qquad
\bT=\begin{pmatrix}1&\bar{\x}\T\\\bzero&\I_k\end{pmatrix}.
\]
\( \bT \) is nonsingular (it is upper triangular with unit diagonal), so
\( [\bone,\X_c]=\X\bT^{-1} \) has rank \( p \) (@prp-mat-rank-product(b)), and hence
\( \X_c \) has rank \( k \). The criterion satisfies \( \norm{\y-\X\bb}^2=\norm{\y-[\bone,\X_c]\bm a}^2 \)
with \( \bm a=\bT\bb \), and \( \bb\mapsto\bT\bb \) is a bijection of \( \Real^p \). So \( \hbeta \)
minimizes the left side iff \( \hat{\bm a}=\bT\hbeta \) minimizes the right side. Because
\( \bone\T\X_c=\bzero\T \), the normal equations for \( [\bone,\X_c] \) are block diagonal,
\[
\begin{pmatrix}n&\bzero\T\\\bzero&\X_c\T\X_c\end{pmatrix}\hat{\bm a}
=\begin{pmatrix}n\bar{y}\\\X_c\T\y\end{pmatrix},
\]
giving \( \hat{a}_0=\bar{y} \) and \( \hat{\bm a}_1=(\X_c\T\X_c)^{-1}\X_c\T\y \). Undoing \( \bT \):
\( \hbeta_1=\hat{\bm a}_1 \) and \( \hat{\beta}_0=\hat{a}_0-\bar{\x}\T\hat{\bm a}_1 \). Finally
\( \X_c\T\bone=\bzero \), so \( \X_c\T\y=\X_c\T(\y-\bar{y}\bone) \).
:::

In words: to find the slopes, centre everything and fit without an intercept; the
intercept then makes the fitted surface pass through the point of means
\( (\bar{\x},\bar{y}) \). The matrix \( \X_c\T\X_c/(n-1) \) is the sample covariance matrix of
the regressors, and \( \X_c\T\y/(n-1) \) their sample covariances with the response. So the
slopes are "covariance divided by variance", as in the straight-line case, with a matrix
inverse in place of the division. This is the simplest case of the Frisch–Waugh–Lovell
theorem of [Section 6.6](../ch06-projections/06-fwl.html), with the intercept playing the
part of the regressors that are "partialled out".

## Indicator variables

The least squares estimates in a model with only indicator columns are simple averages.

::: {#prp-lm-indicator-means}
[Coefficients of indicator variables]

Suppose each case belongs to one of \( c\ge2 \) groups, group \( k \) has \( n_k\ge1 \) cases with
mean response \( \bar{y}_k \), and \( \bz_k \) is the indicator column of group \( k \). Let
\( \X=[\bone,\bz_1,\dots,\bz_{c-1}] \), so that group \( c \) is the reference. Then \( \X \) has
full rank, and
\[
\hat{\beta}_0=\bar{y}_c,\qquad \hat{\beta}_k=\bar{y}_k-\bar{y}_c\quad(k=1,\dots,c-1).
\]
The fitted value of every case is the mean of its group.
:::

::: {.proof}
If \( \X\bb=\bzero \), the entries for a case in group \( c \) give \( b_0=0 \), and those
for a case in group \( k<c \) give \( b_0+b_k=0 \), so \( \bb=\bzero \): the rank is \( p=c \). Let
\( \bb \) be the proposed vector. Its fitted value for a case in group \( k<c \) is
\( \bar{y}_c+(\bar{y}_k-\bar{y}_c)=\bar{y}_k \), and for a case in group \( c \) it is \( \bar{y}_c \).
The residuals \( y_i-\bar{y}_{k(i)} \) therefore sum to zero within each group, so
\( \bz_k\T(\y-\X\bb)=0 \) for every \( k=1,\dots,c \). Every column of \( \X \) is a sum of
group indicators, since \( \bone=\bz_1+\dots+\bz_c \), so \( \X\T(\y-\X\bb)=\bzero \): the
normal equations hold. By @thm-lm-ls-full-rank, \( \bb=\hbeta \).
:::

```{.python .run #cell-rand-dummies-dummies}
import numpy as np
import statsmodels.api as sm
rand = sm.datasets.randhie.load_pandas().data
y = rand["mdvis"].to_numpy(dtype=float)
D = rand[["hlthg", "hlthf", "hlthp"]].to_numpy(dtype=float)   # good, fair, poor
X = np.column_stack([np.ones(len(y)), D])                      # excellent = reference

beta_hat = np.linalg.solve(X.T @ X, X.T @ y)
level = D @ np.array([1, 2, 3])                                # 0 = excellent, ..., 3 = poor
means = np.array([y[level == k].mean() for k in range(4)])
print("coefficients:            ", np.round(beta_hat, 4))
print("group means:             ", np.round(means, 4))
print("means minus reference:   ", np.round(means - means[0], 4))
```

For the health data of @exm-lm-health-coding the group sizes are
11,019 (excellent), 7,309 (good),
1,560 (fair) and \( 302 \) (poor), and the group means
are \( 2.634 \), \( 2.902 \), \( 3.692 \) and
\( 5.795 \) visits. The listing confirms that the coefficients are these
means' differences from the reference.

## A worked example

::: {#exm-lm-money}
[Money demand in Denmark]

Economic theory relates the amount of money people choose to hold to their income,
which drives the volume of transactions, and to interest rates, which measure the return
given up by holding money rather than bonds. A classic data set records, for
\( n=55 \) quarters from 1974 to 1987 in Denmark, the logarithm of real money
holdings (\( y \)), the logarithm of real income, the interest rate on bonds and the
interest rate on deposits (rates as fractions per year). The model has \( p=4 \)
columns: an intercept and the three regressors. Solving the normal equations gives
\[
\hbeta=(4.3945,\ 1.2958,\ -2.6163,\ 0.6186)\T,
\]
with residual sum of squares \( \text{SSE}=0.09253 \) and
\( R^2=0.9262 \). The listing checks that \( \X\T\he \) is zero up to rounding
error, that the two formulas for SSE in @prp-lm-fit-algebra(c) agree, and that the
slopes computed from centred data (@prp-lm-centred) are the same.
[Section 5.6](06-error-variance.html) attaches standard errors to these numbers, and
[Section 5.7](07-interpreting-coefficients.html) interprets them.
:::

```{.python .run #cell-money-demand-data}
import numpy as np
import statsmodels.api as sm
dk = sm.datasets.danish_data.load_pandas().data
y = dk["lrm"].to_numpy()
X = np.column_stack([np.ones(len(y)), dk["lry"], dk["ibo"], dk["ide"]])
n, p = X.shape
print("n =", n, " p =", p)
```

```{.python .run #cell-money-demand-normal}
XtX = X.T @ X
Xty = X.T @ y
beta_hat = np.linalg.solve(XtX, Xty)       # solve X^T X b = X^T y
fitted = X @ beta_hat
resid = y - fitted
SSE = resid @ resid

print("beta_hat =", np.round(beta_hat, 4))
print("X^T e    =", X.T @ resid)           # zero up to rounding
print("SSE two ways:", SSE, y @ y - beta_hat @ Xty)
```

```{.python .run #cell-money-demand-centred}
Z = X[:, 1:] - X[:, 1:].mean(axis=0)       # centred regressors, no intercept column
slopes = np.linalg.solve(Z.T @ Z, Z.T @ (y - y.mean()))
intercept = y.mean() - X[:, 1:].mean(axis=0) @ slopes
print("slopes from centred data:", np.round(slopes, 4))
print("intercept recovered:     ", round(intercept, 4))
```

::: {.warning}
[Formulas are not algorithms]

The formula \( \hbeta=(\X\T\X)^{-1}\X\T\y \) is for reasoning, not for computing. Forming
\( \X\T\X \) squares the condition number of the problem, and inverting it wastes work
and accuracy. The listings solve the normal equations with a linear solver, which is
acceptable for small, well-conditioned problems like these. Good software computes
\( \hbeta \) from an orthogonal factorization of \( \X \) itself.
[Section 6.10](../ch06-projections/10-computation.html) shows the loss of accuracy in an
example, and [Chapter 10](../ch10-computation/index.html) develops the algorithms
(@thm-cmp-cholesky, @thm-cmp-householder).
:::

## Exercises

### A. Check your understanding

::: {#exr-lm-orthogonal-columns}
[A1]

Suppose the columns of \( \X \) are nonzero and mutually orthogonal. Show that
\( \hat{\beta}_j=\x_j\T\y/\x_j\T\x_j \) for each \( j \), so that each coefficient is the same
as in the regression of \( \y \) on \( \x_j \) alone. Show also that deleting a column does not
change the other estimates.
:::

### B. Practice

::: {#exr-lm-add-column}
[B1]

Let \( \X=[\X_1,\x] \) have full rank, and let \( \text{SSE}_1 \) and \( \text{SSE} \) be the residual
sums of squares from regressing \( \y \) on \( \X_1 \) and on \( \X \). Show that
\( \text{SSE}\le\text{SSE}_1 \), directly from @def-lm-least-squares. Deduce that \( R^2 \)
cannot decrease when a column is added to a model with an intercept.
:::

::: {.solution}
Let \( \hbeta_1 \) be the least squares estimate for \( \X_1 \). The vector \( (\hbeta_1\T,0)\T \) is
one candidate in the minimization over \( \bb\in\Real^p \), and it gives the value
\( \norm{\y-\X_1\hbeta_1}^2=\text{SSE}_1 \). The minimum over all candidates is at most this.
With an intercept in \( \X_1 \), \( R^2=1-\text{SSE}/\text{SST} \) and SST does not depend on the
model, so \( R^2 \) cannot decrease. Equality holds iff the extra coefficient in the larger
fit is zero, because by @eq-lm-ls-decomposition the candidate \( (\hbeta_1\T,0)\T \) is optimal only if it
equals \( \hbeta \).
:::

::: {#exr-lm-weighted}
[B2]

Let \( \W=\diag(w_1,\dots,w_n) \) with all \( w_i>0 \), and \( \X \) of full rank. Show that
\( \sum_iw_i(y_i-\x_{(i)}\T\bb)^2 \) has the unique minimizer
\( (\X\T\W\X)^{-1}\X\T\W\y \), by applying @thm-lm-ls-full-rank to the rescaled data
\( \W^{1/2}\y \) and \( \W^{1/2}\X \). Which cases get most influence when the weights are the
reciprocals of the error variances?
:::

::: {.solution}
\( \sum_iw_i(y_i-\x_{(i)}\T\bb)^2=\norm{\W^{1/2}\y-\W^{1/2}\X\bb}^2 \). The matrix \( \W^{1/2}\X \)
has rank \( p \), because \( \W^{1/2} \) is nonsingular (@prp-mat-rank-product(b)). By
@thm-lm-ls-full-rank the unique minimizer is
\( \bigl((\W^{1/2}\X)\T\W^{1/2}\X\bigr)^{-1}(\W^{1/2}\X)\T\W^{1/2}\y=(\X\T\W\X)^{-1}\X\T\W\y \). With
\( w_i=1/\Var(Y_i) \), the least variable cases count most. This is weighted least squares,
the subject of Chapter 21 and, with a general covariance matrix, of
[Section 6.9](../ch06-projections/09-inner-products.html).
:::

::: {#exr-lm-linear-estimator}
[B3]

Show that \( \hbeta=\A\y \) with \( \A=(\X\T\X)^{-1}\X\T \), that \( \A\X=\I_p \), and that
\( \A\A\T=(\X\T\X)^{-1} \). Give an example of a different \( p\times n \) matrix \( \B \) with
\( \B\X=\I_p \) when \( n>p \). ([Chapter 7](../ch07-optimality/index.html) compares all such estimators \( \B\y \).)
:::

### C. Going deeper

::: {#exr-lm-add-case}
[C1]

Let \( \X \) have full rank and let a new case \( (\x_0,y_0) \) be appended. Using
@thm-mat-woodbury, show that the new estimate is
\[
\hbeta_{\text{new}}=\hbeta+\frac{(\X\T\X)^{-1}\x_0\,(y_0-\x_0\T\hbeta)}{1+\x_0\T(\X\T\X)^{-1}\x_0}.
\]
Interpret: the estimate moves in proportion to the new case's prediction error. This is
the basis of *recursive least squares*; [Chapter 10](../ch10-computation/index.html) treats the numerically stable version.
:::
