# The Gauss–Markov theorem

[Section 7.1](01-linear-unbiased.html) turned the search for a best linear unbiased estimator into a
geometric problem. The unbiased coefficient vectors form the flat
\( \mathcal A_{\blambda}=\mathbf{a}_*+\C(\X)\perpc \), and the variance of each is \( \sigma^2 \) times its squared
length. The shortest vector in a flat parallel to \( \C(\X)\perpc \) is the one in \( \C(\X) \), and
that vector is the least squares coefficient vector \( \mathbf{a}_* \). The Gauss–Markov theorem says exactly
this, together with its matrix form and its consequences.

## The theorem

We keep the second-moment assumptions @eq-opt-model, with \( \X \) of any rank \( r \), and write
\( (\X\T\X)\ginv \) for any generalized inverse of \( \X\T\X \).

::: {#thm-opt-gauss-markov}
[Gauss–Markov]

Assume @eq-opt-model. Let \( \blambda\in\C(\X\T) \), say \( \blambda=\X\T\boldsymbol{\uprho} \), and let
\( \hbeta \) be any least squares estimate.

::: {.enumerate options="label=(\alph*)"}
1. \( \blambda\T\hbeta \) is a linear unbiased estimator of \( \blambda\T\bbeta \), with
   \[
\Var(\blambda\T\hbeta)=\sigma^2\blambda\T(\X\T\X)\ginv\blambda .
\]

2. For every linear unbiased estimator \( \mathbf{a}\T\Y \) of \( \blambda\T\bbeta \),
   \[
\Var(\mathbf{a}\T\Y)=\Var(\blambda\T\hbeta)+\sigma^2\norm{(\I-\M)\mathbf{a}}^2 .
\]{#eq-opt-gm-split}

3. Consequently \( \blambda\T\hbeta \) is a BLUE of \( \blambda\T\bbeta \), and it is the only one: an LUE
   \( \mathbf{a}\T\Y \) is a BLUE iff \( \mathbf{a}=\X(\X\T\X)\ginv\blambda \), that is, iff \( \mathbf{a}\T\y=\blambda\T\hbeta \)
   for every \( \y \).

4. *(Matrix form.)* Let \( \bLambda \) be \( p\times q \) with \( \C(\bLambda)\subseteq\C(\X\T) \). For every
   linear unbiased estimator \( \A\Y \) of \( \bLambda\T\bbeta \),
   \[
\Cov(\A\Y)=\Cov(\bLambda\T\hbeta)+\sigma^2\A(\I-\M)\A\T,
\]{#eq-opt-gm-matrix}

   where \( \Cov(\bLambda\T\hbeta)=\sigma^2\bLambda\T(\X\T\X)\ginv\bLambda \). So \( \Cov(\A\Y)-\Cov(\bLambda\T\hbeta) \) is nonnegative definite, and it is zero iff
   \( \A\Y=\bLambda\T\hbeta \) for every \( \y \).
:::
:::

::: {.proof}
Put \( \mathbf{a}_*=\M\boldsymbol{\uprho}=\X(\X\T\X)\ginv\blambda \). By @prp-opt-lue-set(c),
\( \mathbf{a}_*\T\Y=\blambda\T\hbeta \). Since \( \mathbf{a}_*\in\mathcal A_{\blambda} \), it is an LUE by @prp-opt-lue.
By @eq-opt-var-norm its variance is
\[
\sigma^2\norm{\M\boldsymbol{\uprho}}^2=\sigma^2\boldsymbol{\uprho}\T\M\boldsymbol{\uprho}
=\sigma^2\boldsymbol{\uprho}\T\X(\X\T\X)\ginv\X\T\boldsymbol{\uprho}=\sigma^2\blambda\T(\X\T\X)\ginv\blambda ,
\]
using @thm-proj-M-formula. This proves (a).

For (b), let \( \mathbf{a}\T\Y \) be any LUE. Then \( \mathbf{a}\in\mathcal A_{\blambda} \), so \( \M\mathbf{a}=\mathbf{a}_* \) by
@prp-opt-lue-set(b). The decomposition \( \mathbf{a}=\M\mathbf{a}+(\I-\M)\mathbf{a}=\mathbf{a}_*+(\I-\M)\mathbf{a} \) is orthogonal,
so \( \norm{\mathbf{a}}^2=\norm{\mathbf{a}_*}^2+\norm{(\I-\M)\mathbf{a}}^2 \). Multiplying by \( \sigma^2 \) and using
@eq-opt-var-norm twice gives @eq-opt-gm-split.

For (c), the second term of @eq-opt-gm-split is nonnegative, so \( \blambda\T\hbeta \) is a BLUE. Because
\( \sigma^2>0 \), equality holds iff \( (\I-\M)\mathbf{a}=\bzero \), that is, iff \( \mathbf{a}\in\C(\X) \). By
@prp-opt-lue-set(b) the only point of \( \mathcal A_{\blambda} \) in \( \C(\X) \) is \( \mathbf{a}_* \). Two linear
statistics \( \mathbf{a}\T\y \) and \( \mathbf{a}_*\T\y \) agree for every \( \y \) iff \( \mathbf{a}=\mathbf{a}_* \).

For (d), write \( \bLambda=\X\T\mathbf{P} \) for an \( n\times q \) matrix \( \mathbf{P} \), and put
\( \A_*=\mathbf{P}\T\M \), so that \( \A_*\Y=\bLambda\T\hbeta \) by @prp-opt-lue-set(c) applied column by
column. Each row \( \mathbf{a}_k\T \) of an LUE matrix \( \A \) has \( \mathbf{a}_k\in\mathcal A_{\blambda_k} \), where
\( \blambda_k \) is the \( k \)th column of \( \bLambda \). So \( \M\mathbf{a}_k=\M\mathbf{p}_k \), and hence \( \A\M=\A_* \)
(the \( k \)th row of \( \A\M \) is \( (\M\mathbf{a}_k)\T \), since \( \M \) is symmetric). Now split the identity:
\[
\begin{aligned}
\A\A\T&=\A\M\A\T+\A(\I-\M)\A\T\\
&=(\A\M)(\A\M)\T+\A(\I-\M)\A\T=\A_*\A_*\T+\A(\I-\M)\A\T ,
\end{aligned}
\]
using \( \M=\M\M\T \). Multiplying by \( \sigma^2 \) and using @thm-rv-linear gives
\( \Cov(\A\Y)=\Cov(\A_*\Y)+\sigma^2\A(\I-\M)\A\T \). Also
\( \A_*\A_*\T=\mathbf{P}\T\M\mathbf{P}=\bLambda\T(\X\T\X)\ginv\bLambda \), as in (a). The difference
\( \sigma^2\A(\I-\M)\A\T=\sigma^2[\A(\I-\M)][\A(\I-\M)]\T \) is nonnegative definite. It is zero iff
\( \A(\I-\M)=\mathbf{0} \), iff \( \A=\A\M=\A_* \).
:::

Three features of the theorem deserve comment.

**It is distribution-free.** The proof used \( \E(\be)=\bzero \) and \( \Cov(\be)=\sigma^2\I \), and nothing
else. The errors may be skewed, heavy-tailed, discrete, or dependent in any way that leaves them
uncorrelated. They need not even be identically distributed, as long as their variances are equal.

**It holds in every rank.** When \( \X \) is rank deficient, \( \hbeta \) is not unique and neither is the
generalized inverse. But for estimable \( \blambda \), both \( \blambda\T\hbeta \) and
\( \blambda\T(\X\T\X)\ginv\blambda \) are the same for every choice (@thm-proj-invariant-functions and
@lem-proj-ginverse-invariance). A nonestimable \( \blambda\T\bbeta \) has no LUE at all, so the theorem
covers every function for which the question makes sense.

**The proof is a projection.** Equation @eq-opt-gm-split is Pythagoras in \( \Real^n \), applied to the
coefficient vector and not to the data.

::: {.idea}
The least squares coefficient vector \( \mathbf{a}_* \) is the orthogonal projection onto \( \C(\X) \) of the
coefficient vector of *any* linear unbiased estimator. Projecting a coefficient vector onto
\( \C(\X) \) keeps the estimator unbiased and removes the part that only adds noise. That part,
\( (\I-\M)\mathbf{a} \), multiplies the residuals, and the residuals carry no information about
\( \bbeta \).
:::

The decomposition @eq-opt-lue-residual gives the same result in the language of estimators. Every LUE is
\( \blambda\T\hbeta+\bw\T\he \), and the two terms are uncorrelated (@exr-opt-unbiased-zero), so their
variances add. An estimator that uses the residuals only adds uncorrelated noise to least squares.

## What the matrix form buys

The matrix form @eq-opt-gm-matrix is stronger than the scalar statement applied coordinate by coordinate,
and it is the form usually quoted for vector estimators.

::: {#cor-opt-gm-consequences}
[Consequences of the Gauss–Markov theorem]

Assume @eq-opt-model, let \( \C(\bLambda)\subseteq\C(\X\T) \), and let \( \tilde{\boldsymbol{\uptheta}}=\A\Y \) be any LUE of
\( \boldsymbol{\uptheta}=\bLambda\T\bbeta \), with \( \hat{\boldsymbol{\uptheta}}=\bLambda\T\hbeta \).

::: {.enumerate options="label=(\alph*)"}
1. \( \Var(\bb\T\tilde{\boldsymbol{\uptheta}})\ge\Var(\bb\T\hat{\boldsymbol{\uptheta}}) \) for every \( \bb\in\Real^q \); in particular
   each coordinate of \( \hat{\boldsymbol{\uptheta}} \) has the smallest variance.

2. \( \E\norm{\tilde{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}}^2\ge\E\norm{\hat{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}}^2 \).

3. \( \det\Cov(\tilde{\boldsymbol{\uptheta}})\ge\det\Cov(\hat{\boldsymbol{\uptheta}}) \).

4. If \( \rank(\X)=p \), then \( \hbeta=(\X\T\X)^{-1}\X\T\Y \) is the BLUE of \( \bbeta \), and
   \( \Cov(\tilde{\bbeta})-\sigma^2(\X\T\X)^{-1} \) is nonnegative definite for every LUE \( \tilde{\bbeta} \)
   of \( \bbeta \).

5. The fitted vector \( \hY=\M\Y \) is the BLUE of \( \X\bbeta \), and the fitted value \( \hat{y}_i \) is the
   BLUE of \( \x_{(i)}\T\bbeta \), with variance \( \sigma^2h_{ii} \), where \( h_{ii} \) is the leverage.
:::
:::

::: {.proof}
Let \( \mathbf{D}=\Cov(\tilde{\boldsymbol{\uptheta}})-\Cov(\hat{\boldsymbol{\uptheta}})\succeq\mathbf{0} \) (@thm-opt-gauss-markov(d)).
(a) \( \Var(\bb\T\tilde{\boldsymbol{\uptheta}})-\Var(\bb\T\hat{\boldsymbol{\uptheta}})=\bb\T\mathbf{D}\bb\ge0 \) by @def-mat-nnd. The
coordinates are the case \( \bb=\mathbf{e}_k \).
(b) Both estimators are unbiased, so the mean squared errors are the traces of the covariance matrices,
and \( \tr\mathbf{D}=\sum_k\mathbf{e}_k\T\mathbf{D}\mathbf{e}_k\ge0 \).
(c) Write \( \mathbf{C}=\Cov(\hat{\boldsymbol{\uptheta}}) \). If \( \mathbf{C} \) is singular, its determinant is zero and
there is nothing to prove, because covariance matrices are nonnegative definite
(@thm-rv-cov-nnd). Otherwise \( \mathbf{C} \) is positive definite, with a symmetric positive definite
square root (@thm-mat-square-root). The matrix
\( \mathbf{C}^{-1/2}\Cov(\tilde{\boldsymbol{\uptheta}})\mathbf{C}^{-1/2}=\I+\mathbf{C}^{-1/2}\mathbf{D}\mathbf{C}^{-1/2} \) has all
eigenvalues at least \( 1 \): the matrix \( \mathbf{C}^{-1/2}\mathbf{D}\mathbf{C}^{-1/2} \) is nonnegative definite by
@prp-mat-pd-properties(c), so its eigenvalues are nonnegative by @thm-mat-pd-characterizations, and adding
\( \I \) adds \( 1 \) to each of them. Its determinant
\( \det\Cov(\tilde{\boldsymbol{\uptheta}})/\det\mathbf{C} \) is therefore at least \( 1 \).
(d) Take \( \bLambda=\I_p \), which is allowed because \( \C(\X\T)=\Real^p \), and use
\( \Cov(\hbeta)=\sigma^2(\X\T\X)^{-1} \) (@thm-lm-moments).
(e) Take \( \bLambda=\X\T \). Then \( \bLambda\T\hbeta=\X\hbeta=\M\Y \) (@thm-proj-ls-projection),
with covariance \( \sigma^2\M \). Its \( i \)th diagonal entry is \( \sigma^2h_{ii} \) (@def-proj-leverage).
:::

Part (b) says least squares minimizes the total mean squared error among unbiased linear estimators.
Part (c) has a geometric meaning that [Chapter 12](../ch12-intervals-and-bands/index.html) will use. For normal errors, the confidence
ellipsoid for \( \boldsymbol{\uptheta} \) has volume proportional to \( \sqrt{\det\Cov} \), so least squares gives the
smallest ellipsoids.

::: {#exm-opt-gm-simulation}
[The slope estimators under skewed errors]

Return to the three slope estimators of @exm-opt-three-slopes and simulate
\( 200{,}000 \) data sets with \( \beta_0=2 \), \( \beta_1=0.5 \) and errors
\( \varepsilon_i=E_i-1 \), where the \( E_i \) are independent exponential with mean one. These errors
have mean zero and variance one, but they are strongly skewed: most are slightly negative and a few
are large and positive. The average estimates are \( 0.5000 \),
\( 0.4998 \) and \( 0.5001 \), and the simulated variances are
\( 0.00724 \), \( 0.01643 \) and \( 0.01004 \), against the exact
values \( 0.00729 \), \( 0.01653 \) and \( 0.01008 \).
[Figure 7.2.1](#fig-opt-gauss-markov)(a) shows the three sampling distributions. All three are
centred on the truth, and least squares is the most concentrated: its \( \pm2 \) standard deviation bar is
the shortest. The end-point estimator, a difference of just two skewed errors, has a sharp peak but long
tails, so the height of a density is a poor guide to variance. The errors are far from normal and the
theorem does not care.
:::

```{.python .run #cell-gauss-markov-simulate}
import numpy as np

x = np.array([1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 7.0, 8.0, 9.0, 10.5, 12.0])
n = len(x)
X = np.column_stack([np.ones(n), x])
lam = np.array([0.0, 1.0])                      # target: the slope beta_1
a_ls = X @ np.linalg.solve(X.T @ X, lam)        # least squares: a = X (X^T X)^{-1} lambda
a_end = np.zeros(n)
a_end[[0, -1]] = [-1.0, 1.0]
a_end /= x[-1] - x[0]                           # slope through the two end points
low, high = np.arange(n) < n // 2, np.arange(n) >= n // 2
a_grp = (high / high.sum() - low / low.sum()) / (x[high].mean() - x[low].mean())  # group means

estimators = {"ls": a_ls, "end": a_end, "grp": a_grp}
rng = np.random.default_rng(7)
beta, sigma, reps = np.array([2.0, 0.5]), 1.0, 200_000
E = sigma * (rng.exponential(size=(reps, n)) - 1.0)   # skewed errors, mean 0, variance sigma^2
Y = X @ beta + E                                      # one simulated data set per row
sims = {name: Y @ a for name, a in estimators.items()}
for name, est in sims.items():
    print(f"{name:4s} mean {est.mean():.4f}   variance {est.var():.4f}")
```

::: {when-format="html"}
![**Figure 7.2.1.** (a) Sampling distributions of three linear unbiased estimators of a slope, with
skewed errors (200,000 simulated data sets), as smoothed densities, with bars marking the true
value \( \pm2 \) standard deviations. Least squares has the smallest variance, as the
Gauss–Markov theorem requires. (b) In the location model with uniform errors, the midrange,
which is unbiased but not linear, has far smaller variance than the sample mean, which is the
BLUE.](gauss_markov.svg){#fig-opt-gauss-markov width=100%}
:::

::: {when-format="pdf"}
![(a) Sampling distributions of three linear unbiased estimators of a slope, with
skewed errors (200,000 simulated data sets), as smoothed densities, with bars marking the true
value \( \pm2 \) standard deviations. Least squares has the smallest variance, as the
Gauss–Markov theorem requires. (b) In the location model with uniform errors, the midrange,
which is unbiased but not linear, has far smaller variance than the sample mean, which is the
BLUE.](gauss_markov.pdf){width=100%}
:::

## Correlated errors: Aitken's theorem

If \( \Cov(\be)=\sigma^2\V \) with \( \V \) known and positive definite, the theorem applies after a change
of coordinates, and the BLUE becomes the generalized least squares estimate of @def-proj-gls.

::: {#cor-opt-aitken}
[Aitken]

Assume \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\V \), with \( \V \) known and positive definite. For
\( \blambda\in\C(\X\T) \), the BLUE of \( \blambda\T\bbeta \) is \( \blambda\T\hbeta_{\text{GLS}} \), where
\( \hbeta_{\text{GLS}} \) is any solution of \( \X\T\V^{-1}\X\bb=\X\T\V^{-1}\Y \). Its variance is
\( \sigma^2\blambda\T(\X\T\V^{-1}\X)\ginv\blambda \).
:::

::: {.proof}
Let \( \V^{1/2} \) be the symmetric positive definite square root of \( \V \) (@thm-mat-square-root). Put
\( \Z=\V^{-1/2}\Y \) and \( \X_*=\V^{-1/2}\X \). By @thm-rv-linear, \( \E(\Z)=\X_*\bbeta \) and
\( \Cov(\Z)=\sigma^2\I \), so \( \Z \) satisfies @eq-opt-model with model matrix \( \X_* \). Since
\( \V^{-1/2} \) is nonsingular, \( \C(\X_*\T)=\C(\X\T\V^{-1/2})=\C(\X\T) \), so the same functions are
estimable. The map \( \mathbf{a}\T\Y=(\V^{1/2}\mathbf{a})\T\Z \) is a bijection between linear estimators based on
\( \Y \) and those based on \( \Z \). It preserves means and variances, since both sides are the same random
variable. So a BLUE in one model is a BLUE in the other. By @thm-opt-gauss-markov applied to \( \Z \), the
BLUE is \( \blambda\T\bb \) for any least squares solution \( \bb \) of the \( \Z \) problem. That is, \( \bb \)
minimizes \( \norm{\Z-\X_*\bb}^2=(\Y-\X\bb)\T\V^{-1}(\Y-\X\bb) \), which is the generalized least
squares criterion. Its normal equations (@thm-proj-normal-equations) are \( \X_*\T\X_*\bb=\X_*\T\Z \), that is,
\( \X\T\V^{-1}\X\bb=\X\T\V^{-1}\Y \). The variance is part (a) of the theorem with \( \X_* \) in place of \( \X \).
:::

Ordinary least squares remains unbiased when \( \Cov(\Y)=\sigma^2\V \), but in general it is no longer
best. @thm-proj-kruskal characterized the models in which the two estimators coincide, and
@exr-opt-zyskind shows that exactly these models keep ordinary least squares optimal. [Chapter 31](../ch31-general-gauss-markov/index.html)
develops the general theory (@thm-ggm-ols-blue), including singular \( \V \) (@thm-ggm-singular).

The whitening in the proof carries the rest of the least squares toolkit with it. The fitted vector
\( \X\hbeta_{\text{GLS}}=\V^{1/2}\M_*\Z \), where \( \M_* \) is the orthogonal projection onto \( \C(\X_*) \), is the
generalized least squares fit of @def-proj-gls: it is idempotent but not symmetric, an oblique projection
onto \( \C(\X) \) in the sense of @prp-proj-oblique. The residual sum of squares of the whitened problem is
\( \norm{(\I-\M_*)\Z}^2=(\Y-\X\hbeta_{\text{GLS}})\T\V^{-1}(\Y-\X\hbeta_{\text{GLS}}) \), and since
\( \rank(\X_*)=r \), @thm-rv-quadform-mean and @prp-proj-trace-rank give it expectation \( \sigma^2(n-r) \). So
\[
s^2_{\text{GLS}}=\frac{(\Y-\X\hbeta_{\text{GLS}})\T\V^{-1}(\Y-\X\hbeta_{\text{GLS}})}{n-r}
\]
is unbiased for \( \sigma^2 \).

## Design: when is a coefficient estimated best?

The theorem fixes \( \X \) and optimizes over estimators. The variance it leaves,
\( \sigma^2\blambda\T(\X\T\X)\ginv\blambda \), still depends on \( \X \). When the regressor values can be
chosen, one can ask which designs make it small. The Frisch–Waugh–Lovell theorem gives a clean
answer for single coefficients.

::: {#prp-opt-orthogonal-design}
[Orthogonal columns are best]

Let \( \X \) have full column rank, let \( \x_j \) be its \( j \)th column and \( \X_{(j)} \) the other columns,
with projection \( \M_{(j)} \) onto \( \C(\X_{(j)}) \). Under @eq-opt-model,
\[
\Var(\hat{\beta}_j)=\frac{\sigma^2}{\norm{(\I-\M_{(j)})\x_j}^2}\ \ge\ \frac{\sigma^2}{\norm{\x_j}^2},
\]
with equality iff \( \x_j \) is orthogonal to every other column. If all columns are mutually
orthogonal, then \( \hat{\beta}_j=\x_j\T\Y/\norm{\x_j}^2 \), and it is unchanged when any other columns are
deleted from the model.
:::

::: {.proof}
Put \( \tilde{\x}_j=(\I-\M_{(j)})\x_j \). By @thm-proj-fwl(a), with \( \X_1=\x_j \) and \( \X_2=\X_{(j)} \),
\( \hat{\beta}_j=\tilde{\x}_j\T\Y/\norm{\tilde{\x}_j}^2 \). So \( \Var(\hat{\beta}_j)=\sigma^2/\norm{\tilde{\x}_j}^2 \)
by @thm-rv-linear. By @prp-proj-trace-rank(d) applied to \( \I-\M_{(j)} \),
\( \norm{\tilde{\x}_j}\le\norm{\x_j} \), with equality iff \( \M_{(j)}\x_j=\bzero \), that is, iff
\( \x_j\perp\C(\X_{(j)}) \). If all columns are orthogonal, \( \tilde{\x}_j=\x_j \), and the formula for
\( \hat{\beta}_j \) involves no other column, so deleting columns does not change it.
:::

So among designs whose columns have given lengths, orthogonal designs minimize the variance of
every coefficient at once. This is one reason designed experiments are built from orthogonal
contrasts ([Part IV](../ch15-anova-subspaces/index.html)), and why collinearity, meaning near-dependence among columns, inflates variances (@thm-col-variance).
@exr-opt-weighing applies the proposition to a problem of weighing objects on a balance.

## What the theorem does not say

The Gauss–Markov theorem is often quoted as "least squares is optimal". The theorem is narrower
than that. It compares least squares only with estimators that are linear and unbiased, and it
assumes the model's mean and covariance are right. Relax any of these and least squares can lose.

**Nonlinear estimators.** When the errors are not normal, a nonlinear unbiased estimator can
beat least squares by a wide margin.

::: {#exm-opt-midrange}
[Uniform errors and the midrange]

In the location model \( Y_i=\mu+\varepsilon_i \), let the errors be independent and uniform on
\( (-1,1) \), so \( \sigma^2=1/3 \). The BLUE of \( \mu \) is \( \bar{Y} \), with variance \( 1/(3n) \). The
**midrange** \( (Y_{(1)}+Y_{(n)})/2 \), the average of the smallest and largest observation, is not linear.
It is unbiased, because the error distribution is symmetric about zero, and its variance is
\( 2/\{(n+1)(n+2)\} \) (@exr-opt-midrange-variance). For \( n=20 \) the two variances are
\( 0.0167 \) and \( 0.00433 \), a ratio of
\( 3.9 \). For \( n=100 \) they are \( 0.0033 \) and
\( 0.00019 \), a ratio of \( 17.2 \).
The variance of the mean falls like \( 1/n \) and that of the midrange like \( 1/n^2 \)
([Figure 7.2.1](#fig-opt-gauss-markov)(b)).
:::

```{.python .run #cell-gauss-markov-midrange}
ns = np.array([5, 10, 20, 50, 100, 200])
var_mean = 1 / (3 * ns)                          # errors uniform on (-1, 1): variance 1/3
var_mid = 2 / ((ns + 1) * (ns + 2))              # exact variance of the midrange
for m, v1, v2 in zip(ns, var_mean, var_mid):
    print(f"n = {m:3d}   Var(mean) = {v1:.5f}   Var(midrange) = {v2:.6f}   ratio = {v1 / v2:.1f}")
```

Uniform errors have sharp edges, and the extreme observations locate those edges very precisely. A
linear estimator must give every observation a fixed weight whatever the sample looks like, so it
cannot exploit this. Heavy-tailed errors produce the opposite phenomenon: there the extremes are
unreliable, and estimators that downweight them, such as medians and M-estimators, beat least
squares (@exr-opt-laplace-mle; [Chapter 20](../ch20-residuals-leverage-influence/06-robust.html) and [Chapter 45](../ch45-quantile-gamlss/index.html) return to such estimators). [Section 7.4](04-minimum-variance.html) shows that under *normal* errors
neither phenomenon occurs, and least squares is best among all unbiased estimators.

**Biased estimators.** Mean squared error is variance plus squared bias. Accepting a little bias can
reduce variance by more. For an estimator \( \tilde{\bbeta} \) of \( \bbeta \), the natural matrix
criterion is
\[
\E\bigl[(\tilde{\bbeta}-\bbeta)(\tilde{\bbeta}-\bbeta)\T\bigr]
=\Cov(\tilde{\bbeta})+\bigl(\E\tilde{\bbeta}-\bbeta\bigr)\bigl(\E\tilde{\bbeta}-\bbeta\bigr)\T ,
\]
and the Gauss–Markov theorem ranks estimators by it only within the unbiased class. The ridge estimator
\( (\X\T\X+k\I)^{-1}\X\T\Y \) is linear but biased, and some \( k>0 \) always gives a smaller total mean
squared error than least squares (@exr-opt-ridge). Stein (1956) and James and Stein (1961) showed
more: in the model \( \Y\sim\Normal_n(\bbeta,\sigma^2\I) \) with known \( \sigma^2 \) and
\( n\ge3 \), a nonlinear shrinkage estimator has smaller total mean squared error than \( \Y \) itself
*for every* \( \bbeta \) (@thm-shr-james-stein). [Chapter 27](../ch27-shrinkage/index.html) develops shrinkage. [Section 7.6](06-bayes-conjugate.html) derives
a family of biased linear estimators from prior information.

**Wrong assumptions.** If \( \Cov(\be)\ne\sigma^2\I \), @cor-opt-aitken replaces least squares by
generalized least squares. If the mean is misspecified, so that \( \E(\Y)\notin\C(\X) \), then
\( \E(\hY)=\M\E(\Y) \) is the projection of the true mean. Least squares then estimates the best
approximation within the model, not the truth (compare @thm-proj-blp). [Part V](../ch19-theory-of-departures/index.html) is about diagnosing
and repairing these departures.

## Exercises

### A. Check your understanding

::: {#exr-opt-fitted-blue}
[A1]

Show that \( \sum_i\Var(\hat{y}_i)=\sigma^2r \) under @eq-opt-model. Using
@cor-opt-gm-consequences(e), deduce that no linear unbiased estimator of \( \X\bbeta \) has total variance
less than \( \sigma^2r \). Interpret this as "each estimated dimension of the mean costs one \( \sigma^2 \)".
:::

### B. Practice

::: {#exr-opt-ridge}
[B1]

Let \( \X \) have full column rank and, for \( k\ge0 \), let
\( \tilde{\bbeta}_k=(\X\T\X+k\I)^{-1}\X\T\Y \). Under @eq-opt-model, find the bias and covariance matrix of
\( \tilde{\bbeta}_k \). Using a spectral decomposition \( \X\T\X=\Q\bD\Q\T \), show that the total mean
squared error is
\[
f(k)=\sum_{j=1}^p\frac{\sigma^2d_j+k^2\gamma_j^2}{(d_j+k)^2},\qquad \bgamma=\Q\T\bbeta,\quad \bD=\diag(d_1,\dot{s},d_p),
\]
and that \( f'(0)<0 \). Conclude that some \( k>0 \) beats least squares in total mean squared error,
for every \( \bbeta \) and \( \sigma^2 \), although the best \( k \) depends on them.
:::

::: {.solution}
Write \( \bS=\X\T\X \) and \( \W_k=(\bS+k\I)^{-1} \). Then \( \E\tilde{\bbeta}_k=\W_k\bS\bbeta=\bbeta-k\W_k\bbeta \),
so the bias is \( -k\W_k\bbeta \), and \( \Cov(\tilde{\bbeta}_k)=\sigma^2\W_k\bS\W_k \). With
\( \bS=\Q\bD\Q\T \), \( \W_k=\Q(\bD+k\I)^{-1}\Q\T \), so
\( \tr\Cov=\sigma^2\sum_jd_j/(d_j+k)^2 \) and
\( \norm{k\W_k\bbeta}^2=k^2\sum_j\gamma_j^2/(d_j+k)^2 \). Adding gives \( f(k) \). Differentiating,
\[
f'(k)=\sum_j\frac{2d_j(k\gamma_j^2-\sigma^2)}{(d_j+k)^3},
\qquad f'(0)=-2\sigma^2\sum_jd_j^{-2}<0 .
\]
Since \( f \) is differentiable at \( 0 \) with negative derivative, \( f(k)<f(0) \) for all small \( k>0 \), and
\( f(0) \) is the total mean squared error of \( \hbeta \). This does not contradict @thm-opt-gauss-markov,
because \( \tilde{\bbeta}_k \) is biased for \( k>0 \).
:::

::: {#exr-opt-consistency}
[B2]

Suppose the design grows with \( n \) and \( \X\T\X \) has full rank for all large \( n \). Show that
\( \Var(\blambda\T\hbeta)\le\sigma^2\norm{\blambda}^2/\lambda_{\min}(\X\T\X) \). Deduce that if the smallest
eigenvalue of \( \X\T\X \) tends to infinity, then \( \blambda\T\hbeta\to\blambda\T\bbeta \) in mean square,
for every \( \blambda \).
:::

::: {#exr-opt-ols-gls-mean}
[B3]

In the location model with \( \Cov(\Y)=\sigma^2\diag(v_1,\dot{s},v_n) \), compare the variance of the
ordinary mean \( \bar{Y} \) with that of the BLUE given by @cor-opt-aitken. Show that the ratio is at least
one, with equality iff all \( v_i \) are equal.
:::

### C. Going deeper

::: {#exr-opt-weighing}
[C1]

Each of \( k \) objects with unknown weights \( \beta_1,\dot{s},\beta_k \) is weighed in \( n \) weighings on a
two-pan balance. In weighing \( i \) every object is put on the left pan or the right pan, and the reading
\( Y_i \) is the difference of the pan loads plus an error. So \( \X \) is \( n\times k \) with entries
\( \pm1 \), and @eq-opt-model holds.

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \Var(\hat{\beta}_j)\ge\sigma^2/n \) for every \( j \), with equality for all \( j \) iff \( \X\T\X=n\I \).

2. Take \( k=3 \), \( n=4 \), and the two designs whose rows are
   \[
   \X_1:\ (1,1,1),\ (1,-1,-1),\ (-1,1,-1),\ (-1,-1,1);\qquad
   \X_2:\ (1,1,1),\ (1,1,-1),\ (1,-1,1),\ (1,1,1).
   \]
   For each, compute \( \Cov(\hbeta)/\sigma^2 \), its trace and its determinant. Which design is better by each
   criterion, and what do the criteria measure (compare @cor-opt-gm-consequences(b) and (c))?

3. Show that if a \( \pm1 \) matrix with \( \X\T\X=n\I \) exists with \( k=n\ge3 \), then \( n \) is a multiple of \( 4 \).
   (Whether every multiple of \( 4 \) occurs is the Hadamard conjecture, which is open.)
:::
:::

::: {.solution}
(a) Every column has \( \norm{\x_j}^2=n \), so @prp-opt-orthogonal-design gives
\( \Var(\hat{\beta}_j)\ge\sigma^2/n \), with equality iff \( \x_j \) is orthogonal to the other columns. Equality for
all \( j \) means all columns are mutually orthogonal, that is, \( \X\T\X=n\I \).
(b) \( \X_1\T\X_1=4\I \), so \( \Cov(\hbeta)/\sigma^2=\I/4 \), with trace \( 3/4 \) and determinant \( 1/64 \). For
\( \X_2 \), the columns are \( (1,1,1,1) \), \( (1,1,-1,1) \), \( (1,-1,1,1) \), so
\[
\X_2\T\X_2=\begin{pmatrix}4&2&2\\2&4&0\\2&0&4\end{pmatrix},\qquad
(\X_2\T\X_2)^{-1}=\frac1{32}\begin{pmatrix}16&-8&-8\\-8&12&4\\-8&4&12\end{pmatrix},
\]
with trace \( 40/32=5/4 \) and determinant \( 1/32 \). \( \X_1 \) wins on both counts, and on every single
variance (\( 1/4 \) against \( 1/2 \), \( 3/8 \), \( 3/8 \)). The trace is the total mean squared error of \( \hbeta \)
divided by \( \sigma^2 \), and the determinant governs the volume of confidence ellipsoids.
(c) Here \( \X \) is square and \( \X\T\X=n\I \) gives \( \X\X\T=n\I \), so the rows are mutually orthogonal too.
Changing the sign of a column keeps this property, so we may take the first row to be all ones. For the first
three rows, let \( a,b,c,d \) count the columns in which rows \( 2 \) and \( 3 \) have signs \( (+,+) \), \( (+,-) \),
\( (-,+) \), \( (-,-) \). Then \( a+b+c+d=n \), orthogonality of row \( 1 \) to rows \( 2 \) and \( 3 \) gives
\( a+b-c-d=0 \) and \( a-b+c-d=0 \), and orthogonality of rows \( 2 \) and \( 3 \) gives \( a-b-c+d=0 \). Adding the four
equations gives \( 4a=n \).
:::

::: {#exr-opt-midrange-variance}
[C2]

Let \( U_1,\dot{s},U_n \) be independent and uniform on \( (0,1) \). The \( k \)th smallest has a
\( \text{Beta}(k,n-k+1) \) distribution, and \( \Cov(U_{(1)},U_{(n)})=1/\{(n+1)^2(n+2)\} \). Derive the variance
\( 2/\{(n+1)(n+2)\} \) of the midrange in @exm-opt-midrange. Show also that the midrange is unbiased
for \( \mu \) whenever the errors have a distribution symmetric about zero with finite mean.
:::

::: {.solution}
\( \Var(U_{(1)})=\Var(U_{(n)})=n/\{(n+1)^2(n+2)\} \) from the beta variance. Hence
\[
\Var\Bigl(\frac{U_{(1)}+U_{(n)}}2\Bigr)=\frac14\cdot\frac{2n+2}{(n+1)^2(n+2)}=\frac1{2(n+1)(n+2)} .
\]
Errors uniform on \( (-1,1) \) are \( 2U_i-1 \), which multiplies variances by \( 4 \) and gives
\( 2/\{(n+1)(n+2)\} \). For unbiasedness, if \( \be \) and \( -\be \) have the same distribution, then
\( \varepsilon_{(1)}+\varepsilon_{(n)} \) and \( -(\varepsilon_{(n)}+\varepsilon_{(1)}) \) have the same
distribution. So the midrange minus \( \mu \) is symmetric about zero, and it has mean zero because
\( \lvert\varepsilon_{(1)}+\varepsilon_{(n)}\rvert\le2\max_i\lvert\varepsilon_i\rvert \) is integrable.
:::

::: {#exr-opt-zyskind}
[C3]

Let \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\V \) with \( \V \) positive definite. Show that an LUE
\( \mathbf{a}\T\Y \) of \( \blambda\T\bbeta \) is a BLUE iff it is uncorrelated with every linear unbiased estimator of
zero. Deduce that \( \blambda\T\hbeta \) (ordinary least squares) is the BLUE for every estimable
\( \blambda\T\bbeta \) iff \( \C(\V\X)\subseteq\C(\X) \), in agreement with @thm-proj-kruskal.
:::

::: {.solution}
Every LUE is \( \mathbf{a}\T\Y+\bu\T\Y \) with \( \bu\T\Y \) an unbiased estimator of zero, and
\( \Var(\mathbf{a}\T\Y+t\bu\T\Y)=\Var(\mathbf{a}\T\Y)+2t\Cov(\mathbf{a}\T\Y,\bu\T\Y)+t^2\Var(\bu\T\Y) \) for every real
\( t \). If \( \mathbf{a}\T\Y \) is a BLUE, this quadratic in \( t \) is minimized at \( t=0 \), which forces the
covariance to vanish. Conversely, if every such covariance vanishes, the variance of every competitor is
\( \Var(\mathbf{a}\T\Y)+\Var(\bu\T\Y)\ge\Var(\mathbf{a}\T\Y) \). By @exr-opt-unbiased-zero, the unbiased
estimators of zero are \( \bw\T(\I-\M)\Y \). So ordinary least squares \( \boldsymbol{\uprho}\T\M\Y \) is a BLUE for every
\( \boldsymbol{\uprho} \) iff \( \sigma^2\boldsymbol{\uprho}\T\M\V(\I-\M)\bw=0 \) for all \( \boldsymbol{\uprho},\bw \), that is, iff
\( \M\V(\I-\M)=\mathbf{0} \). Transposing, this says \( (\I-\M)\V\M=\mathbf{0} \), so every column of \( \V\M \) lies
in \( \C(\X) \). Since \( \C(\V\M)=\C(\V\X) \), this is \( \C(\V\X)\subseteq\C(\X) \).
:::
