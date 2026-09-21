# Projection in the population: best linear prediction

So far the geometry has lived in observation space \( \Real^n \), and every statement was
about one fixed data vector. There is a parallel geometry one level up, in which the
vectors are *random variables* and the inner product is an expectation. In that
geometry, least squares estimates something even when the linear model is false,
namely the population projection of \( Y \) onto a span of regressors. In this section
\( \mathbf{X} \) is a *random* vector of regressors, not the model matrix \( \X \). The conditional
expectation turns out to be a projection too. This section sets up the correspondence.
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) develops prediction theory from it, and
[Chapter 19](../ch19-theory-of-departures/index.html) uses it to study what least squares estimates under
misspecification.

## Random variables as vectors

Consider random variables with finite second moments on a common probability space.
They form a vector space, and
\[
\inner{U}{V}=\E(UV),\qquad \norm{U}=\bigl(\E U^2\bigr)^{1/2}
\]
is an inner product, provided we identify variables that are equal with probability one,
so that \( \norm{U}=0 \) forces \( U=0 \). This space is usually written \( L^2 \). Two random
variables are orthogonal when \( \E(UV)=0 \). For variables with mean zero, this is the same as
being uncorrelated. The Cauchy–Schwarz inequality @eq-proj-cosine becomes the
statement that correlations lie in \( [-1,1] \). The *cosine* of the angle between two
centred random variables is their correlation.

For a random \( k \)-vector \( \mathbf{X} \) with finite second moments, the set
\( \mathcal L=\{a+\bb\T\mathbf{X}: a\in\Real,\bb\in\Real^k\} \) is a subspace of dimension at most
\( k+1 \). The proofs of @prp-proj-gram-schmidt, @thm-proj-direct-sum and @thm-proj-projection-theorem
used only the axioms of an inner product and a finite spanning set, so they apply
unchanged. Every \( Y \) with \( \E Y^2<\infty \) has a unique nearest point in \( \mathcal L \),
characterized by orthogonality of the error to \( \mathcal L \).

::: {#def-proj-blp}
[Best linear predictor]

The **best linear predictor** of \( Y \) given \( \mathbf{X} \), written \( L(Y\mid\mathbf{X}) \), is the
projection of \( Y \) onto \( \mathcal L \). It is the \( a+\bb\T\mathbf{X} \) that minimizes the mean
squared error \( \E(Y-a-\bb\T\mathbf{X})^2 \).
:::

::: {#thm-proj-blp}
Let \( \bmu_X=\E\mathbf{X} \), \( \mu_Y=\E Y \), \( \bSigma_{XX}=\Cov(\mathbf{X}) \) and
\( \boldsymbol{\upsigma}_{XY}=\Cov(\mathbf{X},Y) \), and suppose \( \bSigma_{XX} \) is positive definite. Then
\[
L(Y\mid\mathbf{X})=\alpha^*+\bbeta^{*\top}\mathbf{X},
\qquad
\bbeta^*=\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY},
\qquad
\alpha^*=\mu_Y-\bbeta^{*\top}\bmu_X .
\]
The prediction error \( U=Y-L(Y\mid\mathbf{X}) \) has mean zero and is uncorrelated with every
component of \( \mathbf{X} \). Moreover
\[
\Var(Y)=\Var\bigl(L(Y\mid\mathbf{X})\bigr)+\E U^2,
\qquad
\frac{\Var\bigl(L(Y\mid\mathbf{X})\bigr)}{\Var(Y)}
=\frac{\boldsymbol{\upsigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY}}{\Var(Y)} .
\]{#eq-proj-population-pythagoras}

:::

::: {.proof}
Orthogonality of \( U=Y-a-\bb\T\mathbf{X} \) to the spanning set \( \{1,X_1,\dots,X_k\} \) of
\( \mathcal L \) gives the *population normal equations*
\[
\E(U)=0,\qquad \E(\mathbf{X} U)=\bzero .
\]
The first gives \( a=\mu_Y-\bb\T\bmu_X \). Substituting into the second and subtracting
\( \bmu_X\E(U)=\bzero \) gives \( \Cov(\mathbf{X},Y)-\Cov(\mathbf{X})\bb=\bzero \), so
\( \bb=\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY} \). By the projection theorem the solution is the
unique minimizer. Since \( U \) has mean zero and is uncorrelated with \( \mathbf{X} \), it is
uncorrelated with \( L(Y\mid\mathbf{X}) \), so the variances add. Finally,
\( \Var(\bbeta^{*\top}\mathbf{X})=\bbeta^{*\top}\bSigma_{XX}\bbeta^*
=\boldsymbol{\upsigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY} \).
:::

The ratio in @eq-proj-population-pythagoras is the population
*squared multiple correlation*, the population counterpart of \( R^2 \). By the same
cosine argument as @thm-proj-r2-cosine, it is the squared correlation between \( Y \)
and \( L(Y\mid\mathbf{X}) \).

## Conditional expectation is also a projection

The linear span \( \mathcal L \) is one subspace of \( L^2 \). A much larger one is
\( \mathcal G=\{g(\mathbf{X}):\E g(\mathbf{X})^2<\infty\} \), all square-integrable functions of
\( \mathbf{X} \). It is infinite-dimensional, but it is closed, and the projection theorem still
holds in it. We accept that fact from measure-theoretic probability and only use its
conclusion.

::: {#prp-proj-conditional-expectation}
Let \( m(\mathbf{X})=\E(Y\mid\mathbf{X}) \) with \( \E Y^2<\infty \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( m(\mathbf{X}) \) is the projection of \( Y \) onto \( \mathcal G \). It minimizes
           \( \E\bigl(Y-g(\mathbf{X})\bigr)^2 \) over \( g \), and \( \E\bigl[(Y-m(\mathbf{X}))g(\mathbf{X})\bigr]=0 \) for all
           \( g(\mathbf{X})\in\mathcal G \);

2. \( L(Y\mid\mathbf{X})=L\bigl(m(\mathbf{X})\mid\mathbf{X}\bigr) \): the best linear predictor of \( Y \)
           equals the best linear predictor of the regression function;

3. \( \E\bigl(Y-L(Y\mid\mathbf{X})\bigr)^2=\E\bigl(Y-m(\mathbf{X})\bigr)^2
        +\E\bigl(m(\mathbf{X})-L(Y\mid\mathbf{X})\bigr)^2 \);

4. if \( m \) is affine, \( m(\mathbf{X})=L(Y\mid\mathbf{X}) \).
:::

:::

::: {.proof}
(a) The orthogonality property follows from the tower rule:
\( \E\bigl[(Y-m(\mathbf{X}))g(\mathbf{X})\bigr]=\E\bigl[g(\mathbf{X})\,\E(Y-m(\mathbf{X})\mid\mathbf{X})\bigr]=0 \).
Orthogonality to the subspace characterizes the nearest point, by the argument of
@thm-proj-projection-theorem, which needs only Pythagoras.
(b) This is @thm-proj-nested(a) in \( L^2 \). Since \( \mathcal L\subseteq\mathcal G \),
projecting onto \( \mathcal L \) can be done in two stages.
(c) \( Y-L=(Y-m)+(m-L) \), and \( m-L\in\mathcal G \) is orthogonal to \( Y-m \) by (a).
(d) If \( m\in\mathcal L \), the second term in (c) is zero.
:::

Part (c) is the population version of [Figure 6.5.1](05-nested.html#fig-proj-nested). The error of the best
linear predictor splits into two orthogonal parts. One is irreducible noise, which no
function of \( \mathbf{X} \) can remove. The other is *approximation error*, the price of
insisting on linearity. Chapters in Part IX reduce the second part by enlarging
\( \mathcal L \) towards \( \mathcal G \).

## Least squares estimates the population projection

Here is the connection with the rest of the chapter. Given data \( (\x_{(i)},y_i) \), put
probability \( 1/n \) on each row. For this *empirical distribution*, the
\( L^2 \) inner product of two variables is \( n^{-1}\sum_i u_iv_i=n^{-1}\bu\T\bv \), which is the
Euclidean inner product of observation space up to the factor \( 1/n \). So the best linear
predictor under the empirical distribution *is* the least squares fit, and the
population normal equations become the normal equations @eq-proj-normal-equations.
Observation-space geometry is \( L^2 \) geometry for the empirical distribution.

As the sample grows, the empirical distribution approaches the population
distribution, and the projections converge with it.

::: {#prp-proj-consistency}
[Consistency for the projection]

Let \( (\mathbf{X}_i,Y_i) \), \( i=1,2,\dots \), be independent copies of \( (\mathbf{X},Y) \) with finite second
moments and \( \bSigma_{XX} \) positive definite. Let \( (\hat{\alpha}_n,\hbeta_n) \) be the least
squares coefficients from the first \( n \) observations, with an intercept. Then
\( (\hat{\alpha}_n,\hbeta_n)\to(\alpha^*,\bbeta^*) \) with probability one. No assumption
about the form of \( \E(Y\mid\mathbf{X}) \), the distribution of the errors, or their variance is
needed.
:::

::: {.proof}
For large \( n \), \( \hbeta_n=\bS_{XX}^{-1}\mathbf{s}_{XY} \), where \( \bS_{XX} \) and \( \mathbf{s}_{XY} \) are
sample covariances with divisor \( n \) ([Section 6.6](06-fwl.html), centring). By the strong law
of large numbers the sample means, second moments and cross-moments converge almost
surely to their population values. So \( \bS_{XX}\to\bSigma_{XX} \) and
\( \mathbf{s}_{XY}\to\boldsymbol{\upsigma}_{XY} \). Matrix inversion is continuous at the nonsingular
\( \bSigma_{XX} \), so eventually \( \bS_{XX} \) is invertible and
\( \hbeta_n\to\bbeta^* \). The intercept follows the same way.
:::

::: {#exm-proj-blp}
[A curved mean]

Let \( X\sim\mathrm{Uniform}(0,2) \) and \( Y=e^X+\varepsilon \) with independent
\( \varepsilon\sim\Normal(0,0.25) \). The regression function is not linear, but
@thm-proj-blp still identifies a best line. A short calculation gives
\( \Cov(X,Y)=1 \) and \( \Var(X)=1/3 \), so the population projection has slope
\( 3.000 \) and intercept \( 0.195 \). The listing fits least
squares lines to samples of increasing size. The estimates settle on these values,
not on any feature of \( e^x \). [Figure 6.11.1](11-population.html#fig-proj-blp) shows the regression function, the
population projection, and a least squares fit.
:::

::: {when-format="html"}
![**Figure 6.11.1.** When the regression function is curved, least squares estimates its best
linear approximation in the \( L^2 \) sense, the population projection. The
approximation error is largest where \( e^x \) bends away from the line.](best_linear_predictor.svg){#fig-proj-blp width=62%}
:::

::: {when-format="pdf"}
![When the regression function is curved, least squares estimates its best
linear approximation in the \( L^2 \) sense, the population projection. The
approximation error is largest where \( e^x \) bends away from the line.](best_linear_predictor.pdf){width=62%}
:::

```{.python .run #cell-best-linear-predictor-population}
import numpy as np
rng = np.random.default_rng(611)

# X ~ Uniform(0, 2), E(Y | X) = exp(X), Var(Y | X) = 0.25.
# Population projection of Y onto span{1, X}: slope Cov(X, Y)/Var(X), intercept E Y - slope * E X.
EX, VarX = 1.0, 4.0 / 12.0
EY = (np.exp(2) - 1) / 2                             # E exp(X)
EXY = (np.exp(2) + 1) / 2                            # E X exp(X)
slope_pop = (EXY - EX * EY) / VarX
intercept_pop = EY - slope_pop * EX

for n in [20, 200, 2000, 20000]:
    x = rng.uniform(0, 2, n)
    y = np.exp(x) + rng.normal(0, 0.5, n)
    b, *_ = np.linalg.lstsq(np.column_stack([np.ones(n), x]), y, rcond=None)
    print(f"n = {n:6d}   intercept {b[0]:.3f}   slope {b[1]:.3f}")
print(f"population    intercept {intercept_pop:.3f}   slope {slope_pop:.3f}")
```

::: {.remark}
[What changes under misspecification]

@prp-proj-consistency is reassuring about the *target*. A least squares
line always estimates a well-defined population quantity. It says nothing about the
*uncertainty* of the estimate. When \( \E(Y\mid\mathbf{X}) \) is not linear, or the variance
is not constant, the approximation error \( m(\mathbf{X})-L(Y\mid\mathbf{X}) \) behaves like extra
noise that depends on \( \mathbf{X} \). The textbook formula \( \sigma^2(\X\T\X)^{-1} \) for the
covariance of \( \hbeta \) is then wrong, even in large samples. The “sandwich”
covariance estimators of [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (@thm-het-sandwich) are designed
for exactly this situation.
:::
