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
\( \bb=\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY} \). By the projection theorem (@thm-proj-projection-theorem) the solution is the
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

## Exercises

### A. Check your understanding

::: {#exr-proj-blp-curved-exact}
[A1]

Verify the numbers in @exm-proj-blp exactly. With \( X\sim\mathrm{Uniform}(0,2) \) and
\( \E(Y\mid X)=e^X \), evaluate \( \E X \), \( \Var(X) \), \( \E Y \) and \( \E(XY) \) in closed form, and show
that the population projection has slope exactly \( 3 \) and intercept exactly \( (e^2-7)/2 \).
Why does neither depend on the error variance?
:::

::: {.solution}
The density is \( 1/2 \) on \( (0,2) \), so \( \E X=1 \) and \( \Var(X)=4/12=1/3 \). Next
\( \E Y=\E e^X=\tfrac12\int_0^2e^x\,dx=(e^2-1)/2 \), and since
\( \int_0^2xe^x\,dx=[(x-1)e^x]_0^2=e^2+1 \), \( \E(XY)=\E(Xe^X)=(e^2+1)/2 \). Hence
\( \Cov(X,Y)=(e^2+1)/2-(e^2-1)/2=1 \) and \( \beta^*=1/(1/3)=3 \), and
\( \alpha^*=(e^2-1)/2-3=(e^2-7)/2 \), the value quoted there. The error enters \( Y \) additively
with mean zero and independently of \( X \), so it changes neither \( \Cov(X,Y) \) nor
\( \E Y-3\E X \); it only enlarges \( \E U^2 \).
:::

::: {#exr-proj-blp-affine}
[A2]

Let \( \A \) be a nonsingular \( k\times k \) matrix and \( \mathbf{d}\in\Real^k \), and let \( a\ne0 \) and
\( b \) be scalars. Show from @def-proj-blp alone, without the formulas of @thm-proj-blp, that
\[
L\bigl(aY+b\mid\A\mathbf{X}+\mathbf{d}\bigr)=a\,L(Y\mid\mathbf{X})+b .
\]
Deduce that the population squared multiple correlation (@prp-rv-multiple-correlation) is
unchanged by a nonsingular affine change of the regressors.
:::

::: {#exr-proj-blp-binary}
[A3]

Let \( X \) take two values only. Show that every function of \( X \) is affine in \( X \), so that
\( \mathcal G=\mathcal L \) and \( L(Y\mid X)=\E(Y\mid X) \): with a binary regressor the best linear
predictor has no approximation error. Extend this to a categorical regressor with \( g \) levels
coded by \( g-1 \) indicators, and say which fitted model of
[Section 6.3](03-orthogonal-projection.html) is its sample counterpart.
:::

::: {.solution}
If \( X\in\{x_0,x_1\} \) with \( x_0\ne x_1 \), then for any function \( \varphi \) the affine
function with \( b=\{\varphi(x_1)-\varphi(x_0)\}/(x_1-x_0) \) and \( a=\varphi(x_0)-bx_0 \)
satisfies \( a+bX=\varphi(X) \) at both values, so \( \mathcal G\subseteq\mathcal L \); the reverse
inclusion is immediate. Since \( m(X)\in\mathcal G \) is then in \( \mathcal L \),
@prp-proj-conditional-expectation(d) gives \( m=L \). With \( g \) levels \( l_1,\dots,l_g \), let
\( Z_k=1\{X=l_k\} \) for \( k=2,\dots,g \). Any function of \( X \) is
\( \varphi(l_1)+\sum_{k\ge2}\{\varphi(l_k)-\varphi(l_1)\}Z_k \), affine in the indicators, so again
\( \mathcal G=\mathcal L \) and \( m=L \). The sample counterpart is the one-way
layout of @exm-proj-oneway-M, whose fit replaces each observation by its group mean and so
reproduces the conditional mean exactly.
:::

### B. Practice

::: {#exr-proj-population-fwl}
[B1]

Split the regressors as \( \mathbf{X}=(X_1,\mathbf{X}_2) \) with \( X_1 \) scalar, write
\( X_{1\cdot2}=X_1-L(X_1\mid\mathbf{X}_2) \) for the part of \( X_1 \) not linearly predictable from the
rest, and assume \( \Var(X_{1\cdot2})>0 \). Show that the coefficient of \( X_1 \) in
\( L(Y\mid\mathbf{X}) \) is
\[
\beta^*_1=\frac{\Cov(X_{1\cdot2},Y)}{\Var(X_{1\cdot2})} .
\]
This is @thm-proj-fwl in the population. What does it say about a regressor that is an exact
linear function of the others?
:::

::: {.solution}
Write \( Y=\alpha^*+\beta^*_1X_1+\bbeta_2^{*\top}\mathbf{X}_2+U \), where \( U \) has mean zero and is
uncorrelated with \( X_1 \) and \( \mathbf{X}_2 \) (@thm-proj-blp). The variable \( X_{1\cdot2} \) is an
affine function of \( (X_1,\mathbf{X}_2) \), so it is uncorrelated with \( U \). By construction it is
uncorrelated with every component of \( \mathbf{X}_2 \) and with the constant, so
\( \Cov(X_{1\cdot2},\bbeta_2^{*\top}\mathbf{X}_2)=0 \). Finally
\( \Cov(X_{1\cdot2},X_1)=\Cov\bigl(X_{1\cdot2},X_{1\cdot2}+L(X_1\mid\mathbf{X}_2)\bigr)=\Var(X_{1\cdot2}) \),
again because the residual of a projection is uncorrelated with the subspace. Collecting terms,
\( \Cov(X_{1\cdot2},Y)=\beta^*_1\Var(X_{1\cdot2}) \). If \( X_1 \) is an exact linear function of
\( \mathbf{X}_2 \), then \( X_{1\cdot2}=0 \), the denominator vanishes, and no coefficient is defined:
this is the population form of rank deficiency, and \( \bSigma_{XX} \) is then singular.
:::

::: {#exr-proj-blp-noise-free}
[B2]

Show that \( \Cov(\mathbf{X},Y)=\Cov\bigl(\mathbf{X},m(\mathbf{X})\bigr) \) and hence
\( \bbeta^*=\bSigma_{XX}^{-1}\Cov\bigl(\mathbf{X},m(\mathbf{X})\bigr) \): the population projection depends on
the conditional distribution of \( Y \) only through the regression function. Which part of
@prp-proj-conditional-expectation is this? Use it to say, in one line, why the slope in
@exm-proj-blp would be \( 3 \) for any error distribution with mean zero, however
heteroscedastic.
:::

::: {.solution}
By the tower rule, \( \E\bigl[\mathbf{X}\{Y-m(\mathbf{X})\}\bigr]=\E\bigl[\mathbf{X}\,\E\{Y-m(\mathbf{X})\mid\mathbf{X}\}\bigr]=\bzero \)
and \( \E\{Y-m(\mathbf{X})\}=0 \), so the two covariance vectors agree; now apply @thm-proj-blp.
This is the content of @prp-proj-conditional-expectation(b), which says that projecting \( Y \)
and projecting \( m(\mathbf{X}) \) onto \( \mathcal L \) give the same answer. In @exm-proj-blp
\( m(x)=e^x \), so \( \Cov(X,Y)=\Cov(X,e^X)=1 \) and \( \beta^*=3 \), whatever the noise does, as long
as \( \E(\varepsilon\mid X)=0 \) and second moments are finite.
:::

::: {#exr-proj-empirical-blp}
[B3]

Make the assertion of this section precise: least squares *is* best linear prediction for the
empirical distribution. Put probability \( 1/n \) on each row \( (\x_{(i)},y_i) \) of a data set with
an intercept, and let \( \mathbf{X} \) and \( Y \) be the resulting random vector and variable.

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \bSigma_{XX}=n^{-1}\tilde{\X}_1\T\tilde{\X}_1 \) and
   \( \boldsymbol{\upsigma}_{XY}=n^{-1}\tilde{\X}_1\T\y \), where \( \tilde{\X}_1 \) is the matrix of centred
   regressors of @prp-proj-leverage-mahalanobis.

2. Deduce from @thm-proj-blp and @prp-lm-centred that \( \bbeta^* \) and \( \alpha^* \) are the
   least squares slopes and intercept, and that the ratio in
   @eq-proj-population-pythagoras is the \( R^2 \) of @thm-proj-r2-cosine.

3. What corresponds to \( \bSigma_{XX} \) being singular?
:::
:::

::: {.solution}
(a) Under the empirical law, \( \E(X_j)=\bar x_j \) and
\( \Cov(X_j,X_k)=n^{-1}\sum_i(x_{ij}-\bar x_j)(x_{ik}-\bar x_k) \), which is the \( (j,k) \) entry of
\( n^{-1}\tilde{\X}_1\T\tilde{\X}_1 \); likewise
\( \Cov(X_j,Y)=n^{-1}\sum_i(x_{ij}-\bar x_j)(y_i-\bar y)=n^{-1}(\tilde{\X}_1\T\y)_j \), since the
centred column is orthogonal to \( \bone \).

(b) The factors \( n^{-1} \) cancel, so
\( \bbeta^*=(\tilde{\X}_1\T\tilde{\X}_1)^{-1}\tilde{\X}_1\T\y \), which is the least squares slope
vector of @prp-lm-centred, and \( \alpha^*=\bar y-\bbeta^{*\top}\bar{\x} \) is its intercept. Then
\( \Var\{L(Y\mid\mathbf{X})\}=n^{-1}\norm{\hY-\bar y\bone}^2 \) and \( \Var(Y)=n^{-1}\norm{\y-\bar y\bone}^2 \),
whose ratio is \( R^2 \).

(c) A singular \( \bSigma_{XX} \) means \( \tilde{\X}_1 \) has linearly dependent columns: the design
is rank deficient after centring, so the coefficients are not identified, exactly
as in @exr-proj-coefficient-set.
:::

### C. Going deeper

::: {#exr-proj-same-moments}
[C1]

The split of the prediction error in @prp-proj-conditional-expectation(c) into noise and
approximation error is invisible to first and second moments. Let \( X \) be uniform on
\( \{-1,0,1\} \) and let \( c>0 \).

::: {.enumerate options="label=(\alph*)"}
1. Let \( Y_1=X+\varepsilon \) with \( \varepsilon \) independent of \( X \), mean zero and variance
   \( c^2\Var(X^2) \). Let \( Y_2=X+c\,(X^2-\tfrac23) \), with no noise at all. Show that
   \( (X,Y_1) \) and \( (X,Y_2) \) have the same means, variances and covariance, hence the same
   best linear predictor and the same \( \E U^2 \).

2. Show that in the first pair the approximation error is zero and all of \( \E U^2 \) is noise,
   while in the second the noise is zero and all of \( \E U^2 \) is approximation error.
   Compute \( \E U^2 \).

3. What does this say about diagnosing a curved mean from a fitted line and its residual
   variance alone? Which later chapter takes the question up, and what extra structure does
   it need?
:::
:::

::: {.solution}
Throughout, \( \E X=0 \), \( \Var(X)=2/3 \), \( \E X^2=2/3 \), \( \E X^3=0 \), \( \E X^4=2/3 \), so
\( \Var(X^2)=2/3-4/9=2/9 \) and \( \Cov(X,X^2)=0 \).

(a) Both responses have mean \( 0 \): \( \E Y_1=0 \) and \( \E Y_2=0+c(2/3-2/3)=0 \). The covariances
agree: \( \Cov(X,Y_1)=\Var(X)=2/3 \) and \( \Cov(X,Y_2)=\Var(X)+c\Cov(X,X^2)=2/3 \). The variances
agree: \( \Var(Y_1)=2/3+c^2\Var(X^2) \) and, since \( X \) and \( X^2-2/3 \) are uncorrelated,
\( \Var(Y_2)=2/3+c^2\Var(X^2) \). By @thm-proj-blp the best linear predictor depends on the joint
law only through these moments, so both are \( L(Y\mid X)=X \), and
\( \E U^2=\Var(Y)-\Var(X)=c^2\Var(X^2)=2c^2/9 \) in both cases.

(b) For \( Y_1 \), \( m_1(X)=X \) is affine, so \( m_1=L \) and the second term of
@prp-proj-conditional-expectation(c) is zero: all of \( \E U^2 \) is
\( \E(Y-m)^2=c^2\Var(X^2) \). For \( Y_2 \), \( Y_2 \) is a function of \( X \), so \( m_2=Y_2 \) and
\( \E(Y-m_2)^2=0 \); the whole of \( \E U^2=\E\{m_2(X)-X\}^2=c^2\E(X^2-2/3)^2=c^2\Var(X^2) \) is
approximation error.

(c) Nothing in the fitted line, the residual variance or the population \( R^2 \) distinguishes
pure noise around a straight line from an exactly curved mean; the two are the same to second
order. Separating them needs replication or an assumed shape for the alternative, which is what
the lack-of-fit tests of [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html)
supply: repeated observations at the same value of \( \mathbf{X} \) estimate the noise on its own, and the rest of
the residual sum of squares is the curvature.
:::

::: {#exr-proj-sign-reversal}
[C2]

A population projection can reverse the sign of a partial effect. Let \( X_2 \) be uniform on
\( \{-1,0,1\} \), let \( Z\sim\Normal(0,1) \) be independent of \( X_2 \), and put
\[
X_1=X_2^2+Z,\qquad Y=m(X_1,X_2)=X_1-11X_2^2 .
\]

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \bSigma_{XX} \) is diagonal with entries \( 11/9 \) and \( 2/3 \), and that
   \( \boldsymbol{\upsigma}_{XY}=(-11/9,\,0)\T \). Deduce
   \( L(Y\mid\mathbf{X})=-6-X_1 \).

2. So the population projection says that \( Y \) falls by one unit for each unit of \( X_1 \),
   while \( \partial m/\partial x_1=1 \) at every point: holding \( x_2 \) fixed, \( Y \) *rises*
   with \( x_1 \). Explain which feature of the design produces the reversal, and check the
   answer against @exr-proj-population-fwl.

3. Add \( X_2^2 \) to the regressor list. Show that the projection then reproduces \( m \)
   exactly, with coefficient \( 1 \) on \( X_1 \).
:::
:::

::: {.solution}
Write \( W=X_2^2 \), so \( \E W=2/3 \), \( \Var(W)=2/9 \) and \( \Cov(W,X_2)=\E X_2^3-\E W\E X_2=0 \).

(a) \( \Var(X_1)=\Var(W)+\Var(Z)=2/9+1=11/9 \), \( \Var(X_2)=2/3 \) and
\( \Cov(X_1,X_2)=\Cov(W,X_2)=0 \), so \( \bSigma_{XX}=\diag(11/9,\,2/3) \), positive definite even
though \( X_1 \) and \( X_2 \) are functionally dependent through \( W \). With \( Y=W+Z-11W=Z-10W \),
\[
\Cov(X_1,Y)=\Cov(W+Z,\,Z-10W)=1-10\Var(W)=1-\tfrac{20}9=-\tfrac{11}9,
\qquad \Cov(X_2,Y)=-10\Cov(X_2,W)=0 .
\]
Hence \( \bbeta^*=(-1,0)\T \), and
\( \alpha^*=\E Y-\bbeta^{*\top}\bmu_X=-\tfrac{20}3+\tfrac23=-6 \).

(b) The regressors are uncorrelated, so \( X_{1\cdot2}=X_1-\E X_1 \) and
@exr-proj-population-fwl gives \( \beta_1^*=\Cov(X_1,Y)/\Var(X_1) \) directly. The trouble is
that the model space contains \( X_2 \) but not \( X_2^2 \), and \( X_2^2 \) is uncorrelated with both
regressors while being a large part of both \( X_1 \) and \( Y \). Its effect therefore cannot be
attributed to \( X_2 \), and lands on \( X_1 \), whose fluctuation it shares. Linear adjustment
adjusts only for what is linearly predictable.

(c) With regressors \( (X_1,X_2,W) \) the function \( m \) is a linear combination of them, so by
@prp-proj-conditional-expectation(d) the projection is \( m \) itself:
\( L(Y\mid X_1,X_2,W)=X_1-11W \), with coefficient \( 1 \) on \( X_1 \) and prediction error zero.
(The three regressors have a nonsingular covariance matrix, since \( Z \) makes \( X_1 \) vary
freely given \( X_2 \).)
:::
