# Prediction: BLUP and Henderson's equations

The twelve laboratories of @exm-mix-proficiency each have a realized effect
\( a_k \). It is not a parameter but the value taken by a random variable, so asking
for it is a prediction problem. The answer is worth the trouble of setting the problem
up properly: it is not the laboratory's own average, and the difference is shrinkage
of exactly the kind [Chapter 27](../ch27-shrinkage/index.html) studied, arriving here
without anybody having decided to be biased.

## The criterion

::: {#def-mix-predictor}
[Mixed functions and linear predictors]

In @def-mix-model, a **mixed function** is a random variable
\( \tau=\blambda\T\bbeta+\mathbf{c}\T\bu \) with known \( \blambda\in\Real^p \) and
\( \mathbf{c}\in\Real^q \). A **linear predictor** of \( \tau \) is a statistic
\( t=\mathbf{a}\T\Y+b \) with known \( \mathbf{a} \) and \( b \); it is **unbiased**
if \( \E(t)=\E(\tau) \) for every \( \bbeta\in\Real^p \). The **best linear unbiased
predictor** (BLUP) of \( \tau \) is the unbiased linear predictor minimizing the
prediction mean squared error \( \E(t-\tau)^2 \).
:::

Since \( \E(t)=\mathbf{a}\T\X\bbeta+b \) and \( \E(\tau)=\blambda\T\bbeta \),
unbiasedness for every \( \bbeta \) forces \( b=0 \) and
\( \X\T\mathbf{a}=\blambda \). So an unbiased linear predictor exists exactly when
\( \blambda\in\C(\X\T) \), that is, when \( \blambda\T\bbeta \) is estimable in the
sense of @def-est-estimable; the random part \( \mathbf{c}\T\bu \) places no
restriction at all, because its mean is zero whatever \( \mathbf{c} \) is. Note also
that the criterion averages over \( \bu \) as well as over \( \be \). This is the
decision that makes the problem solvable, and the warning in the next subsection says
what it costs.

The problem has already been solved.
@thm-cor-blup of
[Section 14.2](../ch14-correlation-lack-of-fit-prediction/02-best-linear-prediction.html)
predicted a new response \( y_0=\x_0\T\bbeta+e_0 \) whose error was correlated with
the training errors, and \( \mathbf{c}\T\bu \) is exactly such an \( e_0 \). What
follows is that theorem in mixed-model notation; its proof is the substitution, plus
the uniqueness and the decomposition (c) that [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) did not need.

::: {#thm-mix-blup}
[Best linear unbiased prediction]

In @def-mix-model let \( \V=\Z\G\Z\T+\R \) be positive definite and let
\( \blambda\in\C(\X\T) \). Put
\[
\hbeta=(\X\T\V^{-1}\X)\ginv\X\T\V^{-1}\Y,\qquad
\hat{\bu}=\G\Z\T\V^{-1}(\Y-\X\hbeta) .
\]{#eq-mix-blup}

::: {.enumerate options="label=(\alph*)"}
1. The BLUP of \( \tau=\blambda\T\bbeta+\mathbf{c}\T\bu \) is
   \( \hat\tau=\blambda\T\hbeta+\mathbf{c}\T\hat{\bu} \). It does not depend on which
   generalized inverse is used, and it is the unique minimizer.

2. Its prediction error variance is
   \[
   \E(\hat\tau-\tau)^2
   =\mathbf{c}\T\bigl(\G-\G\Z\T\V^{-1}\Z\G\bigr)\mathbf{c}
   +\mathbf{d}\T(\X\T\V^{-1}\X)\ginv\mathbf{d},
   \]{#eq-mix-prediction-variance}

   where \( \mathbf{d}=\blambda-\X\T\V^{-1}\Z\G\mathbf{c} \).

3. If \( \bbeta \) were known, the best linear predictor of \( \tau \) would be
   \( \blambda\T\bbeta+\mathbf{c}\T\G\Z\T\V^{-1}(\Y-\X\bbeta) \), with prediction
   error variance the first term of @eq-mix-prediction-variance. The second term is the price of
   estimating \( \bbeta \).
:::

:::

::: {.proof}
Apply @thm-cor-blup to \( \Y=\X\bbeta+(\Z\bu+\be) \), whose error vector has
covariance \( \V \), with \( \sigma^2=1 \), \( y_0=\tau \), \( \x_0=\blambda \) and
\( e_0=\mathbf{c}\T\bu \). Since \( \Cov(\Z\bu+\be,\bu)=\Z\G \), the two covariance
inputs of that theorem are
\[
\bv_0=\Cov(\Z\bu+\be,\,\mathbf{c}\T\bu)=\Z\G\mathbf{c}=:\mathbf{s},\qquad
v_{00}=\Var(\mathbf{c}\T\bu)=\mathbf{c}\T\G\mathbf{c},
\]
and its hypothesis \( \x_0\in\C(\X\T) \) is our \( \blambda\in\C(\X\T) \). So
@eq-cor-blup reads
\( \hat\tau=\blambda\T\hbeta+\mathbf{c}\T\G\Z\T\V^{-1}(\Y-\X\hbeta)
=\blambda\T\hbeta+\mathbf{c}\T\hat{\bu} \), and its mean squared error is
@eq-mix-prediction-variance, because the vector called \( \mathbf{c} \) there is
\( \x_0-\X\T\V^{-1}\bv_0=\blambda-\X\T\V^{-1}\Z\G\mathbf{c}=\mathbf{d} \) here. That
is (a) and (b) apart from the two uniqueness claims.

For those, recall that an unbiased linear predictor is \( t=\mathbf{a}\T\Y \) with
\( \X\T\mathbf{a}=\blambda \), and that \( t-\tau \) then has mean zero, so
\[
\E(t-\tau)^2=\Var(\mathbf{a}\T\Y-\mathbf{c}\T\bu)
=\mathbf{a}\T\V\mathbf{a}-2\mathbf{a}\T\mathbf{s}+\mathbf{c}\T\G\mathbf{c}.
\]
The last term does not involve \( \mathbf{a} \) and \( \V \) is positive definite, so
this is a strictly convex function of \( \mathbf{a} \) on the affine set
\( \{\mathbf{a}:\X\T\mathbf{a}=\blambda\} \), non-empty because
\( \blambda\in\C(\X\T) \), and a strictly convex function has at most one minimizer on
a convex set. So \( \mathbf{a} \), and with it \( \hat\tau \), is unique, by an
argument that nowhere asks \( \X\T \) to have full row rank. The generalized inverse
then cancels, as one can also see directly:
\( \X\hbeta=\X(\X\T\V^{-1}\X)\ginv\X\T\V^{-1}\Y \) is the projection of \( \Y \) onto
\( \C(\X) \) in the \( \V^{-1} \) inner product (@def-proj-gls), and both
\( \blambda\T\hbeta=\mathbf{r}\T\X\hbeta \), for \( \mathbf{r} \) with
\( \blambda=\X\T\mathbf{r} \), and \( \hat{\bu} \) involve \( \hbeta \) only through
\( \X\hbeta \).

For (c), with \( \bbeta \) known the problem is unconstrained prediction of
\( \mathbf{c}\T\bu \) from \( \Y-\X\bbeta \), whose solution is the best linear
predictor of @thm-rv-blp: \( \Cov(\mathbf{c}\T\bu,\Y)\V^{-1}(\Y-\X\bbeta) \), with
error variance
\( \mathbf{c}\T\G\mathbf{c}-\mathbf{s}\T\V^{-1}\mathbf{s} \), which is the first term
of @eq-mix-prediction-variance because \( \mathbf{s}=\Z\G\mathbf{c} \).
:::

Two readings of @eq-mix-blup are worth having. First, \( \hat{\bu} \) is the best
linear predictor of \( \bu \) with \( \bbeta \) replaced by its generalized least
squares estimate; the whole content of "BLUP" is that this substitution is optimal,
not merely convenient. Second, if the model is normal then
\( \E(\bu\mid\Y)=\G\Z\T\V^{-1}(\Y-\X\bbeta) \) by @thm-mvn-conditional, so the BLUP is
the conditional mean with \( \bbeta \) estimated. That is the door into
[Section 32.6](06-bayes.html).

## Shrinkage

Specialize to the balanced one-way model of @def-mix-oneway, with
\( \X=\bone \), \( \Z \) the group indicators, \( \G=\sigma_a^2\I_g \) and
\( \R=\sigma^2\I \). By @eq-mix-balanced-cov,
\( \V^{-1}=\lambda_1^{-1}\M+\sigma^{-2}(\I-\M) \) with
\( \lambda_1=\sigma^2+m\sigma_a^2 \), and since \( \M\Z=\Z \) we have
\( \Z\T\V^{-1}=\lambda_1^{-1}\Z\T \). Also \( \bone\in\C(\Z) \), so the generalized
least squares estimate of \( \mu \) is the ordinary mean \( \bar y
\) (@exm-mix-ols-blue). Hence
\[
\hat a_k=\frac{m\sigma_a^2}{\sigma^2+m\sigma_a^2}\,(\bar y_k-\bar y)
=B\,(\bar y_k-\bar y),
\qquad B=\frac{m\gamma}{1+m\gamma},
\]{#eq-mix-shrinkage}

with \( \gamma=\sigma_a^2/\sigma^2 \). The predicted group mean is
\( \bar y+B(\bar y_k-\bar y) \): the group's own average pulled towards the overall
average by \( B\in(0,1) \). The pull is strong when groups are small or when groups
genuinely differ little, and it vanishes as either \( m \) or \( \gamma \) grows.
Formula @eq-mix-shrinkage is the empirical Bayes shrinkage of
@prp-shr-empirical-bayes with the prior variance supplied by the model, and it is a
ridge estimate of the group effects with penalty \( \sigma^2/\sigma_a^2
\) (@thm-shr-ridge); the next subsection makes that identification exact.

::: {.warning}
Unbiasedness in @def-mix-predictor averages over \( \bu \). Conditionally on the
realized effects it fails, and deliberately so: from @eq-mix-shrinkage,
\( \E(\hat a_k\mid a_1,\dots,a_g)=B\,(a_k-\bar a) \), which is nearer zero than
\( a_k \). Predicted effects are *shrunk*, so their spread understates the spread of
the true effects, and ranking groups by \( \hat a_k \) is not the same as ranking
them by \( a_k \). A league table of hospitals or schools built from BLUPs will move
small units towards the middle, which is right for predicting each unit's next
observation and wrong if the question is "which unit is genuinely worst".
:::

The shrinkage is paid for in mean squared error. Continuing the balanced one-way
case, @eq-mix-prediction-variance with \( \blambda=\bzero \) and \( \mathbf{c}=\mathbf{e}_k \)
gives, after substituting
\( \G-\G\Z\T\V^{-1}\Z\G=\sigma_a^2\sigma^2\lambda_1^{-1}\I \) and
\( \X\T\V^{-1}\X=n/\lambda_1 \),
\[
\E(\hat a_k-a_k)^2
=\frac{\sigma_a^2}{\lambda_1}\Bigl(\sigma^2+\frac{m\sigma_a^2}{g}\Bigr),
\]{#eq-mix-prediction-variance-oneway}

whereas the obvious predictor \( \bar y_k-\bar y \) has mean squared error
\( \sigma^2(1-1/g)/m+\sigma_a^2/g \) (@exr-mix-raw-mse).

::: {#exm-mix-shrinkage}
[Shrinking twelve laboratories]

In @exm-mix-proficiency, \( \hat\sigma^2=0.0620 \) and
\( \hat\sigma_a^2=0.0386 \) give a shrinkage factor
\( B=0.714 \). The laboratory deviations \( \bar y_k-\bar y \) run
from \( -0.450 \) to \( 0.250 \); the
predicted effects run from \( -0.321 \) to
\( 0.178 \). Each has been pulled about
\( 29 \) per cent of the way towards zero.

The gain is real. By @eq-mix-prediction-variance-oneway the prediction error variance is
\( 0.01336 \), a standard error of
\( 0.1156 \), against \( 0.01743 \)
for the raw deviation, a ratio of \( 1.30 \) in mean squared
error. Panel (b) of [Figure 32.3.1](03-blup.html#fig-mix-blup) shows the comparison
over the whole range of intraclass correlations: shrinkage helps most when the
correlation is small, which is exactly when the raw deviations are mostly noise.
:::

::: {when-format="html"}
![**Figure 32.3.1.** (a) The twelve laboratory deviations \( \bar y_k-\bar y \) and
their predicted effects \( \hat a_k \), sorted, with 95% prediction intervals
from @eq-mix-prediction-variance-oneway. (b) Mean squared error of the two predictors of \( a_k \), as a
function of the intraclass correlation, with \( \sigma_a^2+\sigma^2=1 \),
\( g=12 \) and \( m=4 \).](blup_shrinkage.svg){#fig-mix-blup width=100%}
:::

::: {when-format="pdf"}
![(a) The twelve laboratory deviations \( \bar y_k-\bar y \) and their predicted
effects \( \hat a_k \), sorted, with 95% prediction intervals
from @eq-mix-prediction-variance-oneway. (b) Mean squared error of the two predictors of \( a_k \), as a
function of the intraclass correlation, with \( \sigma_a^2+\sigma^2=1 \),
\( g=12 \) and \( m=4 \).](blup_shrinkage.pdf){width=100%}
:::

## Henderson's mixed model equations

Formula @eq-mix-blup asks for \( \V^{-1} \), an \( n\times n \) inverse. Henderson
(1950) found a system of size \( p+q \) that delivers the generalized least squares
estimate and the BLUP together, without ever forming \( \V \).

::: {#thm-mix-henderson}
[Henderson's mixed model equations]

Let \( \G \) and \( \R \) be positive definite and \( \rank(\X)=p \). Define
\[
\mathbf{C}=\begin{pmatrix}
\X\T\R^{-1}\X & \X\T\R^{-1}\Z\\
\Z\T\R^{-1}\X & \Z\T\R^{-1}\Z+\G^{-1}
\end{pmatrix},
\qquad
\mathbf{h}=\begin{pmatrix}\X\T\R^{-1}\y\\ \Z\T\R^{-1}\y\end{pmatrix}.
\]{#eq-mix-henderson}

::: {.enumerate options="label=(\alph*)"}
1. \( \mathbf{C} \) is positive definite, and the unique solution of
   \( \mathbf{C}(\bb\T,\bv\T)\T=\mathbf{h} \) is \( \bb=\hbeta \), the generalized
   least squares estimate, and \( \bv=\hat{\bu} \), the BLUP of @eq-mix-blup.

2. The same pair is the unique minimizer of the penalized least squares criterion
   \[
   Q(\bbeta,\bu)=(\y-\X\bbeta-\Z\bu)\T\R^{-1}(\y-\X\bbeta-\Z\bu)+\bu\T\G^{-1}\bu .
   \]{#eq-mix-penalized}

3. With \( \mathbf{W}=\X\T\V^{-1}\X \), the blocks of \( \mathbf{C}^{-1} \) are
   \[
   (\mathbf{C}^{-1})_{11}=\mathbf{W}^{-1},\qquad
   (\mathbf{C}^{-1})_{12}=-\mathbf{W}^{-1}\X\T\V^{-1}\Z\G,
   \]
   \[
   (\mathbf{C}^{-1})_{22}=\G-\G\Z\T\V^{-1}\Z\G
   +\G\Z\T\V^{-1}\X\mathbf{W}^{-1}\X\T\V^{-1}\Z\G .
   \]
   Its leading block is \( \Cov(\hbeta) \) and its trailing block is
   \( \Cov(\hat{\bu}-\bu) \), the matrix whose quadratic forms are @eq-mix-prediction-variance.
:::

:::

::: {.proof}
We use repeatedly the *push-through identity*
\[
(\G^{-1}+\Z\T\R^{-1}\Z)^{-1}\Z\T\R^{-1}=\G\Z\T\V^{-1},
\]{#eq-mix-push-through}

proved by multiplying on the right by \( \V=\R+\Z\G\Z\T \):
\[
\begin{aligned}
&(\G^{-1}+\Z\T\R^{-1}\Z)^{-1}\bigl(\Z\T+\Z\T\R^{-1}\Z\G\Z\T\bigr)\\
&\qquad=(\G^{-1}+\Z\T\R^{-1}\Z)^{-1}(\G^{-1}+\Z\T\R^{-1}\Z)\G\Z\T\\
&\qquad=\G\Z\T.
\end{aligned}
\]

(a) For \( (\bb\T,\bv\T)\ne\bzero \),
\[
(\bb\T,\bv\T)\mathbf{C}(\bb\T,\bv\T)\T
=\norm{\R^{-1/2}(\X\bb+\Z\bv)}^2+\bv\T\G^{-1}\bv ,
\]
which vanishes only if \( \bv=\bzero \) and \( \X\bb=\bzero \), hence only if
\( \bb=\bzero \) by full column rank. So \( \mathbf{C} \) is positive definite and
the system has a unique solution. The second block row gives
\[
\begin{aligned}
\bv&=(\Z\T\R^{-1}\Z+\G^{-1})^{-1}\Z\T\R^{-1}(\y-\X\bb)\\
&=\G\Z\T\V^{-1}(\y-\X\bb)
\end{aligned}
\]
by @eq-mix-push-through. Substituting into the first block row, and using the
residual identity
\[
\begin{aligned}
\y-\X\bb-\Z\bv&=(\I-\Z\G\Z\T\V^{-1})(\y-\X\bb)\\
&=(\V-\Z\G\Z\T)\V^{-1}(\y-\X\bb)\\
&=\R\V^{-1}(\y-\X\bb),
\end{aligned}
\]
the first block row \( \X\T\R^{-1}(\y-\X\bb-\Z\bv)=\bzero \) becomes
\( \X\T\V^{-1}(\y-\X\bb)=\bzero \), the generalized normal equations
of @def-proj-gls. Hence \( \bb=\hbeta \) and then \( \bv=\hat{\bu} \).

(b) \( Q \) is a quadratic in \( (\bbeta,\bu) \) with matrix \( \mathbf{C} \) and
linear term \( -2\mathbf{h} \), and the display in (a) shows it is strictly convex.
By @prp-mat-quadratic-min its unique minimizer solves
\( \mathbf{C}(\bb\T,\bv\T)\T=\mathbf{h} \), which is @eq-mix-henderson, whose
solution is \( (\hbeta,\hat{\bu}) \) by part (a).

(c) Write \( \mathbf{C}_{11}=\X\T\R^{-1}\X \),
\( \mathbf{C}_{12}=\X\T\R^{-1}\Z \) and
\( \mathbf{C}_{22}=\Z\T\R^{-1}\Z+\G^{-1} \), so that by @eq-mix-v-inverse the Schur
complement of \( \mathbf{C}_{22} \) is
\[
\begin{aligned}
&\mathbf{C}_{11}-\mathbf{C}_{12}\mathbf{C}_{22}^{-1}\mathbf{C}_{12}\T\\
&\qquad=\X\T\bigl(\R^{-1}-\R^{-1}\Z(\G^{-1}+\Z\T\R^{-1}\Z)^{-1}
\Z\T\R^{-1}\bigr)\X\\
&\qquad=\X\T\V^{-1}\X=\mathbf{W}.
\end{aligned}
\]
@thm-mat-partitioned-inverse now gives
\( (\mathbf{C}^{-1})_{11}=\mathbf{W}^{-1} \) and
\( (\mathbf{C}^{-1})_{12}=-\mathbf{W}^{-1}\mathbf{C}_{12}\mathbf{C}_{22}^{-1}
=-\mathbf{W}^{-1}\X\T\V^{-1}\Z\G \), the last step being the transpose
of @eq-mix-push-through. For the trailing block,
\[
\begin{aligned}
(\mathbf{C}^{-1})_{22}&=\mathbf{C}_{22}^{-1}
+\mathbf{C}_{22}^{-1}\mathbf{C}_{12}\T\mathbf{W}^{-1}\mathbf{C}_{12}
\mathbf{C}_{22}^{-1}\\
&=\mathbf{C}_{22}^{-1}+\G\Z\T\V^{-1}\X\mathbf{W}^{-1}\X\T\V^{-1}\Z\G,
\end{aligned}
\]
and \( \mathbf{C}_{22}^{-1}=\G-\G\Z\T\V^{-1}\Z\G \), because
\[
\begin{aligned}
&(\G-\G\Z\T\V^{-1}\Z\G)(\G^{-1}+\Z\T\R^{-1}\Z)\\
&\qquad=\I+\G\Z\T\R^{-1}\Z-\G\Z\T\V^{-1}(\Z\G\Z\T+\R)\R^{-1}\Z=\I,
\end{aligned}
\]
using @eq-mix-push-through once more in the form
\( \G\Z\T\V^{-1}\V=\G\Z\T \). Finally, \( \Cov(\hbeta)=\mathbf{W}^{-1} \) is
@cor-opt-aitken, and the trailing block is @eq-mix-prediction-variance with
\( \blambda=\bzero \) read off entry by entry.
:::

::: {.idea}
Henderson's system says that a mixed model *is* a penalized least squares
problem (@eq-mix-penalized): fit the random effects as though they were free parameters, but
charge \( \bu\T\G^{-1}\bu \) for using them. With \( \R=\sigma^2\I \) and
\( \G=\sigma_u^2\I \) the criterion is
\( \norm{\y-\X\bbeta-\Z\bu}^2+\lambda\norm{\bu}^2 \) with
\( \lambda=\sigma^2/\sigma_u^2 \): ridge regression on the random-effect columns
only (@thm-shr-ridge). The difference from
[Chapter 27](../ch27-shrinkage/index.html) is that the penalty is not chosen by
cross-validation; it is estimated, because the model says what
\( \sigma^2/\sigma_u^2 \) means. [Chapter 43](../ch43-smoothing/index.html) exploits the same equivalence
in the other direction, writing a penalized spline as a mixed model so that its
smoothing parameter can be estimated by REML (@thm-smo-mixed).
:::

```{.python .run #cell-blup-setup}
import numpy as np

g, m = 12, 4
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

X = np.ones((n, 1))                             # one fixed effect, the grand mean
Z = np.repeat(np.eye(g), m, axis=0)             # laboratory indicators

Y = y.reshape(g, m)                             # moment estimates, Section 32.1
lab_mean = Y.mean(axis=1)
ms_between = m * np.sum((lab_mean - Y.mean()) ** 2) / (g - 1)
ms_within = np.sum((Y - lab_mean[:, None]) ** 2) / (n - g)
s2, s2a = ms_within, (ms_between - ms_within) / m

G = s2a * np.eye(g)
R = s2 * np.eye(n)
V = Z @ G @ Z.T + R
```

```{.python .run #cell-blup-blup}
Vinv = np.linalg.inv(V)
beta_gls = np.linalg.solve(X.T @ Vinv @ X, X.T @ Vinv @ y)
u_blup = G @ Z.T @ Vinv @ (y - X @ beta_gls)

# Henderson's mixed model equations: one (p + q) x (p + q) system, no n x n inverse
Rinv = np.linalg.inv(R)
C = np.block([[X.T @ Rinv @ X, X.T @ Rinv @ Z],
              [Z.T @ Rinv @ X, Z.T @ Rinv @ Z + np.linalg.inv(G)]])
rhs = np.concatenate([X.T @ Rinv @ y, Z.T @ Rinv @ y])
solution = np.linalg.solve(C, rhs)

print("GLS estimate   ", beta_gls, solution[:1])
print("first two BLUPs", u_blup[:2], solution[1:3])
```

```{.python .run #cell-blup-shrink}
shrinkage = m * s2a / (s2 + m * s2a)            # the shrinkage factor B
print(f"B = {shrinkage:.4f}")
print("u_hat  ", np.round(u_blup[:4], 4))
print("B(ybar_k - ybar)", np.round(shrinkage * (lab_mean - y.mean())[:4], 4))
```

Two practical points close the section. First, \( \mathbf{C} \) is sparse whenever
\( \Z \) is, which is always: each row of \( \Z \) has one non-zero entry per random
factor, so @eq-mix-henderson can be solved by sparse Cholesky factorization, and a
model with hundreds of thousands of random effects is routine. Second, \( \G \) and
\( \R \) are unknown in practice and estimates are substituted; the result is called
the **empirical** BLUP, and @eq-mix-prediction-variance is then too small, since it
ignores the variability of the estimated covariance. Kackar and Harville (1984)
identified the extra term, which returns in [Section 32.5](05-inference.html) as the
Kenward–Roger correction and in [Section 32.6](06-bayes.html) as the gap between
posterior and plug-in standard errors.

## Exercises

### A. Check your understanding

::: {#exr-mix-blup-limits}
[A1]

In @eq-mix-shrinkage, what happens to \( B \) as (a) \( m\to\infty \),
(b) \( \sigma_a^2\to0 \), (c) \( \sigma^2\to0 \)? Describe in words the predictor in
each limit, and say which of the three limits makes the mixed model agree with a
fixed-effects analysis.
:::

### B. Practice

::: {#exr-mix-raw-mse}
[B1]

Show that the predictor \( \bar Y_k-\bar Y \) of \( a_k \) in the balanced one-way
model has \( \E(\bar Y_k-\bar Y-a_k)^2=\sigma^2(1-1/g)/m+\sigma_a^2/g \), and verify
that @eq-mix-prediction-variance-oneway is smaller for every \( \gamma>0 \). Where does the
difference come from when \( \gamma \) is very small?
:::

::: {.solution}
\( \bar Y_k-\bar Y-a_k=\bar\varepsilon_k-\bar\varepsilon-\bar a \) with
\( \bar a=g^{-1}\sum_ja_j \). The three terms have variances \( \sigma^2/m \),
\( \sigma^2/(gm) \) and \( \sigma_a^2/g \), and
\( \Cov(\bar\varepsilon_k,\bar\varepsilon)=\sigma^2/(gm) \), giving
\( \sigma^2/m-\sigma^2/(gm)+\sigma_a^2/g \) as stated. In terms of
\( \gamma=\sigma_a^2/\sigma^2 \), @eq-mix-prediction-variance-oneway equals
\( \sigma^2\gamma(1+m\gamma/g)/(1+m\gamma) \) while the raw predictor's error is
\( \sigma^2\{(1-1/g)/m+\gamma/g\} \). Putting the two over the denominator
\( 1+m\gamma \) and cancelling, their difference is exactly
\[
\frac{\sigma^2}{m}\Bigl(1-\frac1g\Bigr)\frac{1}{1+m\gamma}
=\frac{\sigma^2}{m}\Bigl(1-\frac1g\Bigr)(1-B),
\]
which is strictly positive for every \( \gamma \). As \( \gamma\to0 \) the BLUP's
error tends to
\( \sigma_a^2 \) while the raw predictor's tends to \( \sigma^2/m \), which is much
larger: the raw predictor pays the full measurement error even when there is nothing
to predict, whereas the BLUP predicts \( 0 \) and pays only \( \sigma_a^2 \).
:::

::: {#exr-mix-blup-of-mean}
[B2]

Take \( \tau=\mu+a_k \), the mean of group \( k \) in the balanced one-way model, so
that \( \blambda=1 \) and \( \mathbf{c}=\mathbf{e}_k \). Use @eq-mix-blup and
@eq-mix-prediction-variance to show that the BLUP is \( \bar y+B(\bar y_k-\bar y) \) with prediction
error variance
\( \sigma^2B/m+(1-B)^2\lambda_1/n \). Show that this equals \( \lambda_1/n \) at
\( B=0 \) and \( \sigma^2/m \) at \( B=1 \); that it is strictly below
\( \sigma^2/m \), the prediction error of the group's own mean; and that it is below
\( \lambda_1/n \) as well exactly when \( m\sigma_a^2>\sigma^2(g-2) \), so that it is
no interpolation between the two — with \( g=3 \), \( m=4 \) and
\( \sigma^2=\sigma_a^2=1 \) it beats both.
:::

::: {.solution}
The predictor is immediate from @eq-mix-blup. For the variance,
\( \mathbf{d}=1-\bone\T\V^{-1}\Z\G\mathbf{e}_k=1-m\sigma_a^2/\lambda_1=1-B \), and
\( \mathbf{W}\ginv=\lambda_1/n \), so the second term of @eq-mix-prediction-variance is
\( (1-B)^2\lambda_1/n \); the first term is
\( \sigma_a^2\sigma^2/\lambda_1=\sigma^2B/m \) since
\( B=m\sigma_a^2/\lambda_1 \). Setting \( B=0 \) and \( B=1 \) gives the two stated
endpoints. For the comparisons, use \( 1-B=\sigma^2/\lambda_1 \) and \( n=gm \) to
collapse the whole expression to
\[
\frac{\sigma^2}{\lambda_1}\Bigl(\sigma_a^2+\frac{\sigma^2}{gm}\Bigr),
\]
and note that \( \sigma^2/m=(\sigma^2/\lambda_1)(\sigma_a^2+\sigma^2/m) \), so the
first comparison holds whenever \( g>1 \). The second asks
\( \sigma^2\sigma_a^2gm+\sigma^4<\lambda_1^2=(\sigma^2+m\sigma_a^2)^2 \), which after
cancelling \( \sigma^4 \) and dividing by \( m\sigma_a^2 \) is
\( \sigma^2g<2\sigma^2+m\sigma_a^2 \). In the numerical case the expression is
\( 13/60 \), below both \( 5/12 \) and \( 1/4 \): the quadratic in \( B \) is not an
interpolation between two strategies but can beat both, and that gain is what the
shrinkage buys.
:::

::: {#exr-mix-henderson-ridge}
[B3]

Take \( \X \) empty (no fixed effects), \( \R=\sigma^2\I \) and \( \G=\sigma_u^2\I \)
in @thm-mix-henderson. Show that \( \hat{\bu} \) is exactly the ridge estimate
\( (\Z\T\Z+\lambda\I)^{-1}\Z\T\y \) with \( \lambda=\sigma^2/\sigma_u^2 \), and that
@eq-mix-prediction-variance becomes \( \sigma_u^2\mathbf{c}\T(\I+\lambda^{-1}\Z\T\Z)^{-1}\mathbf{c} \).
Compare with the mean squared error of ridge in @thm-shr-ridge: which is the same
quantity and which is not?
:::

::: {.solution}
With no \( \X \), @eq-mix-henderson reduces to
\( (\Z\T\Z/\sigma^2+\I/\sigma_u^2)\bv=\Z\T\y/\sigma^2 \), that is
\( \bv=(\Z\T\Z+\lambda\I)^{-1}\Z\T\y \). By @thm-mix-henderson(c) the prediction
error covariance is \( \mathbf{C}^{-1}=\sigma^2(\Z\T\Z+\lambda\I)^{-1}
=\sigma_u^2(\I+\lambda^{-1}\Z\T\Z)^{-1} \). The quantity computed here averages over
\( \bu \); the mean squared error in @thm-shr-ridge is computed for a *fixed*
\( \bbeta \) and is therefore a function of that fixed value. Averaging the ridge
mean squared error over \( \bbeta\sim\Normal(\bzero,\sigma_u^2\I) \) reproduces the
BLUP formula, which is the empirical Bayes reading of @prp-shr-empirical-bayes.
:::

### C. Going deeper

::: {#exr-mix-blup-invariance}
[C1]

Show that the BLUP is invariant to reparameterization of the random effects: if
\( \Z^{*}=\Z\mathbf{T}^{-1} \) and \( \G^{*}=\mathbf{T}\G\mathbf{T}\T \) for a
nonsingular \( \mathbf{T} \), with \( \bu^{*}=\mathbf{T}\bu \), then
\( \V \) is unchanged and \( \hat{\bu}^{*}=\mathbf{T}\hat{\bu} \). What does this say
about predicting a contrast \( \mathbf{c}\T\bu \) directly versus predicting
\( \bu \) and then forming the contrast?
:::

::: {#exr-mix-blup-normal-optimal}
[C2]

Assume the normal mixed model and suppose \( \bbeta \) is known. Show that among
*all* predictors \( t(\Y) \), linear or not, the conditional mean
\( \E(\mathbf{c}\T\bu\mid\Y) \) minimizes \( \E(t-\mathbf{c}\T\bu)^2 \), and that it
coincides with the linear predictor of @thm-mix-blup(c). Explain why the same
argument does not show that \( \hat\tau \) is optimal among all predictors when
\( \bbeta \) is unknown.
:::
