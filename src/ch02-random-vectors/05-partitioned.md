# Moments of partitioned random vectors

Regression with random regressors studies a response vector \( \Y \) (\( p\times1 \)) together with a
regressor vector \( \X \) (\( q\times1 \)) observed on the same unit. In this section \( \X \) is a random
vector, not a model matrix. We stack the two into one vector and read off
what the blocks of its moments mean.

## Blocks of the mean and covariance

::: {#thm-rv-partitioned}
[Partitioned moments]

Let \( \W=(\Y\T,\X\T)\T \) have finite second moments, and write
\[
\E(\W)=\begin{pmatrix}\bmu_Y\\\bmu_X\end{pmatrix},\qquad
\Cov(\W)=\begin{pmatrix}\bSigma_{YY}&\bSigma_{YX}\\\bSigma_{XY}&\bSigma_{XX}\end{pmatrix}.
\]

::: {.enumerate options="label=(\alph*)"}
1. \( \bmu_Y=\E(\Y) \), \( \bmu_X=\E(\X) \), \( \bSigma_{YY}=\Cov(\Y) \), \( \bSigma_{XX}=\Cov(\X) \),
           \( \bSigma_{YX}=\Cov(\Y,\X) \) is \( p\times q \), and \( \bSigma_{XY}=\bSigma_{YX}\T \).

2. For constants \( \A \), \( \B \), \( \mathbf{c} \) of conformable sizes,
           \( \Cov(\A\Y,\B\X)=\A\bSigma_{YX}\B\T \), and
           \[
\Cov(\A\Y+\B\X+\mathbf{c})=\A\bSigma_{YY}\A\T+\A\bSigma_{YX}\B\T+\B\bSigma_{XY}\A\T+\B\bSigma_{XX}\B\T .
\]

3. \( \C(\bSigma_{XY})\subseteq\C(\bSigma_{XX}) \) and \( \C(\bSigma_{YX})\subseteq\C(\bSigma_{YY}) \).
           Equivalently, for every generalized inverse \( \bSigma_{XX}\ginv \),
           \[
\bSigma_{XX}\bSigma_{XX}\ginv\bSigma_{XY}=\bSigma_{XY},\qquad
          \bSigma_{YX}\bSigma_{XX}\ginv\bSigma_{XX}=\bSigma_{YX},
\]
           and similarly with the roles of \( \Y \) and \( \X \) exchanged.
:::

:::

::: {.proof}
(a) The outer product \( (\W-\E\W)(\W-\E\W)\T \) has blocks \( (\Y-\bmu_Y)(\Y-\bmu_Y)\T \),
\( (\Y-\bmu_Y)(\X-\bmu_X)\T \), its transpose, and \( (\X-\bmu_X)(\X-\bmu_X)\T \). Take expectations
blockwise. (b) \( \A\Y+\B\X+\mathbf{c}=[\A,\ \B]\,\W+\mathbf{c} \), and \( \A\Y=[\A,\ \bzero]\W \),
\( \B\X=[\bzero,\ \B]\W \). Apply @thm-rv-linear and multiply out the blocks.

(c) Let \( \mathbf{a}\in\Null(\bSigma_{XX}) \). By @thm-rv-cov-nnd(c), \( \mathbf{a}\T\X \) is constant with
probability one, so its covariance with \( \Y \) vanishes: \( \bSigma_{YX}\mathbf{a}=\Cov(\Y,\mathbf{a}\T\X)=\bzero \).
Thus \( \Null(\bSigma_{XX})\subseteq\Null(\bSigma_{YX}) \). Taking orthogonal complements, and using
\( \Null(\mathbf{B})\perpc=\C(\mathbf{B}\T) \) and the symmetry of \( \bSigma_{XX} \), gives
\( \C(\bSigma_{XY})\subseteq\C(\bSigma_{XX}) \). The same argument with \( \Y \) and \( \X \) exchanged gives
the second inclusion. A column inclusion \( \C(\bSigma_{XY})\subseteq\C(\bSigma_{XX}) \) means
\( \bSigma_{XY}=\bSigma_{XX}\mathbf{K}\T \) for some \( \mathbf{K} \). Then
\( \bSigma_{XX}\bSigma_{XX}\ginv\bSigma_{XY}=\bSigma_{XX}\bSigma_{XX}\ginv\bSigma_{XX}\mathbf{K}\T=\bSigma_{XY} \).
Conversely, that identity exhibits every column of \( \bSigma_{XY} \) as a combination of columns of
\( \bSigma_{XX} \) (@thm-mat-consistency). Transposing gives
\( \bSigma_{YX}=\mathbf{K}\bSigma_{XX} \) and the second identity.
:::

Part (c) looks technical, but it removes every difficulty with singular \( \bSigma_{XX} \) in what
follows. A direction in which \( \X \) does not vary cannot covary with anything.

## Linear prediction and residual covariance

The most important use of the blocks is to split \( \Y \) into a part explained linearly by \( \X \)
and a part uncorrelated with \( \X \).

::: {#thm-rv-blp}
[Best linear predictor]

In the setting of @thm-rv-partitioned, fix a generalized inverse \( \bSigma_{XX}\ginv \) and put
\[
\B_*=\bSigma_{YX}\bSigma_{XX}\ginv,\qquad
\mathbf{E}=\Y-\bmu_Y-\B_*(\X-\bmu_X),\qquad
\bSigma_{YY\cdot X}=\bSigma_{YY}-\bSigma_{YX}\bSigma_{XX}\ginv\bSigma_{XY}.
\]

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\mathbf{E})=\bzero \) and \( \Cov(\mathbf{E},\X)=\bzero \).

2. \( \Cov(\mathbf{E})=\bSigma_{YY\cdot X} \). This matrix is nonnegative definite and does not depend
           on the choice of generalized inverse. Neither does \( \B_*(\X-\bmu_X) \), with probability one.

3. For every constant \( p\times q \) matrix \( \B \) and \( p\times1 \) vector \( \mathbf{c} \),
           \[
\E\bigl[(\Y-\mathbf{c}-\B\X)(\Y-\mathbf{c}-\B\X)\T\bigr]
          =\bSigma_{YY\cdot X}+(\B-\B_*)\bSigma_{XX}(\B-\B_*)\T+\mathbf{d}\mathbf{d}\T,
\]
           where \( \mathbf{d}=\bmu_Y-\mathbf{c}-\B\bmu_X \). So the mean squared error matrix is never smaller,
           in the nonnegative definite ordering, than \( \bSigma_{YY\cdot X} \). Equality holds iff
           \( \mathbf{c}=\bmu_Y-\B\bmu_X \) and \( \B\bSigma_{XX}=\bSigma_{YX} \).
:::

:::

::: {.proof}
By @thm-rv-partitioned(c), write \( \bSigma_{YX}=\mathbf{K}\bSigma_{XX} \). Then
\( \B_*\bSigma_{XX}=\mathbf{K}\bSigma_{XX}\bSigma_{XX}\ginv\bSigma_{XX}=\mathbf{K}\bSigma_{XX}=\bSigma_{YX} \).

(a) The mean is clear. By @thm-rv-partitioned(b),
\( \Cov(\mathbf{E},\X)=\bSigma_{YX}-\B_*\bSigma_{XX}=\bzero \).

(b) Again by @thm-rv-partitioned(b),
\( \Cov(\mathbf{E})=\bSigma_{YY}-\B_*\bSigma_{XY}-\bSigma_{YX}\B_*\T+\B_*\bSigma_{XX}\B_*\T \). Using
\( \B_*\bSigma_{XX}=\bSigma_{YX} \), the last term is \( \bSigma_{YX}\B_*\T \), which cancels the third, so
\( \Cov(\mathbf{E})=\bSigma_{YY}-\B_*\bSigma_{XY}=\bSigma_{YY\cdot X} \). Next,
\( \bSigma_{YX}\bSigma_{XX}\ginv\bSigma_{XY}=\mathbf{K}\bSigma_{XX}\bSigma_{XX}\ginv\bSigma_{XX}\mathbf{K}\T=\mathbf{K}\bSigma_{XX}\mathbf{K}\T \),
which involves no generalized inverse. A covariance matrix is nonnegative definite. Finally,
by @thm-rv-cov-nnd(d), \( \X-\bmu_X=\bSigma_{XX}\mathbf{v} \) for some random \( \mathbf{v} \), with
probability one, and then \( \B_*(\X-\bmu_X)=\bSigma_{YX}\mathbf{v} \) for every choice of generalized inverse.

(c) Decompose \( \Y-\mathbf{c}-\B\X=\mathbf{E}+(\B_*-\B)(\X-\bmu_X)+\mathbf{d} \). The three pieces are,
respectively, mean zero, mean zero and constant, and the first two are uncorrelated by (a).
So the cross products have zero expectation, and the expected outer product is
\( \Cov(\mathbf{E})+(\B_*-\B)\bSigma_{XX}(\B_*-\B)\T+\mathbf{d}\mathbf{d}\T \). Both added matrices are
nonnegative definite. They vanish iff \( \mathbf{d}=\bzero \) and \( (\B-\B_*)\bSigma_{XX}^{1/2}=\bzero \), and the
latter is equivalent to \( (\B-\B_*)\bSigma_{XX}=\bzero \), that is, to \( \B\bSigma_{XX}=\bSigma_{YX} \).
:::

The affine function \( \bmu_Y+\B_*(\X-\bmu_X) \) is the **best linear predictor** of \( \Y \) from \( \X \).
Because the ordering in (c) is the nonnegative definite ordering, it is simultaneously best for
every linear combination \( \mathbf{a}\T\Y \) and for the total squared error. The matrix
\( \bSigma_{YY\cdot X} \) is the **partial covariance matrix** of \( \Y \) given \( \X \). When
\( \bSigma_{XX} \) is nonsingular, \( \B_*=\bSigma_{YX}\bSigma_{XX}^{-1} \) holds the population regression
coefficients, and \( \bSigma_{YY\cdot X} \) is a Schur complement. Two consequences follow from
[Chapter 1](../ch01-matrix-algebra/index.html). If \( \Cov(\W) \) is positive definite, the \( YY \) block of its inverse is
\( \bSigma_{YY\cdot X}^{-1} \) (@thm-mat-partitioned-inverse). And
\( \det\Cov(\W)=\det(\bSigma_{XX})\det(\bSigma_{YY\cdot X}) \) (@thm-mat-block-determinant),
so the generalized variance factorizes into the part of \( \X \) and the part of \( \Y \) that \( \X \) leaves
unexplained.

Nothing here says that \( \E(\Y\mid\X) \) is linear. The best linear predictor is best among
affine functions, whatever the joint distribution. For jointly normal vectors, the conditional
mean happens to be this affine function, and \( \bSigma_{YY\cdot X} \) is the conditional covariance
matrix (@thm-mvn-conditional). The partial correlation of @def-mvn-partial-correlation
is the correlation computed from \( \bSigma_{YY\cdot X} \). [Chapter 6](../ch06-projections/index.html) shows that least
squares estimates \( \B_* \) in large samples, whether or not the regression function is linear.
It also recovers the predictor as an orthogonal projection in the space of square-integrable
random variables. For scalar \( Y \), the coefficient vector used there is \( \B_*\T \).

::: {#prp-rv-multiple-correlation}
[Multiple correlation]

Let \( Y \) be scalar with \( \sigma_{YY}=\Var(Y)>0 \), let \( \bSigma_{XX} \) be positive definite, and write
\( \boldsymbol{\upsigma}_{XY}=\Cov(\X,Y) \). Over all \( \mathbf{a}\ne\bzero \),
\[
\max_{\mathbf{a}}\operatorname{corr}(Y,\mathbf{a}\T\X)^2
=\rho^2_{Y\cdot X}
=\frac{\boldsymbol{\upsigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY}}{\sigma_{YY}},
\]
attained at \( \mathbf{a}=\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY} \), the coefficients of the best linear predictor
(if \( \boldsymbol{\upsigma}_{XY}=\bzero \), every \( \mathbf{a} \) attains the maximum \( 0 \)).
Moreover \( 1-\rho^2_{Y\cdot X}=\sigma_{YY\cdot X}/\sigma_{YY} \).
:::

::: {.proof}
\( \operatorname{corr}(Y,\mathbf{a}\T\X)^2=(\mathbf{a}\T\boldsymbol{\upsigma}_{XY})^2/(\sigma_{YY}\,\mathbf{a}\T\bSigma_{XX}\mathbf{a}) \).
Write \( \mathbf{a}\T\boldsymbol{\upsigma}_{XY}=(\bSigma_{XX}^{1/2}\mathbf{a})\T(\bSigma_{XX}^{-1/2}\boldsymbol{\upsigma}_{XY}) \). The
Cauchy–Schwarz inequality gives
\( (\mathbf{a}\T\boldsymbol{\upsigma}_{XY})^2\le(\mathbf{a}\T\bSigma_{XX}\mathbf{a})(\boldsymbol{\upsigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY}) \),
with equality iff \( \bSigma_{XX}^{1/2}\mathbf{a} \) is proportional to \( \bSigma_{XX}^{-1/2}\boldsymbol{\upsigma}_{XY} \).
Dividing gives the bound and the maximizer. The last identity is the definition of
\( \sigma_{YY\cdot X}=\sigma_{YY}-\boldsymbol{\upsigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY} \), divided by \( \sigma_{YY} \).
:::

\( \rho^2_{Y\cdot X} \) is the population version of \( R^2 \): the fraction of the variance of \( Y \)
that the best linear predictor accounts for. [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) studies
its sample version.

## Conditioning: total expectation and total covariance

When the joint distribution is available, conditioning gives another decomposition of the
covariance matrix. Conditional expectations and covariances are taken entry by entry, as in
@def-rv-expectation, and \( \Cov(\Y\mid\X)=\E\bigl[(\Y-\E(\Y\mid\X))(\Y-\E(\Y\mid\X))\T\bigm|\X\bigr] \).

::: {#prp-rv-total-covariance}
[Laws of total expectation and covariance]

If \( \Y \) has finite second moments, then \( \E(\Y)=\E\bigl[\E(\Y\mid\X)\bigr] \) and
\[
\Cov(\Y)=\E\bigl[\Cov(\Y\mid\X)\bigr]+\Cov\bigl(\E(\Y\mid\X)\bigr).
\]
:::

::: {.proof}
The first identity is the tower property of conditional expectation, applied to each
component. For the second, let \( \mathbf{m}(\X)=\E(\Y\mid\X) \) and write
\( \Y-\bmu=(\Y-\mathbf{m}(\X))+(\mathbf{m}(\X)-\bmu) \). The outer product has two cross terms. The
conditional expectation of \( (\Y-\mathbf{m}(\X))(\mathbf{m}(\X)-\bmu)\T \) given \( \X \) is
\( \bigl[\E(\Y\mid\X)-\mathbf{m}(\X)\bigr](\mathbf{m}(\X)-\bmu)\T=\bzero \), so its expectation is zero, and
likewise for the other. The expectation of \( (\Y-\mathbf{m})(\Y-\mathbf{m})\T \) is
\( \E[\Cov(\Y\mid\X)] \) by the tower property, and that of \( (\mathbf{m}-\bmu)(\mathbf{m}-\bmu)\T \) is
\( \Cov(\mathbf{m}(\X)) \) because \( \E\mathbf{m}(\X)=\bmu \).
:::

::: {#exm-rv-random-intercept}
[A shared random effect]

Measurements \( Y_1,\dots,Y_k \) on the same subject share a subject effect: \( \Y=\alpha\bone+\be \),
where \( \alpha \) has variance \( \tau^2 \), \( \Cov(\be)=\sigma^2\I_k \), and \( \be \) is independent of \( \alpha \)
with mean zero. Given \( \alpha \), the mean is \( \alpha\bone \) and the covariance is \( \sigma^2\I \). By
@prp-rv-total-covariance,
\[
\Cov(\Y)=\sigma^2\I_k+\tau^2\bone\bone\T ,
\]
so every pair of measurements has correlation \( \tau^2/(\tau^2+\sigma^2) \). Independent errors
combined with a shared random component produce exactly this *equicorrelation* structure.
It is the starting point of Chapter 32. Ordinary least squares ignores it,
with consequences that Chapter 33 examines.
:::

## Exercises

### A. Check your understanding

::: {#exr-rv-blp-scalar}
[A1]

For scalar \( Y \) and \( X \) with \( \Var(X)>0 \), show that the best linear predictor has slope
\( \sigma_{YX}/\sigma_{XX} \), and that its mean squared error is \( \sigma_{YY}(1-\rho^2) \).
:::

::: {#exr-rv-nonlinear}
[A2]

Let \( X \) be uniform on \( (-1,1) \) and \( Y=X^2+\varepsilon \), with \( \varepsilon \) independent of \( X \), mean \( 0 \)
and variance \( \sigma^2 \). Find the best linear predictor of \( Y \) from \( X \) and its mean squared error,
and compare with the mean squared error of \( \E(Y\mid X) \).
:::

### B. Practice

::: {#exr-rv-partial-three}
[B1]

Let \( (Y_1,Y_2,Y_3) \) have correlation matrix \( \R \). Compute the partial covariance matrix of
\( (Y_2,Y_3) \) given \( Y_1 \) from @thm-rv-blp, and show that the corresponding correlation is
\[
\frac{\rho_{23}-\rho_{12}\rho_{13}}{\sqrt{(1-\rho_{12}^2)(1-\rho_{13}^2)}} .
\]
:::

::: {#exr-rv-blp-singular}
[B2]

Let \( \X=(X_1,X_1)\T \) with \( \Var(X_1)=1 \), and let \( \Cov(Y,X_1)=c \). Find two different generalized
inverses of \( \bSigma_{XX} \) that give different \( \B_* \), and show that both give the same predictor
\( \B_*(\X-\bmu_X) \).
:::

::: {#exr-rv-mixture}
[B3]

A random label \( K \) takes values \( 1,\dots,g \) with probabilities \( \pi_k \). Given \( K=k \), \( \Y \) has mean
\( \mathbf{m}_k \) and covariance \( \bSigma_k \). Show that
\( \Cov(\Y)=\sum_k\pi_k\bSigma_k+\sum_k\pi_k(\mathbf{m}_k-\bar{\mathbf{m}})(\mathbf{m}_k-\bar{\mathbf{m}})\T \), where
\( \bar{\mathbf{m}}=\sum_k\pi_k\mathbf{m}_k \). Interpret the two terms as within-group and between-group
covariance.
:::

### C. Going deeper

::: {#exr-rv-schur-nnd}
[C1]

Let \( \mathbf{N} \) be the symmetric partitioned matrix
\[
\mathbf{N}=\begin{pmatrix}\bSigma_{YY}&\bSigma_{YX}\\\bSigma_{XY}&\bSigma_{XX}\end{pmatrix}.
\]
Show that \( \mathbf{N} \) is nonnegative definite iff \( \bSigma_{XX} \) is nonnegative definite,
\( \C(\bSigma_{XY})\subseteq\C(\bSigma_{XX}) \), and \( \bSigma_{YY}-\bSigma_{YX}\bSigma_{XX}\ginv\bSigma_{XY} \) is
nonnegative definite for some generalized inverse \( \bSigma_{XX}\ginv \) (equivalently, given the
column-space condition, for every one). Show by a \( 2\times2 \) example that the column-space
condition cannot be dropped.
:::

::: {.solution}
If the matrix is nonnegative definite, it is the covariance matrix
of some \( \W=(\Y\T,\X\T)\T \) (@thm-rv-cov-nnd(b)). Then \( \bSigma_{XX} \) is a covariance matrix, the
column condition is @thm-rv-partitioned(c), and the Schur complement is \( \Cov(\mathbf{E}) \) by
@thm-rv-blp(b). Conversely, the column condition gives \( \bSigma_{YX}=\mathbf{K}\bSigma_{XX} \), and then
\( \bSigma_{YX}\bSigma_{XX}\ginv\bSigma_{XY}=\mathbf{K}\bSigma_{XX}\mathbf{K}\T \). With the nonsingular
\( \bT=\begin{psmallmatrix}\I&-\mathbf{K}\\\bzero&\I\end{psmallmatrix} \),
\[
\bT\begin{pmatrix}\bSigma_{YY}&\bSigma_{YX}\\\bSigma_{XY}&\bSigma_{XX}\end{pmatrix}\bT\T
=\begin{pmatrix}\bSigma_{YY}-\mathbf{K}\bSigma_{XX}\mathbf{K}\T&\bzero\\\bzero&\bSigma_{XX}\end{pmatrix},
\]
which is nonnegative definite. So is the original matrix, since
\( \mathbf{a}\T\mathbf{N}\mathbf{a}=\mathbf{b}\T(\bT\mathbf{N}\bT\T)\mathbf{b} \) with \( \mathbf{b}=(\bT\T)^{-1}\mathbf{a} \). For the example, take
\( \begin{psmallmatrix}1&1\\1&0\end{psmallmatrix} \). Here \( \bSigma_{XX}=0 \) is nonnegative definite, and with the
generalized inverse \( 0 \) the Schur complement is \( 1\ge0 \). But the determinant is \( -1 \), so the matrix is
not nonnegative definite. The column condition fails, because \( \C(1)\not\subseteq\C(0) \). Without
that condition the Schur complement also depends on the choice: every scalar \( t \) is a generalized
inverse of \( 0 \), and the Schur complement is then \( 1-1\cdot t\cdot1=1-t \), which is negative for
\( t>1 \). Under the column condition, \( \bSigma_{YX}\bSigma_{XX}\ginv\bSigma_{XY}=\mathbf{K}\bSigma_{XX}\mathbf{K}\T \) for
every choice, as shown above.
:::
