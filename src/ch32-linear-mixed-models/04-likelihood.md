# Maximum likelihood and REML

Everything so far has treated the variance components as known. They never are. The
moment estimates of [Section 32.1](01-random-effects.html) work in balanced layouts
and lose their footing in unbalanced ones, where no canonical set of quadratic forms
presents itself. Under normality the likelihood has no such difficulty, but it has one
defect, which this section diagnoses and repairs: maximum likelihood underestimates
variance components, because it does not pay for the fixed effects it has fitted.

## The likelihood and its profile

Under the normal linear mixed model, \( \Y\sim\Normal_n(\X\bbeta,\V) \) with
\( \V=\V(\boldsymbol{\uptheta}) \), so by @thm-mvn-density
\[
\ell(\bbeta,\boldsymbol{\uptheta})
=-\tfrac12\bigl\{n\log2\pi+\log\det\V
+(\y-\X\bbeta)\T\V^{-1}(\y-\X\bbeta)\bigr\}.
\]{#eq-mix-loglik}

For fixed \( \boldsymbol{\uptheta} \) only the last term involves \( \bbeta \), and
it is the generalized least squares criterion of @def-proj-gls, minimized by
\( \hbeta(\boldsymbol{\uptheta}) \). Writing
\[
\bP=\V^{-1}-\V^{-1}\X(\X\T\V^{-1}\X)\ginv\X\T\V^{-1},
\]{#eq-mix-P}

the profile log-likelihood of @def-opt-profile is
\[
\ell_P(\boldsymbol{\uptheta})
=-\tfrac12\bigl\{n\log2\pi+\log\det\V+\y\T\bP\y\bigr\}.
\]{#eq-mix-profile}

::: {#lem-mix-P}
[The residual operator P]

Let \( \V \) be positive definite and \( r=\rank(\X) \). The matrix \( \bP \) of
@eq-mix-P does not depend on the choice of generalized inverse, and

::: {.enumerate options="label=(\alph*)"}
1. \( \bP\X=\bzero \), \( \bP\V\bP=\bP \), \( \bP \) is symmetric and
   \( \rank(\bP)=n-r \);

2. \( \bP\y=\V^{-1}(\y-\X\hbeta) \) and
   \( \y\T\bP\y=\min_{\bb}(\y-\X\bb)\T\V^{-1}(\y-\X\bb) \);

3. if \( \V \) depends smoothly on \( \boldsymbol{\uptheta} \) and
   \( \V_j=\partial\V/\partial\theta_j \), then
   \( \partial\bP/\partial\theta_j=-\bP\V_j\bP \).
:::

:::

::: {.proof}
Write \( \X_*=\V^{-1/2}\X \) (@thm-mat-square-root). Then
\( \bP=\V^{-1/2}(\I-\M_*)\V^{-1/2} \), where
\( \M_*=\X_*(\X_*\T\X_*)\ginv\X_*\T \) is the orthogonal projection onto
\( \C(\X_*) \), which by @thm-proj-M-formula is the same for every generalized
inverse and has rank \( r \). Symmetry, \( \bP\X=\V^{-1/2}(\I-\M_*)\X_*=\bzero \),
\( \bP\V\bP=\V^{-1/2}(\I-\M_*)^2\V^{-1/2}=\bP \) and the rank statement follow at
once, which is (a). For (b), \( \bP\y=\V^{-1}\y-\V^{-1}\X\hbeta \) by the definition
of \( \hbeta \), and
\( \y\T\bP\y=\norm{(\I-\M_*)\V^{-1/2}\y}^2 \), the squared distance from
\( \V^{-1/2}\y \) to \( \C(\X_*) \), which is the stated minimum. Part (c) is proved
after @thm-mix-reml, which supplies a shorter route.
:::

## Restricted maximum likelihood

The trouble with @eq-mix-profile is visible in the simplest case. With
\( \V=\sigma^2\I \) it reduces to the profile of @thm-opt-mle, whose maximizer is
\( \hat\sigma^2=\text{SSE}/n \), biased downwards by the factor \( (n-p)/n \). The
likelihood has spent \( p \) dimensions fitting \( \bbeta \) and then measured the
spread of the data as if it had not. The repair is to compute the likelihood of a
part of the data that never saw \( \bbeta \).

::: {#def-mix-reml}
[Error contrasts and restricted maximum likelihood]

Let \( r=\rank(\X) \). An **error contrast** is a random variable
\( \mathbf{k}\T\Y \) with \( \X\T\mathbf{k}=\bzero \); its distribution does not
involve \( \bbeta \). A **maximal set** of error contrasts is
\( \mathbf{K}\T\Y \) for an \( n\times(n-r) \) matrix \( \mathbf{K} \) of full column
rank with \( \C(\mathbf{K})=\C(\X)\perpc \). The **restricted** (or residual)
**log-likelihood** is the log-likelihood of \( \mathbf{K}\T\Y \),
\[
\begin{aligned}
\ell_R(\boldsymbol{\uptheta})=-\tfrac12\bigl\{
&(n-r)\log2\pi+\log\det(\mathbf{K}\T\V\mathbf{K})\\
&+\y\T\mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T\y\bigr\},
\end{aligned}
\]{#eq-mix-reml-loglik}

and the **REML estimate** of \( \boldsymbol{\uptheta} \) maximizes it over the
parameter space.
:::

::: {#thm-mix-reml}
[The restricted likelihood]

In @def-mix-reml, with \( \V \) positive definite:

::: {.enumerate options="label=(\alph*)"}
1. \( \mathbf{K}\T\Y\sim\Normal_{n-r}(\bzero,\mathbf{K}\T\V\mathbf{K}) \) with
   \( \mathbf{K}\T\V\mathbf{K} \) positive definite, so @eq-mix-reml-loglik is its
   log-likelihood.

2. \( \mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T=\bP \) for every
   admissible \( \mathbf{K} \). In particular the quadratic form in
   @eq-mix-reml-loglik equals \( \y\T\bP\y \) and does not depend on
   \( \mathbf{K} \).

3. Two admissible \( \mathbf{K}_1,\mathbf{K}_2 \) give restricted log-likelihoods
   differing by a constant free of \( \boldsymbol{\uptheta} \). Hence the REML
   estimate does not depend on which maximal set of error contrasts is used.

4. If \( \rank(\X)=p \), then for any admissible \( \mathbf{K} \) with
   \( \mathbf{K}\T\mathbf{K}=\I \),
   \[
   \det(\mathbf{K}\T\V\mathbf{K})
   =\frac{\det\V\,\det(\X\T\V^{-1}\X)}{\det(\X\T\X)},
   \]{#eq-mix-reml-determinant}

   so that
   \[
   \begin{aligned}
   \ell_R(\boldsymbol{\uptheta})=\ell_P(\boldsymbol{\uptheta})
   &+\tfrac{p}{2}\log2\pi-\tfrac12\log\det(\X\T\V^{-1}\X)\\
   &+\tfrac12\log\det(\X\T\X).
   \end{aligned}
   \]{#eq-mix-reml-vs-ml}
:::

:::

::: {.proof}
(a) By @thm-mvn-linear, \( \mathbf{K}\T\Y \) is normal with mean
\( \mathbf{K}\T\X\bbeta=\bzero \) and covariance
\( \mathbf{K}\T\V\mathbf{K} \), which is positive definite because
\( \mathbf{K} \) has full column rank and \( \V \) is positive definite. Its density
is @thm-mvn-density.

(b) Fix \( \y \) and consider
\[
\min\{\bz\T\V^{-1}\bz:\ \mathbf{K}\T\bz=\mathbf{K}\T\y\}.
\]
The feasible set is \( \y+\C(\X) \): indeed \( \mathbf{K}\T(\bz-\y)=\bzero \) says
\( \bz-\y\in\C(\mathbf{K})\perpc=\C(\X) \). So the minimum equals
\( \min_{\bb}(\y-\X\bb)\T\V^{-1}(\y-\X\bb)=\y\T\bP\y \) by @lem-mix-P(b). On the
other hand the constrained problem can be solved directly. The criterion is strictly
convex, so by @prp-mat-lagrange the minimizer satisfies
\( \V^{-1}\bz=\mathbf{K}\boldsymbol{\uptheta} \), that is
\( \bz=\V\mathbf{K}\boldsymbol{\uptheta} \); the constraint gives
\( \mathbf{K}\T\V\mathbf{K}\boldsymbol{\uptheta}=\mathbf{K}\T\y \), and the minimum
value is
\( \boldsymbol{\uptheta}\T\mathbf{K}\T\V\mathbf{K}\boldsymbol{\uptheta}
=\y\T\mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T\y \). The two expressions
agree for every \( \y \), and both are quadratic forms with symmetric matrices, so
the matrices coincide.

(c) Since \( \C(\mathbf{K}_1)=\C(\mathbf{K}_2) \) and both have \( n-r \) columns of
full rank, \( \mathbf{K}_2=\mathbf{K}_1\A \) for a nonsingular \( (n-r)\times(n-r) \)
matrix \( \A \). Then
\( \det(\mathbf{K}_2\T\V\mathbf{K}_2)=(\det\A)^2\det(\mathbf{K}_1\T\V\mathbf{K}_1) \),
while the quadratic forms are equal by (b). So the two log-likelihoods differ by
\( -\log\lvert\det\A\rvert \), which does not involve \( \boldsymbol{\uptheta} \).

(d) Let \( \mathbf{T}=[\X,\mathbf{K}] \), which is nonsingular. Because
\( \X\T\mathbf{K}=\bzero \) and \( \mathbf{K}\T\mathbf{K}=\I \),
\[
\det(\mathbf{T}\T\mathbf{T})
=\det\begin{pmatrix}\X\T\X&\bzero\\\bzero&\I\end{pmatrix}=\det(\X\T\X),
\]
so \( (\det\mathbf{T})^2=\det(\X\T\X) \) and
\( \det(\mathbf{T}\T\V\mathbf{T})=\det(\X\T\X)\det\V \). Next,
\( \mathbf{T}^{-1} \) has rows \( (\X\T\X)^{-1}\X\T \) stacked on
\( \mathbf{K}\T \), as one checks by multiplying. Hence the leading \( p\times p \)
block of \( (\mathbf{T}\T\V\mathbf{T})^{-1}=\mathbf{T}^{-1}\V^{-1}
\mathbf{T}^{-\top} \) is \( (\X\T\X)^{-1}\X\T\V^{-1}\X(\X\T\X)^{-1} \). By
@thm-mat-partitioned-inverse that block is the inverse of the Schur complement of the
trailing block \( \mathbf{K}\T\V\mathbf{K} \), so by
@thm-mat-block-determinant(b),
\[
\det(\mathbf{T}\T\V\mathbf{T})
=\det(\mathbf{K}\T\V\mathbf{K})\cdot
\frac{\det(\X\T\X)^2}{\det(\X\T\V^{-1}\X)} .
\]
Equating the two expressions for \( \det(\mathbf{T}\T\V\mathbf{T}) \) gives
@eq-mix-reml-determinant, and substituting it together with (b) into
@eq-mix-reml-loglik gives @eq-mix-reml-vs-ml.
:::

The proof of @lem-mix-P(c) is now immediate: differentiating
\( \bP=\mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T \) with \( \mathbf{K} \)
fixed gives
\( \partial\bP/\partial\theta_j
=-\mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T\V_j
\mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T=-\bP\V_j\bP \).

::: {.idea}
REML fits \( \boldsymbol{\uptheta} \) to the residual directions only. It is the
general form of the familiar divisor \( n-p \): with \( \V=\sigma^2\I \),
\( \rank(\X)=p \) and an orthonormal \( \mathbf{K} \) (that is
\( \mathbf{K}\T\mathbf{K}=\I \), as in @thm-mix-reml(d)), @eq-mix-reml-loglik is the
likelihood of \( \mathbf{K}\T\Y\sim\Normal_{n-p}(\bzero,\sigma^2\I) \), maximized at
\( \hat\sigma^2=\norm{\mathbf{K}\T\y}^2/(n-p)=\text{SSE}/(n-p)=s^2 \); by
@thm-mix-reml(c) any other admissible \( \mathbf{K} \) gives the same answer. Ordinary
maximum likelihood gives \( \text{SSE}/n \). So \( s^2 \), the estimator every
earlier chapter used, was the REML estimator all along. Its bias in a genuine mixed
model is computed exactly in @prp-mix-oneway-reml.
:::

::: {.warning}
Since \( \mathbf{K} \) is determined by \( \C(\X) \), two models with *different
fixed effects* have restricted likelihoods of different data, and comparing them is
meaningless; [Section 32.5](05-inference.html) makes the rule explicit. Worse,
@eq-mix-reml-vs-ml shows that the constant
\( \tfrac12\log\det(\X\T\X) \) is part of the definition, and software often omits
it; then the reported restricted log-likelihood changes if a column of \( \X \) is
measured in thousands rather than units, by
\( \log\lvert\det\mathbf{S}\rvert \) for the rescaling matrix \( \mathbf{S} \).
Comparisons of covariance structures with the *same* \( \X \) are unaffected, since
the constant cancels.
:::

## Why maximum likelihood is biased, exactly

In the balanced one-way model both likelihoods can be maximized in closed form, and
the bias can be computed rather than asserted.

::: {#prp-mix-oneway-reml}
[Balanced one-way: ML and REML]

In the balanced normal one-way random effects model, write
\( \lambda_1=\sigma^2+m\sigma_a^2 \) and \( \lambda_2=\sigma^2 \), and treat
\( (\lambda_1,\lambda_2) \) with \( \lambda_1\ge\lambda_2>0 \) as the parameters.
Then, up to constants free of the parameters,
\[
-2\ell_P=g\log\lambda_1+(n-g)\log\lambda_2
+\frac{\text{SSB}}{\lambda_1}+\frac{\text{SSE}}{\lambda_2},
\]{#eq-mix-oneway-ml}

\[
-2\ell_R=(g-1)\log\lambda_1+(n-g)\log\lambda_2
+\frac{\text{SSB}}{\lambda_1}+\frac{\text{SSE}}{\lambda_2}.
\]{#eq-mix-oneway-reml}

The two criteria differ in one coefficient, but leave the interior of the parameter
space at different moments.

::: {.enumerate options="label=(\alph*)"}
1. When \( \text{MSB}\ge\text{MSE} \), REML gives \( \hat\lambda_1=\text{MSB} \),
   \( \hat\lambda_2=\text{MSE} \), so \( \hat\sigma^2=\text{MSE} \) and
   \( \hat\sigma_a^2=(\text{MSB}-\text{MSE})/m \): the moment estimates of
   @prp-mix-oneway(c), which are unbiased.

2. When \( \text{MSB}\ge g\,\text{MSE}/(g-1) \), maximum likelihood gives
   \( \hat\lambda_1=\text{SSB}/g \), \( \hat\lambda_2=\text{MSE} \), so
   \( \hat\sigma^2_{\text{ML}}=\text{MSE} \) and
   \[
   \hat\sigma^2_{a,\text{ML}}
   =\frac{1}{m}\Bigl(\frac{g-1}{g}\text{MSB}-\text{MSE}\Bigr).
   \]

   The right-hand side is the *unconstrained* maximizer whatever the data, it is
   nonnegative exactly under the stated condition, and
   \[
   \E\Bigl\{\frac{1}{m}\Bigl(\frac{g-1}{g}\text{MSB}-\text{MSE}\Bigr)\Bigr\}
   -\sigma_a^2=-\frac{\lambda_1}{gm}.
   \]{#eq-mix-ml-bias}

3. Otherwise the maximizer lies on the boundary \( \lambda_1=\lambda_2 \) and the
   method returns \( \hat\sigma_a^2=0 \) with \( \hat\sigma^2 \) the pooled mean
   square: for REML when \( \text{MSB}<\text{MSE} \), for maximum likelihood on the
   wider range \( \text{MSB}<g\,\text{MSE}/(g-1) \). Between the two thresholds the
   REML estimate of \( \sigma_a^2 \) is positive and the ML estimate is zero.
:::

:::

::: {.proof}
By @eq-mix-balanced-cov, \( \V=\lambda_1\M+\lambda_2(\I-\M) \), so
\( \det\V=\lambda_1^g\lambda_2^{n-g} \) and
\( \V^{-1}=\lambda_1^{-1}\M+\lambda_2^{-1}(\I-\M) \). Since \( \bone\in\C(\Z) \),
\( \V^{-1}\bone=\lambda_1^{-1}\bone \), so the generalized least squares estimate of
\( \mu \) is \( \bar y \) and the residual is \( \y-\bar y\bone \). Splitting it by
@thm-proj-nested into its parts in \( \C(\Z)\cap\bone\perpc \) and
\( \C(\Z)\perpc \),
\[
\y\T\bP\y=(\y-\bar y\bone)\T\V^{-1}(\y-\bar y\bone)
=\frac{\text{SSB}}{\lambda_1}+\frac{\text{SSE}}{\lambda_2}.
\]
Substituting into @eq-mix-profile gives @eq-mix-oneway-ml. For the restricted
version, \( \X\T\V^{-1}\X=n/\lambda_1 \) and \( \X\T\X=n \), so the correction in
@eq-mix-reml-vs-ml is \( -\tfrac12\log(n/\lambda_1)+\tfrac12\log n
=\tfrac12\log\lambda_1 \), which converts the coefficient \( g \) of
\( \log\lambda_1 \) into \( g-1 \) after multiplying by \( -2 \). That
is @eq-mix-oneway-reml.

Both right-hand sides have the form \( c_1\log\lambda_1+S_1/\lambda_1 \) plus
\( c_2\log\lambda_2+S_2/\lambda_2 \), a sum of two functions of separate arguments,
each strictly convex in \( 1/\lambda \) and minimized at \( \lambda=S/c \). So
without the constraint \( \lambda_1\ge\lambda_2 \), REML gives
\( \hat\lambda_1=\text{SSB}/(g-1)=\text{MSB} \) and
\( \hat\lambda_2=\text{SSE}/(n-g)=\text{MSE} \), while maximum likelihood gives
\( \hat\lambda_1=\text{SSB}/g \) and the same \( \hat\lambda_2 \). These are feasible
exactly when \( \hat\lambda_1\ge\hat\lambda_2 \), that is when
\( \text{MSB}\ge\text{MSE} \) for REML and when
\( \text{SSB}/g=(g-1)\text{MSB}/g\ge\text{MSE} \), i.e.
\( \text{MSB}\ge g\,\text{MSE}/(g-1) \), for maximum likelihood. Translating back
through \( \sigma^2=\lambda_2 \) and \( \sigma_a^2=(\lambda_1-\lambda_2)/m \) gives
(a) and (b), except for @eq-mix-ml-bias. For that, use
\( \E\,\text{MSB}=\lambda_1 \) and \( \E\,\text{MSE}=\sigma^2 \):
\[
\E\Bigl\{\frac{1}{m}\Bigl(\frac{g-1}{g}\text{MSB}-\text{MSE}\Bigr)\Bigr\}
=\frac{1}{m}\Bigl(\frac{g-1}{g}\lambda_1-\sigma^2\Bigr)
=\sigma_a^2-\frac{\lambda_1}{gm}.
\]
For (c), suppose the unconstrained minimizer of one of the two objectives violates
\( \lambda_1\ge\lambda_2 \), which by the previous paragraph happens exactly on the
stated range. Each objective is a sum of two terms, one
in each argument, each strictly unimodal on \( (0,\infty) \) and tending to
\( +\infty \) at both ends, so a constrained minimum exists; were it interior to
\( \{\lambda_1>\lambda_2\} \) the gradient would vanish there, making it the
infeasible unconstrained minimizer. Hence it lies on the boundary
\( \lambda_1=\lambda_2=\lambda \), where the objective becomes
\( c\log\lambda+(\text{SSB}+\text{SSE})/\lambda \) with \( c=n \) for ML and
\( c=n-1 \) for REML, minimized at the corresponding pooled mean square.
:::

The bias \( -\lambda_1/(gm) \) is \( \lambda_1 \) divided by the total sample size:
the price of the single degree of freedom used up by the grand mean, charged entirely
to the between-group component. In the proficiency study,
\( \lambda_1=\sigma^2+m\sigma_a^2=0.25^2+4(0.20)^2=0.2225 \) and
\( gm=48 \), so the bias is
\( -0.0046 \), which is
\( 11.6 \) per cent of \( \sigma_a^2 \). A simulation of
\( 40\,000 \) studies reproduces it to the last recorded digit. Both the formula
@eq-mix-ml-bias and that simulation describe the *interior* expression of
@prp-mix-oneway-reml(b); truncating it at zero, as the estimator itself does, can only
raise the estimate, so the estimator's own bias is slightly smaller in modulus. Here
the expression is negative in \( 2.0 \) per cent of the simulated
studies and the truncated bias is \( -0.0045 \).
With \( p \) fixed effects instead of one the same calculation gives a bias of order
\( p/(gm) \), which is why REML matters most when the fixed-effects model is large
relative to the data.

## Score equations and algorithms

::: {#prp-mix-scores}
[Score equations]

Let \( \V_j=\partial\V/\partial\theta_j \) and \( \bP \) be as in @eq-mix-P. Then
\[
\begin{aligned}
\frac{\partial\ell_P}{\partial\theta_j}
&=-\tfrac12\bigl\{\tr(\V^{-1}\V_j)-\y\T\bP\V_j\bP\y\bigr\},\\
\frac{\partial\ell_R}{\partial\theta_j}
&=-\tfrac12\bigl\{\tr(\bP\V_j)-\y\T\bP\V_j\bP\y\bigr\},
\end{aligned}
\]{#eq-mix-scores}

The restricted score exceeds the profile score by
\( \tfrac12\tr\{(\X\T\V^{-1}\X)\ginv\X\T\V^{-1}\V_j\V^{-1}\X\} \), which is
nonnegative whenever \( \V_j \) is nonnegative definite, as it is for a
variance-component parameterization \( \V=\sum_j\theta_j\Z_j\Z_j\T \).
:::

::: {.proof}
For the profile, \( \partial\log\det\V/\partial\theta_j=\tr(\V^{-1}\V_j) \) and
\( \partial(\y\T\bP\y)/\partial\theta_j=-\y\T\bP\V_j\bP\y \) by @lem-mix-P(c). For
the restricted version, write \( \ell_R \) through @eq-mix-reml-loglik with a fixed
\( \mathbf{K} \): then
\( \partial\log\det(\mathbf{K}\T\V\mathbf{K})/\partial\theta_j
=\tr\{(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T\V_j\mathbf{K}\}=\tr(\bP\V_j) \) by
@thm-mix-reml(b) and the cyclic property of the trace (@thm-mat-trace-cyclic), and the
quadratic term differentiates as before. Subtracting, and using @eq-mix-P for
\( \bP \), gives the stated
difference. Its value is the same for every generalized inverse, because
\( \C(\X\T\V^{-1}\V_j\V^{-1}\X)\subset\C(\X\T)=\C(\X\T\V^{-1}\X) \), so take the
Moore–Penrose inverse, which is nonnegative definite; when \( \V_j \) is nonnegative
definite too, the difference is a trace of a product of two nonnegative definite
matrices, hence nonnegative.
:::

Setting @eq-mix-scores to zero gives equations that are not solvable in closed form
outside balanced designs, and three families of algorithms are used.

**Fisher scoring.** The expected information for the restricted likelihood has
entries \( \tfrac12\tr(\bP\V_j\bP\V_k) \), so a Newton-type step costs one trace per
pair of components. Convergence is fast near the maximum and unreliable far from it,
and nothing keeps the iterates inside the parameter space.

**Average information REML.** Gilmour, Thompson and Cullis (1995) noticed that the
average of the observed and expected information has entries
\( \tfrac12\y\T\bP\V_j\bP\V_k\bP\y \), which need no traces at all. Average
information REML is the workhorse of large animal-breeding and genomic analyses.

**EM.** Treat \( \bu \) as missing data (Dempster, Laird and Rubin, 1977). With
\( \R=\sigma^2\I \) and \( \G=\sigma_u^2\I_q \), the complete-data maximum likelihood
estimates are averages of squares, so the E-step and M-step give the maximum
likelihood update
\[
\sigma_u^{2,\text{new}}=\frac{1}{q}\bigl\{\norm{\hat{\bu}}^2+\tr\mathbf{S}\bigr\},
\]
\[
\sigma^{2,\text{new}}=\frac{1}{n}\bigl\{\norm{\y-\X\hbeta-\Z\hat{\bu}}^2
+\tr(\Z\mathbf{S}\Z\T)\bigr\},
\]
because \( \E(\bu\T\bu\mid\y)=\norm{\hat{\bu}}^2+\tr\mathbf{S} \) and
\( \E(\be\T\be\mid\y)=\norm{\y-\X\hbeta-\Z\hat{\bu}}^2+\tr(\Z\mathbf{S}\Z\T) \),
where
\[
\mathbf{S}=\Cov(\bu\mid\y)=\bigl(\G^{-1}+\Z\T\R^{-1}\Z\bigr)^{-1}
=\G-\G\Z\T\V^{-1}\Z\G
\]
is the inverse of the trailing block of \( \mathbf{C} \) in @eq-mix-henderson, the
matrix the push-through identity @eq-mix-push-through was applied to. It is *not* the
trailing block of \( \mathbf{C}^{-1} \) in @thm-mix-henderson(c), which carries a
third term because \( \bbeta \) has been estimated; that matrix gives an iteration
maximizing neither likelihood. Each step is cheap, increases
the likelihood and never leaves the parameter space; the price is slow convergence. A
restricted version runs the same argument on the error contrasts. Laird
and Ware (1982) introduced this algorithm to longitudinal data analysis, and it
remains the most reliable starting procedure.

::: {.warning}
One starting value is not enough. Fit the Grunfeld panel of @exm-mix-grunfeld with a
random slope as well as a random intercept, on the logarithmic scale for the variances
that keeps the iterates positive: a search begun with a negligible slope variance
crawls along the boundary and stops at the random-intercept fit, with restricted
log-likelihood \( -1176.53 \) against the maximum
\( -1159.01 \), a gap of \( 17.52 \)
reported as successful convergence. Here thirty random starts all reach the same
maximum, which is some evidence that it is the only one; but such likelihoods can
also be genuinely multimodal,
increasingly so as the covariance structure grows (Searle, Casella and McCulloch,
1992).
:::

## The proficiency study, computed

::: {#exm-mix-reml-lab}
[REML and ML for the twelve laboratories]

For the data of @exm-mix-proficiency, @prp-mix-oneway-reml gives
\( \hat\sigma^2=0.0620 \) by both methods, and
\[
\hat\sigma_{a,\text{REML}}^2=0.0386,\qquad
\hat\sigma_{a,\text{ML}}^2=0.0341,
\]
the maximum likelihood estimate being \( 0.883 \) of the
restricted one, because it replaces \( \text{MSB} \) by
\( (g-1)\text{MSB}/g \) before subtracting \( \text{MSE} \). The
restricted log-likelihood at its maximum is
\( -8.228 \) and the ordinary log-likelihood is
\( -8.360 \); the two are not comparable, being
likelihoods of \( 47 \) and \( 48 \) numbers
respectively. Both agree with the general matrix formulas coded from
@eq-mix-profile and @eq-mix-reml-loglik, and with statsmodels' `MixedLM`.

[Figure 32.4.1](04-likelihood.html#fig-mix-reml) plots the two profiles against
\( \sigma_a^2 \), each with \( \sigma^2 \) maximized out: the same shape in different
places, the restricted one flatter on the right, which is why REML intervals are
wider. Cutting the restricted profile \( 1.92 \) units below its maximum gives
\( \sigma_a^2\in(0.0086,\,0.1275) \), an interval comparable
to the exact one of @cor-mix-icc-interval and available in models where no exact
interval exists.
:::

::: {when-format="html"}
![**Figure 32.4.1.** Profile log-likelihoods of \( \sigma_a^2 \) for the proficiency
study, with \( \sigma^2 \) maximized out and each curve shifted to have maximum zero.
The horizontal line is \( 1.92 \) below the maximum, the cut that gives an
approximate 95% interval.](reml_profile.svg){#fig-mix-reml width=68%}
:::

::: {when-format="pdf"}
![Profile log-likelihoods of \( \sigma_a^2 \) for the proficiency study, with
\( \sigma^2 \) maximized out and each curve shifted to have maximum zero. The
horizontal line is \( 1.92 \) below the maximum, the cut that gives an approximate 95%
interval.](reml_profile.pdf){width=68%}
:::

```{.python .run #cell-reml-likelihood}
import numpy as np

g, m = 12, 4
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

X = np.ones((n, 1))
Z = np.repeat(np.eye(g), m, axis=0)

def loglik(s2, s2a, reml=True):
    """Log-likelihood (reml=False) or restricted log-likelihood of the variances."""
    V = s2a * (Z @ Z.T) + s2 * np.eye(n)
    L = np.linalg.cholesky(V)
    Vi_X = np.linalg.solve(V, X)
    W = X.T @ Vi_X                              # X' V^{-1} X
    r = y - X @ np.linalg.solve(W, Vi_X.T @ y)  # GLS residual
    quad = r @ np.linalg.solve(V, r)            # = y' P y
    logdetV = 2 * np.sum(np.log(np.diag(L)))
    out = -0.5 * (n * np.log(2 * np.pi) + logdetV + quad)
    if reml:
        out += -0.5 * (np.linalg.slogdet(W)[1] - np.linalg.slogdet(X.T @ X)[1]
                       - X.shape[1] * np.log(2 * np.pi))
    return out
```

```{.python .run #cell-reml-closed}
Y = y.reshape(g, m)
lab_mean = Y.mean(axis=1)
ms_between = m * np.sum((lab_mean - Y.mean()) ** 2) / (g - 1)
ms_within = np.sum((Y - lab_mean[:, None]) ** 2) / (n - g)

reml_s2, reml_s2a = ms_within, (ms_between - ms_within) / m
ml_s2 = ms_within
ml_s2a = ((g - 1) / g * ms_between - ms_within) / m

print(f"REML  sigma^2 {reml_s2:.4f}  sigma_a^2 {reml_s2a:.4f}")
print(f"ML    sigma^2 {ml_s2:.4f}  sigma_a^2 {ml_s2a:.4f}")
print("restricted log-likelihood", round(loglik(reml_s2, reml_s2a), 4))
print("log-likelihood           ", round(loglik(ml_s2, ml_s2a, reml=False), 4))
```

Two hundred EM steps from a deliberately poor start reach those maximum likelihood
values to nine decimal places, which is what a correct \( \mathbf{S} \) buys; the
trailing block of \( \mathbf{C}^{-1} \) in its place converges instead to
\( (0.0644,\,0.0380) \), maximizing neither likelihood.

```{.python .run #cell-reml-em}
ZZt = Z @ Z.T
s2, s2a = 0.1, 0.1                              # a deliberately poor starting value
for _ in range(200):
    V = s2a * ZZt + s2 * np.eye(n)
    Vi = np.linalg.inv(V)
    W = X.T @ Vi @ X
    beta_hat = np.linalg.solve(W, X.T @ Vi @ y)
    u_hat = s2a * Z.T @ Vi @ (y - X @ beta_hat)
    S = np.linalg.inv(np.eye(g) / s2a + Z.T @ Z / s2)   # Cov(u | y), not (C^{-1})_22
    resid = y - X @ beta_hat - Z @ u_hat
    s2a = (u_hat @ u_hat + np.trace(S)) / g
    s2 = (resid @ resid + np.trace(Z @ S @ Z.T)) / n

print(f"EM    sigma^2 {s2:.4f}  sigma_a^2 {s2a:.4f}")
```

## Exercises

### A. Check your understanding

::: {#exr-mix-reml-ordinary}
[A1]

Take \( \V=\sigma^2\I \) and \( \rank(\X)=p \). Show from @eq-mix-reml-vs-ml that
\( -2\ell_R=(n-p)\log\sigma^2+\text{SSE}/\sigma^2 \) up to constants, and deduce that
the REML estimate of \( \sigma^2 \) is \( s^2 \). Then use @prp-opt-divisor to find
which multiple of SSE has the smallest mean squared error, and show that the maximum
likelihood estimator beats \( s^2 \) for some \( (n,p) \) and loses to it for others.
:::

::: {.solution}
With \( \V=\sigma^2\I \): \( \log\det\V=n\log\sigma^2 \),
\( \y\T\bP\y=\text{SSE}/\sigma^2 \) and
\( \det(\X\T\V^{-1}\X)=\sigma^{-2p}\det(\X\T\X) \), so
\( \log\det(\X\T\V^{-1}\X)-\log\det(\X\T\X)=-p\log\sigma^2 \); adding this to
\( -2\ell_P=n\log\sigma^2+\text{SSE}/\sigma^2 \) leaves
\( (n-p)\log\sigma^2+\text{SSE}/\sigma^2 \). Its minimizer is
\( \text{SSE}/(n-p)=s^2 \). By @prp-opt-divisor the best multiple of SSE is
\( \text{SSE}/(n-p+2) \), which neither estimator in the question attains; writing
\( k=n-p \), the same proposition gives relative mean squared errors \( 2/k \) for
\( s^2 \) and \( (2k+p^2)/n^2 \) for \( \text{SSE}/n \). Multiplying out
\( n^2=(k+p)^2 \), the inequality \( (2k+p^2)/n^2\le2/k \) reduces to
\( k(p-4)\le2p \). So maximum likelihood wins for every \( p\le4 \), whatever
\( n \), and loses once \( p \) is moderately large: with \( n=100 \) and
\( p=10 \) it gives \( (180+100)/100^2=0.0280 \) against
\( 2/90=0.0222 \). Fitting many fixed effects and then dividing by
\( n \) is the bad case, which is the case REML exists for.
:::

::: {#exr-mix-reml-scale}
[A2]

A colleague reports two restricted log-likelihoods for the same model fitted to the
same data, differing by \( 6.91 \), the only difference between the runs being that
one measured a regressor in units and the other in thousands. Explain the discrepancy
using @eq-mix-reml-vs-ml, and confirm the value \( 6.91 \).
:::

### B. Practice

::: {#exr-mix-reml-contrasts}
[B1]

Verify @thm-mix-reml(b) directly in the simplest case: \( n=3 \),
\( \X=\bone_3 \), \( \V=\diag(v_1,v_2,v_3) \). Choose a basis \( \mathbf{K} \) of
\( \bone\perpc \), compute
\( \mathbf{K}(\mathbf{K}\T\V\mathbf{K})^{-1}\mathbf{K}\T \) and \( \bP \), and check
they agree. Then change \( \mathbf{K} \) to a different basis and confirm that
\( \det(\mathbf{K}\T\V\mathbf{K}) \) changes by a factor free of \( \V \).
:::

::: {#exr-mix-reml-two-groups}
[B2]

Two groups of sizes \( n_1,n_2 \) have independent normal responses with a common
mean and variances \( \sigma^2 \) and \( c\sigma^2 \) with \( c \) known. Write down
\( \ell_P \) and \( \ell_R \) for \( \sigma^2 \), find both maximizers, and confirm
that REML divides by \( n_1+n_2-1 \) and ML by \( n_1+n_2 \).
:::

::: {#exr-mix-reml-score-check}
[B3]

Specialize the score equations @eq-mix-scores of @prp-mix-scores to the balanced
one-way model with
\( \boldsymbol{\uptheta}=(\lambda_1,\lambda_2) \), for which
\( \V=\lambda_1\M+\lambda_2(\I-\M) \), \( \V_1=\M \) and \( \V_2=\I-\M \). Show that
the restricted score equations are \( \hat\lambda_1=\text{MSB} \) and
\( \hat\lambda_2=\text{MSE} \), recovering @prp-mix-oneway-reml(a) without any
optimization.
:::

::: {.solution}
Here \( \bP=\lambda_1^{-1}(\M-\bP_1)+\lambda_2^{-1}(\I-\M) \), so
\( \bP\V_1=\lambda_1^{-1}(\M-\bP_1) \) with trace \( (g-1)/\lambda_1 \), and
\( \y\T\bP\V_1\bP\y=\text{SSB}/\lambda_1^2 \). The first equation in
@eq-mix-scores therefore reads \( (g-1)/\lambda_1=\text{SSB}/\lambda_1^2 \), that is
\( \lambda_1=\text{SSB}/(g-1) \). Likewise \( \tr(\bP\V_2)=(n-g)/\lambda_2 \) and
\( \y\T\bP\V_2\bP\y=\text{SSE}/\lambda_2^2 \), giving
\( \lambda_2=\text{SSE}/(n-g) \). Using \( \V^{-1} \) instead of \( \bP \) replaces
\( g-1 \) by \( \tr(\V^{-1}\V_1)\lambda_1=g \), which is the ML answer.
:::

### C. Going deeper

::: {#exr-mix-reml-bayes-preview}
[C1]

Show that
\[
\int\exp\{\ell(\bbeta,\boldsymbol{\uptheta})\}\,d\bbeta
=\exp\{\ell_R(\boldsymbol{\uptheta})\}
\]
when \( \rank(\X)=p \) and the constant
\( \tfrac12\log\det(\X\T\X) \) is dropped from @eq-mix-reml-loglik. (Hint: complete
the square in \( \bbeta \) and use the normal integral.) This identity is the subject
of [Section 32.6](06-bayes.html); state in words what it says about REML.
:::

::: {.solution}
Write \( (\y-\X\bbeta)\T\V^{-1}(\y-\X\bbeta)
=\y\T\bP\y+(\bbeta-\hbeta)\T\X\T\V^{-1}\X(\bbeta-\hbeta) \), which is the
decomposition behind @lem-mix-P(b). Integrating
\( \exp\{-\tfrac12(\bbeta-\hbeta)\T\X\T\V^{-1}\X(\bbeta-\hbeta)\} \) over
\( \Real^p \) gives \( (2\pi)^{p/2}\det(\X\T\V^{-1}\X)^{-1/2} \). Hence
\[
\begin{aligned}
\int e^{\ell}\,d\bbeta
&=e^{-\frac12\{n\log2\pi+\log\det\V+\y\T\bP\y\}}\\
&\qquad\times(2\pi)^{p/2}\det(\X\T\V^{-1}\X)^{-1/2},
\end{aligned}
\]
whose logarithm is @eq-mix-reml-vs-ml without the \( \det(\X\T\X) \) term. In words:
the restricted likelihood is the marginal likelihood of
\( \boldsymbol{\uptheta} \) after integrating the fixed effects out with a flat prior.
:::
