# The additive model

Suppose a response depends on five continuous covariates and nothing is known about the
shape of that dependence. The honest model is \( \E(Y\mid\mathbf{z})=f(z_1,\dots,z_5) \),
\( f \) unknown and smooth. [Chapter 43](../ch43-smoothing/index.html) estimated such an
\( f \) from one covariate by averaging the responses near the point of interest; in five
dimensions that fails, for a reason that has nothing to do with the estimator.

## Why additivity

Take \( n \) points spread over the unit cube in \( q \) dimensions. To collect a fixed
fraction \( c \) of them, a cube-shaped neighbourhood must have side \( c^{1/q} \): with
\( c=0.05 \) that is \( 0.05 \) in one dimension, about \( 0.55 \) in five and about
\( 0.74 \) in ten. An average over three quarters of the range of every covariate is not local, and no
estimator repairs that: the data are not dense enough. This is the *curse of dimensionality*,
and the rates make it quantitative. For a twice-differentiable curve @thm-smo-local-poly
gives a mean squared error of order \( n^{-4/5} \); in \( q \) dimensions the same
expansion gives \( n^{-4/(4+q)} \), and no estimator does better for twice-differentiable
regression functions (Stone, 1980, whose conditions and proof are not reproduced here). At
\( q=5 \) a twelve-fold increase in sample size buys only what a four-fold increase buys in
one dimension.

The response is to restrict the shape of \( f \) — not to a linear function, usually too
severe a restriction, but to a sum of functions of one variable each.

::: {#def-add-model}
[Additive model]

Let \( (y_i,\x_{(i)},z_{i1},\dots,z_{iq}) \), \( i=1,\dots,n \), be observed, where
\( \x_{(i)}\T \) is the \( i \)th row of an \( n\times p \) matrix \( \X \) whose first
column is \( \bone \). The **additive model** is
\[
y_i=\x_{(i)}\T\bbeta+\sum_{j=1}^{q}f_j(z_{ij})+\varepsilon_i ,
\qquad \E(\varepsilon_i)=0,\ \Var(\varepsilon_i)=\sigma^2 ,
\]{#eq-add-model}

with the \( \varepsilon_i \) uncorrelated, the \( f_j \) unknown smooth functions, and the
**identifiability constraints**
\[
\sum_{i=1}^{n}f_j(z_{ij})=0,\qquad j=1,\dots,q .
\]{#eq-add-centring}

The vector \( \boldsymbol{\upeta}=\X\bbeta+\sum_j\mathbf{f}_j \), where
\( \mathbf{f}_j=\{f_j(z_{1j}),\dots,f_j(z_{nj})\}\T \), is the **additive predictor**. If
\( \X=\bone \) the model is *purely additive*; otherwise it is a **partial linear** or
semiparametric model.
:::

Two comments at once. The model contains no interactions: the effect of \( z_1 \) is the
same curve whatever \( z_2 \) does, and [Section 44.3](03-varying-coefficients.html) buys
some of that flexibility back. And @eq-add-centring is not a restriction on the model but a
choice of representation, because the additive predictor fixes the pieces only up to
constants.

::: {#prp-add-identifiability}
[Identifiability of the components]

Write \( \mathcal{V}_j \) for the linear space of vectors
\( \{g(z_{1j}),\dots,g(z_{nj})\}\T \) as \( g \) ranges over the functions the \( j \)th
term is allowed to take, and assume \( \bone\in\mathcal{V}_j \) for each \( j \).

::: {.enumerate options="label=(\alph*)"}
1. If \( (\bbeta,f_1,\dots,f_q) \) and \( (\tilde{\bbeta},\tilde{f}_1,\dots,\tilde{f}_q) \)
   give the same additive predictor and both satisfy @eq-add-centring, then every
   difference \( \mathbf{f}_j-\tilde{\mathbf{f}}_j \) is orthogonal to \( \bone \), and so
   is \( \X(\bbeta-\tilde{\bbeta})=-\sum_j(\mathbf{f}_j-\tilde{\mathbf{f}}_j) \).
   If in addition the centred spaces
   \( \C(\X)\cap\bone\perpc,\ \mathcal{V}_1\cap\bone\perpc,\dots,\mathcal{V}_q\cap\bone\perpc \)
   admit no nontrivial relation — no *exact concurvity*, the condition studied in
   [Section 44.2](02-backfitting.html) — then \( \X\bbeta=\X\tilde{\bbeta} \) and
   \( \mathbf{f}_j=\tilde{\mathbf{f}}_j \) for every \( j \).

2. Without @eq-add-centring the components are determined only up to constants that sum
   to zero: replacing \( f_j \) by \( f_j+c_j \) with \( \sum_jc_j=0 \) leaves
   \( \boldsymbol{\upeta} \) unchanged.

3. If \( \mathbf{f}_j=\Z_j\bgamma_j \) for a design block \( \Z_j \), the constraint
   \( \bone\T\Z_j\bgamma_j=0 \) is imposed exactly by writing \( \bgamma_j=\bU_j\boldsymbol{\upalpha}_j \),
   where the columns of \( \bU_j \) are an orthonormal basis of
   \( \{\bgamma:\bone\T\Z_j\bgamma=0\} \): minimizing any criterion over \( \boldsymbol{\upalpha}_j \)
   with design \( \Z_j\bU_j \) and penalty \( \bU_j\T\bP_j\bU_j \) is the same as
   minimizing it over \( \bgamma_j \) subject to the constraint.
:::

:::

::: {.proof}
(b) is immediate from \( \sum_j(f_j+c_j)=\sum_jf_j+\sum_jc_j \). For (a), equality of the
predictors gives \( \X(\bbeta-\tilde{\bbeta})=-\sum_j(\mathbf{f}_j-\tilde{\mathbf{f}}_j) \).
Each difference \( \mathbf{f}_j-\tilde{\mathbf{f}}_j \) lies in
\( \mathcal{V}_j\cap\bone\perpc \), because both terms satisfy @eq-add-centring, so the
left-hand side lies in \( \C(\X)\cap\bone\perpc \) and the displayed equation is a relation
among the centred spaces. That alone does not separate the pieces. Under the additional
hypothesis the relation must be trivial, so \( \X(\bbeta-\tilde{\bbeta})=\bzero \) and
every difference vanishes.

For (c), \( \{\bgamma:\bone\T\Z_j\bgamma=0\} \) is the null space of the row vector
\( \mathbf{c}\T=\bone\T\Z_j \), of dimension \( d_j-1 \) when \( \mathbf{c}\ne\bzero \),
and \( \bU_j \) is a bijection from \( \Real^{d_j-1} \) onto it. At
\( \bgamma_j=\bU_j\boldsymbol{\upalpha}_j \) the design term is
\( (\Z_j\bU_j)\boldsymbol{\upalpha}_j \) and the penalty
\( \boldsymbol{\upalpha}_j\T(\bU_j\T\bP_j\bU_j)\boldsymbol{\upalpha}_j \), so the
two minimizations agree at corresponding points and the correspondence is onto the feasible
set.
:::

A QR factorization of \( \Z_j\T\bone \) produces \( \bU_j \) in one line: the first
column of \( \Q \) spans \( \mathbf{c} \), the rest its orthogonal complement. Absorbing
the constraint this way rather than centring the columns of \( \Z_j \) matters, because a
basis with the partition-of-unity property of @prp-smo-bspline has
\( \Z_j\bone=\bone \): subtracting column means leaves \( \Z_j \) with a null vector and
the penalized cross-product singular (@exr-add-column-centring).

## Terms, designs and penalties

Every term in this chapter is represented in a basis. Choose \( d_j \) basis functions for
the \( j \)th component and write
\[
f_j(z)=\sum_{l=1}^{d_j}\gamma_{jl}B_{jl}(z),\qquad
\mathbf{f}_j=\Z_j\bgamma_j,\qquad (\Z_j)_{il}=B_{jl}(z_{ij}),
\]{#eq-add-basis}

after which the constraint of @prp-add-identifiability(c) has been absorbed, so that
\( \Z_j \) has \( d_j \) columns summing to zero. With the B-spline basis of
@def-smo-bspline and the difference penalty of @def-smo-pspline this is a P-spline block,
the default term throughout the chapter; but nothing below uses that particular basis. All
that is used is that the roughness of the \( j \)th component is a quadratic form
\( \bgamma_j\T\bP_j\bgamma_j \) with \( \bP_j \) symmetric nonnegative
definite.

The estimate is then the minimizer of a penalized least squares criterion in the sense of
@def-reg-penalized, with one penalty per term:
\[
Q(\bbeta,\bgamma_1,\dots,\bgamma_q)
=\Bigl\lVert\y-\X\bbeta-\sum_{j=1}^{q}\Z_j\bgamma_j\Bigr\rVert^{2}
+\sum_{j=1}^{q}\lambda_j\,\bgamma_j\T\bP_j\bgamma_j .
\]{#eq-add-criterion}

Collect the pieces: let \( \Z=[\X,\Z_1,\dots,\Z_q] \) with \( d=p+\sum_jd_j \) columns,
let \( \bgamma=(\bbeta\T,\bgamma_1\T,\dots,\bgamma_q\T)\T \), and let
\[
\bP=\begin{pmatrix}
\bzero & & & \\
 & \lambda_1\bP_1 & & \\
 & & \ddots & \\
 & & & \lambda_q\bP_q
\end{pmatrix}
\]{#eq-add-blockpenalty}

be the block-diagonal penalty, with a zero block for the unpenalized columns. Then
@eq-add-criterion is \( \norm{\y-\Z\bgamma}^2+\bgamma\T\bP\bgamma \), a ridge-type
criterion with a structured penalty, and differentiating gives the **penalized normal
equations**
\[
\A\bgamma=\Z\T\y,\qquad \A=\Z\T\Z+\bP .
\]{#eq-add-normal}

::: {#exm-add-constraint}
[Two P-spline blocks]

Take \( q=2 \), a cubic B-spline basis of \( 15 \) functions on equally spaced knots for
each covariate, and the second-difference penalty of @def-smo-pspline. Unconstrained, each
block has \( 15 \) columns summing to \( \bone \) row by row and a penalty of rank
\( 13 \); after @prp-add-identifiability(c) it has \( 14 \) columns summing to zero and a
penalty of the same rank, so \( \Z=[\bone,\Z_1,\Z_2] \) has \( 29 \) columns. Of the two
unpenalized directions per block, constant and linear trend, the constraint removes one: the
term still fits a straight line in \( z_j \) free of charge, but not an arbitrary level.
:::

## The fit is a linear smoother

Because @eq-add-criterion is quadratic, the fit is linear in \( \y \) — a linear smoother
in the sense of @prp-smo-linear-smoother — and every quantity of interest is a trace,
exactly as for ridge regression.

::: {#prp-add-df}
[Smoother matrix and degrees of freedom]

Suppose \( \A=\Z\T\Z+\bP \) is nonsingular, and let
\( \mathbf{C}=\A^{-1}\Z\T \), so that \( \hat{\bgamma}=\mathbf{C}\y \). Write
\( \mathbf{C}_j=\Z_j\mathbf{C}_{[j]} \), where \( \mathbf{C}_{[j]} \) holds the rows of
\( \mathbf{C} \) belonging to the \( j \)th block.

::: {.enumerate options="label=(\alph*)"}
1. The fitted vector is \( \hat{\y}=\bS\y \) with \( \bS=\Z\A^{-1}\Z\T \), a symmetric
   matrix whose eigenvalues lie in \( [0,1] \); the fitted \( j \)th component is
   \( \hat{\mathbf{f}}_j=\mathbf{C}_j\y \), and
   \( \bS=\X\mathbf{C}_{[0]}+\sum_j\mathbf{C}_j \).

2. The **effective degrees of freedom** of the fit and of its terms are
   \[
   \mathrm{df}=\tr\bS=\tr(\A^{-1}\Z\T\Z),\qquad
   \mathrm{df}_j=\tr\mathbf{C}_j=d_j-\lambda_j\tr\bigl\{(\A^{-1})_{jj}\bP_j\bigr\},
   \]{#eq-add-edf}

   and they add: \( \mathrm{df}=p+\sum_{j}\mathrm{df}_j \). Each
   \( \mathrm{df}_j \) decreases from \( d_j \) at \( \lambda_j=0 \) towards
   \( \dim\Null(\bP_j) \) as \( \lambda_j\to\infty \).

3. If \( \Cov(\Y)=\sigma^2\I \), then
   \( \Cov(\hat{\bgamma})=\sigma^2\A^{-1}\Z\T\Z\A^{-1} \) and
   \( \Cov(\hat{\mathbf{f}}_j)=\sigma^2\mathbf{C}_j\mathbf{C}_j\T \).
:::

:::

::: {.proof}
(a) Solving @eq-add-normal gives \( \hat{\bgamma}=\A^{-1}\Z\T\y \) and
\( \hat{\y}=\Z\hat{\bgamma}=\bS\y \). Symmetry is clear, and \( \bv\T\bS\bv\ge0 \) because \( \A^{-1} \) is positive
definite. For the upper bound, fix \( \bv \) and put \( \bu=\A^{-1}\Z\T\bv \), so
that \( \A\bu=\Z\T\bv \). Then
\[
\bv\T\bS\bv=\bu\T\A\bu=\norm{\Z\bu}^2+\bu\T\bP\bu\ \ge\ \norm{\Z\bu}^2 ,
\]
while also \( \bv\T\bS\bv=\bv\T(\Z\bu)\le\norm{\bv}\norm{\Z\bu}
\le\norm{\bv}\sqrt{\bv\T\bS\bv} \) by Cauchy--Schwarz. Dividing by
\( \sqrt{\bv\T\bS\bv} \) when it is nonzero gives
\( \bv\T\bS\bv\le\norm{\bv}^2 \), so the eigenvalues lie in \( [0,1] \).
Splitting \( \Z\hat{\bgamma} \) by blocks gives the last identity.
(b) \( \tr\bS=\tr(\Z\A^{-1}\Z\T)=\tr(\A^{-1}\Z\T\Z) \) by the cyclic property, and
\( \A^{-1}\Z\T\Z=\A^{-1}(\A-\bP)=\I-\A^{-1}\bP \). The \( j \)th diagonal
block of \( \A^{-1}\bP \) is \( \lambda_j(\A^{-1})_{jj}\bP_j \), and the block
belonging to the unpenalized columns is zero; taking traces block by block gives both
formulas in @eq-add-edf and their sum. Monotonicity in \( \lambda_j \) is @exr-add-df-monotone,
and as \( \lambda_j\to\infty \) the term is forced into \( \Null(\bP_j) \), which
costs \( \dim\Null(\bP_j) \) dimensions.
(c) Both are immediate from
\( \Cov(\mathbf{M}\Y)=\mathbf{M}\Cov(\Y)\mathbf{M}\T \) (@thm-rv-linear).
:::

Part (b) is the additive-model version of @thm-shr-ridge(b): each term is charged for the
dimensions it uses, and the charges add up. The split makes a printed summary readable —
"living area: 5.3 degrees of freedom" says at a glance how far from a straight line the
curve is — and it is what [Chapter 29](../ch29-model-selection/index.html) needs, since
@def-sel-gcv and @prp-sel-aic-bic apply verbatim with \( \mathrm{df} \) in place of the
number of parameters.

::: {.warning}
The degrees of freedom of a *term* are not the degrees of freedom of a test: they are a
trace computed at a fixed \( \boldsymbol{\uplambda} \), and nothing so far says that
\( \mathrm{df}_j \) is the right divisor for a sum of squares or the right parameter for a
reference distribution. [Section 44.6](06-generalized-additive-models.html) returns to it.
:::

## Exercises

### A. Check your understanding

::: {#exr-add-curse}
[A1]

Repeat the neighbourhood calculation of [Section 44.1](01-additive-models.html): what side
length does a cube-shaped neighbourhood need in order to contain \( 1\% \) of \( n \)
points uniformly spread over the unit cube in \( q \) dimensions, for
\( q=1,3,10 \)? At what \( q \) does the required side exceed \( 0.5 \)?
:::

### B. Practice

::: {#exr-add-column-centring}
[B1]

Let \( \Z_j \) be a B-spline block with the partition-of-unity property, so
\( \Z_j\bone=\bone \), and let \( \bP_j \) be a difference penalty, so
\( \bP_j\bone=\bzero \). Show that centring the *columns*, that is replacing
\( \Z_j \) by \( (\I-n^{-1}\bone\bone\T)\Z_j \), makes \( \Z_j\T\Z_j+\lambda_j\bP_j \)
singular, while the reparameterization of @prp-add-identifiability(c) does not.
:::

::: {.solution}
Write \( \tilde{\Z}_j=(\I-n^{-1}\bone\bone\T)\Z_j \). Then
\( \tilde{\Z}_j\bone=\bzero \), so
\( \bone\T(\tilde{\Z}_j\T\tilde{\Z}_j+\lambda_j\bP_j)\bone=0 \) and the matrix
is singular. Under (c) the matrix is
\( \bU_j\T(\Z_j\T\Z_j+\lambda_j\bP_j)\bU_j \); if it killed some
\( \boldsymbol{\upalpha}\ne\bzero \), then \( \bgamma=\bU_j\boldsymbol{\upalpha} \)
would satisfy \( \Z_j\bgamma=\bzero \) and \( \bP_j\bgamma=\bzero \). The second
puts \( \bgamma \) in the span of the constant, or of the constant and linear sequences for
a second-difference penalty, and then \( \Z_j\bgamma=\bzero \) forces
\( \bgamma=\bzero \), since the B-spline basis reproduces those functions exactly.
:::

::: {#exr-add-df-monotone}
[B2]

Show that \( \mathrm{df}_j \) in @eq-add-edf is nonincreasing in \( \lambda_j \). *Hint:*
work with one term and no parametric part first, where
\( \mathrm{df}(\lambda)=\sum_k\theta_k/(\theta_k+\lambda) \) for the eigenvalues
\( \theta_k \) of \( \Z_j\T\Z_j \) in the metric of \( \bP_j \).
:::

::: {.solution}
Simultaneously diagonalize \( \Z_j\T\Z_j \) and \( \bP_j \): some nonsingular
\( \bT \) gives \( \bT\T\Z_j\T\Z_j\bT=\I \) and
\( \bT\T\bP_j\bT=\diag(\nu_k) \), \( \nu_k\ge0 \), by
@cor-mat-generalized-rayleigh(b). Then
\( \mathrm{df}(\lambda)=\sum_k1/(1+\lambda\nu_k) \), decreasing in \( \lambda \) except
where \( \nu_k=0 \), which contributes \( \dim\Null(\bP_j) \) whatever
\( \lambda \) is. With several terms, argue the same way on the \( j \)th block of
\( \A^{-1}\bP \), treating the others as part of \( \Z\T\Z \).
:::

::: {#exr-add-partial-linear}
[B3]

In @def-add-model take \( q=1 \) and let \( \X \) be a general parametric design. Show
that the penalized estimate of \( \bbeta \) is
\( \hbeta=\{\X\T(\I-\bS_1)\X\}^{-1}\X\T(\I-\bS_1)\y \) where
\( \bS_1=\Z_1(\Z_1\T\Z_1+\lambda_1\bP_1)^{-1}\Z_1\T \) — the Frisch–Waugh–Lovell
form of @thm-proj-fwl with the smoother in place of a projection. What changes in the
proof?
:::

### C. Going deeper

::: {#exr-add-interaction}
[C1]

An additive model has no interactions. Show that if the true mean is
\( f(z_1,z_2)=z_1z_2 \) on \( [0,1]^2 \) with \( z_1,z_2 \) independent and uniform, the
best additive approximation in mean square,
\( \argmin_{g_1,g_2}\E\{(Z_1Z_2-g_1(Z_1)-g_2(Z_2))^2\} \) subject to
\( \E g_1(Z_1)=\E g_2(Z_2)=0 \), is
\( g_1(z)=(z-\tfrac12)/2 \), \( g_2(z)=(z-\tfrac12)/2 \) plus the constant \( 1/4 \), and
compute the mean squared approximation error.
:::

