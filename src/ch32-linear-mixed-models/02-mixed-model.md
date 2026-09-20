# The linear mixed model

The one-way random effects model of [Section 32.1](01-random-effects.html) has one
random factor and no regressors. Real studies have both: a trial has fixed treatment
effects and random patient effects, a panel has fixed economic regressors and random
firm effects, a field trial has fixed varieties and random blocks. The model that
holds them all keeps the usual fixed part and adds a second design matrix whose
coefficients are random.

## Definition and the marginal model

::: {#def-mix-model}
[Linear mixed model]

Let \( \X \) be a known \( n\times p \) matrix and \( \Z \) a known \( n\times q \)
matrix. The **linear mixed model** is
\[
\Y=\X\bbeta+\Z\bu+\be ,
\]{#eq-mix-model}

where \( \bbeta\in\Real^p \) is a vector of unknown constants, the **fixed
effects**; \( \bu \) is a \( q\times1 \) random vector, the **random effects**, with
\( \E(\bu)=\bzero \) and \( \Cov(\bu)=\G \); \( \be \) is an \( n\times1 \) random
vector with \( \E(\be)=\bzero \) and \( \Cov(\be)=\R \); and
\( \Cov(\bu,\be)=\bzero \). The matrices \( \G \) and \( \R \) are nonnegative
definite and are known functions \( \G(\boldsymbol{\uptheta}) \),
\( \R(\boldsymbol{\uptheta}) \) of an unknown parameter
\( \boldsymbol{\uptheta} \), the **variance components**. The **normal** linear
mixed model adds \( \bu\sim\Normal_q(\bzero,\G) \) independent of
\( \be\sim\Normal_n(\bzero,\R) \).
:::

Taking expectations in @eq-mix-model and using @thm-rv-linear gives the **marginal
model**
\[
\E(\Y)=\X\bbeta,\qquad
\Cov(\Y)=\V(\boldsymbol{\uptheta})=\Z\G\Z\T+\R ,
\]{#eq-mix-marginal}

and under normality \( \Y\sim\Normal_n(\X\bbeta,\V) \). Everything Part VII knows
about the general Gauss–Markov model applies to @eq-mix-marginal at once. What
distinguishes a mixed model from the general model of
[Chapter 31](../ch31-general-gauss-markov/index.html) is that \( \V \) is not an
arbitrary unknown matrix but is built from a handful of parameters by a structure the
experiment dictates. An \( n\times n \) covariance has \( n(n+1)/2 \) free entries; a
random intercept model has two.

::: {.remark}
[Conditional and marginal]

Two means live in @eq-mix-model. Conditionally on the random effects,
\( \E(\Y\mid\bu)=\X\bbeta+\Z\bu \): this is what the subjects in *this* sample did.
Marginally, \( \E(\Y)=\X\bbeta \): this is what the population does on average. In
the linear model with normal errors the two lead to the same \( \bbeta \), and the
distinction looks pedantic. It is not pedantic at all once the link is nonlinear.
For a generalized linear mixed model (Chapter 40) the marginal mean of \( g^{-1}
(\X\bbeta+\Z\bu) \) is not \( g^{-1}(\X\bbeta) \), and the subject-specific and
population-averaged coefficients genuinely differ.
:::

::: {.warning}
The pair \( (\G,\R) \) is not identified by the marginal distribution; only
\( \V=\Z\G\Z\T+\R \) is. In the one-way layout with \( \Z \) the group indicators,
the model with \( \G=\sigma_a^2\I_g \), \( \R=\sigma^2\I_n \) has exactly the same
marginal distribution as the model with no random effects at all and
\( \R=\sigma^2\I+\sigma_a^2\Z\Z\T \), that is, with equicorrelated errors within
groups. Every likelihood-based procedure in this chapter estimates \( \V \), so it
cannot tell the two apart. What the random-effects reading adds is a *target for
prediction*: the vector \( \bu \), whose realized value
[Section 32.3](03-blup.html) predicts. It also adds a constraint, since
\( \sigma_a^2\ge0 \) forces a nonnegative intraclass correlation, whereas the
equicorrelated model allows \( \rho \) down to \( -1/(m-1) \).
:::

## Inverting the covariance matrix

Formulas in this chapter are written with \( \V^{-1} \), which is \( n\times n \).
When \( q \) is much smaller than \( n \), or when the data fall into independent
clusters, \( \V^{-1} \) never has to be formed.

::: {#prp-mix-vinverse}
[The inverse and determinant of V]

Let \( \R \) and \( \G \) be positive definite in @def-mix-model. Then \( \V \) is
positive definite, and
\[
\V^{-1}=\R^{-1}-\R^{-1}\Z\bigl(\G^{-1}+\Z\T\R^{-1}\Z\bigr)^{-1}\Z\T\R^{-1},
\]{#eq-mix-v-inverse}

\[
\det\V=\det\R\,\det\G\,\det\bigl(\G^{-1}+\Z\T\R^{-1}\Z\bigr).
\]{#eq-mix-v-determinant}

If in addition the observations fall into \( g \) clusters with
\( \R=\bigoplus_k\R_k \), \( \Z=\bigoplus_k\Z_k \) and \( \G=\bigoplus_k\G_k \) block
diagonal, the last saying that different clusters' random effects are independent,
then
\( \V=\bigoplus_k(\Z_k\G_k\Z_k\T+\R_k) \) is block diagonal too, and everything
decomposes cluster by cluster.
:::

::: {.proof}
\( \V=\R+\Z\G\Z\T \) is a sum of a positive definite and a nonnegative definite
matrix, hence positive definite. Now use @thm-mat-woodbury with the roles
\( \B\mapsto\R \), \( \bU\mapsto\Z \), \( \mathbf{C}\mapsto\G \) and (in the
notation of that theorem) the second outer factor also \( \Z \). The matrix
\( \G^{-1}+\Z\T\R^{-1}\Z \) is positive definite, hence nonsingular, so the
hypothesis holds and @eq-mat-woodbury is @eq-mix-v-inverse. For the determinant,
@thm-mat-block-determinant(c) with \( \B=\R \), \( \bU=\Z\G \) and the second factor
\( \Z\T \) gives
\[
\begin{aligned}
\det\V&=\det\R\,\det(\I_q+\Z\T\R^{-1}\Z\G)\\
&=\det\R\,\det\G\,\det(\G^{-1}+\Z\T\R^{-1}\Z),
\end{aligned}
\]
the last step by \( \det(\I_q+\A\G)=\det(\G^{-1}+\A)\det\G \). The block diagonal
statement is immediate from the definitions.
:::

A study of \( n=10^5 \) pupils in \( q=500 \) schools therefore needs a
\( 500\times500 \) inverse, not a \( 10^5\times10^5 \) one, and if \( \R \) is
diagonal the rest of @eq-mix-v-inverse costs nothing.
[Section 32.3](03-blup.html) improves on this again.

## Four examples

::: {#exm-mix-random-intercept}
[Random intercepts]

Observations are grouped, group \( k \) contributing \( n_k \) of them. Take
\( \Z \) to be the \( n\times g \) matrix of group indicators,
\( \G=\sigma_a^2\I_g \) and \( \R=\sigma^2\I_n \), so that
\( \V=\sigma^2\I+\sigma_a^2\Z\Z\T \) as in @prp-mix-oneway(a). The fixed part
\( \X\bbeta \) may hold anything: an intercept alone, recovering
@def-mix-oneway; a treatment indicator constant within group, as in a
cluster-randomized trial; or regressors that vary within a group, as in a
longitudinal study. Each group's block of \( \V \) is
\( \sigma^2\I_{n_k}+\sigma_a^2\bone\bone\T \), the **compound symmetric** or
equicorrelated form.
:::

::: {#exm-mix-random-slope}
[Random intercepts and slopes]

Subject \( i \) is measured at times \( t_{i1},\dots,t_{in_i} \), and each subject has
its own intercept and its own slope in time:
\[
Y_{ij}=\x_{(ij)}\T\bbeta+u_{0i}+u_{1i}t_{ij}+\varepsilon_{ij},
\qquad
\begin{pmatrix}u_{0i}\\u_{1i}\end{pmatrix}\sim\Normal_2(\bzero,\bD),
\]
independently across subjects, with \( \be_i\sim\Normal_{n_i}(\bzero,\sigma^2\I) \).
Here \( \Z_i=[\bone,\mathbf{t}_i] \), \( \G=\I_N\otimes\bD \) for \( N \) subjects,
and subject \( i \)'s covariance block is
\( \V_i=\Z_i\bD\Z_i\T+\sigma^2\I \), with entries
\[
\Cov(Y_{ij},Y_{ij'})=d_{11}+d_{12}(t_{ij}+t_{ij'})+d_{22}t_{ij}t_{ij'}
+\sigma^2\mathbf{1}\{j=j'\} .
\]
Two features deserve notice. The variance of a response is quadratic in time, so a
random slope model necessarily says that subjects spread apart, or, if
\( d_{12}<0 \), first converge and then spread. And the correlation between two
occasions depends on both times, not only on their separation: it is not stationary.
Whether this is the right shape is a modelling question that
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) takes up.
:::

::: {#exm-mix-nested-crossed}
[Nested and crossed random factors]

With pupils inside classrooms inside schools, write
\( \Z=[\Z_{\text{sch}},\Z_{\text{cls}}] \) and
\( \G=\begin{psmallmatrix}\sigma_S^2\I&\bzero\\\bzero&\sigma_C^2\I\end{psmallmatrix} \).
Because every classroom belongs to exactly one school, \( \V \) is block diagonal by
school, and the analysis splits school by school; this is the covariance
behind @exr-tw-random-classrooms. With subjects responding to items, each subject seeing
every item, the two factors are **crossed**: \( \Z_{\text{sub}} \) and
\( \Z_{\text{item}} \) generate overlapping patterns and \( \V \) has no block
diagonal form at all. Nothing in the theory changes, but the computation does, and
crossed random effects are the reason modern software stores \( \Z \) sparsely.
:::

::: {#exm-mix-grunfeld}
[Grunfeld's panel with a random firm effect]

[Section 6.6](../ch06-projections/06-fwl.html) fitted gross investment on market
value and capital stock for \( 11 \) US firms over
\( 20 \) years, once with a separate intercept per firm (the
within estimator) and once with a single intercept. The mixed model puts the firm
effect in between: \( \X \) holds the intercept, value and capital,
\( \Z \) holds the firm indicators, \( \G=\sigma_f^2\I_{11} \) and
\( \R=\sigma^2\I \). Fitted by restricted maximum likelihood
([Section 32.4](04-likelihood.html)), it gives
\( \hat\sigma^2=2527.6 \) and
\( \hat\sigma_f^2=6740.9 \), a firm standard deviation of
\( 82.1 \) and an intraclass correlation of
\( 0.727 \): most of the variation in investment is between firms.
The coefficients are \( 0.1094 \) for value and
\( 0.3082 \) for capital, against
\( 0.1101 \) and \( 0.3100 \) for the
within estimator and \( 0.1145 \) and
\( 0.2275 \) with a single intercept. The mixed estimate has
landed essentially on the within estimate, and
[Section 32.3](03-blup.html) explains why: with \( 20 \) years per
firm, the weight the model gives to between-firm comparisons is
\( 1-0.9816 \), less than two per cent.
:::

```{.python .run #cell-grunfeld-mixed-data}
import numpy as np
import pandas as pd
import statsmodels.api as sm

df = sm.datasets.grunfeld.load_pandas().data.sort_values(["firm", "year"])
firms = df["firm"].unique()
y = df["invest"].to_numpy()
X = np.column_stack([np.ones(len(y)), df["value"], df["capital"]])
blocks = [np.where(df["firm"].to_numpy() == f)[0] for f in firms]
n, p, g = len(y), X.shape[1], len(firms)
print(f"{n} observations, {g} firms, {n // g} years each")
```

## When ordinary least squares is already best

The marginal model @eq-mix-marginal is a general Gauss–Markov model, so
@cor-opt-aitken says the best linear unbiased estimator of an estimable
\( \blambda\T\bbeta \) is the generalized least squares estimate built from
\( \V^{-1} \), and @thm-ggm-ols-blue says that ordinary least squares achieves this
for every \( \y \) exactly when \( \C(\V\X)\subseteq\C(\X) \). Mixed models supply
the standard examples on both sides of that line.

::: {#exm-mix-ols-blue}
[Cluster-constant regressors]

Take @exm-mix-random-intercept with a balanced layout, \( n_k=m \), and suppose every
column of \( \X \) is constant within groups, so that \( \C(\X)\subseteq\C(\Z) \).
By @eq-mix-balanced-cov, \( \V=\lambda_1\M+\sigma^2(\I-\M) \) with \( \M \) the
projection onto \( \C(\Z) \), so \( \V\X=\lambda_1\X \) and
\( \C(\V\X)=\C(\X) \). Ordinary least squares is therefore the best linear unbiased
estimator, whatever the variance components are, and the analysis of a balanced
cluster-randomized trial can be carried out on the cluster means without loss. The
covariance of the estimate is *not* \( \sigma^2(\X\T\X)^{-1} \), however: it is
\( \lambda_1(\X\T\X)^{-1} \), larger by the factor \( 1+m\gamma \). Getting the point
estimate right and the standard error wrong is the classic error in clustered data.
:::

::: {#exm-mix-ols-not-blue}
[Regressors that vary within a cluster]

Keep the same \( \V \), but let \( \X=[\bone,\x] \) with \( \x \) varying within
groups. Decomposing \( \x \) into its between and within parts,
\[
\V\x=\lambda_1\M\x+\sigma^2(\I-\M)\x ,
\]
which lies in \( \C(\X) \) only if one of the two parts vanishes or
\( \lambda_1=\sigma^2 \). Ordinary least squares is then not efficient, and the
reason is visible: generalized least squares downweights the between-group part of
\( \x \) by \( \sigma^2/\lambda_1=1/(1+m\gamma) \) relative to the within-group
part, because between-group comparisons carry the extra variance \( \sigma_a^2 \).
The mixed estimate is a compromise between the within estimator of
[Section 6.6](../ch06-projections/06-fwl.html), which uses only \( (\I-\M)\x \), and
the between estimator, which uses only \( \M\x \); as \( m\gamma\to\infty \) the
weight on the between part goes to zero. @exm-mix-grunfeld is that limit in
practice.
:::

::: {.idea}
A random effect is a statement about *where the information is*. Declaring the
grouping factor random says that comparisons between groups carry an extra variance
\( \sigma_a^2 \) and should be downweighted accordingly, but not discarded. Fixed
effects for the groups discard them entirely; ignoring the grouping keeps them at
full weight. The mixed model interpolates, and the interpolation weight is estimated
from the data.
:::

The remaining difficulty is that \( \V \) depends on the unknown
\( \boldsymbol{\uptheta} \). The estimator actually used is the feasible generalized
least squares estimator of @def-ggm-feasible,
\[
\hbeta(\hat{\boldsymbol{\uptheta}})
=\bigl(\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\X\bigr)^{-1}
\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\Y ,
\]{#eq-mix-feasible}

with \( \hat{\boldsymbol{\uptheta}} \) from
[Section 32.4](04-likelihood.html). @thm-ggm-feasible gives the large-sample
justification; [Section 32.5](05-inference.html) asks what it costs in small samples,
and the answer is: more than people usually allow for.

## Exercises

### A. Check your understanding

::: {#exr-mix-write-down-v}
[A1]

A study measures \( 3 \) subjects on \( 2 \) occasions each, with a random subject
effect of variance \( \sigma_s^2 \) and independent errors of variance
\( \sigma^2 \). Write \( \Z \), \( \G \), \( \R \) and \( \V \) explicitly as
\( 6\times6 \) or \( 6\times3 \) matrices, and check @eq-mix-v-inverse numerically for
\( \sigma_s^2=2 \), \( \sigma^2=1 \).
:::

::: {#exr-mix-identifiability}
[A2]

A colleague argues that in @exm-mix-random-intercept the two variance components
cannot be separated, since raising \( \sigma_a^2 \) and lowering \( \sigma^2 \) both
change the diagonal of \( \V \) in the same direction. Is the argument right? What
happens if every group has \( n_k=1 \)?
:::

### B. Practice

::: {#exr-mix-vinverse-oneway}
[B1]

Specialize @eq-mix-v-inverse to @exm-mix-random-intercept with group sizes \( n_k \),
and show that the \( k \)th block of \( \V^{-1} \) is
\[
\frac{1}{\sigma^2}\Bigl(\I_{n_k}
-\frac{\sigma_a^2}{\sigma^2+n_k\sigma_a^2}\bone\bone\T\Bigr).
\]
Deduce that the generalized least squares criterion is
\( \sigma^{-2}\sum_k\{\sum_j r_{kj}^2-w_kn_k\bar r_k^2\} \) with
\( w_k=n_k\sigma_a^2/(\sigma^2+n_k\sigma_a^2) \), and interpret \( w_k \).
:::

::: {.solution}
Each block of \( \V \) is \( \sigma^2\I+\sigma_a^2\bone\bone\T \); apply
@eq-mat-sherman-morrison with \( \B=\sigma^2\I \), \( \bu=\bv=\sigma_a\bone \):
the correction is
\( \sigma_a^2\sigma^{-4}\bone\bone\T/(1+n_k\sigma_a^2/\sigma^2) \), which rearranges
to the stated form. The criterion \( \br\T\V^{-1}\br \) then equals
\( \sigma^{-2}\sum_k\{\sum_jr_{kj}^2-n_k\bar r_k^2\,w_k\} \) after writing
\( \bone\T\br_k=n_k\bar r_k \). The weight \( w_k\in(0,1) \) increases with
\( n_k \): the model removes a fraction \( w_k \) of the group's mean residual from
the criterion, because that part of the residual is explained by the group's random
effect rather than counted against the fit.
:::

::: {#exr-mix-ols-blue-blocks}
[B2]

Let \( \V=\lambda_1\M+\sigma^2(\I-\M) \) as in @eq-mix-balanced-cov. Show that
\( \C(\V\X)\subseteq\C(\X) \) holds if and only if \( \C(\X) \) decomposes as
\( (\C(\X)\cap\C(\Z))\dirsum(\C(\X)\cap\C(\Z)\perpc) \), that is, iff \( \M \) maps
\( \C(\X) \) into itself. Give an \( \X \) with both between and within columns for
which ordinary least squares *is* best.
:::

::: {.solution}
\( \V\X=\sigma^2\X+(\lambda_1-\sigma^2)\M\X \), so (assuming
\( \lambda_1\ne\sigma^2 \)) \( \C(\V\X)\subseteq\C(\X) \) iff
\( \C(\M\X)\subseteq\C(\X) \), that is, \( \M \) maps \( \C(\X) \) into itself. A
subspace invariant under an orthogonal projection splits into its parts inside
\( \C(\Z) \) and inside \( \C(\Z)\perpc \). Example: \( \X=[\bone,\x] \) with
\( \bone \) constant within groups and \( \x \) having zero mean *within every
group*, so \( \M\x=\bzero \); then \( \M\X=[\bone,\bzero] \) lies in \( \C(\X) \). A
balanced split-plot design has exactly this structure, which is why
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) can analyse it
stratum by stratum.
:::

### C. Going deeper

::: {#exr-mix-crossed-rank}
[C1]

Two crossed random factors with \( a \) and \( b \) levels and one observation per
cell give \( \Z=[\Z_A,\Z_B] \) with \( n=ab \), \( q=a+b \). Show that
\( \rank(\Z)=a+b-1 \), so \( \Z\T\Z \) is singular, and explain why this causes no
difficulty for @eq-mix-marginal but would if one tried to treat \( \bu \) as a
parameter vector to be estimated. Compute \( \V \) and show it is *not* block
diagonal under any ordering of the observations when \( a,b\ge2 \).
:::
