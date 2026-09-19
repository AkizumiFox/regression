# Testable hypotheses

Scientific questions are usually phrased in terms of coefficients rather than subspaces: "the
coefficients of poverty and single parenthood are equal", "the four regions have the same
intercept", "the elasticity is one". Such statements have the form
\[
H:\ \bLambda\T\bbeta=\mathbf{d} ,
\]
with a known \( p\times m \) matrix \( \bLambda \) and a known \( \mathbf{d}\in\Real^m \). This is the
**general linear hypothesis**. This section decides when such a statement is a hypothesis about the
distribution of \( \Y \) at all, gives its \( F \) statistic in terms of the unrestricted fit alone,
shows that the statistic agrees with the reduced-model statistic of
[Section 11.1](01-reduced-models.html), and describes least squares under the hypothesis. Nothing
here assumes that \( \X \) has full column rank.

## Which constraints are hypotheses

In the full model the distribution of \( \Y \) depends on \( \bbeta \) only through \( \X\bbeta \). A
statement about \( \bbeta \) that cannot be read off from \( \X\bbeta \) is therefore not a statement about
the data-generating distribution, and no data can test it.

::: {#def-glh-testable}
[Testable hypothesis]

Let \( \bLambda \) be a \( p\times m \) matrix and \( \mathbf{d}\in\Real^m \). The hypothesis
\( H:\bLambda\T\bbeta=\mathbf{d} \) is **testable** if

::: {.enumerate options="label=(\roman*)"}
1. \( \bLambda\T\bbeta \) is estimable, that is, \( \C(\bLambda)\subseteq\C(\X\T) \) (@def-est-estimable and
   @thm-est-characterization), so that \( \bLambda=\X\T\bT \) for some \( n\times m \) matrix \( \bT \); and

2. the equations \( \bLambda\T\bb=\mathbf{d} \) are consistent, that is, \( \mathbf{d}\in\C(\bLambda\T) \)
   (@thm-mat-consistency).
:::

Its number of degrees of freedom is \( q=\rank(\bLambda) \). We assume \( q\ge1 \).
:::

Condition (ii) is automatic when \( \bLambda \) has full column rank, for then \( \C(\bLambda\T)=\Real^m \).
It matters when the columns of \( \bLambda \) are redundant, as when all six pairwise differences of four
group effects are set to zero. Redundant constraints are harmless if they agree with each other. If
they do not, as in \( \alpha_1-\alpha_2=1 \), \( \alpha_2-\alpha_3=1 \), \( \alpha_1-\alpha_3=0 \), the "hypothesis" is
false for every \( \bbeta \), and no test is possible.

Condition (i) is the substantive one. Without it the constraint cannot be settled by any amount of
data.

::: {#prp-glh-nonestimable}
[Nonestimable constraints cannot be tested]

Let \( \blambda\notin\C(\X\T) \) and \( c\in\Real \).

::: {.enumerate options="label=(\alph*)"}
1. For every \( \bbeta\in\Real^p \) there is \( \bbeta' \) with \( \X\bbeta'=\X\bbeta \) and
   \( \blambda\T\bbeta'=c \).

2. Let \( \phi \) be any test of \( H:\blambda\T\bbeta=c \) in the model @eq-opt-normal-model, with rejection
   probability \( \pi(\bbeta,\sigma^2)=\E_{\bbeta,\sigma^2}\phi(\Y) \). For every point \( (\bbeta,\sigma^2) \) of the
   alternative there is a point \( (\bbeta',\sigma^2) \) of the hypothesis with
   \( \pi(\bbeta',\sigma^2)=\pi(\bbeta,\sigma^2) \). So no test has power exceeding its size anywhere.
:::

:::

::: {.proof}
(a) By @thm-est-characterization(c), \( \blambda\notin\C(\X\T) \) means that \( \blambda \) is not orthogonal to
\( \Null(\X) \): there is \( \bv \) with \( \X\bv=\bzero \) and \( \blambda\T\bv\ne0 \). Put
\( \bbeta'=\bbeta+t\bv \) with \( t=(c-\blambda\T\bbeta)/\blambda\T\bv \). Then \( \X\bbeta'=\X\bbeta \) and
\( \blambda\T\bbeta'=c \). (b) The law of \( \Y \), and hence \( \pi \), depends on \( (\bbeta,\sigma^2) \) only through
\( (\X\bbeta,\sigma^2) \). Apply (a).
:::

For a matrix \( \bLambda \), @prp-ss-effective-df gives the general statement: the constraint
\( \bLambda\T\bb=\bzero \) restricts the mean vector only through its estimable part
\( \C(\bLambda)\cap\C(\X\T) \), and the extra sum of squares has \( \dim(\C(\bLambda)\cap\C(\X\T)) \) degrees of
freedom. In the one-way layout \( \E(Y_{ij})=\mu+\alpha_i \), the constraint \( \alpha_1=0 \) is a *side
condition* (@thm-est-side-conditions). It selects one of the equivalent parameterizations and says
nothing about the data. Programs that report a "test" of \( \alpha_1=0 \) are testing something else,
usually a contrast such as \( \alpha_1-\alpha_g \) that the program's coding makes estimable.

## The F statistic of a testable hypothesis

The sum of squares of a testable hypothesis was computed in @thm-ss-hypothesis for \( \bLambda \) of full
column rank. The next theorem removes that restriction, identifies the reduced model, and supplies
the distribution.

::: {#thm-glh-general-f}
[The \( F \) test of a testable hypothesis]

Assume @eq-opt-normal-model with \( r<n \), and let \( H:\bLambda\T\bbeta=\mathbf{d} \) be testable with
\( q=\rank(\bLambda) \). Let \( \G \) be any generalized inverse of \( \X\T\X \), \( \hbeta \) any least squares
estimate, \( \W=\bLambda\T\G\bLambda \) and \( \W\ginv \) any generalized inverse of \( \W \). Put
\[
\text{SS}_H=(\bLambda\T\hbeta-\mathbf{d})\T\,\W\ginv\,(\bLambda\T\hbeta-\mathbf{d}),
\qquad
F_H=\frac{\text{SS}_H/q}{\text{SSE}/(n-r)} .
\]{#eq-glh-general-F}

::: {.enumerate options="label=(\alph*)"}
1. \( \text{SS}_H \) does not depend on the choices of \( \hbeta \), \( \G \) and \( \W\ginv \).

2. Let \( \bb_0 \) be any solution of \( \bLambda\T\bb_0=\mathbf{d} \) and
   \( \mathcal S_0=\{\X\bb:\bLambda\T\bb=\bzero\} \). Then \( \mathcal S_0 \) is a subspace of \( \C(\X) \) of
   dimension \( r-q \), the set of mean vectors allowed by \( H \) is \( \X\bb_0+\mathcal S_0 \), and
   \[
\text{SS}_H=\min_{\bLambda\T\bb=\mathbf{d}}\norm{\y-\X\bb}^2-\text{SSE}
   =\norm{(\M-\bP_{\mathcal S_0})(\y-\X\bb_0)}^2 .
\]

3. \( F_H\sim F(q,\,n-r,\,\gamma) \) with
   \[
\gamma=\frac{(\bLambda\T\bbeta-\mathbf{d})\T\,\W\ginv\,(\bLambda\T\bbeta-\mathbf{d})}{\sigma^2},
\]{#eq-glh-general-noncentrality}

   and \( \gamma=0 \) iff \( \bLambda\T\bbeta=\mathbf{d} \). In particular \( F_H \) is the \( F \) statistic of
   @thm-glh-f-test for the data \( \y-\X\bb_0 \) and the reduced model space \( \mathcal S_0 \).
:::

:::

::: {.proof}
*Reduction to full column rank.* Choose \( q \) linearly independent columns of \( \bLambda \), forming a
\( p\times q \) matrix \( \bLambda_1 \). Every column of \( \bLambda \) is a combination of them, so
\( \bLambda=\bLambda_1\mathbf{K} \) for a \( q\times m \) matrix \( \mathbf{K} \), which has rank \( q \) because
\( \rank\bLambda=q \). Hence \( \mathbf{K}\mathbf{K}\T \) is nonsingular. Put \( \mathbf{d}_1=\bLambda_1\T\bb_0 \), so that
\( \mathbf{d}=\mathbf{K}\T\mathbf{d}_1 \). Since \( \mathbf{K}\T \) has full column rank it is one-to-one, so \( \mathbf{d}_1 \) is
determined by \( \mathbf{d} \) and does not depend on the choice of \( \bb_0 \). Moreover
\( \bLambda\T\bb=\mathbf{d} \iff \mathbf{K}\T(\bLambda_1\T\bb-\mathbf{d}_1)=\bzero\iff\bLambda_1\T\bb=\mathbf{d}_1 \). The same
argument with \( \mathbf{d} \) replaced by \( \bzero \) shows that \( \mathcal S_0 \) is also
\( \{\X\bb:\bLambda_1\T\bb=\bzero\} \). The columns of \( \bLambda_1 \) are columns of \( \bLambda \), so
\( \bLambda_1\T\bbeta \) is estimable.

*(a).* By @thm-ss-hypothesis(a), \( \bu_1=\bLambda_1\T\hbeta-\mathbf{d}_1 \) and the positive definite matrix
\( \W_1=\bLambda_1\T\G\bLambda_1 \) do not depend on \( \hbeta \) and \( \G \). Now
\( \bu=\bLambda\T\hbeta-\mathbf{d}=\mathbf{K}\T\bu_1 \) and \( \W=\mathbf{K}\T\W_1\mathbf{K} \). Let
\( \bv=\mathbf{K}\T(\mathbf{K}\mathbf{K}\T)^{-1}\W_1^{-1}\bu_1 \), so that \( \mathbf{K}\bv=\W_1^{-1}\bu_1 \). Then
\( \W\bv=\mathbf{K}\T\W_1\mathbf{K}\bv=\mathbf{K}\T\bu_1=\bu \), and for every generalized inverse,
\[
\bu\T\W\ginv\bu=\bv\T\W\W\ginv\W\bv=\bv\T\W\bv=(\mathbf{K}\bv)\T\W_1(\mathbf{K}\bv)=\bu_1\T\W_1^{-1}\bu_1 .
\]
The right side involves none of the choices.

*(b).* \( \mathcal S_0 \) is the image of the subspace \( \Null(\bLambda\T) \) under \( \bb\mapsto\X\bb \), so it is a
subspace of \( \C(\X) \). Its dimension is \( r-\dim(\C(\bLambda)\cap\C(\X\T))=r-q \) by @prp-ss-effective-df,
because \( \C(\bLambda)\subseteq\C(\X\T) \). The solutions of \( \bLambda\T\bb=\mathbf{d} \) are \( \bb_0+\Null(\bLambda\T) \),
so the allowed means are \( \X\bb_0+\mathcal S_0 \). For the minimum, write \( \y_*=\y-\X\bb_0 \). As \( \bb \)
ranges over the constraint set, \( \X\bb-\X\bb_0 \) ranges over \( \mathcal S_0 \), so
\[
\min_{\bLambda\T\bb=\mathbf{d}}\norm{\y-\X\bb}^2=\min_{\mathbf{m}\in\mathcal S_0}\norm{\y_*-\mathbf{m}}^2=\norm{(\I-\bP_{\mathcal S_0})\y_*}^2
\]
by @thm-proj-projection-theorem. Also \( \text{SSE}=\norm{(\I-\M)\y}^2=\norm{(\I-\M)\y_*}^2 \), since
\( \X\bb_0\in\C(\X) \). Because \( \mathcal S_0\subseteq\C(\X) \), @eq-proj-extra-ss turns the difference into
\( \norm{(\M-\bP_{\mathcal S_0})\y_*}^2 \). By @thm-ss-hypothesis(c) applied to \( \bLambda_1 \) and
\( \mathbf{d}_1 \), the same difference equals \( \bu_1\T\W_1^{-1}\bu_1 \), which is \( \text{SS}_H \) by (a).

*(c).* \( \Y-\X\bb_0\sim\Normal_n(\X(\bbeta-\bb_0),\sigma^2\I) \) is a normal linear model with mean in
\( \C(\X) \), and by (b) \( F_H \) is its \( F \) statistic for the reduced model space \( \mathcal S_0 \), which has
dimension \( r-q \). @thm-glh-f-test gives \( F_H\sim F(q,n-r,\gamma) \) with
\( \sigma^2\gamma=\norm{(\M-\bP_{\mathcal S_0})(\X\bbeta-\X\bb_0)}^2 \). To evaluate this, note that the identity in
(b) holds for every data vector. Apply it to the vector \( \X\bbeta \), for which \( \bbeta \) itself is a least
squares estimate: \( \sigma^2\gamma=(\bLambda\T\bbeta-\mathbf{d})\T\W\ginv(\bLambda\T\bbeta-\mathbf{d}) \). Finally, by (a)
this equals \( \bu_1(\bbeta)\T\W_1^{-1}\bu_1(\bbeta) \) with \( \bu_1(\bbeta)=\bLambda_1\T\bbeta-\mathbf{d}_1 \). It is zero
iff \( \bLambda_1\T\bbeta=\mathbf{d}_1 \), that is, iff \( \bLambda\T\bbeta=\mathbf{d} \).
:::

Formula @eq-glh-general-F needs only the unrestricted fit. Its ingredients have a statistical
meaning: \( \bLambda\T\hbeta-\mathbf{d} \) is the estimated departure from the hypothesis, and
\( \sigma^2\W \) is its covariance matrix (@thm-opt-sampling(a)). So \( \text{SS}_H/\sigma^2 \) is the
squared Mahalanobis length of the estimated departure, and the noncentrality @eq-glh-general-noncentrality
is the same length for the true departure. Under \( H \) the quadratic form \( \text{SS}_H/\sigma^2 \) is
\( \chi^2(q) \), which is @cor-opt-quadratic. The theorem adds the reduced model, the distribution under
alternatives, and the fact that a redundant \( \bLambda \) can be used as it stands, with any generalized
inverse of \( \W \).

The test space of the hypothesis is the orthogonal complement of \( \mathcal S_0 \) within \( \C(\X) \). By
@thm-proj-constraint-space it is \( \C(\M\bT) \), where \( \bLambda=\X\T\bT \): the hypothesis constrains the
mean vector to be orthogonal to the fitted values of the columns of \( \bT \).

## Least squares under the hypothesis

The minimum in @thm-glh-general-f(b) is attained, and its minimizers have a closed form. The
following extends @thm-ss-restricted, which assumed full rank, to every model matrix.

::: {#prp-glh-restricted-ls}
[Least squares under a testable hypothesis]

Let \( H:\bLambda\T\bbeta=\mathbf{d} \) be testable with \( \bLambda \) of full column rank \( q \), and write
\( \bLambda=\X\T\bT \). Let \( \G \) be a generalized inverse of \( \X\T\X \), \( \hbeta=\G\X\T\y \),
\( \W=\bLambda\T\G\bLambda \) and \( \bu=\bLambda\T\hbeta-\mathbf{d} \). Then \( \W \) is positive definite, and:

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta_H=\hbeta-\G\bLambda\W^{-1}\bu \) satisfies \( \bLambda\T\hbeta_H=\mathbf{d} \) and minimizes
   \( \norm{\y-\X\bb}^2 \) subject to \( \bLambda\T\bb=\mathbf{d} \). The set of all minimizers is
   \( \hbeta_H+\bigl(\Null(\X)\cap\Null(\bLambda\T)\bigr) \).

2. The restricted fitted vector is unique:
   \( \X\hbeta_H=\M\y-\M\bT\W^{-1}\bu=\X\bb_0+\bP_{\mathcal S_0}(\y-\X\bb_0) \), and
   \( \norm{\y-\X\hbeta_H}^2=\text{SSE}+\text{SS}_H \), with \( \text{SS}_H=\bu\T\W^{-1}\bu \).

3. With \( \hat{\boldsymbol{\theta}}=\W^{-1}\bu \), the restricted estimate solves the constrained normal equations (compare @thm-proj-normal-equations)
   \[
\X\T\X\hbeta_H+\bLambda\hat{\boldsymbol{\theta}}=\X\T\y,\qquad \bLambda\T\hbeta_H=\mathbf{d},
\]

   and \( \text{SS}_H=\hat{\boldsymbol{\theta}}\T\W\hat{\boldsymbol{\theta}} \).
:::

:::

::: {.proof}
\( \W=\bT\T\X\G\X\T\bT=\bT\T\M\bT \) by @thm-proj-M-formula. It is symmetric, and positive definite by
@thm-ss-hypothesis(a). Also \( \X\G\bLambda=\X\G\X\T\bT=\M\bT \).

(a) \( \bLambda\T\hbeta_H=\bLambda\T\hbeta-\W\W^{-1}\bu=\mathbf{d} \). For any feasible \( \bb \), @eq-proj-distance-split with
\( \mathbf{u}=\X\bb \) gives \( \norm{\y-\X\bb}^2=\text{SSE}+\norm{\X\hbeta-\X\bb}^2 \). Write \( \bb=\hbeta_H+\bv \) with
\( \bLambda\T\bv=\bzero \). Then \( \X\hbeta-\X\bb=\M\bT\W^{-1}\bu-\X\bv \), and the cross term vanishes:
\[
(\M\bT\W^{-1}\bu)\T\X\bv=\bu\T\W^{-1}\bT\T\X\bv=\bu\T\W^{-1}\bLambda\T\bv=0,
\]
using \( \M\X=\X \). Hence
\( \norm{\y-\X\bb}^2=\text{SSE}+\norm{\M\bT\W^{-1}\bu}^2+\norm{\X\bv}^2 \), which is smallest exactly when
\( \X\bv=\bzero \). The minimizers are \( \hbeta_H+\bv \) with \( \bv\in\Null(\X)\cap\Null(\bLambda\T) \).

(b) \( \X\hbeta_H=\X\hbeta-\X\G\bLambda\W^{-1}\bu=\M\y-\M\bT\W^{-1}\bu \), and both terms are unique by
@thm-ss-hypothesis(a). The minimum value is
\( \text{SSE}+\bu\T\W^{-1}\bT\T\M\bT\W^{-1}\bu=\text{SSE}+\bu\T\W^{-1}\bu \), and \( \bu\T\W^{-1}\bu=\text{SS}_H \)
by @thm-glh-general-f(a). The fitted vector \( \X\hbeta_H \) is the point of \( \X\bb_0+\mathcal S_0 \) nearest to
\( \y \), which is \( \X\bb_0+\bP_{\mathcal S_0}(\y-\X\bb_0) \) by the argument in the proof
of @thm-glh-general-f(b).

(c) \( \X\T\X\hbeta_H=\X\T\X\hbeta-\X\T\X\G\X\T\bT\W^{-1}\bu=\X\T\y-\X\T\M\bT\hat{\boldsymbol{\theta}}=\X\T\y-\bLambda\hat{\boldsymbol{\theta}} \),
because \( \X\T\X\hbeta=\X\T\y \) and \( \X\T\M=\X\T \). The last identity is
\( \bu\T\W^{-1}\bu=\hat{\boldsymbol{\theta}}\T\W\hat{\boldsymbol{\theta}} \).
:::

Part (c) identifies \( \hat{\boldsymbol{\theta}} \) as a vector of Lagrange multipliers (compare @prp-mat-lagrange).
Each multiplier measures how hard the data push against its constraint, and the hypothesis sum of
squares is their squared Mahalanobis length. This is the "Lagrange multiplier" form of the test met
again in [Section 11.4](04-likelihood-ratio.html).

## Nonzero right-hand sides: offsets

When \( \mathbf{d}\ne\bzero \) the hypothesis restricts the mean to the *affine* set \( \X\bb_0+\mathcal S_0 \), not
to a subspace. @thm-glh-general-f(c) shows how to handle it: subtract a known point of the set from the
data. The reduced model becomes
\[
\E(\Y-\X\bb_0)\in\mathcal S_0 ,
\]
a linear model for the adjusted response \( \Y-\X\bb_0 \) with a **known offset** \( \X\bb_0 \). Any solution
\( \bb_0 \) of \( \bLambda\T\bb_0=\mathbf{d} \) may be used. Two solutions differ by an element of \( \Null(\bLambda\T) \),
so their offsets differ by an element of \( \mathcal S_0 \), which the reduced model absorbs. In practice
one fits \( \y-\X\bb_0 \) on a basis of \( \mathcal S_0 \), for example the columns of \( \X\mathbf{U} \) where the columns of
\( \mathbf{U} \) span \( \Null(\bLambda\T) \), and takes the difference of residual sums of squares. The formula
@eq-glh-general-F gives the same number with no second fit.

::: {#exm-glh-regions-rank-deficient}
[Regional hypotheses with an over-parameterized model]

Return to the state murder data, now with an intercept *and* an indicator for each of the four
regions, so that
\( \bbeta=(\mu,\alpha_N,\alpha_M,\alpha_S,\alpha_W,\beta_{\text{pov}},\beta_{\text{sin}},\beta_{\text{urb}})\T \). The
model matrix is \( 50\times8 \) of rank \( 7 \). No region effect is estimable on its own, but every
contrast among them is.

*All four regions equal.* Write the hypothesis with all six pairwise differences
\( \alpha_i-\alpha_j \) as the columns of \( \bLambda \). They are consistent and have rank \( q=3 \). With the
Moore–Penrose inverses of \( \X\T\X \) and of \( \W \), @eq-glh-general-F gives \( \text{SS}_H=14.83 \) and
\( F=2.440 \), exactly the reduced-model test of @exm-glh-region, as @thm-glh-general-f(b)
requires. A second generalized inverse, made by deleting the row and column of \( \alpha_W \) from
\( \X\T\X \) before inverting, gives a different \( \hbeta \) but the same \( \bLambda\T\hbeta \), \( \W \) and \( F \).

*Midwest, South and West equal.* With \( \bLambda \) made of \( \alpha_M-\alpha_S \) and \( \alpha_S-\alpha_W \), the sum of
squares is \( 3.881 \) on \( 2 \) degrees of freedom, \( F=0.958 \), \( p=0.392 \). The
restricted fit of @prp-glh-restricted-ls coincides with the fit of the model in which the three
regions are merged. Under the constraint, the Northeast lies \( 1.379 \) below the common
level of the other three, and the coefficient of single parenthood moves from \( 0.428 \) to
\( 0.441 \).

*A nonzero right-hand side.* The contrast \( \alpha_N-(\alpha_M+\alpha_S+\alpha_W)/3 \), the Northeast against the average
of the other regions, is estimated as \( \hat{\psi}=-1.426 \) with standard error \( 0.595 \). To test
whether it equals \( -1 \), take the offset \( \bb_0=(0,-1,0,\dots,0)\T \) and fit
\( \y-\X\bb_0 \) on a basis of \( \mathcal S_0 \). The drop in residual sum of squares is \( 1.040 \), equal
to \( (\hat{\psi}+1)^2/\blambda\T\G\blambda \) as @eq-glh-general-F predicts, and \( F=0.513 \) with
\( p=0.478 \). Any other solution of \( \blambda\T\bb_0=-1 \) gives the same answer.

*A constraint that is not a hypothesis.* The constraint \( \alpha_N=0 \) is not estimable. Minimizing the
residual sum of squares subject to it gives exactly \( \text{SSE} \): the constrained model spans the whole
of \( \C(\X) \), in agreement with @prp-glh-nonestimable.
:::

```{.python .run #cell-testable-general}
import numpy as np
import statsmodels.api as sm
from scipy import linalg, stats

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
D = np.column_stack([(region == g).astype(float) for g in "NMSW"])   # all four regions
X = np.column_stack([np.ones(n), D, data["poverty"], data["single"], data["urban"]])
r = np.linalg.matrix_rank(X)                                          # 7, not 8

G = np.linalg.pinv(X.T @ X)                     # one generalized inverse of X'X
b = G @ X.T @ y                                 # one least squares solution
s2 = np.sum((y - X @ b) ** 2) / (n - r)

def glh_test(Lam, d):
    """F test of the testable hypothesis Lam' beta = d, via any g-inverse of Lam' G Lam."""
    u = Lam.T @ b - d
    W = Lam.T @ G @ Lam
    q = np.linalg.matrix_rank(Lam)
    ss_h = u @ np.linalg.pinv(W) @ u
    F = ss_h / q / s2
    return ss_h, q, F, stats.f.sf(F, q, n - r)

def contrast(pairs):
    """Columns e_i - e_j on the region coefficients (positions 1..4 = N, M, S, W)."""
    L = np.zeros((8, len(pairs)))
    for k, (i, j) in enumerate(pairs):
        L[1 + i, k], L[1 + j, k] = 1.0, -1.0
    return L

all_pairs = contrast([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])   # rank 3, redundant
three_equal = contrast([(1, 2), (2, 3)])                                 # M = S = W
for name, L in (("four regions equal", all_pairs), ("M = S = W", three_equal)):
    ss_h, q, F, p = glh_test(L, np.zeros(L.shape[1]))
    print(f"{name:20s} SS_H = {ss_h:.3f}  q = {q}  F = {F:.3f}  p = {p:.4f}")
```

The offset computation fits the reduced model directly, as a check on the formula.

```{.python .run #cell-testable-offset}
def sse(Z, target):
    coef, *_ = np.linalg.lstsq(Z, target, rcond=None)
    return np.sum((target - Z @ coef) ** 2)

sse_full = sse(X, y)
lam = np.array([0, 1, -1 / 3, -1 / 3, -1 / 3, 0, 0, 0])   # Northeast minus the other three
d = -1.0                                                  # hypothesized difference
b0 = np.array([0, -1.0, 0, 0, 0, 0, 0, 0])                 # one solution of lam' b0 = d
U = linalg.null_space(lam[None, :])                       # beta - b0 ranges over C(U)
sse_H = sse(X @ U, y - X @ b0)                            # reduced model, offset X b0
est = lam @ b
se = np.sqrt(s2 * lam @ G @ lam)
print(f"estimate {est:.3f} (se {se:.3f}); SSE_H - SSE = {sse_H - sse_full:.3f}")
print(f"formula: {(est - d) ** 2 / (lam @ G @ lam):.3f}")
```

The restricted estimates, from two different generalized inverses:

```{.python .run #cell-testable-restricted}
keep = [0, 1, 2, 3, 5, 6, 7]
G2 = np.zeros((8, 8))
G2[np.ix_(keep, keep)] = np.linalg.inv((X.T @ X)[np.ix_(keep, keep)])

def restricted(Lam, d, Ginv):
    """Least squares under the testable constraint Lam' b = d (Lam of full column rank)."""
    bh = Ginv @ X.T @ y
    W = Lam.T @ Ginv @ Lam
    return bh - Ginv @ Lam @ np.linalg.solve(W, Lam.T @ bh - d)

bH = restricted(three_equal, np.zeros(2), G)
bH2 = restricted(three_equal, np.zeros(2), G2)
print("constraint holds:", np.allclose(three_equal.T @ bH, 0))
print("fitted values agree:", np.allclose(X @ bH, X @ bH2))
```

::: {.warning}
[Hypotheses chosen after seeing the data]

The same Northeast contrast, tested against zero, gives \( t=-2.398 \) and \( p=0.021 \). That
p-value would be valid had the contrast been specified before the data were examined. It was not: it
was suggested by the estimates, after the joint test of all regional differences had failed to reject.
Among the many contrasts one could have tested, the one that looks largest is, by construction, the
one most likely to look significant by chance. [Chapter 13](../ch13-multiplicity/index.html) gives
procedures, notably Scheffé's, that remain valid when the contrast is chosen by looking at the data.
:::

## Exercises

### A. Check your understanding

::: {#exr-glh-testable-oneway}
[A1]

In the one-way layout \( \E(Y_{ij})=\mu+\alpha_i \), \( i=1,2,3 \), with all groups nonempty, decide which of the
following are testable hypotheses: (i) \( \alpha_1=\alpha_2 \); (ii) \( \alpha_1=0 \); (iii)
\( \alpha_1+\alpha_2+\alpha_3=0 \); (iv) \( \mu+\alpha_1=5 \); (v) \( \alpha_1-\alpha_2=1 \), \( \alpha_2-\alpha_3=2 \),
\( \alpha_1-\alpha_3=3 \); (vi) \( \alpha_1-\alpha_2=1 \), \( \alpha_2-\alpha_3=2 \), \( \alpha_1-\alpha_3=4 \). Give \( q \) for those
that are.
:::

::: {.solution}
The estimable functions are the combinations \( c_0\mu+\sum_ic_i\alpha_i \) with \( c_0=\sum_ic_i \) (@exm-proj-ginverse-numeric and [Section 8.2](../ch08-estimability/02-estimable-functions.html)).
(i) is testable, \( q=1 \). (ii) and (iii) are not estimable: their coefficient of \( \mu \) is \( 0 \) while the
\( \alpha \) coefficients sum to \( 1 \) and \( 3 \). They are side conditions (@thm-est-side-conditions). (iv) is testable, \( q=1 \). (v) is
testable with \( q=2 \): the three differences are estimable, have rank \( 2 \), and the right-hand sides
are consistent because \( 1+2=3 \). (vi) is estimable but inconsistent, since
\( (\alpha_1-\alpha_2)+(\alpha_2-\alpha_3)=\alpha_1-\alpha_3 \) would need \( 1+2=4 \). No parameter satisfies it.
:::

### B. Practice

::: {#exr-glh-reformulate}
[B1]

Let \( \bLambda \) be \( p\times m \) and \( \mathbf{A} \) a nonsingular \( m\times m \) matrix. Show that the hypotheses \( \bLambda\T\bbeta=\mathbf{d} \) and
\( (\bLambda\mathbf{A})\T\bbeta=\mathbf{A}\T\mathbf{d} \) are the same, and that @eq-glh-general-F gives the same \( \text{SS}_H \) for
both. Conclude that the test does not depend on how a set of constraints is written.
:::

::: {.solution}
\( (\bLambda\mathbf{A})\T\bb=\mathbf{A}\T\mathbf{d}\iff\mathbf{A}\T(\bLambda\T\bb-\mathbf{d})=\bzero\iff\bLambda\T\bb=\mathbf{d} \), since \( \mathbf{A}\T \)
is nonsingular. So the constraint sets, and by @thm-glh-general-f(b) the sums of squares, agree.
Directly: with \( \bu=\bLambda\T\hbeta-\mathbf{d} \), the new departure is \( \mathbf{A}\T\bu \) and the new middle matrix is
\( \mathbf{A}\T\W\mathbf{A} \), and one generalized inverse of it is \( \mathbf{A}^{-1}\W\ginv\mathbf{A}^{-\top} \), which gives
\( \bu\T\mathbf{A}\mathbf{A}^{-1}\W\ginv\mathbf{A}^{-\top}\mathbf{A}\T\bu=\bu\T\W\ginv\bu \). The choice of generalized inverse does not
matter by @thm-glh-general-f(a).
:::

::: {#exr-glh-restricted-moments}
[B2]

In @prp-glh-restricted-ls, show that \( \E(\X\hbeta_H)=\X\bbeta-\M\bT\W^{-1}(\bLambda\T\bbeta-\mathbf{d}) \) and
\( \Cov(\X\hbeta_H)=\sigma^2\bP_{\mathcal S_0} \), where \( \bP_{\mathcal S_0}=\M-\M\bT\W^{-1}\bT\T\M \). Conclude that
restricted fitted values are unbiased when \( H \) holds, have smaller variance than \( \M\Y \) in every
direction, and are biased when \( H \) fails.
:::

::: {.solution}
\( \X\hbeta_H=\M\Y-\M\bT\W^{-1}(\bT\T\M\Y-\mathbf{d}) \), using \( \bLambda\T\hbeta=\bT\T\M\Y \). Its mean is
\( \X\bbeta-\M\bT\W^{-1}(\bT\T\X\bbeta-\mathbf{d}) \), and \( \bT\T\X\bbeta=\bLambda\T\bbeta \). The random part is
\( (\M-\M\bT\W^{-1}\bT\T\M)\Y \). The matrix in parentheses is the projection onto \( \C(\X) \) minus the
projection onto \( \C(\M\bT) \) (@thm-proj-M-formula, since \( \W=\bT\T\M\bT \)), which is \( \bP_{\mathcal S_0} \) by @thm-proj-constraint-space. It is symmetric idempotent, so the covariance is \( \sigma^2\bP_{\mathcal S_0} \).
The difference \( \sigma^2(\M-\bP_{\mathcal S_0}) \) from \( \Cov(\M\Y) \) is a projection times \( \sigma^2 \), so it is
nonnegative definite. The bias is \( -\M\bT\W^{-1}(\bLambda\T\bbeta-\mathbf{d}) \), zero iff \( H \) holds, because
\( \M\bT \) has full column rank.
:::

::: {#exr-glh-slope-one}
[B3]

In simple regression \( \E(Y_i)=\beta_0+\beta_1x_i \), derive the \( F \) statistic for \( \beta_1=1 \) in two
ways: from @eq-glh-general-F, and by regressing \( y_i-x_i \) on an intercept and comparing with the full
fit. Show that both give \( F=(\hat{\beta}_1-1)^2S_{xx}/s^2 \).
:::

::: {#exr-glh-mean-specified}
[B4]

Show that the hypothesis \( \X\bbeta=\X\bbeta_0 \), with \( \bbeta_0 \) known, is testable, that \( q=r \), and that
its \( F \) statistic is \( \{\norm{\M(\y-\X\bbeta_0)}^2/r\}/s^2 \). What is the reduced model?
:::

### C. Going deeper

::: {#exr-glh-nonestimable-formula}
[C1]

In the one-way layout with two groups of sizes \( 2 \) and \( 3 \), take \( \bLambda=(0,1,0)\T \), so that the
"hypothesis" is \( \alpha_1=0 \). Evaluate \( (\bLambda\T\hbeta)^2/(\bLambda\T\G\bLambda) \) for two different
generalized inverses \( \G \) of \( \X\T\X \) and show that the values differ. Explain which step of the proof of
@thm-glh-general-f fails.
:::

::: {#exr-glh-eliminate}
[C2]

Let \( \X \) have full column rank and write \( \bLambda\T=[\bLambda_1\T,\bLambda_2\T] \) with \( \bLambda_2\T \) a nonsingular
\( q\times q \) block, after permuting the coefficients so that \( \bbeta=(\bbeta_1\T,\bbeta_2\T)\T \) conformably.
Show that \( \bLambda\T\bbeta=\bzero \) is equivalent to \( \bbeta_2=-\bLambda_2^{-\top}\bLambda_1\T\bbeta_1 \), that the reduced
model matrix is \( \X_1-\X_2\bLambda_2^{-\top}\bLambda_1\T \), and that this matrix has full column rank \( p-q \).
This is how the reduced model is fitted when no offset or null-space basis is at hand.
:::
