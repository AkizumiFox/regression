# When ordinary least squares is best

The analysis in @exm-ggm-elnino left a puzzle: the generalized and the ordinary
estimates of the trend agreed to three decimals and the seasonal coefficients to
four, yet their standard errors disagreed by a factor of more than two in one
direction and by half as much again in the other. Both halves have one
explanation, the subject of this section: for a large class of designs ordinary
least squares is *already* the best linear unbiased estimator, or very nearly so,
while its reported standard error is wrong all the same.

Throughout, \( \V \) is nonnegative definite and may be singular, \( \M \) is the
orthogonal projection onto \( \C(\X) \), and "\( \A\Y \) is a BLUE of \( \X\bbeta \)" means
that \( \boldsymbol{\uprho}\T\A\Y \) is best linear unbiased for \( \boldsymbol{\uprho}\T\X\bbeta \) at every
\( \boldsymbol{\uprho}\in\Real^n \).

## A criterion for a best estimator

Everything in this section, and most of the next, rests on one small lemma, the
covariance version of @exr-opt-zyskind: a linear estimator is
best exactly when it is uncorrelated with every linear statistic that estimates
zero. What is new is that \( \V \) may be singular, that the criterion is
read off one vector \( \mathbf{c} \) at a time, and that uniqueness holds on
\( \C(\X:\V) \), not merely with probability one.

::: {#lem-ggm-blue-criterion}
[Criterion for a BLUE]

Assume @eq-ggm-model with \( \V \) nonnegative definite. Let \( \mathbf{c}\in\Real^n \) and
\( \blambda=\X\T\mathbf{c} \), so that \( \mathbf{c}\T\Y \) is a linear unbiased estimator of
\( \blambda\T\bbeta \).

::: {.enumerate options="label=(\alph*)"}
1. \( \mathbf{c}\T\Y \) is a BLUE of \( \blambda\T\bbeta \) if and only if \( \V\mathbf{c}\in\C(\X) \).

2. Any two BLUEs of \( \blambda\T\bbeta \) are equal with probability one. More
   precisely, if \( \mathbf{c}\T\Y \) and \( \mathbf{b}\T\Y \) are both BLUEs then
   \( \mathbf{c}\T\y=\mathbf{b}\T\y \) for every \( \y\in\C(\X:\V) \).

3. A matrix estimator \( \A\Y \) with \( \A\X=\X \) is a BLUE of \( \X\bbeta \) if and only
   if \( \A\V(\I-\M)=\mathbf{0} \), equivalently \( \C(\V\A\T)\subseteq\C(\X) \).
:::
:::

::: {.proof}
(a) Every linear unbiased estimator of \( \blambda\T\bbeta \) is \( (\mathbf{c}+\mathbf{d})\T\Y \) with
\( \X\T\mathbf{d}=\bzero \), that is, with \( \mathbf{d}\in\C(\X)\perpc \). For real \( t \),
\[
\Var\bigl((\mathbf{c}+t\mathbf{d})\T\Y\bigr)=\sigma^2\bigl(\mathbf{c}\T\V\mathbf{c}+2t\,\mathbf{c}\T\V\mathbf{d}+t^2\mathbf{d}\T\V\mathbf{d}\bigr),
\]
a quadratic in \( t \) with nonnegative leading coefficient. If \( \mathbf{c}\T\Y \) is best,
this quadratic is minimized at \( t=0 \) for every such \( \mathbf{d} \), which forces
\( \mathbf{c}\T\V\mathbf{d}=0 \). Conversely, if \( \mathbf{c}\T\V\mathbf{d}=0 \) for every \( \mathbf{d}\perp\C(\X) \), then
\( \Var((\mathbf{c}+\mathbf{d})\T\Y)=\Var(\mathbf{c}\T\Y)+\sigma^2\mathbf{d}\T\V\mathbf{d}\ge\Var(\mathbf{c}\T\Y) \). So
\( \mathbf{c}\T\Y \) is best iff \( \V\mathbf{c}\perp\C(\X)\perpc \), which is \( \V\mathbf{c}\in\C(\X) \).

(b) Let \( \mathbf{c}\T\Y \) and \( \mathbf{b}\T\Y \) be BLUEs, with common variance \( v \) and
covariance \( \kappa \). The average \( \frac12(\mathbf{c}+\mathbf{b})\T\Y \) is unbiased, so
\( v\le\frac14(v+v+2\kappa) \), that is \( v\le\kappa \). Therefore
\[
\Var\bigl((\mathbf{c}-\mathbf{b})\T\Y\bigr)=2v-2\kappa\le0,
\]
so it is zero, and \( (\mathbf{c}-\mathbf{b})\T\Y \) has mean zero and zero variance: it vanishes
with probability one. For the deterministic statement, zero variance means
\( (\mathbf{c}-\mathbf{b})\T\V(\mathbf{c}-\mathbf{b})=0 \), hence \( \V(\mathbf{c}-\mathbf{b})=\bzero \) because \( \V \) is
nonnegative definite, so \( \mathbf{c} \) and \( \mathbf{b} \) agree on \( \C(\V) \); and
\( \X\T(\mathbf{c}-\mathbf{b})=\bzero \) makes them agree on \( \C(\X) \). They therefore agree on
\( \C(\X)+\C(\V)=\C(\X:\V) \).

(c) Apply (a) to \( \mathbf{c}=\A\T\boldsymbol{\uprho} \), which satisfies
\( \X\T\A\T\boldsymbol{\uprho}=(\A\X)\T\boldsymbol{\uprho}=\X\T\boldsymbol{\uprho} \), so \( \mathbf{c}\T\Y \) is unbiased for
\( \boldsymbol{\uprho}\T\X\bbeta \). The condition \( \V\A\T\boldsymbol{\uprho}\in\C(\X) \) for every \( \boldsymbol{\uprho} \) is
\( \C(\V\A\T)\subseteq\C(\X) \), which says \( (\I-\M)\V\A\T=\mathbf{0} \); transposing gives
\( \A\V(\I-\M)=\mathbf{0} \).
:::

## The characterization

::: {#thm-ggm-ols-blue}
[When ordinary least squares is best]

Assume @eq-ggm-model with \( \V \) nonnegative definite, and let \( \W \) be any
\( n\times(n-r) \) matrix with \( \C(\W)=\C(\X)\perpc \).

::: {.enumerate options="label=(\alph*)"}
1. The following are equivalent.

   ::: {.enumerate options="label=(\roman*)"}
   1. \( \M\Y \) is a BLUE of \( \X\bbeta \).
   2. \( \M\V(\I-\M)=\mathbf{0} \), equivalently \( \M\V=\V\M \).
   3. \( \C(\V\X)\subseteq\C(\X) \).
   4. \( \V=\X\A\X\T+\W\B\W\T \) for some nonnegative definite \( \A \) and \( \B \).
   5. \( \C(\X) \) has an orthonormal basis of eigenvectors of \( \V \).
   :::

2. For \( \boldsymbol{\uprho}\in\Real^n \), the estimator \( \boldsymbol{\uprho}\T\M\Y \) of \( \boldsymbol{\uprho}\T\X\bbeta \) is a BLUE
   if and only if \( \V\M\boldsymbol{\uprho}\in\C(\X) \). When \( \rank(\X)=p \) and
   \( \blambda\in\Real^p \), this reads \( \V\X(\X\T\X)^{-1}\blambda\in\C(\X) \).

3. If \( \V \) is positive definite, the conditions of (a) are also equivalent to
   \( \M\y=\A_{\V}\y \) for every \( \y \), and then ordinary and generalized least
   squares give the same fit, the same estimate of every estimable function, and
   the same \( \text{SSE}_{\V} \).
:::
:::

::: {.proof}
For positive definite \( \V \) the chain (i) \( \Leftrightarrow \) (ii) \( \Leftrightarrow \) (iii) is
@exr-opt-zyskind; the argument below needs no invertibility.

(a) (i) \( \Leftrightarrow \) (ii). Take \( \A=\M \) in @lem-ggm-blue-criterion(c), which
applies because \( \M\X=\X \). For the second form of (ii): if \( \M\V(\I-\M)=\mathbf{0} \)
then \( \M\V=\M\V\M \), and \( \M\V\M \) is symmetric, so \( \M\V \) is symmetric and equals
\( (\M\V)\T=\V\M \). Conversely \( \M\V=\V\M \) gives \( \M\V(\I-\M)=\V\M-\V\M=\mathbf{0} \).

(ii) \( \Leftrightarrow \) (iii). If \( \M\V=\V\M \) then \( \V\X=\V\M\X=\M\V\X\in\C(\M)=\C(\X) \).
Conversely, if \( \C(\V\X)\subseteq\C(\X) \) then \( \M\V\X=\V\X \). Since
\( \C(\M)=\C(\X) \), every column of \( \M \) is a combination of columns of \( \X \), so
\( \M\V\M=\V\M \). The left side is symmetric, so \( \V\M=(\V\M)\T=\M\V \).

(ii) \( \Rightarrow \) (iv). Split the identity as \( \I=\M+(\I-\M) \):
\[
\begin{aligned}
\V&=\M\V\M+\M\V(\I-\M)\\
&\qquad{}+(\I-\M)\V\M+(\I-\M)\V(\I-\M)\\
&=\M\V\M+(\I-\M)\V(\I-\M),
\end{aligned}
\]
the cross terms vanishing by (ii). Now \( \M=\X(\X\T\X)\ginv\X\T \) by
@thm-proj-M-formula, so \( \M\V\M=\X\A\X\T \) with
\( \A=(\X\T\X)\ginv\X\T\V\X(\X\T\X)\ginv \), which is nonnegative definite because
\( \V \) is. Likewise \( \I-\M=\W(\W\T\W)^{-1}\W\T \), since \( \W \) has full column rank
and spans \( \C(\X)\perpc \), so \( (\I-\M)\V(\I-\M)=\W\B\W\T \) with
\( \B=(\W\T\W)^{-1}\W\T\V\W(\W\T\W)^{-1} \), also nonnegative definite.

(iv) \( \Rightarrow \) (iii). \( \W\T\X=\mathbf{0} \), so \( \V\X=\X\A\X\T\X\in\C(\X) \).

(ii) \( \Leftrightarrow \) (v). Condition (iii) says \( \C(\X) \) is invariant under \( \V \).
If it is, the restriction of \( \V \) to \( \C(\X) \) is a symmetric operator on that
subspace, so by @thm-mat-spectral it has an orthonormal eigenbasis there, and
these are eigenvectors of \( \V \). Conversely, if \( \C(\X) \) is spanned by
eigenvectors of \( \V \) then \( \V \) maps each of them, and hence \( \C(\X) \), into
\( \C(\X) \).

(b) By @lem-ggm-blue-criterion(a) with \( \mathbf{c}=\M\boldsymbol{\uprho} \), which satisfies
\( \X\T\M\boldsymbol{\uprho}=\X\T\boldsymbol{\uprho} \), so that \( \boldsymbol{\uprho}\T\M\Y \) is unbiased for
\( \boldsymbol{\uprho}\T\X\bbeta \). For the full-rank form, \( \blambda\T\hbeta=\boldsymbol{\uprho}\T\M\Y \) with
\( \boldsymbol{\uprho}=\X(\X\T\X)^{-1}\blambda \), and then \( \M\boldsymbol{\uprho}=\boldsymbol{\uprho} \).

(c) With \( \V \) positive definite, (iii) is Kruskal's condition, and
@thm-proj-kruskal says it holds iff \( \M\y=\A_{\V}\y \) for every \( \y \). Equal fits
give equal estimates of every estimable function by
@thm-proj-invariant-functions, and equal residuals, hence equal
\( \text{SSE}_{\V} \).
:::

Form (iv) is Rao's **simple covariance structure**: the error splits as
\( \be=\be_1+\be_2 \) with \( \be_1\in\C(\X) \), \( \be_2\perp\C(\X) \) and the two parts
uncorrelated (@exr-ggm-simple-structure). Ordinary least squares is best exactly
when the noise splits along the model space in this way, for then the residuals
carry no usable information about \( \bbeta \). Form (v) is the one to check when
the eigenstructure of \( \V \) is known.

## Balanced layouts and a common shock

The most important case in practice is a balanced layout whose errors share a
shock at some level of the design.

::: {#exm-ggm-balanced}
[Equicorrelation in a balanced layout]

Let the \( n=Gm \) observations fall into \( G \) groups of equal size \( m \), let \( \Z \)
be the \( n\times G \) matrix of group indicators, and let
\( \V=(1-\rho)\I+\rho\Z\Z\T \), the covariance produced by a group-level shock.
Since \( \Z\Z\T \) is \( m \) times the orthogonal projection onto the
group-constant vectors, condition (iii) holds iff \( \C(\X) \) is closed under
replacing a vector by its group means — which happens whenever \( \C(\X) \) is
spanned by group-constant vectors together with vectors whose group means
vanish, exactly what a balanced design gives through its between-group and
within-group contrasts.

In the El Niño data of @exm-ggm-elnino the groups are years. If the seasonal
part of the model is the full set of eleven month effects, every column of \( \X \)
splits into a year-level part and a part whose year means vanish, both again in
\( \C(\X) \), so the condition holds *exactly*: the two fits agree to nine
decimals, whatever \( \rho \). Their standard errors do not. At \( \rho=0.5 \) the
root of the expected reported variance \( \E(s^2)\{(\X\T\X)^{-1}\}_{22} \) is only
\( 0.391 \) of the true standard error of the trend.
:::

The same argument covers randomized complete block designs with random blocks,
and any design whose every regressor is constant within a block or sums to zero
in each block. The analysis of variance tables of
[Part IV](../ch15-anova-subspaces/index.html) therefore survive random block
effects as *point estimates*; only their error terms need rethinking, and
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) develops the
split-plot case, where condition (iv) becomes a decomposition into error strata.
Incomplete block designs break the condition: between-block comparisons then
carry treatment information that a within-block analysis discards, and recovering
it — Yates's recovery of interblock information — is generalized least squares.

## When some contrasts are best and others are not

Part (b) of @thm-ggm-ols-blue is finer than part (a): even when \( \M\Y \) is not
best for \( \X\bbeta \) as a whole, particular functions still can be. The
mechanism is easiest to see when the balance of @exm-ggm-balanced is broken.

::: {#exm-ggm-clusters}
[Unequal clusters]

Take \( G \) clusters of sizes \( m_1,\dots,m_G \), a covariate \( \bw \) centred inside
every cluster, and the model \( \X=[\bone,\bw] \) with
\( \V=(1-\rho)\I+\rho\Z\Z\T \). Because \( \bone\T\bw=0 \), the vector attached to the
slope is \( \mathbf{c}=\bw/\bw\T\bw \), and \( \Z\Z\T\bw=\bzero \) because \( \bw \) has zero mean
in each cluster; so \( \V\mathbf{c}=(1-\rho)\mathbf{c}\in\C(\X) \) and the within-cluster slope is
estimated best by ordinary least squares, for every \( \rho \) and every set of
sizes.

The vector attached to the level is \( \mathbf{c}=\bone/n \), and
\( \Z\Z\T\bone \) is the vector whose entry for an observation in cluster \( k \) is
\( m_k \). That step function lies in \( \C(\X)=\spn\{\bone,\bw\} \) only if all the
\( m_k \) are equal, since \( \bw \) is not constant within clusters. So with unequal
clusters the overall level is *not* estimated best by ordinary least squares;
@eq-ggm-cluster-weights gives the estimator that is.

With sizes \( 2, 3, 6 \) and \( \rho=0.5 \) the weights of
@eq-ggm-cluster-weights are \( 1.333 \), \( 1.500 \) and \( 1.714 \), summing to
\( 4.5476 \), so the generalized estimate of the level has variance
\( 0.2199\,\sigma^2 \) against \( 0.2479\,\sigma^2 \) for the sample mean: an
efficiency of \( 0.887 \), well above the Kantorovich bound \( 0.437 \)
that @prp-ggm-ols-loss(c) allows for the eigenvalues \( 3.500 \) and
\( 0.500 \) of this \( \V \). The slope loses nothing. Meanwhile the usual
standard errors are wrong for both: the root of the expected reported variance is
\( 0.563 \) of the true standard error for the level and \( 1.314 \) of it
for the slope.
:::

::: {when-format="html"}
![**Figure 31.2.1.** Efficiency of the ordinary least squares estimate of the overall
level, for three cluster layouts with equicorrelation \( \rho \) inside clusters (solid).
The balanced layout loses nothing at any \( \rho \). The dotted curves are the
Kantorovich bound of @thm-dep-efficiency for the same \( \V \), far below the
efficiency the design actually attains.](unequal_clusters.svg){#fig-ggm-clusters width=62%}
:::

::: {when-format="pdf"}
![Efficiency of the ordinary least squares estimate of the overall
level, for three cluster layouts with equicorrelation \( \rho \) inside clusters (solid).
The balanced layout loses nothing at any \( \rho \). The dotted curves are the
Kantorovich bound of @thm-dep-efficiency for the same \( \V \), far below the
efficiency the design actually attains.](unequal_clusters.pdf){width=62%}
:::

```{.python .run #cell-unequal-clusters-setup}
import numpy as np

sizes = [2, 3, 6]
w = np.array([-1.0, 1.0, -1.0, 0.0, 1.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
n = len(w)
Z = np.zeros((n, len(sizes)))                       # cluster indicators
Z[np.arange(n), np.repeat(np.arange(len(sizes)), sizes)] = 1.0
X = np.column_stack([np.ones(n), w])
rho = 0.5
V = (1 - rho) * np.eye(n) + rho * Z @ Z.T

Ginv = np.linalg.inv(X.T @ X)
Cov_ols = Ginv @ X.T @ V @ X @ Ginv                 # sigma^2 = 1 throughout
Cov_gls = np.linalg.inv(X.T @ np.linalg.solve(V, X))
print("efficiency of OLS:  level %.4f   slope %.4f"
      % (Cov_gls[0, 0] / Cov_ols[0, 0], Cov_gls[1, 1] / Cov_ols[1, 1]))
```

## How much is lost, and what is actually broken

::: {#prp-ggm-ols-loss}
[The cost of ordinary least squares]

Assume @eq-ggm-model with \( \V \) positive definite and \( \rank(\X)=p \). Write
\( \hbeta=(\X\T\X)^{-1}\X\T\Y \) and let \( \lambda_1\ge\dots\ge\lambda_n>0 \) be the
eigenvalues of \( \V \). For \( \blambda\ne\bzero \) define the efficiency
\( e(\blambda)=\Var(\blambda\T\hbeta_{\V})/\Var(\blambda\T\hbeta) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \Cov(\hbeta)-\Cov(\hbeta_{\V}) \) is nonnegative definite, so
   \( 0<e(\blambda)\le1 \) for every \( \blambda \).

2. \( e(\blambda)=1 \) if and only if \( \blambda\T\hbeta=\blambda\T\hbeta_{\V} \) with
   probability one, which by @thm-ggm-ols-blue(b) happens iff
   \( \V\X(\X\T\X)^{-1}\blambda\in\C(\X) \).

3. \( e(\blambda)\ge4\lambda_1\lambda_n/(\lambda_1+\lambda_n)^2 \) for every
   \( \blambda \) and every \( \X \), and the bound is attained.
:::
:::

::: {.proof}
(a) Both \( \blambda\T\hbeta \) and \( \blambda\T\hbeta_{\V} \) are linear unbiased
estimators of \( \blambda\T\bbeta \), and the latter is best by @cor-opt-aitken, so
\( \Var(\blambda\T\hbeta)\ge\Var(\blambda\T\hbeta_{\V})>0 \). Taking
\( \blambda \) over \( \Real^p \) gives \( \blambda\T\bigl(\Cov(\hbeta)-\Cov(\hbeta_{\V})\bigr)\blambda\ge0 \).

(b) If \( e(\blambda)=1 \) then \( \blambda\T\hbeta \) attains the minimum variance and
is therefore a BLUE, so it equals \( \blambda\T\hbeta_{\V} \) with probability one by
@lem-ggm-blue-criterion(b); the converse is immediate. The condition is part (b)
of @thm-ggm-ols-blue.

(c) This is @thm-dep-efficiency, proved in
[Section 19.2](../ch19-theory-of-departures/02-wrong-covariance.html) by a matrix form
of Kantorovich's inequality, together with the design given there that attains it.
:::

Part (c) is the worst case over all designs and all \( \blambda \). A determinant
version is also known: if \( p\le n/2 \), then
\[
\frac{\det\Cov(\hbeta_{\V})}{\det\Cov(\hbeta)}\ \ge\ \prod_{j=1}^{p}\frac{4\lambda_j\lambda_{n-j+1}}{(\lambda_j+\lambda_{n-j+1})^2},
\]
with equality for a design built from the corresponding pairs of eigenvectors.
This is the Bloomfield–Watson–Knott inequality, proved independently by
Bloomfield and Watson (1975) and by Knott (1975) in back-to-back papers. We state
it without proof; the argument needs more matrix theory than this book develops.

In practice the loss is far smaller, because real designs are not adversarial.
The El Niño model with three harmonics is typical: Kruskal's condition fails,
with \( \max_{s,t}\lvert\{(\I-\M)\V\X\}_{st}\rvert=0.126 \), yet the efficiency of
the trend coefficient is \( 0.9996 \) and that of every seasonal coefficient is
\( 1.0000 \) to four decimals. The intercept is exactly best, because
\( \bone \) is orthogonal to every other column, so \( \mathbf{c}=\bone/n \) and
\( \V\bone=\{1+11\rho\}\bone \) is a multiple of \( \bone \). What is badly wrong is
the standard error: the root of the expected value of
\( s^2\{(\X\T\X)^{-1}\}_{jj} \) is \( 0.390 \) of the true standard error for
the trend and \( 1.406 \) of it for the first cosine, and
[Figure 31.2.2](#fig-ggm-efficiency) shows both ratios running away as \( \rho \)
grows while the efficiency never falls below \( 0.99 \).

::: {when-format="html"}
![**Figure 31.2.2.** The El Niño model with three harmonics, as a function of the
intraclass correlation \( \rho \). The efficiency of the ordinary least squares
trend estimate stays above \( 0.99 \), while the ratio of the reported to the true
standard error falls to a third for the trend and rises past three for a seasonal
coefficient.](elnino_efficiency.svg){#fig-ggm-efficiency width=62%}
:::

::: {when-format="pdf"}
![The El Niño model with three harmonics, as a function of the
intraclass correlation \( \rho \). The efficiency of the ordinary least squares
trend estimate stays above \( 0.99 \), while the ratio of the reported to the true
standard error falls to a third for the trend and rises past three for a seasonal
coefficient.](elnino_efficiency.pdf){width=62%}
:::

::: {.idea}
Correlated errors damage the *standard error* far more than the *estimate*.
@thm-dep-efficiency bounds the efficiency below by a quantity depending only
on the extreme eigenvalues of \( \V \), and for ordinary designs the loss is a few
per cent. Nothing bounds the error in \( s^2(\X\T\X)^{-1} \), which can be too small
or too large by any factor. If only one of the two problems can be fixed, fix the
standard error.
:::

That is what the sandwich estimator of
[Section 21.4](../ch21-nonnormality-heteroscedasticity-serial/04-sandwich-estimators.html)
does. It leaves \( \hbeta \) alone and replaces \( s^2(\X\T\X)^{-1} \) by
\( (\X\T\X)^{-1}\bigl(\sum_k\X_k\T\he_k\he_k\T\X_k\bigr)(\X\T\X)^{-1} \), summed over
units that can be treated as independent; for clustered data the unit is a
cluster, and the argument of @thm-het-sandwich applies cluster by cluster
rather than observation by observation. We do not prove that extension:
@thm-het-sandwich assumes independent observations, consistency of the cluster
form as the number of clusters grows is due to Liang and Zeger (1986), and
@exr-bs-wild-cluster derives the same matrix as a bootstrap covariance. It
recovers no efficiency and needs many independent units, as
[Section 31.4](04-feasible-gls.html) shows; the two repairs are complementary.

## Exercises

### A. Check your understanding

::: {#exr-ggm-eigen-columns}
[A1]

Show directly from @thm-ggm-ols-blue(a)(v) that ordinary least squares is best
when the columns of \( \X \) are eigenvectors of \( \V \), and use this to recover the
equicorrelation result of @exr-proj-equicorrelation.
:::

::: {.solution}
If each column of \( \X \) is an eigenvector of \( \V \), then \( \C(\X) \) is spanned by
eigenvectors, which is (v). For \( \V=(1-\rho)\I+\rho\bone\bone\T \) the vector
\( \bone \) is an eigenvector with eigenvalue \( 1+(n-1)\rho \) and every vector
orthogonal to \( \bone \) is an eigenvector with eigenvalue \( 1-\rho \). A model
containing an intercept has \( \C(\X)=\spn\{\bone\}\dirsum(\C(\X)\cap\bone\perpc) \),
a span of eigenvectors, so (v) holds.
:::

::: {#exr-ggm-multiple-of-i}
[A2]

Show that \( \M\V=\V\M \) for *every* \( \X \) if and only if \( \V \) is a multiple of
\( \I \).
:::

### B. Practice

::: {#exr-ggm-simple-structure}
[B1]

Assume the simple covariance structure \( \V=\X\A\X\T+\W\B\W\T \) of
@thm-ggm-ols-blue(a)(iv). Show that \( \Cov(\M\Y,(\I-\M)\Y)=\mathbf{0} \), and that the
error can be written \( \be=\be_1+\be_2 \) with \( \be_1\in\C(\X) \),
\( \be_2\in\C(\X)\perpc \) and \( \Cov(\be_1,\be_2)=\mathbf{0} \) with probability one.
:::

::: {.solution}
\( \Cov(\M\Y,(\I-\M)\Y)=\sigma^2\M\V(\I-\M)=\mathbf{0} \) by (ii). Put \( \be_1=\M\be \) and
\( \be_2=(\I-\M)\be \); these lie in the stated spaces, add to \( \be \), and have
cross-covariance \( \sigma^2\M\V(\I-\M)=\mathbf{0} \). Their covariance matrices are
\( \sigma^2\X\A\X\T \) and \( \sigma^2\W\B\W\T \).
:::

::: {#exr-ggm-gls-unknown-theta}
[B2]

Suppose \( \Cov(\Y)=\bSigma(\boldsymbol{\uptheta}) \) with \( \boldsymbol{\uptheta} \) unknown, and that
\( \C(\bSigma(\boldsymbol{\uptheta})\X)\subseteq\C(\X) \) for every \( \boldsymbol{\uptheta} \) in the parameter set.
Show that \( \M\Y \) is a BLUE for every \( \boldsymbol{\uptheta} \), so that the best estimator is known
even though the covariance is not. Give an example with two unknown variance
components.
:::

::: {#exr-ggm-cluster-efficiency}
[B3]

In @exm-ggm-clusters with two clusters of sizes \( m_1 \) and \( m_2 \), show that the
efficiency of the sample mean as an estimate of the level is
\[
e=\frac{(m_1+m_2)^2}{\bigl\{(1-\rho)(m_1+m_2)+\rho(m_1^2+m_2^2)\bigr\}(w_1+w_2)},
\qquad w_k=\frac{m_k}{1-\rho+m_k\rho}.
\]
Show that \( e=1 \) when \( m_1=m_2 \), and find the limit as \( \rho\to1 \) with
\( m_1=1 \) and \( m_2=m \).
:::

::: {.solution}
\( \Var(\bar Y)=\sigma^2\bone\T\V\bone/n^2 \) with
\( \bone\T\V\bone=(1-\rho)n+\rho\sum_km_k^2 \), and
\( \Var(\hat\mu_{\V})=\sigma^2/\sum_kw_k \) by @eq-ggm-cluster-weights; divide. If
\( m_1=m_2=m \) both variances are \( \sigma^2(1-\rho+m\rho)/(2m) \), so \( e=1 \). With
\( m_1=1 \), \( m_2=m \) and \( \rho\to1 \), \( w_k\to1 \) and \( e\to(1+m)^2/\{2(1+m^2)\} \),
which falls to \( 1/2 \) as \( m \) grows.
:::

### C. Going deeper

::: {#exr-ggm-blue-subspace}
[C1]

Let \( \rank(\X)=p \) and put
\( \mathcal B=\{\blambda\in\Real^p:\V\X(\X\T\X)^{-1}\blambda\in\C(\X)\} \), the set of
coefficient functions that ordinary least squares estimates best. Show that
\( \mathcal B \) is a subspace, that \( \mathcal B=\Real^p \) iff
\( \C(\V\X)\subseteq\C(\X) \), and that
\( \dim\mathcal B=p-\rank\bigl((\I-\M)\V\X(\X\T\X)^{-1}\bigr) \). Evaluate it
for @exm-ggm-clusters.
:::

::: {.solution}
A vector \( \blambda \) lies in \( \mathcal B \) iff
\( (\I-\M)\V\X(\X\T\X)^{-1}\blambda=\bzero \), so \( \mathcal B \) is the null space of
that matrix, hence a subspace, and @thm-mat-rank-nullity gives its dimension. It
is all of \( \Real^p \) iff the matrix is zero, which is
\( \C(\V\X)\subseteq\C(\X) \) because \( (\X\T\X)^{-1} \) is nonsingular. In
@exm-ggm-clusters, \( p=2 \) and the map kills the slope direction but not the
level direction, so \( \dim\mathcal B=1 \): the functions estimated best are exactly
the multiples of the slope.
:::

::: {#exr-ggm-two-eigenvalues}
[C2]

Suppose \( \V \) has only two distinct eigenvalues, \( \lambda_1 \) with multiplicity
\( k \) and \( \lambda_2 \) with multiplicity \( n-k \), and let \( \mathcal E_1 \) be the
eigenspace of \( \lambda_1 \). Show that \( \M\V=\V\M \) iff
\( \C(\X)=(\C(\X)\cap\mathcal E_1)\dirsum(\C(\X)\cap\mathcal E_1\perpc) \). Deduce that
when this splitting holds every \( \blambda\T\hbeta \) has efficiency \( 1 \), and that
when it fails there is a design attaining the bound
\( 4\lambda_1\lambda_2/(\lambda_1+\lambda_2)^2 \).
:::
