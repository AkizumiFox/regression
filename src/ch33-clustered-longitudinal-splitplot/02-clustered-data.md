# Clustered data and compound symmetry

Whole plots are one instance of a pattern that is everywhere in applied work:
observations arrive in groups — households, schools, hospitals, litters, survey
clusters — and within a group resemble one another more than they resemble
observations from elsewhere. The simplest description gives every pair within a
group the same correlation and every pair across groups none. What that does to
least squares is sharp and exactly computable: it inflates the variance of some
coefficients and *deflates* that of others.

## The equicorrelated cluster model

::: {#def-cls-cluster}
[Clustered errors with compound symmetry]

The \( n \) observations fall into \( G \) **clusters**, cluster \( g \) containing
\( m_g \) of them, with \( n=\sum_gm_g \). The model is \( \Y=\X\bbeta+\be \) with
\( \E\be=\bzero \) and
\[
\Cov(\be)=\V=\diag(\bSigma_1,\dots,\bSigma_G),\qquad
\bSigma_g=\sigma^2\bigl\{(1-\rho)\I_{m_g}+\rho\bone\bone\T\bigr\},
\]{#eq-cls-compound}

a block diagonal matrix with one block per cluster. This covariance is called
**compound symmetry** or **equicorrelation**, and \( \rho \) is the **intraclass
correlation**. The matrix is positive definite iff
\( -1/(\max_gm_g-1)<\rho<1 \). The clusters are **balanced** if \( m_g=m \) for all
\( g \).
:::

Compound symmetry is what a single shared random effect produces: if
\( \varepsilon_i=u_{g(i)}+e_i \) with \( u_g \) and \( e_i \) uncorrelated, of variances
\( \sigma_u^2 \) and \( \sigma_e^2 \), then @eq-cls-compound holds with
\( \sigma^2=\sigma_u^2+\sigma_e^2 \) and \( \rho=\sigma_u^2/\sigma^2 \) — the one-way
random effects model of @prp-mix-oneway and @exr-ss-random-oneway, whose power to
inflate an \( F \) statistic was met at @exr-aov-intraclass. So read, \( \rho \)
is the fraction of variance attributable to the cluster and is nonnegative;
negative intraclass correlations occur — a fixed quota within a cluster,
competition among littermates — but not from a random-effect representation.

The results below take the clusters balanced, \( m_g=m \); the clean statements
fail otherwise (@exr-cls-unequal-clusters), balance being what makes the two
strata of @lem-cls-strata eigenspaces of \( \V \).

## Which coefficients suffer

Call \( \bv \) **cluster-constant** if it is constant on every cluster and
**cluster-centred** if its entries sum to zero within every cluster. With
\( \bP_W \) the orthogonal projection replacing each entry by its cluster mean, the
two conditions are \( \bP_W\bv=\bv \) and \( \bP_W\bv=\bzero \), and every vector is
the sum of one of each.

::: {#prp-cls-compound}
[Compound symmetry and the design effect]

Let @eq-cls-compound hold with \( m_g=m \) for all \( g \), let
\( \lambda_w=\sigma^2\{1+(m-1)\rho\} \) and \( \lambda_s=\sigma^2(1-\rho) \), and let
\( \rank(\X)=p \).

::: {.enumerate options="label=(\alph*)"}
1. \( \V=\lambda_w\bP_W+\lambda_s(\I-\bP_W) \). Its eigenvalues are \( \lambda_w \) with
   multiplicity \( G \) and \( \lambda_s \) with multiplicity \( n-G \), and
   \[
\V^{-1}=\frac{1}{\lambda_s}\Bigl(\I-\frac{m\rho}{1+(m-1)\rho}\bP_W\Bigr).
\]

2. \( \C(\V\X)\subseteq\C(\X) \) iff \( \C(\bP_W\X)\subseteq\C(\X) \). In particular, if
   \( \X=[\X_b,\X_w] \) with every column of \( \X_b \) cluster-constant and every column
   of \( \X_w \) cluster-centred, then least squares is generalized least squares and
   every least squares estimate of an estimable function is best linear unbiased.

3. Let \( \blambda\T\bbeta \) be estimable and let
   \( \mathbf{c}=\X(\X\T\X)\ginv\blambda \), so that
   \( \blambda\T\hbeta=\mathbf{c}\T\y \). If \( \mathbf{c} \) is cluster-constant then
   \[
\Var(\blambda\T\hbeta)=\{1+(m-1)\rho\}\cdot\sigma^2\norm{\mathbf{c}}^2,
\]
   and if \( \mathbf{c} \) is cluster-centred then
   \( \Var(\blambda\T\hbeta)=(1-\rho)\,\sigma^2\norm{\mathbf{c}}^2 \). In both cases
   \( \sigma^2\norm{\mathbf{c}}^2 \) is the variance the estimate would have under
   \( \Cov(\be)=\sigma^2\I \); the multiplier is the **design effect**. In either
   case \( \blambda\T\hbeta \) is best linear unbiased even if (b) fails.

4. Suppose \( \C(\X)=\mathcal S_b\dirsum\mathcal S_w \) with
   \( \mathcal S_b\subseteq\C(\bP_W) \) of dimension \( p_b \) and
   \( \mathcal S_w\subseteq\C(\bP_W)\perpc \) of dimension \( p_w=p-p_b \). Then the usual
   residual mean square has
   \[
\E(s^2)=\frac{\lambda_w(G-p_b)+\lambda_s(n-G-p_w)}{n-p},
\]{#eq-cls-naive-s2}

   and for any \( \mathbf{c} \) as in (c) the usual variance estimate
   \( v=s^2\blambda\T(\X\T\X)\ginv\blambda \) satisfies
   \( \E(v)/\Var(\blambda\T\hbeta)=\E(s^2)/\lambda \), with \( \lambda=\lambda_w \) or
   \( \lambda_s \) according to the stratum of \( \mathbf{c} \).
:::

:::

::: {.proof}
(a) Within one cluster, \( \bone\bone\T=m\bar{\mathbf{J}}_m \), so the block is
\( \sigma^2\{(1-\rho)\I_m+m\rho\bar{\mathbf{J}}_m\}
 =\lambda_s(\I_m-\bar{\mathbf{J}}_m)+\lambda_w\bar{\mathbf{J}}_m \), using
\( \lambda_s+m\rho\sigma^2=\lambda_w \). Assembling the blocks gives the stated form,
since \( \bP_W \) is block diagonal with blocks \( \bar{\mathbf{J}}_m \). A matrix
\( \lambda_w\bP+\lambda_s(\I-\bP) \) with \( \bP \) an orthogonal projection of rank
\( G \) has the two eigenvalues with those multiplicities and inverse
\( \lambda_w^{-1}\bP+\lambda_s^{-1}(\I-\bP)
 =\lambda_s^{-1}\{\I-(1-\lambda_s/\lambda_w)\bP\} \), and
\( 1-\lambda_s/\lambda_w=m\rho/\{1+(m-1)\rho\} \).

(b) By (a), \( \V\X=\lambda_s\X+(\lambda_w-\lambda_s)\bP_W\X \), so
\( \C(\V\X)\subseteq\C(\X) \) iff \( \C(\bP_W\X)\subseteq\C(\X) \) (the case
\( \lambda_w=\lambda_s \), that is \( \rho=0 \), being trivial). Under the stated
partition, \( \bP_W\X_b=\X_b \) and \( \bP_W\X_w=\bzero \), so
\( \C(\bP_W\X)=\C(\X_b)\subseteq\C(\X) \). The conclusion is @thm-proj-kruskal
together with @cor-opt-aitken, as in @thm-cls-splitplot(a); the general statement
is @thm-ggm-ols-blue.

(c) \( \Var(\mathbf{c}\T\Y)=\mathbf{c}\T\V\mathbf{c} \). If \( \bP_W\mathbf{c}=\mathbf{c} \) then
\( \V\mathbf{c}=\lambda_w\mathbf{c} \) and the variance is \( \lambda_w\norm{\mathbf{c}}^2 \); if
\( \bP_W\mathbf{c}=\bzero \) then it is \( \lambda_s\norm{\mathbf{c}}^2 \). Under
\( \Cov(\be)=\sigma^2\I \) it would be \( \sigma^2\norm{\mathbf{c}}^2 \), and dividing gives
the two multipliers. For the last claim, write \( \lambda \) for \( \lambda_w \) or
\( \lambda_s \) according to the stratum of \( \mathbf{c} \), so that
\( \V\mathbf{c}=\lambda\mathbf{c} \) in either case, and let \( \mathbf{a}\T\Y \) be
any linear unbiased estimate of \( \blambda\T\bbeta \). Then
\( \X\T(\mathbf{a}-\mathbf{c})=\bzero \), so \( \mathbf{a}=\mathbf{c}+\mathbf{d} \)
with \( \mathbf{d}\perp\C(\X) \) and in particular
\( \mathbf{c}\T\mathbf{d}=0 \). The cross term
\( \mathbf{c}\T\V\mathbf{d}=\lambda\mathbf{c}\T\mathbf{d} \) vanishes and
\[
\mathbf{a}\T\V\mathbf{a}=\lambda\norm{\mathbf{c}}^2+\mathbf{d}\T\V\mathbf{d}
 \ \ge\ \lambda\norm{\mathbf{c}}^2,
\]
with equality iff \( \mathbf{d}\T\V\mathbf{d}=0 \), that is \( \mathbf{d}=\bzero \),
\( \V \) being positive definite. So \( \mathbf{c}\T\Y \) is the unique best linear
unbiased estimate. This is the criterion behind @thm-ggm-ols-blue: a linear
estimate \( \mathbf{c}\T\Y \) of an estimable function is best as soon as
\( \V\mathbf{c}\in\C(\X) \).

(d) By @thm-dep-covariance(b), whose part (b) uses \( \X \) only through \( \M \)
and \( \rank(\X)=p \) and so holds verbatim here,
\( \E(s^2)=\tr\{(\I-\M)\V\}/(n-p) \) with \( \M \) the
projection onto \( \C(\X) \). Now
\( \tr\{(\I-\M)\V\}=\lambda_w\tr\{(\I-\M)\bP_W\}+\lambda_s\tr\{(\I-\M)(\I-\bP_W)\} \).
Because \( \M=\bP_{\mathcal S_b}+\bP_{\mathcal S_w} \) with
\( \bP_{\mathcal S_b}\bP_W=\bP_{\mathcal S_b} \) and \( \bP_{\mathcal S_w}\bP_W=\bzero \),
\[
\tr\{(\I-\M)\bP_W\}=\tr(\bP_W)-\tr(\bP_{\mathcal S_b})=G-p_b,
\]
and similarly \( \tr\{(\I-\M)(\I-\bP_W)\}=(n-G)-p_w \). This
gives @eq-cls-naive-s2. Finally \( v=s^2\norm{\mathbf{c}}^2 \) because
\( \norm{\mathbf{c}}^2=\blambda\T(\X\T\X)\ginv\X\T\X(\X\T\X)\ginv\blambda
 =\blambda\T(\X\T\X)\ginv\blambda \) for estimable \( \blambda \), so
\( \E(v)/\Var(\mathbf{c}\T\Y)=\E(s^2)\norm{\mathbf{c}}^2/(\lambda\norm{\mathbf{c}}^2) \).
:::

The design effect cuts both ways. A regressor measured at the cluster level — the
school's curriculum, the village's water supply — is cluster-constant, and its
variance is multiplied by \( 1+(m-1)\rho \): with \( m=50 \) and \( \rho=0.05 \)
the factor is \( 3.45 \), so \( 5000 \) pupils in \( 100 \) schools carry the
information of about \( 1450 \) independent pupils, the **effective sample size**
\( n/\{1+(m-1)\rho\} \). Since the factor grows linearly in \( m \), extra
clusters buy far more than extra members per cluster (@exr-cls-effective-n). A
cluster-centred regressor is in the opposite position: its design effect
\( 1-\rho \) is less than one, the shared component cancelling out of
within-cluster comparisons — the subplot stratum again.

::: {.idea}
Clustering does not simply lose information; it moves it from between-cluster to
within-cluster comparisons. A design that can ask its question within clusters
gains, one that must ask it between clusters loses, in proportion to the cluster
size.
:::

Part (d) says what the usual output does with all this: \( s^2 \) averages \( \V \)
over the residual space while the variance of an estimate reads \( \V \) along one
direction, so reported and true variance stand in the ratio \( \E(s^2)/\lambda \),
the same for every direction in a stratum.

::: {#exm-cls-clusters}
[Forty clusters of five]

Take \( G=40 \) clusters of \( m=5 \), with \( \rho=0.4 \), \( \sigma^2=1 \), one
cluster-level and one within-cluster regressor besides the intercept, so
\( p_b=2 \) and \( p_w=1 \). Then \( \lambda_w=2.6 \) and \( \lambda_s=0.6 \), and the
design effects are exactly \( 2.60 \) and \( 0.60 \), as the listing
confirms from the full \( 200\times200 \) covariance matrix. By @eq-cls-naive-s2,
\( \E(s^2)=\{2.6(38)+0.6(159)\}/197=0.9858 \), so the reported variance of the
cluster-level coefficient is \( 0.379 \) times the truth and that of the
within-cluster coefficient \( 1.643 \) times the truth. A nominal
\( 95\% \) interval for the first is about \( 1.62 \) times too short, and
the simulation below finds its coverage to be \( 0.759 \).
:::

```{.python .run #cell-design-effect-structure}
import numpy as np

m, rho, sigma2 = 5, 0.4, 1.0                         # cluster size, intraclass correlation
S = sigma2 * ((1 - rho) * np.eye(m) + rho * np.ones((m, m)))
w = np.linalg.eigvalsh(S)                            # one large, m-1 equal small eigenvalues
print("eigenvalues of the cluster covariance:", w.round(4))
print("large:", sigma2 * (1 + (m - 1) * rho), "  small:", sigma2 * (1 - rho))
```

```{.python .run #cell-design-effect-deff}
def design(G, m, seed):
    """Cluster indicator, a cluster-level regressor and a within-cluster regressor."""
    rng = np.random.default_rng(seed)
    g = np.repeat(np.arange(G), m)
    xb = np.repeat(rng.normal(0, 1, G), m)           # constant inside each cluster
    xw = rng.normal(0, 1, G * m)
    xw = xw - np.repeat(xw.reshape(G, m).mean(axis=1), m)   # sums to zero inside each cluster
    X = np.column_stack([np.ones(G * m), xb, xw])
    return g, X

G = 40
g, X = design(G, m, 2024)
V = np.kron(np.eye(G), S)                            # block diagonal, one block per cluster
XtXinv = np.linalg.inv(X.T @ X)
cov_true = XtXinv @ X.T @ V @ X @ XtXinv             # the real covariance of the OLS estimate
cov_iid = sigma2 * XtXinv                            # what the usual formula reports
deff = np.diag(cov_true) / np.diag(cov_iid)
print("design effects (intercept, between, within):", deff.round(4))
print("predicted:", [1 + (m - 1) * rho, 1 + (m - 1) * rho, 1 - rho])
```

## Standard errors that survive the model

Generalized least squares with the right \( \V \) fixes everything, and by
@prp-cls-compound(b) is not even needed for the point estimate when the regressors
split into cluster-constant and cluster-centred parts. But \( \rho \) is rarely
known and compound symmetry is itself a guess. With many clusters the alternative
is to keep ordinary least squares and estimate its covariance from the clusters
themselves: the sandwich idea of
[Section 21.4](../ch21-nonnormality-heteroscedasticity-serial/04-sandwich-estimators.html)
with the cluster in place of the observation. Writing \( \X_g \) and \( \he_g \) for
the rows and least squares residuals of cluster \( g \), the **cluster-robust**
covariance estimator is, as in @exr-bs-wild-cluster,
\[
\hat{\bU}_{\text{CR}}=(\X\T\X)^{-1}\Bigl(\sum_{g=1}^G\X_g\T\he_g\he_g\T\X_g\Bigr)(\X\T\X)^{-1}.
\]{#eq-cls-cluster-sandwich}

Nothing about \( \V \) is assumed except that different clusters are uncorrelated:
within a cluster the covariance may be arbitrary, and may differ from cluster to
cluster — far weaker than @eq-cls-compound.

::: {.remark}
[Consistency of the cluster-robust estimator]

Let the clusters be independent, cluster \( g \) contributing rows \( \X_g \) and
errors \( \be_g \) of mean zero and covariance \( \bSigma_g \), with \( p \) fixed,
and assume (K1) \( m_g\le m_{\max} \) for all \( g \), and
\( \E\norm{\be_g}^4\le K \); (K2) \( \norm{\x_{(i)}}\le C \) for all \( i \); (K3)
\( n^{-1}\X\T\X\to\A \), positive definite. Write
\( \bU=\Cov(\hbeta)=(\X\T\X)^{-1}(\sum_g\X_g\T\bSigma_g\X_g)(\X\T\X)^{-1} \). Then
\( n(\hat{\bU}_{\text{CR}}-\bU)\to\mathbf{0} \) in probability as \( G\to\infty \);
and if also \( n^{-1}\sum_g\X_g\T\bSigma_g\X_g\to\B \), positive definite, then
\( \sqrt n(\hbeta-\bbeta)\to\Normal_p(\bzero,\A^{-1}\B\A^{-1}) \) and the Wald
statistic of a hypothesis \( \bLambda\T\bbeta \) of rank \( q \), formed with
\( \hat{\bU}_{\text{CR}} \), tends to \( \chi^2(q) \).

The proof is that of @thm-het-sandwich with the cluster in place of the
observation: the cluster score \( \mathbf{s}_g=\X_g\T\be_g \) replaces
\( \x_{(i)}e_i \), every sum has \( G \) independent terms instead of \( n \), and
each bound goes through with \( m_{\max}C \) in place of \( C \) (@exr-cls-crve-steps).
:::

Note where the asymptotics live: \( G\to\infty \), with cluster sizes bounded; the
number of *observations* is irrelevant. A study of a million patients in eight
hospitals has eight independent pieces of information about anything varying at
the hospital level, and no covariance estimator manufactures more.

```{.python .run #cell-design-effect-robust}
def cluster_sandwich(X, resid, g):
    """The cluster-robust covariance estimate: one outer product per cluster."""
    XtXinv = np.linalg.inv(X.T @ X)
    meat = np.zeros((X.shape[1], X.shape[1]))
    for k in np.unique(g):
        s = X[g == k].T @ resid[g == k]               # the cluster's score contribution
        meat += np.outer(s, s)
    return XtXinv @ meat @ XtXinv

rng_one = np.random.default_rng(55)
n, p = X.shape
e_one = (np.linalg.cholesky(S) @ rng_one.standard_normal((m, G))).T.reshape(-1)
y_one = X[:, 1] + e_one                              # the between coefficient is 1
b_one = XtXinv @ X.T @ y_one
res_one = y_one - X @ b_one
s2_one = res_one @ res_one / (n - p)
se_naive_one = np.sqrt(s2_one * XtXinv[1, 1])
se_cr_one = np.sqrt(cluster_sandwich(X, res_one, g)[1, 1] * G / (G - 1))
se_true_one = np.sqrt(cov_true[1, 1])
print(f"between coefficient {b_one[1]:.4f}: naive se {se_naive_one:.4f}, "
      f"cluster-robust {se_cr_one:.4f}, true sd {se_true_one:.4f}")
```

On the data generated here the cluster-level coefficient is \( 1.0376 \)
against a true value of \( 1 \), with a usual standard error of \( 0.0856 \),
a cluster-robust one of \( 0.1600 \) (the listing scales the estimator by the
simple factor \( G/(G-1) \)) and an exact standard deviation of
\( 0.1239 \): the naive interval is far too short, the robust one
conservative on this sample.

[Figure 33.2.1](#fig-cls-design-effect) shows the pattern over repeated samples,
with the fuller correction \( G/(G-1)\cdot(n-1)/(n-p) \). The naive interval
covers about three quarters of the time — \( 0.757 \) to
\( 0.781 \) across the six designs — whatever the number of clusters,
its error being a bias rather than a variance. The cluster-robust interval is
badly short with \( 8 \) clusters (coverage \( 0.874 \)), still short
with \( 20 \) (\( 0.917 \)) and acceptable by \( 160 \)
(\( 0.951 \)); generalized least squares with the correct covariance is
right throughout.

::: {when-format="html"}
![**Figure 33.2.1.** (a) The design effect of @prp-cls-compound(c): inflation
\( 1+(m-1)\rho \) for a cluster-level regressor, deflation \( 1-\rho \) for a
within-cluster one. (b) Simulated coverage of nominal \( 95\% \) intervals for a
cluster-level coefficient, with \( m=5 \) and \( \rho=0.4 \).](design_effect.svg){#fig-cls-design-effect width=100%}
:::

::: {when-format="pdf"}
![(a) The design effect of @prp-cls-compound(c): inflation
\( 1+(m-1)\rho \) for a cluster-level regressor, deflation \( 1-\rho \) for a
within-cluster one. (b) Simulated coverage of nominal \( 95\% \) intervals for a
cluster-level coefficient, with \( m=5 \) and \( \rho=0.4 \).](design_effect.pdf){width=100%}
:::

::: {.warning}
With few clusters the cluster-robust standard error is not a safe default: its
bias is downward and substantial, and referring the statistic to \( t(G-1) \), the
usual remedy and the one used in the figure, is not enough. Better are a bias
correction of the residuals, a bootstrap over whole
clusters (@exr-bs-wild-cluster), or a covariance model fitted by restricted
maximum likelihood, which does not need \( G \) large. Twenty clusters is a common
rule of thumb; it is a floor, not a target.
:::

## Exercises

### A. Check your understanding

::: {#exr-cls-icc-range}
[A1]

Show from @prp-cls-compound(a) that the compound symmetric matrix is positive
definite iff \( -1/(m-1)<\rho<1 \), and that the random-effects representation
\( \varepsilon_i=u_{g(i)}+e_i \) forces \( \rho\ge0 \). Give a mechanism that produces
\( \rho<0 \).
:::

::: {#exr-cls-effective-n}
[A2]

A trial can afford \( n=1200 \) pupils and must randomize schools, not pupils, to
the two arms. The intraclass correlation of the outcome is \( \rho=0.10 \). Compare
the effective sample size for \( 12 \) schools of \( 100 \) pupils, \( 40 \) schools of
\( 30 \), and \( 120 \) schools of \( 10 \).
:::

::: {.solution}
The effective sample size is \( n/\{1+(m-1)\rho\} \): \( 1200/10.9=110 \),
\( 1200/3.9=308 \) and \( 1200/1.9=632 \) pupils. As \( m \) grows the effective size
tends to \( n/(m\rho)=G/\rho \), which depends only on the number of schools. For a
school-level treatment the number of schools is nearly the whole story, and
pupils within a school are cheap but nearly redundant.
:::

### B. Practice

::: {#exr-cls-two-strata-decomposition}
[B1]

Let \( \X=[\X_b,\X_w] \) as in @prp-cls-compound(b), with \( \bone\in\C(\X_b) \). Show
that \( \M=\M_b+\M_w \) with \( \M_b \) and \( \M_w \) the projections onto \( \C(\X_b) \)
and \( \C(\X_w) \), that the residual sum of squares splits into a between-cluster
and a within-cluster part, and that the two corresponding mean squares are
unbiased for \( \lambda_w \) and \( \lambda_s \).
:::

::: {.solution}
Any column of \( \X_b \) is fixed by \( \bP_W \) and any column of \( \X_w \) is
annihilated by it, so \( \C(\X_b)\perp\C(\X_w) \) and @thm-proj-sum gives
\( \M=\M_b+\M_w \). Then \( \I-\M=(\bP_W-\M_b)+(\I-\bP_W-\M_w) \), a sum of two
orthogonal projections of ranks \( G-p_b \) and \( n-G-p_w \). By
@prp-cls-compound(a), \( \V \) acts as \( \lambda_w \) on the range of the first and
\( \lambda_s \) on the range of the second, so @thm-rv-quadform-mean gives the
expectations \( \lambda_w(G-p_b) \) and \( \lambda_s(n-G-p_w) \). Dividing by the
ranks gives the unbiased estimates, of which @eq-cls-naive-s2 is the pooled
average.
:::

::: {#exr-cls-unequal-clusters}
[B2]

Let the clusters have unequal sizes \( m_g \). Show that \( \V \) of @eq-cls-compound
has eigenvalues \( \sigma^2\{1+(m_g-1)\rho\} \), one per cluster, together with
\( \sigma^2(1-\rho) \) with multiplicity \( n-G \). Deduce that \( \bone \) is an
eigenvector only if all clusters have the same size, so that the intercept alone
breaks @prp-cls-compound(b).
:::

### C. Going deeper

::: {#exr-cls-crve-steps}
[C1]

Carry out in detail the three steps of the cluster-sandwich argument above,
with the constants: show that (K1)–(K3) give \( \norm{\X_g}\le m_{\max}C \) and
\( \E\norm{\X_g\T\be_g}^2\le m_{\max}^2C^2K^{1/2} \); that
\( \Var(\mathbf{a}\T\X_g\T(\be_g\be_g\T-\bSigma_g)\X_g\mathbf{a}) \) is bounded
uniformly in \( g \); and that
\( n^{-1}\sum_g\norm{\mathbf{t}_g}^2\le m_{\max}^2C^2(\hbeta-\bbeta)\T\A_n(\hbeta-\bbeta) \).
Where exactly does the number of clusters, rather than of observations, enter?
:::

::: {#exr-cls-crve-bias}
[C2]

Show that \( \E(\X_g\T\he_g\he_g\T\X_g)\neq\X_g\T\bSigma_g\X_g \) in general, by
writing \( \he_g=\be_g-\X_g(\X\T\X)^{-1}\X\T\be \) and expanding, and show that the
leading correction is of order \( p/G \). Why is the bias downward, and what does
this suggest about the factor \( G/(G-1)\cdot(n-1)/(n-p) \) used in the coverage
simulation behind [Figure 33.2.1](#fig-cls-design-effect)?
:::
