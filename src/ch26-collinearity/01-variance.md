# What collinearity does to variance

The columns of a model matrix are **collinear** when a nontrivial linear combination of them is
zero, and **nearly collinear** when a combination with coefficients of reasonable size is a short
vector. The exact case is rank deficiency, settled in [Chapter 8](../ch08-estimability/index.html).
Here \( \X \) has full column rank and every linear function of \( \bbeta \) is estimable, yet some are
estimated so poorly that in practice they have not been estimated at all; "collinearity" below
means this near case.

[Section 19.6](../ch19-theory-of-departures/06-collinearity.html) wrote the variance of
\( \mathbf{a}\T\hbeta \) as a sum over the eigenvectors of \( \X\T\X \) (@prp-dep-collinear-directions). This section
says which functions collinearity damages, links them to estimability, and separates the statistical problem from the numerical one of
[Chapter 10](../ch10-computation/index.html).

## Canonical coordinates

Throughout the section \( \X \) is \( n\times p \) of rank \( p \), \( \E(\Y)=\X\bbeta \) and
\( \Cov(\Y)=\sigma^2\I \). Write the spectral decomposition (@thm-mat-spectral) as
\[
\X\T\X=\sum_{\ell=1}^p\lambda_\ell\bv_\ell\bv_\ell\T,\qquad
\lambda_1\ge\lambda_2\ge\dots\ge\lambda_p>0,
\]
with orthonormal eigenvectors \( \bv_1,\dots,\bv_p \), and put \( \gamma_\ell=\bv_\ell\T\bbeta \). The numbers
\( \gamma_\ell \) are the coordinates of \( \bbeta \) in the eigenvector basis, the **canonical
coordinates** of the regression. Since \( \X\bv_\ell \) has squared length \( \lambda_\ell \), the
mean vector is
\[
\X\bbeta=\sum_{\ell=1}^p\gamma_\ell\,\X\bv_\ell=\sum_{\ell=1}^p\sqrt{\lambda_\ell}\,\gamma_\ell\,\bu_\ell,
\qquad \bu_\ell=\X\bv_\ell/\sqrt{\lambda_\ell},
\]{#eq-col-canonical}

and the vectors \( \bu_\ell \) are the left singular vectors of \( \X \) (@thm-mat-svd). A small
\( \lambda_\ell \) means that a unit change in \( \gamma_\ell \) moves the mean vector only by
\( \sqrt{\lambda_\ell} \). The data then carry little information about \( \gamma_\ell \), and the
following theorem makes this precise.

::: {#thm-col-variance}
[Variances in canonical coordinates]

Under the assumptions above, let \( \hat{\gamma}_\ell=\bv_\ell\T\hbeta \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\gamma}_\ell=\bu_\ell\T\Y/\sqrt{\lambda_\ell} \), and \( \hbeta=\sum_\ell\hat{\gamma}_\ell\bv_\ell \). The
   \( \hat{\gamma}_\ell \) are uncorrelated, with \( \Var(\hat{\gamma}_\ell)=\sigma^2/\lambda_\ell \).

2. For every nonzero \( \mathbf{a}\in\Real^p \),
   \[
   \frac{\sigma^2\norm{\mathbf{a}}^2}{\lambda_1}\le\Var(\mathbf{a}\T\hbeta)\le\frac{\sigma^2\norm{\mathbf{a}}^2}{\lambda_p},
   \]
   with equality on the left at \( \mathbf{a}=\bv_1 \) and on the right at \( \mathbf{a}=\bv_p \). The ratio of
   the largest to the smallest variance of \( \mathbf{a}\T\hbeta \) over unit vectors \( \mathbf{a} \) is \( \kappa(\X)^2 \).

3. If \( \mathbf{a}\in\spn\{\bv_1,\dots,\bv_k\} \), then \( \Var(\mathbf{a}\T\hbeta)\le\sigma^2\norm{\mathbf{a}}^2/\lambda_k \):
   the eigenvalues \( \lambda_{k+1},\dots,\lambda_p \) do not enter. For every \( \mathbf{a} \),
   \( \Var(\mathbf{a}\T\hbeta)\ge\sigma^2(\mathbf{a}\T\bv_p)^2/\lambda_p \).

4. \( \E\norm{\hbeta-\bbeta}^2=\sigma^2\sum_\ell1/\lambda_\ell\ge\sigma^2/\lambda_p \), while, by
   @prp-dep-collinear-directions(b), \( \E\norm{\X\hbeta-\X\bbeta}^2=\sigma^2p \) whatever the eigenvalues.

5. For each column, \( \norm{(\I-\M_{(j)})\x_j}^2\ge\lambda_p \), where \( \M_{(j)} \) projects onto the span of
   the other columns. Hence \( \Var(\hat{\beta}_j)\le\sigma^2/\lambda_p \).
:::

:::

::: {.proof}
(a) By @thm-lm-moments, \( \hbeta=(\X\T\X)^{-1}\X\T\Y \) with covariance
\( \sigma^2(\X\T\X)^{-1}=\sigma^2\sum_\ell\lambda_\ell^{-1}\bv_\ell\bv_\ell\T \). Hence
\( \bv_\ell\T\hbeta=\lambda_\ell^{-1}(\X\bv_\ell)\T\Y=\bu_\ell\T\Y/\sqrt{\lambda_\ell} \), and
\( \Cov(\bv_\ell\T\hbeta,\bv_m\T\hbeta)=\sigma^2\bv_\ell\T(\X\T\X)^{-1}\bv_m \), which is \( \sigma^2/\lambda_\ell \) if
\( \ell=m \) and zero otherwise. The expansion of \( \hbeta \) holds because the \( \bv_\ell \) form an
orthonormal basis.

(b) Write \( c_\ell=\mathbf{a}\T\bv_\ell \), so that \( \sum_\ell c_\ell^2=\norm{\mathbf{a}}^2 \). By
@prp-dep-collinear-directions(a), \( \Var(\mathbf{a}\T\hbeta)=\sigma^2\sum_\ell c_\ell^2/\lambda_\ell \), which lies between
\( \sigma^2\norm{\mathbf{a}}^2/\lambda_1 \) and \( \sigma^2\norm{\mathbf{a}}^2/\lambda_p \), with the stated equality cases. The
singular values of \( \X \) are \( \sqrt{\lambda_\ell} \) (@thm-mat-svd), so the ratio of the extremes is
\( \lambda_1/\lambda_p=\kappa(\X)^2 \) by @def-mat-condition-number.

(c) If \( \mathbf{a}\in\spn\{\bv_1,\dots,\bv_k\} \), then \( c_\ell=0 \) for \( \ell>k \) and
\( \sum_{\ell\le k}c_\ell^2/\lambda_\ell\le\norm{\mathbf{a}}^2/\lambda_k \). The lower bound keeps only the term
\( \ell=p \).

(d) \( \E\norm{\hbeta-\bbeta}^2=\tr\Cov(\hbeta)=\sigma^2\tr(\X\T\X)^{-1} \), and the eigenvalues of
\( (\X\T\X)^{-1} \) are the \( 1/\lambda_\ell \). The statement about fitted values is
@prp-dep-collinear-directions(b) multiplied by \( n \).

(e) \( (\I-\M_{(j)})\x_j=\x_j-\X_{(j)}\mathbf{c} \) for some \( \mathbf{c} \), which is \( \X\bb \) for a vector
\( \bb \) with \( b_j=1 \). By @thm-mat-extremal-rayleigh, \( \norm{\X\bb}^2\ge\lambda_p\norm{\bb}^2\ge\lambda_p \).
The variance bound follows from @prp-opt-orthogonal-design.
:::

Part (b) is the precise sense in which \( \kappa(\X)^2 \) measures collinearity: it is the ratio of the
worst to the best precision among linear functions of the same length. Part (e) is the Frisch–Waugh–Lovell picture of
@thm-proj-fwl: a coefficient is imprecise when its column's part orthogonal to the others is short,
and no such part is shorter than \( \sqrt{\lambda_p} \).

::: {.idea}
Collinearity does not make least squares imprecise in general. It makes it imprecise in a few
directions of coefficient space, those of the eigenvectors with small eigenvalues, and precise
in the others. Whether it matters depends on whether the question asked points along a weak
direction.
:::

::: {#exm-col-two-regressors}
[Two regressors that move together]

Draw \( n=30 \) values \( x_{i1} \) from a standard normal distribution and set
\( x_{i2}=x_{i1}+0.15z_i \) with independent standard normal \( z_i \). The sample correlation is
\( 0.9951 \). Centre both regressors and fit \( \E(Y)=\beta_0+\beta_1x_1+\beta_2x_2 \), so that the
intercept is orthogonal to the slopes and the slope block of \( \X\T\X \) carries all the collinearity.
Its eigenvalues are \( 70.05 \) and \( 0.1718 \). The weak eigenvector is close to
\( (1,-1)\T/\sqrt2 \) (its first coordinate is \( 0.718 \) in absolute value), so the weak canonical
coordinate is essentially the difference \( \beta_1-\beta_2 \). In units of \( \sigma \), the standard
deviations are

| Function | Standard deviation |
|:---|:---:|
| \( \hat{\beta}_1 \) | \( 1.734 \) |
| \( \hat{\beta}_2 \) | \( 1.682 \) |
| \( \hat{\beta}_1+\hat{\beta}_2 \) | \( 0.177 \) |
| \( \hat{\beta}_1-\hat{\beta}_2 \) | \( 3.412 \) |

Had \( x_1 \) been orthogonal to \( x_2 \) with the same spread, the first entry would have been \( 0.171 \), smaller
by the factor \( 10.1=\sqrt{102.4} \), the root of the variance inflation factor of Section 26.2. The sum of the slopes is estimated about as well as a slope in an orthogonal design, the
difference worst of all.
The ratio of the largest to the smallest variance over unit vectors is \( 407.8 \), the squared condition
number of the centred slope columns (@thm-col-variance(b)).

[Figure 26.1.1](#fig-col-canonical) shows what this means for a single data set. With the design fixed and
true slopes \( \beta_1=\beta_2=1 \), the estimates from \( 2000 \) replicate responses
lie along the line \( \beta_1+\beta_2=2 \). The estimate of \( \beta_2 \) ranges from \( -4.535 \) to \( 7.359 \) and is
negative in \( 26.1 \) per cent of the replicates, while \( \hat{\beta}_1+\hat{\beta}_2 \) stays between \( 1.440 \) and \( 2.639 \).
Every data set says clearly that the two regressors together matter, yet a quarter of them give
\( \beta_2 \) the wrong sign, and \( 54.1 \) per cent give at least one of the two slopes the wrong sign.
:::

::: {when-format="html"}
![**Figure 26.1.1.** (a) The design of @exm-col-two-regressors. (b) Estimates \( (\hat{\beta}_1,\hat{\beta}_2) \) from \( 2000 \)
replicate responses, true value marked: long in the weak direction \( (1,-1) \), thin along \( (1,1) \); dashed,
\( \beta_1+\beta_2=2 \).](canonical_directions.svg){#fig-col-canonical width=100%}
:::

::: {when-format="pdf"}
![(a) The design of @exm-col-two-regressors. (b) Estimates \( (\hat{\beta}_1,\hat{\beta}_2) \) from \( 2000 \)
replicate responses, true value marked: long in the weak direction \( (1,-1) \), thin along \( (1,1) \); dashed,
\( \beta_1+\beta_2=2 \).](canonical_directions.pdf){width=100%}
:::

```{.python .run #cell-variance-design}
import numpy as np
rng = np.random.default_rng(2601)
n, sigma = 30, 1.0
x1 = rng.normal(size=n)
x2 = x1 + 0.15 * rng.normal(size=n)                  # x2 nearly repeats x1
x1c, x2c = x1 - x1.mean(), x2 - x2.mean()             # centring makes the intercept orthogonal
X = np.column_stack([np.ones(n), x1c, x2c])
G = np.linalg.inv(X.T @ X)                            # Cov(beta_hat) / sigma^2

lam, V = np.linalg.eigh(X[:, 1:].T @ X[:, 1:])        # canonical coordinates of the slopes
print("correlation of x1 and x2:", round(np.corrcoef(x1, x2)[0, 1], 4))
print("eigenvalues:", np.round(lam, 3), " weak direction:", np.round(V[:, 0], 3))

def sd(a):
    """Standard deviation of a' beta_hat, in units of sigma."""
    return np.sqrt(a @ G @ a)

print("sd of b1, b2   :", round(sd(np.array([0, 1, 0])), 3), round(sd(np.array([0, 0, 1])), 3))
print("sd of b1 + b2  :", round(sd(np.array([0, 1, 1])), 3))
print("sd of b1 - b2  :", round(sd(np.array([0, 1, -1])), 3))
```

```{.python .run #cell-variance-simulate}
beta = np.array([0.0, 1.0, 1.0])                      # equal true slopes
reps = 2000
Y = X @ beta + sigma * rng.normal(size=(reps, n))     # 2000 replicate responses, same design
B = np.linalg.solve(X.T @ X, X.T @ Y.T).T             # one row of estimates per replicate
print("share of replicates with b2 < 0:", np.mean(B[:, 2] < 0))
print("share with b1 < 0 or b2 < 0:", np.mean((B[:, 1:] < 0).any(axis=1)))
print("range of b1 + b2:", np.round([B[:, 1:].sum(1).min(), B[:, 1:].sum(1).max()], 3))
```

## The symptoms

The classical symptoms all follow from @thm-col-variance: large standard errors for single
coefficients alongside small ones for some combinations; coefficients of implausible size or sign;
instability under deletion of a few cases; and a significant \( F \) test for a group in which no single
\( t \) test is significant, because the group statistic is governed by the best combination of its
members (@prp-glh-max-t). None of them means the model is wrong.

## Near-estimability

When the near-dependence is a small perturbation of an exact one, the functions that stay precise
are exactly those estimable in the limiting, rank-deficient design.

::: {#prp-col-near-estimable}
[Near-estimability]

Let \( \X_1 \) be \( n\times(p-1) \) of full column rank, \( \mathbf{c}\in\Real^{p-1} \), and \( \bz \) a unit vector
orthogonal to \( \C(\X_1) \). For \( \delta>0 \) let \( \X_\delta=[\X_1,\ \X_1\mathbf{c}+\delta\bz] \), with least squares
estimate \( \hbeta_\delta=(\hbeta_1\T,\hat{\beta}_2)\T \), and let \( \X_0=[\X_1,\ \X_1\mathbf{c}] \), which has rank
\( p-1 \). Write \( \mathbf{a}=(\mathbf{a}_1\T,a_2)\T \).

::: {.enumerate options="label=(\alph*)"}
1. \( \Var(\mathbf{a}\T\hbeta_\delta)=\sigma^2\mathbf{a}_1\T(\X_1\T\X_1)^{-1}\mathbf{a}_1+\sigma^2(a_2-\mathbf{c}\T\mathbf{a}_1)^2/\delta^2 \).

2. \( \mathbf{a}\T\bbeta \) is estimable in the model with matrix \( \X_0 \) iff \( a_2=\mathbf{c}\T\mathbf{a}_1 \). For such
   \( \mathbf{a} \), \( \Var(\mathbf{a}\T\hbeta_\delta) \) does not depend on \( \delta \) and equals the variance of the least
   squares estimator of \( \mathbf{a}\T\bbeta \) in the limiting model. For every other \( \mathbf{a} \) it grows like
   \( \delta^{-2} \) as \( \delta\to0 \).
:::

:::

::: {.proof}
(a) Apply @thm-proj-fwl with \( \X_2=\X_1 \) and the last column in the role of \( \X_1 \) there. The last
column residualized on \( \C(\X_1) \) is \( (\I-\M_1)(\X_1\mathbf{c}+\delta\bz)=\delta\bz \), so
\( \hat{\beta}_2=\bz\T\Y/\delta \). By part (c) of the same theorem and \( \X_1\T\bz=\bzero \),
\[
\hbeta_1=(\X_1\T\X_1)^{-1}\X_1\T\bigl(\Y-(\X_1\mathbf{c}+\delta\bz)\hat{\beta}_2\bigr)=\hbeta_1^{\mathrm{S}}-\mathbf{c}\,\hat{\beta}_2,
\qquad \hbeta_1^{\mathrm{S}}=(\X_1\T\X_1)^{-1}\X_1\T\Y .
\]
Hence \( \mathbf{a}\T\hbeta_\delta=\mathbf{a}_1\T\hbeta_1^{\mathrm{S}}+(a_2-\mathbf{c}\T\mathbf{a}_1)\hat{\beta}_2 \). The two terms are
uncorrelated, because \( \Cov(\hbeta_1^{\mathrm{S}},\hat{\beta}_2)=\sigma^2(\X_1\T\X_1)^{-1}\X_1\T\bz/\delta=\bzero \), and
\( \Var(\hat{\beta}_2)=\sigma^2/\delta^2 \).

(b) By @thm-est-characterization, \( \mathbf{a}\T\bbeta \) is estimable under \( \X_0 \) iff
\( \mathbf{a}\in\C(\X_0\T) \). Since \( \X_1 \) has full column rank, \( \C(\X_1\T)=\Real^{p-1} \), so
\( \C(\X_0\T)=\{(\mathbf{a}_1\T,\mathbf{c}\T\mathbf{a}_1)\T:\mathbf{a}_1\in\Real^{p-1}\} \). For such \( \mathbf{a} \) the second term in (a)
vanishes. In the limiting model \( \X_0\bbeta=\X_1(\bbeta_1+\mathbf{c}\beta_2) \), so
\( \mathbf{a}\T\bbeta=\mathbf{a}_1\T(\bbeta_1+\mathbf{c}\beta_2) \) is estimated by \( \mathbf{a}_1\T\hbeta_1^{\mathrm{S}} \), with variance
\( \sigma^2\mathbf{a}_1\T(\X_1\T\X_1)^{-1}\mathbf{a}_1 \). If \( a_2\ne\mathbf{c}\T\mathbf{a}_1 \), the second term is a positive
multiple of \( \delta^{-2} \).
:::

The orthogonality of \( \bz \) to \( \C(\X_1) \) matters only for the size of the bounded term. A general
perturbation \( \bw \) splits as \( \M_1\bw+(\I-\M_1)\bw \), and the first part, being \( \X_1\mathbf{d} \) for some
\( \mathbf{d} \), can be absorbed into \( \mathbf{c} \). The dichotomy survives: a function estimable in the limit keeps a
variance that does not depend on \( \delta \), and every other function's variance grows like \( \delta^{-2} \). But the
bounded value can exceed the limiting-model variance (@exr-col-near-general).
[Figure 8.3.1](../ch08-estimability/03-characterizations.html#fig-est-near) showed the \( \delta^{-2} \) growth
numerically; the proposition identifies the rate and the exact set of functions that escape it.

## Statistical, not numerical

[Chapter 10](../ch10-computation/index.html) met ill-conditioning as a numerical problem: a large
\( \kappa(\X) \) amplifies rounding errors (@thm-cmp-perturbation). The statistical cost is different. It is
\( \sigma^2/\lambda_\ell \) along each weak direction, so it involves the error variance and the absolute size of
the eigenvalues, not only their ratio. Rescaling a column changes \( \kappa(\X) \) (@exr-mat-cond-scaling)
but no estimand's variance, so diagnostics must be computed for a stated scaling
([Section 26.3](03-condition-indices.html)) or be scale-free ([Section 26.2](02-vif.html)). And the
eigenvalues grow with \( n \): if the rows are drawn independently from a distribution with nonsingular
second-moment matrix \( \bSigma \), then \( \lambda_p/n\to\lambda_{\min}(\bSigma) \) by the strong law of large numbers and @prp-cmp-weyl, so
collinearity in the population is paid for by a larger \( n \). Longley's regression shows how far
apart the two problems can be: QR computes its coefficients to about eleven digits (@exm-cmp-longley),
yet its variance inflation factors reach \( 1788.5 \) (@exm-dep-longley).

Measurement error in the regressors sits between the two. It perturbs \( \X \) as rounding does, by an
amount that more digits cannot reduce, and a perturbation comparable to \( \sqrt{\lambda_p} \) can turn the
weak eigenvector and change which combinations appear well estimated (@thm-cmp-perturbation, with
\( \kappa(\X) \) as the amplifier). The directions that collinearity leaves ill-determined are thus the least
robust to errors in \( \X \); [Chapter 24](../ch24-errors-in-variables/index.html) treats the resulting bias (@thm-eiv-attenuation).

## Exercises

### A. Check your understanding

::: {#exr-col-check-bounds}
[A1]

In @exm-col-two-regressors, the difference \( \beta_1-\beta_2 \) corresponds to \( \mathbf{a}=(0,1,-1)\T \), of squared
length \( 2 \). Check that the variance quoted for \( \hat{\beta}_1-\hat{\beta}_2 \) lies between the two bounds of
@thm-col-variance(b), applied to the slope block, and explain why it is close to the upper one.
:::

::: {.solution}
Computed without rounding, the variance is \( 11.640\,\sigma^2 \), and the bounds are \( 2\sigma^2/70.05=0.0286\,\sigma^2 \) and
\( 2\sigma^2/\lambda_2=11.643\,\sigma^2 \). (The check needs more digits than the table gives: \( 3.412^2 \) and
\( 2/0.1718 \) agree to four figures.) The variance is \( \sigma^2\sum_\ell(\mathbf{a}\T\bv_\ell)^2/\lambda_\ell \), and
\( (\mathbf{a}\T\bv_2)^2=1.9995 \) of the total \( \norm{\mathbf{a}}^2=2 \), because \( \mathbf{a} \) points almost exactly along
the weak eigenvector. Nearly all of its weight sits on the term with the small eigenvalue.
:::

::: {#exr-col-sign-probability}
[A2]

Under normal errors, \( \hat{\beta}_2\sim\Normal(\beta_2,\Var(\hat{\beta}_2)) \). With \( \beta_2=1 \) and the standard
deviation of @exm-col-two-regressors, compute \( \Pr(\hat{\beta}_2<0) \) and compare it with the share of
replicates in which \( \hat{\beta}_2 \) was negative.
:::

::: {#exr-col-units-direction}
[A3]

Let \( \X=[\x_1,\x_2] \) have orthogonal columns with \( \norm{\x_1}=1 \) and \( \norm{\x_2}=10 \). Find the canonical
coordinates and the weak direction. Now measure the second regressor in units ten times larger, so
that the column becomes \( \x_2/10 \). What happens to the weak direction and to \( \kappa(\X) \)? Which
variances of estimable quantities have changed?
:::

### B. Practice

::: {#exr-col-max-coefficient}
[B1]

Show that \( \max_j\Var(\hat{\beta}_j)\ge\sigma^2/(p\lambda_p) \), so that a small \( \lambda_p \) always hurts at
least one coefficient. Give an example with \( p=3 \) in which \( \lambda_p \) is tiny but only one
coefficient has a large variance.
:::

::: {.solution}
By @thm-col-variance(d), \( \sum_j\Var(\hat{\beta}_j)=\sigma^2\sum_\ell1/\lambda_\ell\ge\sigma^2/\lambda_p \), and the largest of \( p \)
numbers is at least their average. For the example take \( \X \) with orthonormal columns and multiply the
third column by \( \epsilon \). Then \( \X\T\X=\diag(1,1,\epsilon^2) \), \( \lambda_p=\epsilon^2 \), and only
\( \Var(\hat{\beta}_3)=\sigma^2/\epsilon^2 \) is large: a short column, not a near-dependence, which is why
[Section 26.3](03-condition-indices.html) scales the columns.
:::

::: {#exr-col-near-general}
[B2]

In @prp-col-near-estimable, replace \( \bz \) by an arbitrary vector \( \bw\notin\C(\X_1) \). Write
\( \bw=\X_1\mathbf{d}+\norm{(\I-\M_1)\bw}\,\bz \) with \( \bz \) a unit vector orthogonal to \( \C(\X_1) \), and show that
\[
\Var(\mathbf{a}\T\hbeta_\delta)=\sigma^2\mathbf{a}_1\T(\X_1\T\X_1)^{-1}\mathbf{a}_1
+\sigma^2\frac{\bigl(a_2-(\mathbf{c}+\delta\mathbf{d})\T\mathbf{a}_1\bigr)^2}{\delta^2\norm{(\I-\M_1)\bw}^2}.
\]
Deduce that if \( a_2=\mathbf{c}\T\mathbf{a}_1 \) the variance does not depend on \( \delta \) and equals the
limiting-model variance plus \( \sigma^2(\mathbf{d}\T\mathbf{a}_1)^2/\norm{(\I-\M_1)\bw}^2 \), which vanishes iff
\( \mathbf{d}\T\mathbf{a}_1=0 \), and that for every other \( \mathbf{a} \) it still grows like \( \delta^{-2} \).
:::

::: {.solution}
\( \X_1\mathbf{c}+\delta\bw=\X_1(\mathbf{c}+\delta\mathbf{d})+\delta\norm{(\I-\M_1)\bw}\,\bz \), which is the setting of the
proposition with \( \mathbf{c} \) replaced by \( \mathbf{c}+\delta\mathbf{d} \) and \( \delta \) by \( \delta\norm{(\I-\M_1)\bw} \). Part (a)
gives the formula. If \( a_2=\mathbf{c}\T\mathbf{a}_1 \), the numerator of the second term is \( \delta^2(\mathbf{d}\T\mathbf{a}_1)^2 \), the
\( \delta^2 \) cancels, and the second term is the constant \( \sigma^2(\mathbf{d}\T\mathbf{a}_1)^2/\norm{(\I-\M_1)\bw}^2 \); the first
is the limiting-model variance, by part (b) of the proposition. Otherwise the numerator tends to \( (a_2-\mathbf{c}\T\mathbf{a}_1)^2>0 \) and the term grows like
\( \delta^{-2} \).
:::

::: {#exr-col-ch08-sum}
[B3]

[Section 8.3](../ch08-estimability/03-characterizations.html) used twelve equally spaced
\( t_i \) on \( [-1,1] \), \( \X_1=[\bone,\mathbf{t}] \), \( \mathbf{c}=(0,1)\T \) and \( \bw=(\cos\pi t_i) \), and reported
\( \Var(\hat{\beta}_1+\hat{\beta}_2)=0.212\,\sigma^2 \) whatever \( \delta \). Use @exr-col-near-general to explain
the number: show that \( \mathbf{d}=(d_1,0)\T \) because \( \bw \) is even and \( \mathbf{t} \) is odd, and that the variance equals
\( \sigma^2/\sum_it_i^2 \).
:::

::: {.solution}
Here \( \mathbf{a}=(0,1,1)\T \), so \( \mathbf{a}_1=(0,1)\T \) and \( a_2=\mathbf{c}\T\mathbf{a}_1=1 \): the function is estimable in the
limit. The coefficients \( \mathbf{d} \) of the projection of \( \bw \) on \( \C(\bone,\mathbf{t}) \) are
\( d_1=\bar w \) and \( d_2=\mathbf{t}\T\bw/\mathbf{t}\T\mathbf{t}=0 \), because the \( t_i \) are symmetric about zero and \( \cos \) is even.
So \( \mathbf{d}\T\mathbf{a}_1=d_2=0 \) and the second term vanishes. The first is
\( \sigma^2[(\X_1\T\X_1)^{-1}]_{22}=\sigma^2/\sum_it_i^2 \), since \( \bone\perp\mathbf{t} \). With
\( t_i=-1+2i/11 \), \( i=0,\dots,11 \), \( \sum_it_i^2=(4/121)\sum_{i=0}^{11}(i-5.5)^2=(4/121)\cdot143=4.727 \), and
\( 1/4.727=0.2115 \), which rounds to \( 0.212 \).
:::
