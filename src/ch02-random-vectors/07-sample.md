# Sample moments and their computation

The population quantities of this chapter have sample analogues, and the operator rules
give their expectations at once. Let \( \y_1,\dots,\y_n \) be \( p\times1 \) observations, and stack them as
the rows of the \( n\times p \) data matrix \( \boldsymbol{\mathcal{Y}}=[\y_1,\dots,\y_n]\T \). The script letter
keeps the data matrix apart from a single random vector \( \Y \). The rows keep the lowercase
\( \y_i \) of observed data, but in @prp-rv-sample-moments they are treated as random. With
\( \mathbf{C}=\I_n-n^{-1}\bone\bone\T \), the sample mean vector and the **sample covariance matrix** are
\[
\bar{\y}=\frac1n\sum_{i=1}^n\y_i=\frac1n\boldsymbol{\mathcal{Y}}\T\bone,
\qquad
\bS=\frac{1}{n-1}\sum_{i=1}^n(\y_i-\bar{\y})(\y_i-\bar{\y})\T=\frac{1}{n-1}\boldsymbol{\mathcal{Y}}\T\mathbf{C}\boldsymbol{\mathcal{Y}} .
\]
The second form of \( \bS \) holds because \( \mathbf{C}\boldsymbol{\mathcal{Y}} \) is the data matrix with its column means
subtracted, and \( \mathbf{C} \) is symmetric and idempotent. The sample correlation matrix is
\( \hat{\R}=\hat{\bD}^{-1}\bS\hat{\bD}^{-1} \), where \( \hat{\bD} \) holds the sample standard deviations.

::: {#prp-rv-sample-moments}
[Moments of sample moments]

Suppose the rows \( \y_i \) are pairwise uncorrelated random vectors with common mean \( \bmu \) and
common covariance matrix \( \bSigma \). Then

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\bar{\y})=\bmu \) and \( \Cov(\bar{\y})=\bSigma/n \);

2. \( \E(\bS)=\bSigma \);

3. \( \bS \) is nonnegative definite, and \( \rank(\bS)\le\min(p,\,n-1) \) for every realization.
:::

:::

::: {.proof}
(a) is @prp-rv-sums(c) with \( c_i=1/n \). For (b), let \( \mathbf{u}_j \) denote the \( j \)th column of
\( \boldsymbol{\mathcal{Y}} \), the \( n \) observations of variable \( j \). The \( (j,k) \) entry of \( \boldsymbol{\mathcal{Y}}\T\mathbf{C}\boldsymbol{\mathcal{Y}} \) is the bilinear
form \( \mathbf{u}_j\T\mathbf{C}\mathbf{u}_k \). Because different rows are uncorrelated,
\( \Cov(\mathbf{u}_k,\mathbf{u}_j)=\sigma_{kj}\I_n \), and \( \E(\mathbf{u}_j)=\mu_j\bone \). By
@cor-rv-quadform-special(d),
\[
\E(\mathbf{u}_j\T\mathbf{C}\mathbf{u}_k)=\sigma_{kj}\tr(\mathbf{C})+\mu_j\mu_k\,\bone\T\mathbf{C}\bone=(n-1)\sigma_{jk},
\]
since \( \tr(\mathbf{C})=n-1 \) and \( \mathbf{C}\bone=\bzero \). (c) For any \( \mathbf{a} \),
\( \mathbf{a}\T\bS\mathbf{a}=\norm{\mathbf{C}\boldsymbol{\mathcal{Y}}\mathbf{a}}^2/(n-1)\ge0 \). The rank of \( \boldsymbol{\mathcal{Y}}\T\mathbf{C}\boldsymbol{\mathcal{Y}}=(\mathbf{C}\boldsymbol{\mathcal{Y}})\T(\mathbf{C}\boldsymbol{\mathcal{Y}}) \) equals
\( \rank(\mathbf{C}\boldsymbol{\mathcal{Y}}) \), which is at most \( \rank(\mathbf{C})=n-1 \) and at most \( p \).
:::

The divisor \( n-1 \) is the trace of \( \mathbf{C} \), the dimension of the space the centred data live in,
exactly as \( n-r \) is for regression residuals in @exm-rv-rss-bias. Part (c) has a practical
side. With \( n\le p \), as in genomics or finance, the sample covariance matrix is always
singular, even when \( \bSigma \) is positive definite. It has a linear combination with zero sample
variance whether or not one exists in the population. [Chapter 28](../ch28-high-dimensional/index.html) returns to
this.

**Computation.**  In practice \( \mathbf{C} \) is never formed: one subtracts column means. The
eigenvalues of \( \bS \) are the squared singular values of \( \mathbf{C}\boldsymbol{\mathcal{Y}}/\sqrt{n-1} \)
(@thm-mat-svd), which can be computed without forming a cross-product matrix
([Chapter 10](../ch10-computation/index.html)).

::: {#exm-rv-longley}
[Macroeconomic series]

Longley's data give \( n=16 \) annual observations (1947–1962) on
\( p=7 \) US series: total employment, the GNP deflator, GNP, unemployment, armed
forces, population and the year. The eigenvalues of \( \bS \) run from \( 9.94\times 10^{9} \) down
to \( 1.07\times 10^{-2} \), a ratio of \( 9.3\times 10^{11} \). That spread says little about
the structure of the data. It mostly reflects units: the standard deviation of GNP is
\( 99395 \), that of the deflator \( 10.79 \).

The correlation matrix removes the units (@prp-rv-correlation(d)). Its eigenvalues sum to
\( p=7 \). The largest is \( 5.533 \), so one direction carries about
\( 79 \)% of the standardized variation, and the smallest is
\( 2.59\times 10^{-4} \). By @thm-mat-extremal-rayleigh, the unit vector \( \mathbf{v} \) belonging
to the smallest eigenvalue minimizes the sample variance of \( \mathbf{v}\T\bz \) over standardized
combinations. Its two largest weights are \( 0.681 \) on year and
\( -0.674 \) on GNP, and the combination has standard deviation
\( 0.016 \), against \( 1 \) for each standardized variable. In these sixteen years GNP
is almost an exact linear function of time (their correlation is
\( 0.995 \)). This is the sample version of @thm-rv-cov-nnd(c): a nearly
zero eigenvalue marks a nearly constant linear combination. For regression it signals
collinearity, which [Chapter 26](../ch26-collinearity/index.html) studies. Finally, the first
\( 5 \) rows alone give a sample covariance matrix of rank
\( 4 \), as @prp-rv-sample-moments(c) requires.
:::

```{.python .run #cell-longley-moments-sample}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.longley.load_pandas().data
names = list(data.columns)

Y = data.to_numpy(dtype=float)                    # n x p data matrix, one row per year
n, p = Y.shape
C = np.eye(n) - np.ones((n, n)) / n               # centring matrix I - J/n
S = Y.T @ C @ Y / (n - 1)                         # sample covariance matrix
d = np.sqrt(np.diag(S))
R = S / np.outer(d, d)                            # D^{-1} S D^{-1}, D = diag(d)

lam_S = np.linalg.eigvalsh(S)
lam_R, V = np.linalg.eigh(R)                      # ascending eigenvalues
print("eigenvalues of S:", lam_S)
print("eigenvalues of R:", lam_R)
v = V[:, 0]                                       # direction of least variance
Z = (Y - Y.mean(0)) / d                           # standardized variables
print("sample variance of Z v:", np.var(Z @ v, ddof=1))
```

## Exercises

### B. Practice

::: {#exr-rv-sample-serial}
[B1]

Suppose the rows of \( \boldsymbol{\mathcal{Y}} \) are observed in time order, with common mean \( \bmu \) and
covariance \( \bSigma \), and that neighbouring rows are correlated: \( \Cov(\y_i,\y_{i+1})=\rho\bSigma \),
while \( \Cov(\y_i,\y_j)=\bzero \) for \( \lvert i-j\rvert\ge2 \). Show that \( \E(\bS)=(1-2\rho/n)\bSigma \) and
\( \Cov(\bar{\y})=n^{-2}[n+2(n-1)\rho]\bSigma \). Which of the two effects disappears as \( n\to\infty \)?
:::

::: {.solution}
As in the proof of @prp-rv-sample-moments,
\( \Cov(\mathbf{u}_k,\mathbf{u}_j)=\sigma_{kj}\bT \), where \( \bT \) is the \( n\times n \) tridiagonal matrix with ones on
the diagonal and \( \rho \) next to it. Since \( \mathbf{C}\bone=\bzero \), @cor-rv-quadform-special(d) gives
\( \E(\mathbf{u}_j\T\mathbf{C}\mathbf{u}_k)=\sigma_{kj}\tr(\mathbf{C}\bT)=\sigma_{kj}\bigl[\tr(\bT)-n^{-1}\bone\T\bT\bone\bigr] \). Here
\( \tr(\bT)=n \) and \( \bone\T\bT\bone=n+2(n-1)\rho \), so the bracket is \( (n-1)-2(n-1)\rho/n \), and dividing by
\( n-1 \) gives \( \E(\bS)=(1-2\rho/n)\bSigma \). For the mean, \( \bar{\y}=n^{-1}\sum_i\y_i \) and
@prp-rv-sums gives \( n^{-2} \) times the sum of all \( \Cov(\y_i,\y_j) \): \( n \) diagonal terms
\( \bSigma \) and \( 2(n-1) \) neighbouring terms \( \rho\bSigma \). The bias of \( \bS \) is \( O(1/n) \) and disappears,
but \( n\Cov(\bar{\y})\to(1+2\rho)\bSigma \), not \( \bSigma \). For \( \rho>0 \) the naive standard error
\( \sqrt{s_{jj}/n} \) stays too small by a fixed factor however large the sample.
:::

### C. Going deeper

::: {#exr-rv-two-observations}
[C1]

::: {.enumerate options="label=(\alph*)"}
1. With \( n=2 \), show that \( \bS=\tfrac12\mathbf{d}\mathbf{d}\T \) with \( \mathbf{d}=\y_1-\y_2 \). Deduce
           that \( \rank(\bS)=1 \) whenever \( \y_1\ne\y_2 \), and that every sample correlation is \( \pm1 \),
           whatever \( \bSigma \) is.

2. With \( n=3 \) and \( p=3 \), show that \( \det(\bS)=0 \), hence that the three sample correlations
           satisfy \( 1+2\hat{r}_{12}\hat{r}_{13}\hat{r}_{23}-\hat{r}_{12}^2-\hat{r}_{13}^2-\hat{r}_{23}^2=0 \)
           when the sample variances are positive. Conclude that \( \hat{r}_{23} \) is always one of the
           two *endpoints* of the interval allowed by @prp-rv-three-correlations, never an interior
           point.

3. What does this say about reading a sample correlation matrix computed from few
           observations, and about the nearly singular \( \hat{\R} \) of @exm-rv-longley?
:::

:::

::: {.solution}
(a) Here \( \bar{\y}=\tfrac12(\y_1+\y_2) \), so
\( \y_1-\bar{\y}=\tfrac12\mathbf{d} \) and \( \y_2-\bar{\y}=-\tfrac12\mathbf{d} \). Summing the two outer
products and dividing by \( n-1=1 \) gives \( \bS=\tfrac12\mathbf{d}\mathbf{d}\T \), of rank one. Then
\( s_{jk}=d_jd_k/2 \) and \( \hat{r}_{jk}=d_jd_k/\sqrt{d_j^2d_k^2}=\operatorname{sign}(d_jd_k) \),
provided \( d_j,d_k\ne0 \). The two observations determine one direction, and every pair of
coordinates moves along it.

(b) By @prp-rv-sample-moments(c), \( \rank(\bS)\le n-1=2<3 \), so \( \det(\bS)=0 \). Writing
\( \bS=\hat{\bD}\hat{\R}\hat{\bD} \) with \( \hat{\bD} \) the diagonal matrix of sample standard
deviations, which is nonsingular when all sample variances are positive, gives
\( \det(\hat{\R})=0 \). Expanding the \( 3\times3 \) determinant gives the stated identity, which is a
quadratic in \( \hat{r}_{23} \) with roots
\( \hat{r}_{12}\hat{r}_{13}\pm\sqrt{(1-\hat{r}_{12}^2)(1-\hat{r}_{13}^2)} \), the two endpoints
in @prp-rv-three-correlations. A population correlation matrix hits these endpoints only when
one variable is an exact linear function of the other two; three observations force it.

(c) With \( n\le p \) the sample correlation matrix is singular whatever the population is, so
its smallest eigenvalues and any statement derived from them describe the sample size as
much as the data. Longley's \( \hat{\R} \) has \( n=16 \) and \( p=7 \), so its near singularity is not
forced by the rank bound; it is a genuine feature of the series.
:::

::: {#exr-rv-smallest-eigenvalue}
[C2]

Let \( \bS \) be positive definite and let \( \lambda_{\min} \) be its smallest eigenvalue.

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \lambda_{\min}=\min\{\mathbf{a}\T\bS\mathbf{a}:\norm{\mathbf{a}}=1\} \)
           (@thm-mat-extremal-rayleigh), while
           \( 1/(\bS^{-1})_{jj}=\min\{\mathbf{a}\T\bS\mathbf{a}:a_j=1\} \) (@thm-mat-partitioned-inverse).

2. Deduce that \( \lambda_{\min}\le1/(\bS^{-1})_{jj} \) for every \( j \), with equality iff \( \mathbf{e}_j \) is
           an eigenvector of \( \bS \) belonging to \( \lambda_{\min} \).

3. For centred data, \( \mathbf{a}\T\bS\mathbf{a}=\norm{\mathbf{C}\boldsymbol{\mathcal{Y}}\mathbf{a}}^2/(n-1) \), so
           \( 1/(\bS^{-1})_{jj} \) is the sample variance, with divisor \( n-1 \), of the residuals of
           the least squares regression of variable \( j \) on the other \( p-1 \) variables, while
           \( \lambda_{\min} \) minimizes over *all* normalized combinations. Using the smallest eigenvalue \( 2.59\times 10^{-4} \) of the Longley
           correlation matrix in @exm-rv-longley, what can you say about how well any one
           standardized series is predicted by the other six?
:::

:::

::: {.solution}
(a) The first statement is the Rayleigh quotient characterization. For the
second, order the variables so that \( j=1 \) and partition \( \mathbf{a}=(1,\bb\T)\T \). Then
\( \mathbf{a}\T\bS\mathbf{a}=s_{11}+2\bb\T\bS_{21}+\bb\T\bS_{22}\bb \), a positive definite quadratic in
\( \bb \), minimized at \( \bb=-\bS_{22}^{-1}\bS_{21} \) with value the Schur complement
\( s_{11}-\bS_{12}\bS_{22}^{-1}\bS_{21} \), which by @thm-mat-partitioned-inverse is
\( 1/(\bS^{-1})_{11} \).

(b) Any \( \mathbf{a} \) with \( a_j=1 \) has \( \norm{\mathbf{a}}\ge1 \), so
\( \mathbf{a}\T\bS\mathbf{a}\ge\lambda_{\min}\norm{\mathbf{a}}^2\ge\lambda_{\min} \); minimizing over such \( \mathbf{a} \)
gives the inequality. Equality forces a minimizer with \( \norm{\mathbf{a}}=1 \) and \( a_j=1 \), hence
\( \mathbf{a}=\mathbf{e}_j \), and \( \mathbf{e}_j \) attains the Rayleigh minimum only if it is an eigenvector for
\( \lambda_{\min} \). Conversely, if it is, both sides equal \( \lambda_{\min} \).

(c) The residual sample variance of any one standardized series regressed on the other six is
at least \( 2.59\times 10^{-4} \), so by @prp-rv-multiple-correlation its multiple correlation
with the others satisfies \( 1-R^2\ge2.59\times 10^{-4} \). Collinearity as measured by the
smallest eigenvalue is therefore always at least as severe as collinearity measured one
variable at a time: the direction of least variance need not be close to any coordinate
direction, and in Longley's data it is a contrast between year and GNP rather than a single
series.
:::
