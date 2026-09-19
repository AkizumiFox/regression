# Sample moments and their computation

The population quantities of this chapter have sample analogues, and the operator rules
give their expectations at once. Let \( \y_1,\dots,\y_n \) be \( p\times1 \) observations, and stack them as
the rows of the \( n\times p \) data matrix \( \bm{\mathcal{Y}}=[\y_1,\dots,\y_n]\T \). The script letter
keeps the data matrix apart from a single random vector \( \Y \). The rows keep the lowercase
\( \y_i \) of observed data, but in @prp-rv-sample-moments they are treated as random. With
\( \bm C=\I_n-n^{-1}\bone\bone\T \), the sample mean vector and the **sample covariance matrix** are
\[
\bar{\y}=\frac1n\sum_{i=1}^n\y_i=\frac1n\bm{\mathcal{Y}}\T\bone,
\qquad
\bS=\frac{1}{n-1}\sum_{i=1}^n(\y_i-\bar{\y})(\y_i-\bar{\y})\T=\frac{1}{n-1}\bm{\mathcal{Y}}\T\bm C\bm{\mathcal{Y}} .
\]
The second form of \( \bS \) holds because \( \bm C\bm{\mathcal{Y}} \) is the data matrix with its column means
subtracted, and \( \bm C \) is symmetric and idempotent. The sample correlation matrix is
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
(a) is @prp-rv-sums(c) with \( c_i=1/n \). For (b), let \( \bm u_j \) denote the \( j \)th column of
\( \bm{\mathcal{Y}} \), the \( n \) observations of variable \( j \). The \( (j,k) \) entry of \( \bm{\mathcal{Y}}\T\bm C\bm{\mathcal{Y}} \) is the bilinear
form \( \bm u_j\T\bm C\bm u_k \). Because different rows are uncorrelated,
\( \Cov(\bm u_k,\bm u_j)=\sigma_{kj}\I_n \), and \( \E(\bm u_j)=\mu_j\bone \). By
@cor-rv-quadform-special(d),
\[
\E(\bm u_j\T\bm C\bm u_k)=\sigma_{kj}\tr(\bm C)+\mu_j\mu_k\,\bone\T\bm C\bone=(n-1)\sigma_{jk},
\]
since \( \tr(\bm C)=n-1 \) and \( \bm C\bone=\bzero \). (c) For any \( \bm a \),
\( \bm a\T\bS\bm a=\norm{\bm C\bm{\mathcal{Y}}\bm a}^2/(n-1)\ge0 \). The rank of \( \bm{\mathcal{Y}}\T\bm C\bm{\mathcal{Y}}=(\bm C\bm{\mathcal{Y}})\T(\bm C\bm{\mathcal{Y}}) \) equals
\( \rank(\bm C\bm{\mathcal{Y}}) \), which is at most \( \rank(\bm C)=n-1 \) and at most \( p \).
:::

The divisor \( n-1 \) is the trace of \( \bm C \), the dimension of the space the centred data live in,
exactly as \( n-r \) is for regression residuals in @exm-rv-rss-bias. Part (c) has a practical
side. With \( n\le p \), as in genomics or finance, the sample covariance matrix is always
singular, even when \( \bSigma \) is positive definite. It has a linear combination with zero sample
variance whether or not one exists in the population. Chapter 28 returns to
this.

**Computation.**  In practice \( \bm C \) is never formed: one subtracts column means. The
eigenvalues of \( \bS \) are the squared singular values of \( \bm C\bm{\mathcal{Y}}/\sqrt{n-1} \)
(@thm-mat-svd), which can be computed without forming a cross-product matrix
(Chapter 10).

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
\( 2.59\times 10^{-4} \). By @thm-mat-extremal-rayleigh, the unit vector \( \bm v \) belonging
to the smallest eigenvalue minimizes the sample variance of \( \bm v\T\bz \) over standardized
combinations. Its two largest weights are \( 0.681 \) on year and
\( -0.674 \) on GNP, and the combination has standard deviation
\( 0.016 \), against \( 1 \) for each standardized variable. In these sixteen years GNP
is almost an exact linear function of time (their correlation is
\( 0.995 \)). This is the sample version of @thm-rv-cov-nnd(c): a nearly
zero eigenvalue marks a nearly constant linear combination. For regression it signals
collinearity, which Chapter 26 studies. Finally, the first
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

Suppose the rows of \( \bm{\mathcal{Y}} \) are observed in time order, with common mean \( \bmu \) and
covariance \( \bSigma \), and that neighbouring rows are correlated: \( \Cov(\y_i,\y_{i+1})=\rho\bSigma \),
while \( \Cov(\y_i,\y_j)=\bzero \) for \( \lvert i-j\rvert\ge2 \). Show that \( \E(\bS)=(1-2\rho/n)\bSigma \) and
\( \Cov(\bar{\y})=n^{-2}[n+2(n-1)\rho]\bSigma \). Which of the two effects disappears as \( n\to\infty \)?
:::

::: {.solution}
As in the proof of @prp-rv-sample-moments,
\( \Cov(\bm u_k,\bm u_j)=\sigma_{kj}\bT \), where \( \bT \) is the \( n\times n \) tridiagonal matrix with ones on
the diagonal and \( \rho \) next to it. Since \( \bm C\bone=\bzero \), @cor-rv-quadform-special(d) gives
\( \E(\bm u_j\T\bm C\bm u_k)=\sigma_{kj}\tr(\bm C\bT)=\sigma_{kj}\bigl[\tr(\bT)-n^{-1}\bone\T\bT\bone\bigr] \). Here
\( \tr(\bT)=n \) and \( \bone\T\bT\bone=n+2(n-1)\rho \), so the bracket is \( (n-1)-2(n-1)\rho/n \), and dividing by
\( n-1 \) gives \( \E(\bS)=(1-2\rho/n)\bSigma \). For the mean, \( \bar{\y}=n^{-1}\sum_i\y_i \) and
@prp-rv-sums gives \( n^{-2} \) times the sum of all \( \Cov(\y_i,\y_j) \): \( n \) diagonal terms
\( \bSigma \) and \( 2(n-1) \) neighbouring terms \( \rho\bSigma \). The bias of \( \bS \) is \( O(1/n) \) and disappears,
but \( n\Cov(\bar{\y})\to(1+2\rho)\bSigma \), not \( \bSigma \). For \( \rho>0 \) the naive standard error
\( \sqrt{s_{jj}/n} \) stays too small by a fixed factor however large the sample.
:::
