# Moment generating function and density

## Moment generating functions of random vectors

The moment generating function of a random vector \( \Y\in\Real^n \) is
\( M_{\Y}(\mathbf{t})=\E\exp(\mathbf{t}\T\Y) \), and @thm-rv-mgf proves the two
properties we need when it is finite on a neighbourhood of \( \bzero \).
*Uniqueness*: two random vectors whose moment generating functions are finite
and equal near \( \bzero \) have the same distribution. *Factorization*: two blocks
\( \Y_1 \), \( \Y_2 \) are independent iff \( M_{\Y}(\mathbf{t})=M_{\Y_1}(\mathbf{t}_1)M_{\Y_2}(\mathbf{t}_2) \) near
\( \bzero \). The second extends to \( k \) blocks by induction. If \( M_{\Y} \) is the product of the
\( k \) block functions, setting \( \mathbf{t}_1=\bzero \) shows that \( (\Y_2,\dots,\Y_k) \) has the product
of the last \( k-1 \), so \( \Y_1 \) is independent of \( (\Y_2,\dots,\Y_k) \), whose blocks are
mutually independent by the induction hypothesis. The converse is immediate.

::: {#thm-mvn-mgf}
[Moment generating function]

A random vector \( \Y \) is \( \Normal_n(\bmu,\bSigma) \) iff
\[
M_{\Y}(\mathbf{t})=\exp\bigl(\mathbf{t}\T\bmu+\tfrac12\mathbf{t}\T\bSigma\mathbf{t}\bigr)
\qquad\text{for every }\mathbf{t}\in\Real^n .
\]{#eq-mvn-mgf}

Consequently the distribution of a normal vector is determined by its mean vector
and covariance matrix.
:::

::: {.proof}
Suppose \( \Y \) is normal with mean \( \bmu \) and covariance \( \bSigma \). By
@eq-mvn-linear-combination, \( \mathbf{t}\T\Y\sim\Normal(\mathbf{t}\T\bmu,\mathbf{t}\T\bSigma\mathbf{t}) \),
and \( M_{\Y}(\mathbf{t}) \) is the moment generating function of this scalar evaluated at
\( s=1 \). @lem-mvn-univariate(a) gives @eq-mvn-mgf.

Conversely, suppose @eq-mvn-mgf holds for a symmetric \( \bSigma \). Fix
\( \mathbf{a} \). The moment generating function of \( \mathbf{a}\T\Y \) at \( s \) is
\( M_{\Y}(s\mathbf{a})=\exp\{(\mathbf{a}\T\bmu)s+\tfrac12(\mathbf{a}\T\bSigma\mathbf{a})s^2\} \). The
coefficient \( \mathbf{a}\T\bSigma\mathbf{a} \) is nonnegative, since otherwise this function
would tend to zero as \( s\to\pm\infty \), whereas every random variable \( X \) has
\( \E e^{sX}+\E e^{-sX}\ge\E e^{s\lvert X\rvert}\ge1 \). By
@lem-mvn-univariate(b), \( \mathbf{a}\T\Y\sim\Normal(\mathbf{a}\T\bmu,\mathbf{a}\T\bSigma\mathbf{a}) \).
So \( \Y \) is normal, and matching means and variances of all linear combinations
shows \( \E\Y=\bmu \) and \( \Cov(\Y)=\bSigma \).

Finally, two normal vectors with the same \( \bmu \) and \( \bSigma \) have the same moment
generating function @eq-mvn-mgf, which is finite everywhere, so they have
the same distribution by the uniqueness property.
:::

The last sentence is what makes \( \Normal_n(\bmu,\bSigma) \) a well-defined
distribution. Combined with @thm-mvn-existence, it says that
\( \Y\sim\Normal_n(\bmu,\bSigma) \) holds iff \( \Y \) has the same distribution as
\( \bmu+\A\Z \) for some, and then every, matrix \( \A \) with \( \A\A\T=\bSigma \). The
number of columns of \( \A \), and so the length of \( \Z \), does not matter. Some
books take this representation as the definition
(Christensen 2020; Seber and Lee 2003). The proof above shows the two
definitions agree.

The moment generating function also encodes moments. Differentiating
@eq-mvn-mgf at \( \mathbf{t}=\bzero \) gives \( \partial M_{\Y}/\partial\mathbf{t}=\bmu \)
and \( \partial^2M_{\Y}/\partial\mathbf{t}\,\partial\mathbf{t}\T=\bSigma+\bmu\bmu\T=\E(\Y\Y\T) \).
Higher derivatives give the fourth moments needed for variances of quadratic
forms in [Chapter 4](../ch04-quadratic-forms/index.html) (@exr-mvn-fourth-moments).

## The density

::: {#thm-mvn-density}
[Density]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \).

::: {.enumerate options="label=(\alph*)"}
1. If \( \bSigma \) is positive definite, \( \Y \) has the density
     \[
f(\y)=(2\pi)^{-n/2}\det(\bSigma)^{-1/2}
      \exp\Bigl\{-\tfrac12(\y-\bmu)\T\bSigma^{-1}(\y-\bmu)\Bigr\},
    \qquad \y\in\Real^n .
\]{#eq-mvn-density}

2. Conversely, a random vector with density @eq-mvn-density, for some
           positive definite \( \bSigma \), is \( \Normal_n(\bmu,\bSigma) \).

3. If \( \bSigma \) is singular, \( \Y \) has no density on \( \Real^n \).
:::

:::

::: {.proof}
(a) By @thm-mvn-mgf, \( \Y \) has the same distribution as
\( \bmu+\bSigma^{1/2}\Z \) with \( \Z\sim\Normal_n(\bzero,\I) \), so it suffices to find the
density of the latter. The entries of \( \Z \) are independent standard normals,
so \( \Z \) has density
\( \prod_i(2\pi)^{-1/2}e^{-z_i^2/2}=(2\pi)^{-n/2}\exp(-\tfrac12\bz\T\bz) \).
Since \( \bSigma \) is positive definite, so is \( \bSigma^{1/2} \) (its eigenvalues are
the positive square roots of those of \( \bSigma \)), and the map
\( \bz\mapsto\y=\bmu+\bSigma^{1/2}\bz \) is an affine bijection of \( \Real^n \) with inverse
\( \bz=\bSigma^{-1/2}(\y-\bmu) \). The change-of-variables formula multiplies the
density of \( \Z \) at \( \bz(\y) \) by \( \lvert\det\bSigma^{-1/2}\rvert=\det(\bSigma)^{-1/2} \).
Finally \( \bz\T\bz=(\y-\bmu)\T\bSigma^{-1/2}\bSigma^{-1/2}(\y-\bmu)
=(\y-\bmu)\T\bSigma^{-1}(\y-\bmu) \).
(b) A density determines a distribution, and by (a) the \( \Normal_n(\bmu,\bSigma) \)
distribution has this density.
(c) \( \Y \) lies with probability one in \( \bmu+\C(\bSigma) \), an affine subspace of
dimension less than \( n \). It has Lebesgue measure zero, so it would receive
probability zero under any distribution with a density.
:::

A singular normal vector does have a density with respect to the natural
measure on the affine subspace where it lives
(@prp-mvn-rank-representation and @exr-mvn-singular-density), but that is
rarely needed. In this book, arguments are made with the moment generating
function or with linear transformations, which work in every case, and densities
are used only for pictures and likelihoods.

::: {#exm-mvn-bivariate}
[The bivariate normal]

For \( n=2 \), write \( \sigma_1^2,\sigma_2^2>0 \) for the variances and \( \rho \) for the
correlation. Then
\[
\bSigma=\begin{pmatrix}\sigma_1^2&\rho\sigma_1\sigma_2\\\rho\sigma_1\sigma_2&\sigma_2^2\end{pmatrix},
\qquad
\det\bSigma=\sigma_1^2\sigma_2^2(1-\rho^2),
\]
so \( \bSigma \) is positive definite iff \( \lvert\rho\rvert<1 \). With standardized
coordinates \( u_i=(y_i-\mu_i)/\sigma_i \), the density @eq-mvn-density is
\[
f(y_1,y_2)=\frac{1}{2\pi\sigma_1\sigma_2\sqrt{1-\rho^2}}
\exp\Bigl\{-\frac{u_1^2-2\rho u_1u_2+u_2^2}{2(1-\rho^2)}\Bigr\}.
\]
When \( \rho=0 \) it factors into the two univariate normal densities. When
\( \rho=\pm1 \) there is no density: \( \Var(u_2\mp u_1)=0 \) in standardized units, so
all the probability lies on a line.
:::

## Contours and the Mahalanobis distance

For positive definite \( \bSigma \), the density @eq-mvn-density depends on
\( \y \) only through
\[
\Delta^2(\y)=(\y-\bmu)\T\bSigma^{-1}(\y-\bmu),
\]
the squared Mahalanobis distance from \( \y \) to \( \bmu \) introduced in [Section 2.2](../ch02-random-vectors/02-linear.html). The density is a
decreasing function of \( \Delta^2 \), so its contours are the ellipsoids
\( \{\y: \Delta^2(\y)=c^2\} \). Their shape comes from the spectral decomposition
\( \bSigma=\sum_i\lambda_i\bu_i\bu_i\T \) (@thm-mat-spectral). In the
coordinates \( w_i=\bu_i\T(\y-\bmu) \) we have \( \Delta^2=\sum_iw_i^2/\lambda_i \), so the
ellipsoid has principal axes along the eigenvectors \( \bu_i \), with half-lengths
\( c\sqrt{\lambda_i} \). Large-variance directions are long axes. The Mahalanobis
distance measures distance from the centre in units of standard deviation along
each principal direction. The volume of \( \{\y:\Delta^2(\y)\le c^2\} \) is proportional to
\( \prod_i\sqrt{\lambda_i}=\det(\bSigma)^{1/2} \), so the generalized variance \( \det\bSigma \) ([Section 2.1](../ch02-random-vectors/01-mean-cov.html))
measures the size of the region the distribution occupies. The factor \( \det(\bSigma)^{-1/2} \) in
@eq-mvn-density spreads the same total probability over that region: a larger generalized
variance means a flatter density.

::: {#prp-mvn-mahalanobis}
If \( \Y\sim\Normal_n(\bmu,\bSigma) \) with \( \bSigma \) positive definite, then
\( \Delta^2(\Y) \) has the chi-squared distribution with \( n \) degrees of freedom,
\( \chi^2(n) \), which is \( \chi^2(n,0) \) in the notation of @def-qf-noncentral-chisq. So the ellipsoid
\( \{\y:\Delta^2(\y)\le q\} \), where \( q \) is the \( 1-\alpha \) quantile of \( \chi^2(n) \), has
probability exactly \( 1-\alpha \).
:::

::: {.proof}
Put \( \Z=\bSigma^{-1/2}(\Y-\bmu) \). For every \( \mathbf{a} \),
\( \mathbf{a}\T\Z=(\bSigma^{-1/2}\mathbf{a})\T\Y-\mathbf{a}\T\bSigma^{-1/2}\bmu \) is a linear combination of \( \Y \)
minus a constant, hence normal by @def-mvn, so \( \Z \) is normal. Its mean is zero and its
covariance is \( \I \) (@prp-rv-whitening(a)). By
@thm-mvn-mgf its moment generating function is
\( \exp(\tfrac12\mathbf{t}\T\mathbf{t})=\prod_i\exp(\tfrac12t_i^2) \), which factorizes, so its entries
are independent standard normals.
\( \Delta^2(\Y)=\Z\T\Z=\sum_iZ_i^2 \) is a sum of \( n \) independent squared standard normals,
which is the definition of \( \chi^2(n) \).
:::

[Chapter 4](../ch04-quadratic-forms/index.html) generalizes this considerably, to \( \Y\T\A\Y \) with a
nonzero mean and a general \( \A \) (@thm-qf-chisq). Because the density is
largest inside these ellipsoids, they are the smallest regions with their
probability content, which is why they reappear as confidence regions for
\( \bbeta \) in [Chapter 12](../ch12-intervals-and-bands/index.html) (@thm-ci-ellipsoid).

**Computing with the density.**  Neither \( \bSigma^{-1} \) nor \( \det\bSigma \)
should be formed explicitly. With the Cholesky factorization \( \bSigma=\bL\bL\T \),
we have \( \Delta^2(\y)=\norm{\bL^{-1}(\y-\bmu)}^2 \), computed by a triangular solve, and
\( \log\det\bSigma=2\sum_i\log l_{ii} \). The listing evaluates the log density this
way, checks it against `scipy`, and counts how many of the simulated vectors
fall inside the \( 95\% \) ellipsoid. The fraction is \( 0.9500 \). The
average of \( \Delta^2 \) over the draws is \( 3.004 \) and its variance is
\( 5.994 \), against the \( \chi^2(3) \) values \( 3 \) and \( 6 \).

```{.python .run #cell-sampling-logpdf}
import numpy as np
from scipy import stats
rng = np.random.default_rng(301)
mu = np.array([1.0, -2.0, 0.5])
Sigma = np.array([[4.0, 1.2, -0.8],
                  [1.2, 2.0, 0.6],
                  [-0.8, 0.6, 1.5]])
N = 400_000
a = np.array([0.3, -1.0, 2.0])
q = [0.05, 0.25, 0.5, 0.75, 0.95]
exact = stats.norm.ppf(q, loc=a @ mu, scale=np.sqrt(a @ Sigma @ a))
cov_err = max(np.abs(np.cov(Y, rowvar=False) - Sigma).max() for Y in (Y_chol, Y_root))
def draw_normal(mu, A, size, rng):
    """Rows are independent copies of mu + A Z with Z ~ N(0, I)."""
    Z = rng.standard_normal((size, A.shape[1]))
    return mu + Z @ A.T

L = np.linalg.cholesky(Sigma)                  # lower triangular, L L^T = Sigma
lam, U = np.linalg.eigh(Sigma)
root = U @ np.diag(np.sqrt(lam)) @ U.T         # symmetric square root

Y_chol = draw_normal(mu, L, N, rng)
Y_root = draw_normal(mu, root, N, rng)
for Y in (Y_chol, Y_root):
    print(np.round(np.cov(Y, rowvar=False) - Sigma, 3))

def normal_logpdf(y, mu, Sigma):
    """log density of N(mu, Sigma) at the rows of y, Sigma positive definite."""
    L = np.linalg.cholesky(Sigma)
    w = np.linalg.solve(L, (y - mu).T)         # w = L^{-1}(y - mu), so w^T w = Delta^2
    D2 = np.sum(w**2, axis=0)
    logdet = 2.0 * np.sum(np.log(np.diag(L)))
    return -0.5 * (len(mu) * np.log(2 * np.pi) + logdet + D2), D2

logf, D2 = normal_logpdf(Y_chol, mu, Sigma)
inside = np.mean(D2 <= stats.chi2.ppf(0.95, df=3))
print(f"fraction inside the 95% ellipsoid: {inside:.4f}")
```

## Exercises

### A. Check your understanding

::: {#exr-mvn-identify-density}
[A1]

Show that \( f(y_1,y_2)=k\exp\{-\tfrac12(2y_1^2+2y_1y_2+y_2^2-4y_1-2y_2)\} \) is a bivariate
normal density for a suitable constant \( k \). Find \( k \), \( \bmu \), \( \bSigma \) and the
correlation.
:::

::: {.solution}
The quadratic part is \( \y\T\mathbf{P}\y \) with
\( \mathbf{P}=\begin{psmallmatrix}2&1\\1&1\end{psmallmatrix} \), which is positive definite with
\( \det\mathbf{P}=1 \). Matching \( (\y-\bmu)\T\mathbf{P}(\y-\bmu)=\y\T\mathbf{P}\y-2\bmu\T\mathbf{P}\y+\bmu\T\mathbf{P}\bmu \)
with the linear terms requires \( \mathbf{P}\bmu=(2,1)\T \), so \( \bmu=(1,0)\T \) and
\( \bmu\T\mathbf{P}\bmu=2 \). The exponent is \( -\tfrac12(\y-\bmu)\T\mathbf{P}(\y-\bmu)+1 \). By
@thm-mvn-density(b) this is the \( \Normal_2(\bmu,\bSigma) \) density with
\( \bSigma=\mathbf{P}^{-1}=\begin{psmallmatrix}1&-1\\-1&2\end{psmallmatrix} \) provided
\( k e=(2\pi)^{-1}(\det\bSigma)^{-1/2}=(2\pi)^{-1} \). So \( k=e^{-1}/(2\pi) \), and the correlation
is \( -1/\sqrt2 \).
:::

::: {#exr-mvn-ellipse}
[A2]

Let \( \bSigma \) have unit variances and correlation \( \rho \), with \( \lvert\rho\rvert<1 \). Find the
principal axes and half-lengths of the ellipse \( \Delta^2=c^2 \), show that its area is
\( \pi c^2\sqrt{1-\rho^2} \), and show that it contains probability \( 1-e^{-c^2/2} \).
:::

### B. Practice

::: {#exr-mvn-fourth-moments}
[B1]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) and \( \mathbf{X}=\Y-\bmu \). Expand the moment generating
function of \( \mathbf{X} \) in a power series to show that all odd-order moments of \( \mathbf{X} \) vanish and
\[
\E(X_iX_jX_kX_l)=\sigma_{ij}\sigma_{kl}+\sigma_{ik}\sigma_{jl}+\sigma_{il}\sigma_{jk}.
\]
Deduce that \( \Var(\mathbf{X}\T\A\mathbf{X})=2\tr(\A\bSigma\A\bSigma) \) for symmetric \( \A \).
:::

::: {.solution}
By @thm-mvn-mgf,
\( M_{\mathbf{X}}(\mathbf{t})=\exp(\tfrac12\mathbf{t}\T\bSigma\mathbf{t})=\sum_{m\ge0}(\mathbf{t}\T\bSigma\mathbf{t})^m/(2^mm!) \), a
series containing only terms of even degree in \( \mathbf{t} \). Mixed moments of order \( d \) are the
\( d \)th partial derivatives at \( \bzero \), which come from the degree-\( d \) terms, so odd-order moments
vanish. The degree-four term is \( p(\mathbf{t})=\tfrac18(\mathbf{t}\T\bSigma\mathbf{t})^2=T(\mathbf{t},\mathbf{t},\mathbf{t},\mathbf{t}) \), where
\( T(\bu,\bv,\bw,\mathbf{x})=\tfrac18(\bu\T\bSigma\bv)(\bw\T\bSigma\mathbf{x}) \) is \( 4 \)-linear. Differentiating
\( T(\mathbf{t},\mathbf{t},\mathbf{t},\mathbf{t}) \) once in each of the directions \( \mathbf{e}_i,\mathbf{e}_j,\mathbf{e}_k,\mathbf{e}_l \) gives
the sum of \( T \) over the \( 24 \) ways of placing these four vectors in the four slots. This
holds whether or not the indices are distinct, since it is the product rule applied to a
multilinear form. A placement contributes \( \tfrac18\sigma_{ab}\sigma_{cd} \), where the first two
slots receive the indices \( a,b \) and the last two receive \( c,d \). Each of the three ways of
splitting the four *directions* into two pairs, \( \{ij,kl\} \), \( \{ik,jl\} \), \( \{il,jk\} \), arises
from \( 8 \) placements (choose which pair goes first, and the order within each pair), which
gives the formula.
Then, for symmetric \( \A \),
\[
\E(\mathbf{X}\T\A\mathbf{X})^2=\sum_{i,j,k,l}a_{ij}a_{kl}\E(X_iX_jX_kX_l)
=(\tr\A\bSigma)^2+2\tr(\A\bSigma\A\bSigma),
\]
and subtracting \( (\E\mathbf{X}\T\A\mathbf{X})^2=(\tr\A\bSigma)^2 \) gives the variance.
:::

::: {#exr-mvn-singular-density}
[B2]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) with \( \rank(\bSigma)=r<n \), and use the notation of @prp-mvn-rank-representation.
Show that \( \mathbf{V}=\Q\T(\Y-\bmu) \) has a density on \( \Real^r \),
and that this density can be written as
\( (2\pi)^{-r/2}(\lambda_1\cdots\lambda_r)^{-1/2}\exp\{-\tfrac12(\y-\bmu)\T\bSigma^{+}(\y-\bmu)\} \)
with \( \y=\bmu+\Q\mathbf{v} \), where \( \bSigma^+=\Q\bLambda^{-1}\Q\T \) is the Moore–Penrose inverse.
:::

::: {#exr-mvn-noncentral-mahalanobis}
[B3]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) with \( \bSigma \) positive definite and let \( \boldsymbol{\upnu} \) be a fixed
vector. Show that \( (\Y-\boldsymbol{\upnu})\T\bSigma^{-1}(\Y-\boldsymbol{\upnu}) \) has the noncentral chi-squared
distribution \( \chi^2(n,\gamma) \) with
\( \gamma=(\bmu-\boldsymbol{\upnu})\T\bSigma^{-1}(\bmu-\boldsymbol{\upnu}) \), in the convention of
@def-qf-noncentral-chisq, and find its mean.
:::

### C. Going deeper

::: {#exr-mvn-smallest-region}
[C1]

Let \( \Y \) have density @eq-mvn-density, and let \( E=\{\y:\Delta^2(\y)\le q\} \). Show that every
Borel set \( A \) with \( \Pr(\Y\in A)\ge\Pr(\Y\in E) \) has volume at least that of \( E \).
*Hint:* compare \( \int_{A\setminus E}f \) with \( \int_{E\setminus A}f \).
:::

::: {.solution}
The density is a decreasing function of \( \Delta^2 \), so
with \( c=(2\pi)^{-n/2}\det(\bSigma)^{-1/2}e^{-q/2}>0 \) we have \( f\ge c \) on \( E \) and \( f<c \) off \( E \). We
may assume \( A \) has finite volume. Subtracting \( \Pr(\Y\in A\cap E) \) from both sides of
\( \Pr(\Y\in A)\ge\Pr(\Y\in E) \) gives \( \int_{A\setminus E}f\ge\int_{E\setminus A}f \). Hence
\[
c\operatorname{vol}(A\setminus E)\ge\int_{A\setminus E}f\ge\int_{E\setminus A}f\ge c\operatorname{vol}(E\setminus A),
\]
so \( \operatorname{vol}(A\setminus E)\ge\operatorname{vol}(E\setminus A) \). Adding
\( \operatorname{vol}(A\cap E) \) to both sides gives \( \operatorname{vol}(A)\ge\operatorname{vol}(E) \). The
argument uses only that \( E \) is a set on which \( f \) is at least as large as anywhere outside
it, so it applies to the level sets of any density.
:::
