# Definition through linear combinations

Nearly every exact distributional statement in the first half of this book has
the same form. The response is \( \Y=\X\bbeta+\be \) with normal errors, and the
object of interest is a linear function of \( \Y \) (a coefficient estimate, a
fitted value, a residual) or a quadratic function of \( \Y \) (a sum of squares).
A usable notion of normality for random vectors must therefore do two things.
It must survive every linear map, including maps that collapse dimensions. And it
must say when two linear functions of the same vector are independent.

The residual vector shows why the first requirement has to include singular
cases. In a model with \( p \) linearly independent columns its covariance matrix has
rank \( n-p \) ([Chapter 6](../ch06-projections/index.html)), so the residuals have no density on
\( \Real^n \). Defining multivariate normality through a density would
exclude it. This chapter uses a definition that does not mention densities at
all: a random vector is normal when every linear combination of its entries is
normal. The density, the moment generating function and the representation
\( \bmu+\A\Z \) all follow as theorems.

## The univariate normal, including the degenerate case

For \( \mu\in\Real \) and \( \sigma^2>0 \), a random variable \( X \) has the
\( \Normal(\mu,\sigma^2) \) distribution when it has density
\( (2\pi\sigma^2)^{-1/2}\exp\{-(x-\mu)^2/(2\sigma^2)\} \). We also allow
\( \sigma^2=0 \), and take \( \Normal(\mu,0) \) to be the point mass at \( \mu \). Treating
constants as normal variables with zero variance costs nothing and removes a
great many special cases later.

::: {#lem-mvn-univariate}
Let \( X\sim\Normal(\mu,\sigma^2) \) with \( \sigma^2\ge0 \).

::: {.enumerate options="label=(\alph*)"}
1. \( \E e^{sX}=\exp(\mu s+\tfrac12\sigma^2s^2) \) for every real \( s \).

2. A random variable whose moment generating function is
           \( \exp(\mu s+\tfrac12\sigma^2s^2) \) on an interval around \( 0 \) has the
           \( \Normal(\mu,\sigma^2) \) distribution.

3. \( aX+b\sim\Normal(a\mu+b,a^2\sigma^2) \), and if \( X_1\sim\Normal(\mu_1,\sigma_1^2) \)
           and \( X_2\sim\Normal(\mu_2,\sigma_2^2) \) are independent, then
           \( X_1+X_2\sim\Normal(\mu_1+\mu_2,\sigma_1^2+\sigma_2^2) \).
:::

:::

::: {.proof}
(a) For \( \sigma^2=0 \) this is \( e^{\mu s} \). Otherwise, completing the square,
\( sx-(x-\mu)^2/(2\sigma^2)=\mu s+\sigma^2s^2/2-(x-\mu-\sigma^2s)^2/(2\sigma^2) \),
and the remaining Gaussian integral equals \( (2\pi\sigma^2)^{1/2} \).
(b) is the uniqueness theorem for moment generating functions
(@thm-rv-mgf(a)), in its scalar form.
(c) The moment generating function of \( aX+b \) is
\( e^{bs}\E e^{(as)X}=\exp\{(a\mu+b)s+\tfrac12a^2\sigma^2s^2\} \). For the sum,
independence gives \( \E e^{s(X_1+X_2)}=\E e^{sX_1}\,\E e^{sX_2} \), which is the
normal moment generating function with the added parameters. Apply (b).
:::

## The definition

::: {#def-mvn}
[Multivariate normal distribution]

A random vector \( \Y=(Y_1,\dots,Y_n)\T \) has a **multivariate normal
distribution** if \( \mathbf{a}\T\Y \) has a univariate normal distribution (possibly
degenerate) for every \( \mathbf{a}\in\Real^n \). Such a \( \Y \) is also said to be
**jointly normal**, and its entries are said to be jointly normal random
variables.
:::

Taking \( \mathbf{a} \) to be a coordinate vector shows that each \( Y_i \) is normal, so
\( \E Y_i^2<\infty \). By the Cauchy–Schwarz inequality every product \( Y_iY_j \) is
integrable, so the mean vector \( \bmu=\E\Y \) and the covariance matrix
\( \bSigma=\Cov(\Y) \) of @def-rv-mean-cov exist. By @thm-rv-linear,
\( \E(\mathbf{a}\T\Y)=\mathbf{a}\T\bmu \) and \( \Var(\mathbf{a}\T\Y)=\mathbf{a}\T\bSigma\mathbf{a} \), and a normal
variable is determined by its mean and variance. Hence
\[
\mathbf{a}\T\Y\sim\Normal\bigl(\mathbf{a}\T\bmu,\ \mathbf{a}\T\bSigma\mathbf{a}\bigr)
\qquad\text{for every }\mathbf{a}\in\Real^n .
\]{#eq-mvn-linear-combination}

We write \( \Y\sim\Normal_n(\bmu,\bSigma) \). The notation claims that the joint
distribution of \( \Y \) depends only on \( \bmu \) and \( \bSigma \). That is true, but it
is not obvious: the definition fixes the law of each one-dimensional
projection, and it is not yet clear that these laws pin down the joint law. The
claim is proved in @thm-mvn-mgf. Until then, \( \Normal_n(\bmu,\bSigma) \)
should be read as “normal, with mean \( \bmu \) and covariance \( \bSigma \)”.

## Existence

::: {#thm-mvn-existence}
[Existence and the representation \( \bmu+\A\Z \)]

Let \( Z_1,\dots,Z_m \) be independent \( \Normal(0,1) \) variables and
\( \Z=(Z_1,\dots,Z_m)\T \).

::: {.enumerate options="label=(\alph*)"}
1. \( \Z\sim\Normal_m(\bzero,\I_m) \).

2. For any \( n\times m \) matrix \( \A \) and \( \bmu\in\Real^n \), \( \Y=\bmu+\A\Z \) is
           multivariate normal with mean \( \bmu \) and covariance \( \A\A\T \).

3. For every \( \bmu\in\Real^n \) and every nonnegative definite \( n\times n \)
           matrix \( \bSigma \) there is a random vector with the
           \( \Normal_n(\bmu,\bSigma) \) distribution, for example
           \( \bmu+\bSigma^{1/2}\Z \) with \( m=n \).
:::

No other matrices occur: the covariance matrix of a normal vector is
nonnegative definite.
:::

::: {.proof}
(a) \( \mathbf{a}\T\Z=\sum_i a_iZ_i \) is a sum of independent normal variables
\( a_iZ_i\sim\Normal(0,a_i^2) \), which is normal by @lem-mvn-univariate(c).
(b) \( \mathbf{a}\T\Y=\mathbf{a}\T\bmu+(\A\T\mathbf{a})\T\Z \) is a constant plus a linear
combination of \( \Z \), hence normal by (a). Its moments come from
@thm-rv-linear: \( \E\Y=\bmu \) and \( \Cov(\Y)=\A\I\A\T=\A\A\T \).
(c) By @thm-mat-square-root a nonnegative definite \( \bSigma \) has a symmetric
nonnegative definite square root with \( \bSigma^{1/2}\bSigma^{1/2}=\bSigma \). Take
\( \A=\bSigma^{1/2} \) in (b). The last sentence is @thm-rv-cov-nnd.
:::

Any factorization \( \bSigma=\A\A\T \) serves in (b). The symmetric square root is
convenient in proofs. The Cholesky factor, a lower triangular \( \bL \) with
\( \bL\bL\T=\bSigma \), is cheaper to compute and is what software uses. A
factor with fewer columns than rows is natural when \( \bSigma \) is singular.
Once @thm-mvn-mgf is proved, all of these give the same distribution.
The first listing checks this by simulation. It draws
\( 400{,}000 \) vectors from a three-dimensional normal distribution with each of
two factors. Both sample covariance matrices agree with \( \bSigma \) to within
\( 0.016 \) in every entry, and the quantiles of a linear
combination match the exact normal quantiles.

```{.python .run #cell-sampling-draw}
import numpy as np
from scipy import stats
rng = np.random.default_rng(301)
mu = np.array([1.0, -2.0, 0.5])
Sigma = np.array([[4.0, 1.2, -0.8],
                  [1.2, 2.0, 0.6],
                  [-0.8, 0.6, 1.5]])
N = 400_000

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
```

## Singular normal distributions

When \( \bSigma \) is singular, @def-mvn still applies, and the distribution
is called **singular** or degenerate. If \( \bSigma\bv=\bzero \), then
\( \Var(\bv\T\Y)=\bv\T\bSigma\bv=0 \), so \( \bv\T\Y=\bv\T\bmu \) with probability
one. Every vector in \( \Null(\bSigma) \) gives such a linear constraint, and
together they say that \( \Y \) lies in the affine subspace \( \bmu+\C(\bSigma) \) with
probability one (@thm-rv-cov-nnd). A subspace of dimension
\( \rank(\bSigma)<n \) has Lebesgue measure zero, so no density on \( \Real^n \) is
possible.

::: {#exm-mvn-overlapping-sums}
[Three overlapping sums]

Let \( \Z\sim\Normal_3(\bzero,\I) \) and \( \bU=(Z_1+Z_2,\ Z_2-Z_3,\ Z_1+Z_3)\T \). Then
\( \bU=\A\Z \) with rows \( (1,1,0) \), \( (0,1,-1) \) and \( (1,0,1) \), so \( \bU \) is normal with
mean zero and covariance
\[
\A\A\T=\begin{pmatrix}2&1&1\\1&2&-1\\1&-1&2\end{pmatrix}.
\]
The vector \( (1,-1,-1)\T \) lies in \( \Null(\A\A\T) \), and the rank is \( 2 \).
Correspondingly \( U_1-U_2-U_3=0 \) for every outcome, not merely with
probability one. Each \( U_i \) is \( \Normal(0,2) \), each pair \( (U_i,U_j) \) has a
nonsingular covariance matrix (determinant \( 3 \)) and a bivariate density, and yet
the triple lives on a plane.
:::

In regression the singular case is the rule. With \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \),
the fitted values and the residuals are normal vectors of ranks \( \rank(\X) \) and
\( n-\rank(\X) \) (@exm-mvn-ls-estimator). The last part of the listing draws a
normal vector in \( \Real^3 \) from a \( 3\times2 \) factor. Every draw lies on the plane
\( \bmu+\C(\A) \) to within \( 2.7\times 10^{-15} \), which is rounding error.

## Why linear combinations suffice

@def-mvn describes a joint distribution entirely through its
one-dimensional projections. That this loses no information is a general fact
about random vectors, the Cramér–Wold device proved with characteristic
functions in @prp-rv-cramer-wold: if \( \mathbf{a}\T\Y \) and \( \mathbf{a}\T\W \) have the same
distribution for every \( \mathbf{a}\in\Real^n \), then \( \Y \) and \( \W \) have the same
distribution. We do not need it in this generality. For normal vectors the
moment generating function is finite everywhere, and the next section runs the
same argument with it. The device explains why the definition is natural. It also
warns that all directions \( \mathbf{a} \) matter. Normality of the \( n \) coordinate
variables, which are only \( n \) directions, is not enough
([Section 3.3](03-linear.html#sec-mvn-normal-marginals) and @exr-mvn-sign-flip-directions).

## Exercises

### A. Check your understanding

::: {#exr-mvn-independent-entries}
[A1]

Let \( Y_1,\dots,Y_n \) be independent with \( Y_i\sim\Normal(\mu_i,\sigma_i^2) \). Show directly from
@def-mvn that \( \Y=(Y_1,\dots,Y_n)\T \) is multivariate normal, and find its mean and
covariance. Show also that a constant vector \( \mathbf{c} \) is \( \Normal_n(\mathbf{c},\bzero) \).
:::

::: {#exr-mvn-cholesky}
[A2]

For the bivariate normal of @exm-mvn-bivariate, find the Cholesky factor \( \bL \) of
\( \bSigma \) and write \( Y_1 \) and \( Y_2 \) in terms of two independent standard normals. Do the
same with the symmetric square root when \( \sigma_1=\sigma_2=1 \). Why do the two
representations define the same distribution?
:::

### B. Practice

::: {#exr-mvn-rank-one}
[B1]

Let \( Z\sim\Normal(0,1) \) and \( \Y=(Z,\ 2Z-1,\ 3-Z)\T \). Show that \( \Y \) is normal, find \( \bmu \)
and \( \bSigma \), and find the rank of \( \bSigma \). Describe the set in which \( \Y \) lies, and write
it as the solution set of two linear equations.
:::

::: {.solution}
\( \Y=\bmu+\A Z \) with \( \bmu=(0,-1,3)\T \) and
\( \A=(1,2,-1)\T \), so \( \Y \) is normal by @thm-mvn-existence with
\[
\bSigma=\A\A\T=\begin{pmatrix}1&2&-1\\2&4&-2\\-1&-2&1\end{pmatrix},\qquad\rank(\bSigma)=1 .
\]
\( \Y \) lies on the line \( \{\bmu+s\A:s\in\Real\} \). The null space of \( \bSigma \) is
\( \A\perpc \), spanned by \( (2,-1,0)\T \) and \( (1,0,1)\T \). The corresponding constraints
\( \bv\T\Y=\bv\T\bmu \) are \( 2Y_1-Y_2=1 \) and \( Y_1+Y_3=3 \), and both hold identically.
:::

### C. Going deeper

::: {#exr-mvn-factor-uniqueness}
[C1]

The text says that every factorization \( \bSigma=\A\A\T \) gives the same distribution, and proves
it later from the moment generating function. Here is the algebraic reason.

::: {.enumerate options="label=(\alph*)"}
1. Let \( \A \) and \( \B \) be \( n\times m \) with \( \A\A\T=\B\B\T \). Show that \( \B=\A\Q \) for some
           \( m\times m \) orthogonal \( \Q \). *Hint:* \( \norm{\A\T\mathbf{x}}=\norm{\B\T\mathbf{x}} \) for every
           \( \mathbf{x} \), so \( \A\T\mathbf{x}\mapsto\B\T\mathbf{x} \) is a well-defined linear isometry of
           \( \C(\A\T) \) onto \( \C(\B\T) \); extend it to \( \Real^m \).

2. Deduce that \( \bmu+\A\Z \) and \( \bmu+\B\Z \) have the same distribution, with no appeal
           to @thm-mvn-mgf. *Hint:* the density of \( \Z \) depends on \( \bz \) only through
           \( \norm{\bz} \). Conclude for the Cholesky factor and the symmetric square root of a
           positive definite \( \bSigma \).

3. Show that \( \rank(\A)=\rank(\bSigma) \) for every such factor, so that no representation
           \( \bmu+\A\Z \) uses fewer than \( \rank(\bSigma) \) standard normal variables, and that
           \( \rank(\bSigma) \) of them suffice.
:::

:::

::: {.solution}
(a) For every \( \mathbf{x} \),
\( \norm{\A\T\mathbf{x}}^2=\mathbf{x}\T\A\A\T\mathbf{x}=\mathbf{x}\T\B\B\T\mathbf{x}=\norm{\B\T\mathbf{x}}^2 \). If
\( \A\T\mathbf{x}=\A\T\mathbf{x}' \) then \( \norm{\B\T(\mathbf{x}-\mathbf{x}')}=\norm{\A\T(\mathbf{x}-\mathbf{x}')}=0 \), so the
assignment \( \Q_1:\A\T\mathbf{x}\mapsto\B\T\mathbf{x} \) is a well-defined linear map from \( \C(\A\T) \) onto
\( \C(\B\T) \), and it preserves norms, hence inner products by polarization. The two column
spaces have the same dimension \( r=\rank(\A)=\rank(\A\A\T)=\rank(\B\B\T)=\rank(\B) \) by
@prp-mat-rank-product, so their orthogonal complements in \( \Real^m \) have dimension
\( m-r \); extend \( \Q_1 \) by mapping an orthonormal basis of \( \C(\A\T)\perpc \) onto one of
\( \C(\B\T)\perpc \). The extension is an orthogonal matrix with \( \Q_1\A\T=\B\T \), so
\( \B=\A\Q_1\T \) and \( \Q=\Q_1\T \) works.

(b) By (a), \( \B\Z=\A(\Q\Z) \). The density of \( \Z \) is
\( (2\pi)^{-m/2}\exp(-\norm{\bz}^2/2) \), and the substitution \( \bz=\Q\T\mathbf{w} \) has Jacobian
\( \lvert\det\Q\T\rvert=1 \) and leaves \( \norm{\bz} \) unchanged, so \( \Q\Z \) has that same density.
Hence \( \B\Z \) and \( \A\Z \) have the same distribution, and so do \( \bmu+\A\Z \) and
\( \bmu+\B\Z \). With \( \bSigma=\bL\bL\T=\bSigma^{1/2}\bSigma^{1/2} \) both factors are square, so the
Cholesky and square-root constructions differ by an orthogonal factor and produce the same
law.

(c) \( \rank(\A)=\rank(\A\A\T)=\rank(\bSigma) \) by @prp-mat-rank-product, and an \( n\times m \)
matrix has rank at most \( m \), so \( m\ge\rank(\bSigma) \). The factor
\( \A=\boldsymbol{\Gamma}_r\bLambda_r^{1/2} \) built from the \( r=\rank(\bSigma) \) positive eigenvalues
attains it.
:::

::: {#exr-mvn-normal-limits}
[C2]

Suppose \( \Y_k\sim\Normal_n(\bmu_k,\bSigma_k) \) converges in distribution to a random vector \( \Y \).

::: {.enumerate options="label=(\alph*)"}
1. Treat first \( n=1 \): if \( \Normal(m_k,v_k) \) converges in distribution, show that \( m_k\to m \)
           and \( v_k\to v \) for finite limits with \( v\ge0 \), and that the limit law is
           \( \Normal(m,v) \), degenerate when \( v=0 \). *Hint:* the characteristic functions are
           \( \exp(itm_k-t^2v_k/2) \); a limit that is continuous at \( 0 \) keeps \( v_k \) from escaping
           to \( \infty \), and \( e^{itm_k} \) can converge for every \( t \) only if \( m_k \) does.

2. Deduce for general \( n \), using @eq-mvn-linear-combination and @def-mvn, that \( \bmu_k\to\bmu \),
           \( \bSigma_k\to\bSigma \) entrywise and \( \Y\sim\Normal_n(\bmu,\bSigma) \). So the normal family,
           singular members included, is closed under limits in distribution.

3. Show conversely that every singular \( \Normal_n(\bmu,\bSigma) \) is a limit of nonsingular
           ones. The degenerate distributions this section insists on admitting are exactly the
           ones that would otherwise be missing.
:::

:::

::: {.solution}
(a) Let \( \phi \) be the characteristic function of the limit, so that
\( \exp(itm_k-t^2v_k/2)\to\phi(t) \) for every \( t \), and \( \phi \) is continuous at \( 0 \) with
\( \phi(0)=1 \). Taking absolute values, \( e^{-t^2v_k/2}\to\lvert\phi(t)\rvert \). Choose \( t_0>0 \) with
\( \lvert\phi(t_0)\rvert>0 \); then \( v_k\to v=-2\log\lvert\phi(t_0)\rvert/t_0^2 \), finite and
nonnegative. Hence \( e^{itm_k}\to\psi(t)=\phi(t)e^{t^2v/2} \) for every \( t \), with
\( \lvert\psi\rvert\equiv1 \). If \( \lvert m_k\rvert\to\infty \) along a subsequence, then
\( \int_0^T e^{itm_k}dt=(e^{iTm_k}-1)/(im_k)\to0 \) for every \( T \), while dominated convergence
makes the limit \( \int_0^T\psi \); so \( \psi=0 \) almost everywhere, contradicting
\( \lvert\psi\rvert=1 \). So \( (m_k) \) is bounded, and if one subsequence converged to \( m \) and
another to \( m' \), then \( e^{itm}=e^{itm'} \) for all \( t \), forcing \( m=m' \). Thus \( m_k\to m \) and
\( \phi(t)=\exp(itm-t^2v/2) \), the characteristic function of \( \Normal(m,v) \).

(b) For each \( \mathbf{a} \), the continuous mapping theorem gives
\( \mathbf{a}\T\Y_k\to\mathbf{a}\T\Y \) in distribution, and
\( \mathbf{a}\T\Y_k\sim\Normal(\mathbf{a}\T\bmu_k,\mathbf{a}\T\bSigma_k\mathbf{a}) \) by @eq-mvn-linear-combination. By
(a), \( \mathbf{a}\T\bmu_k \) and \( \mathbf{a}\T\bSigma_k\mathbf{a} \) converge and \( \mathbf{a}\T\Y \) is normal. Taking
\( \mathbf{a}=\mathbf{e}_i \) gives \( \bmu_k\to\bmu \) entrywise and \( \sigma_{k,ii}\to\sigma_{ii} \); taking
\( \mathbf{a}=\mathbf{e}_i+\mathbf{e}_j \) and subtracting gives \( \sigma_{k,ij}\to\sigma_{ij} \). The limit \( \bSigma \)
is nonnegative definite as a limit of nonnegative definite matrices, and every \( \mathbf{a}\T\Y \) is
normal with mean \( \mathbf{a}\T\bmu \) and variance \( \mathbf{a}\T\bSigma\mathbf{a} \), which is @def-mvn.

(c) \( \bSigma+\varepsilon\I \) is positive definite for \( \varepsilon>0 \), and the characteristic
functions \( \exp(i\mathbf{t}\T\bmu-\tfrac12\mathbf{t}\T\bSigma\mathbf{t}-\tfrac12\varepsilon\norm{\mathbf{t}}^2) \)
converge pointwise as \( \varepsilon\downarrow0 \) to that of \( \Normal_n(\bmu,\bSigma) \).
:::
