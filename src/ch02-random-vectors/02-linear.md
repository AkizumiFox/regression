# Linear and affine transformations

Almost every statistic in a linear model is a linear function of the response: least
squares coefficients, fitted values, residuals, contrasts. Their means and covariance
matrices follow from one theorem.

## The transformation rules

::: {#thm-rv-linear}
[Moments of affine transformations]

Let \( \Y \) be a \( p\times1 \) random vector with finite second moments, \( \E(\Y)=\bmu \) and
\( \Cov(\Y)=\bSigma \). Let \( \A \) be \( k\times p \), \( \B \) be \( m\times p \), and \( \bb \) be \( k\times1 \),
all constant. Then

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\A\Y+\bb)=\A\bmu+\bb \);

2. \( \Cov(\A\Y+\bb)=\A\bSigma\A\T \);

3. \( \Cov(\A\Y,\B\Y)=\A\bSigma\B\T \).
:::

More generally, if \( \bU \) and \( \V \) have finite second moments and \( \A \), \( \B \), \( \bb \),
\( \bm c \) are constants of conformable sizes, then
\( \Cov(\A\bU+\bb,\,\B\V+\bm c)=\A\,\Cov(\bU,\V)\,\B\T \).
:::

::: {.proof}
(a) is @prp-rv-expectation-linear(a). For the general statement, the centred
vectors are \( \A\bU+\bb-\E(\A\bU+\bb)=\A(\bU-\E\bU) \) and similarly \( \B(\V-\E\V) \). Hence
\[
\Cov(\A\bU+\bb,\B\V+\bm c)
=\E\bigl[\A(\bU-\E\bU)(\V-\E\V)\T\B\T\bigr]
=\A\,\E\bigl[(\bU-\E\bU)(\V-\E\V)\T\bigr]\B\T ,
\]
again by @prp-rv-expectation-linear(a). Parts (c) and (b) are the cases
\( \bU=\V=\Y \), and \( \B=\A \) in (b).
:::

For a single linear combination, \( \A=\bm a\T \) is a row vector, and (b) and (c) give
\[
\Var(\bm a\T\Y)=\bm a\T\bSigma\bm a=\sum_{i=1}^p\sum_{j=1}^p a_ia_j\sigma_{ij},
\qquad
\Cov(\bm a\T\Y,\bb\T\Y)=\bm a\T\bSigma\bb .
\]{#eq-rv-linear-combination}

Every variance and covariance enters \( \Var(\bm a\T\Y) \), weighted by the products of the
coefficients. Conversely, \( \bSigma \) is determined by the variances of linear
combinations, because \( \sigma_{ij}=\tfrac12[\Var(Y_i+Y_j)-\Var(Y_i)-\Var(Y_j)] \).

::: {#prp-rv-sums}
[Covariance of sums]

Let \( \bU,\bU_2,\V \) have finite second moments, with \( \bU \) and \( \bU_2 \) of the same size in
(a), and \( \bU \) and \( \V \) of the same size in (b). Then

::: {.enumerate options="label=(\alph*)"}
1. \( \Cov(\bU+\bU_2,\V)=\Cov(\bU,\V)+\Cov(\bU_2,\V) \);

2. \( \Cov(\bU+\V)=\Cov(\bU)+\Cov(\bU,\V)+\Cov(\V,\bU)+\Cov(\V) \);

3. if \( \Y_1,\dots,\Y_n \) are pairwise uncorrelated \( p\times1 \) vectors and
           \( c_1,\dots,c_n \) are constants, then \( \Cov\bigl(\sum_ic_i\Y_i\bigr)=\sum_ic_i^2\Cov(\Y_i) \).
           In particular, if every \( \Cov(\Y_i)=\bSigma \), then \( \Cov(\bar{\Y})=\bSigma/n \).
:::

:::

::: {.proof}
(a) Centre each vector and use linearity of the outer product and of \( \E \). (b) Apply (a)
in each argument; note that the two cross terms are transposes of each other, not equal.
(c) By (a) and (b), applied repeatedly,
\( \Cov\bigl(\sum_ic_i\Y_i\bigr)=\sum_i\sum_jc_ic_j\Cov(\Y_i,\Y_j) \), using
\( \Cov(c\bU,\V)=c\Cov(\bU,\V) \) from @thm-rv-linear. The terms with \( i\ne j \) vanish.
:::

::: {#exm-rv-centering}
[Centring a sample]

Let \( \Y \) be \( n\times1 \) with \( \Cov(\Y)=\sigma^2\I_n \), and write \( \bm J=\bone\bone\T \). The mean is
\( \bar{Y}=n^{-1}\bone\T\Y \), and the vector of deviations is \( \bm d=(\I-n^{-1}\bm J)\Y \). The
matrix \( \I-n^{-1}\bm J \) is symmetric and idempotent, so @thm-rv-linear gives
\[
\Cov(\bm d)=\sigma^2(\I-n^{-1}\bm J),
\qquad
\Cov(\bar{Y},\bm d)=\sigma^2 n^{-1}\bone\T(\I-n^{-1}\bm J)=\bzero\T .
\]
The deviations have variances \( \sigma^2(1-1/n) \) and are negatively correlated, with
correlation \( -1/(n-1) \) between each pair, because they are forced to sum to zero. And each
deviation is uncorrelated with the mean. Under normality, uncorrelated becomes
independent (@exm-mvn-mean-deviations), which is the first half of
@cor-qf-sample-variance.
:::

## Covariance matrices are nonnegative definite

A variance cannot be negative, and by @eq-rv-linear-combination this simple
fact constrains the whole matrix. Recall from @def-mat-nnd that a symmetric
\( \bSigma \) is nonnegative definite if \( \bm a\T\bSigma\bm a\ge0 \) for all \( \bm a \), and positive
definite if in addition \( \bm a\T\bSigma\bm a>0 \) for \( \bm a\ne\bzero \).

::: {#thm-rv-cov-nnd}
[Covariance matrices]

Let \( \Y \) have finite second moments, mean \( \bmu \) and covariance matrix \( \bSigma \).

::: {.enumerate options="label=(\alph*)"}
1. \( \bSigma \) is nonnegative definite.

2. Conversely, every \( p\times p \) nonnegative definite matrix is the covariance matrix
           of some random vector.

3. For a constant \( \bm a \), \( \bSigma\bm a=\bzero \) iff \( \bm a\T\Y=\bm a\T\bmu \) with
           probability one.

4. \( \Pr(\Y-\bmu\in\C(\bSigma))=1 \). If \( \mathcal S \) is any subspace with
           \( \Pr(\Y-\bmu\in\mathcal S)=1 \), then \( \C(\bSigma)\subseteq\mathcal S \).
:::

In particular, \( \bSigma \) is singular iff some nontrivial linear combination of \( \Y \) is
constant with probability one. Then \( \Y \) lies, with probability one, in the affine subspace
\( \bmu+\C(\bSigma) \) of dimension \( \rank(\bSigma)<p \), and this is the smallest such flat.
:::

::: {.proof}
(a) \( \bSigma \) is symmetric, and \( \bm a\T\bSigma\bm a=\Var(\bm a\T\Y)\ge0 \) by
@eq-rv-linear-combination.

(b) Let \( \bSigma \) be nonnegative definite and let \( \bSigma^{1/2} \) be its symmetric
nonnegative definite square root (@thm-mat-square-root). Let \( \Z \) have independent
components with mean \( 0 \) and variance \( 1 \), for instance independent standard normal
variables, so that \( \Cov(\Z)=\I \). Then
\( \Cov(\bSigma^{1/2}\Z)=\bSigma^{1/2}\I\bSigma^{1/2}=\bSigma \).

(c) If \( \bSigma\bm a=\bzero \), the random variable \( \bm a\T(\Y-\bmu) \) has mean zero and
variance \( \bm a\T\bSigma\bm a=0 \), so it is zero with probability one. Conversely, if
\( \bm a\T(\Y-\bmu)=0 \) with probability one, then
\( \bSigma\bm a=\E\bigl[(\Y-\bmu)\,(\Y-\bmu)\T\bm a\bigr]=\bzero \).

(d) Let \( \bm a_1,\dots,\bm a_s \) be a basis of \( \Null(\bSigma) \). By (c), each event
\( \{\bm a_j\T(\Y-\bmu)=0\} \) has probability one, and so does their intersection. On that
intersection \( \Y-\bmu \) is orthogonal to \( \Null(\bSigma) \). Because \( \bSigma \) is symmetric,
\( \Null(\bSigma) \) is the orthogonal complement of \( \C(\bSigma) \), and so \( \Y-\bmu\in\C(\bSigma) \).
Now let \( \mathcal S \) be a subspace with \( \Pr(\Y-\bmu\in\mathcal S)=1 \). For every
\( \bm a\in\mathcal S\perpc \) we have \( \bm a\T(\Y-\bmu)=0 \) with probability one, hence
\( \bm a\in\Null(\bSigma)=\C(\bSigma)\perpc \) by (c). So
\( \mathcal S\perpc\subseteq\C(\bSigma)\perpc \), and taking orthogonal complements gives
\( \C(\bSigma)\subseteq\mathcal S \). The final statements follow, with
\( \dim\C(\bSigma)=\rank(\bSigma) \). For minimality among all flats, suppose
\( \Pr(\Y\in\bm c+\mathcal S)=1 \) for a vector \( \bm c \) and a subspace \( \mathcal S \). Every
\( \bm a\in\mathcal S\perpc \) has \( \bm a\T(\Y-\bm c)=0 \) with probability one, and taking
expectations gives \( \bm a\T(\bmu-\bm c)=0 \). Hence \( \bmu-\bm c\in\mathcal S \), so
\( \bm c+\mathcal S=\bmu+\mathcal S \), and the subspace case applies.
:::

A covariance matrix is therefore positive definite exactly when no linear combination of
the components is degenerate. Singular covariance matrices are not pathological. They
arise whenever the components satisfy an exact linear constraint, as with the deviations of
@exm-rv-centering, whose covariance matrix has rank \( n-1 \) because \( \bone\T\bm d=0 \). The
residual vector of a regression is another case: it is orthogonal to every column of \( \X \),
and when \( \Cov(\Y)=\sigma^2\I \) its covariance matrix has rank \( n-\rank(\X) \)
([Section 2.2](02-linear.html#sec-rv-linear-model), [Chapter 6](../ch06-projections/index.html)).

::: {#exm-rv-multinomial}
[Multinomial counts]

In \( m \) independent trials, each outcome falls into one of \( k \) categories with
probabilities \( \bm\pi=(\pi_1,\dots,\pi_k)\T \), where \( \sum_j\pi_j=1 \). Let \( \bm Z_t \) be the indicator
vector of trial \( t \), with a single \( 1 \) in the position of its category. Then
\( \E(\bm Z_t)=\bm\pi \), and \( \bm Z_t\bm Z_t\T=\diag(\bm Z_t) \) because the indicator has only one
nonzero entry. So \( \Cov(\bm Z_t)=\diag(\bm\pi)-\bm\pi\bm\pi\T \). The count vector
\( \Y=\sum_t\bm Z_t \) is a sum of independent vectors, and @prp-rv-sums(c) gives
\[
\bSigma=\Cov(\Y)=m\bigl(\diag(\bm\pi)-\bm\pi\bm\pi\T\bigr).
\]
Because \( \bSigma\bone=m(\bm\pi-\bm\pi\,\bm\pi\T\bone)=\bzero \), the matrix is singular, and
@thm-rv-cov-nnd(c) recovers the constraint \( \bone\T\Y=m \). For \( m=20 \) and
\( \bm\pi=(0.5,0.3,0.2)\T \) the eigenvalues are \( 7.762 \), \( 4.638 \)
and \( 0 \). All \( 5000 \) simulated count vectors lie in the plane
\( \bmu+\C(\bSigma)=\{\y:\bone\T\y=m\} \).
:::

Linear maps cannot increase the rank of a covariance matrix, since
\( \rank(\A\bSigma\A\T)\le\rank(\bSigma) \) (@exr-rv-cov-rank). They can, however, discard
degenerate directions. If \( \Y=(Z,Z)\T \) with \( \Var(Z)=1 \), then \( \bSigma \) has all entries equal
to \( 1 \) and is singular, while \( \A=(1,0) \) gives \( \Var(\A\Y)=1>0 \). By @thm-rv-linear, the
covariance matrix of \( \A\Y \) is \( \A\bSigma\A\T \). If \( \bSigma \) is positive definite and \( \A \) has
full row rank \( k\le p \), then for \( \bm b\neq\bzero \) the vector \( \A\T\bm b \) is nonzero and
\( \bm b\T\A\bSigma\A\T\bm b>0 \). So \( \A\bSigma\A\T \) is positive definite. If \( k>p \), then
\( \A\bSigma\A\T \) has rank at most \( p \) and is singular. A vector of more than \( p \) linear
combinations of \( p \) variables always satisfies a linear constraint.

## Standardization and Mahalanobis distance

A scalar with mean \( \mu \) and variance \( \sigma^2>0 \) is standardized by \( (Y-\mu)/\sigma \). The
vector version must remove correlation as well as scale.

::: {#prp-rv-whitening}
[Whitening]

Let \( \bSigma=\Cov(\Y) \) be positive definite, with \( \E(\Y)=\bmu \).

::: {.enumerate options="label=(\alph*)"}
1. \( \Z=\bSigma^{-1/2}(\Y-\bmu) \) has \( \E(\Z)=\bzero \) and \( \Cov(\Z)=\I_p \).

2. A \( p\times p \) matrix \( \W \) satisfies \( \W\bSigma\W\T=\I_p \) iff \( \W=\Q\bSigma^{-1/2} \) for
           some orthogonal \( \Q \). In particular, if \( \bSigma=\bL\bL\T \) is the Cholesky
           factorization, \( \W=\bL^{-1} \) is such a matrix.

3. For every \( \W \) as in (b),
           \( \norm{\W(\Y-\bmu)}^2=(\Y-\bmu)\T\bSigma^{-1}(\Y-\bmu) \).
:::

:::

::: {.proof}
(a) By @thm-rv-linear,
\( \Cov(\Z)=\bSigma^{-1/2}\bSigma\bSigma^{-1/2}=\I \), since \( \bSigma^{-1/2} \) is symmetric and
commutes with \( \bSigma \). (b) If \( \W=\Q\bSigma^{-1/2} \), then \( \W\bSigma\W\T=\Q\Q\T=\I \).
Conversely, if \( \W\bSigma\W\T=\I \), put \( \Q=\W\bSigma^{1/2} \). Then \( \Q\Q\T=\W\bSigma\W\T=\I \), so the
square matrix \( \Q \) is orthogonal, and \( \W=\Q\bSigma^{-1/2} \). For the Cholesky factor,
\( \bL^{-1}\bSigma(\bL\T)^{-1}=\bL^{-1}\bL\bL\T(\bL\T)^{-1}=\I \). (c) \( \W\T\W=\bSigma^{-1/2}\Q\T\Q\bSigma^{-1/2}=\bSigma^{-1} \).
:::

The quantity in (c),
\[
\Delta^2(\Y)=(\Y-\bmu)\T\bSigma^{-1}(\Y-\bmu),
\]{#eq-rv-mahalanobis}

is the squared **Mahalanobis distance** of \( \Y \) from \( \bmu \). It measures distance in
units of the covariance structure. Its level sets are ellipsoids with axes along the
eigenvectors of \( \bSigma \), with half-lengths proportional to the square roots of the
eigenvalues. Whitening maps these ellipsoids to spheres. Part (b) says that whitening is
unique only up to a rotation or reflection. Different choices give different coordinates,
but the same distances.

::: {#exm-rv-whitening}
[Two whitenings]

Take \( \bmu=(3,1)\T \) and
\( \bSigma=\begin{psmallmatrix}4&2.4\\2.4&2.25\end{psmallmatrix} \), so the standard deviations
are \( 2 \) and \( 1.5 \) and the correlation is \( 0.8 \). The eigenvalues are \( 5.680 \) and
\( 0.570 \). The symmetric and Cholesky whitening matrices are
\[
\bSigma^{-1/2}=\begin{pmatrix}0.717&-0.425\\-0.425&1.027\end{pmatrix},
\qquad
\bL^{-1}=\begin{pmatrix}0.500&0\\-0.667&1.111\end{pmatrix}.
\]
The Cholesky version standardizes \( Y_1 \) alone and then removes from \( Y_2 \) its linear
dependence on \( Y_1 \), in the spirit of Gram–Schmidt. The symmetric version treats both
coordinates alike and moves points as little as possible
(@exr-rv-whitening-closest). [Figure 2.2.1](02-linear.html#fig-rv-whitening) applies both to
\( 400 \) simulated points. The two whitened clouds differ by a rotation through
\( 22.5^\circ \), and every point keeps its Mahalanobis distance.
:::

::: {when-format="html"}
![**Figure 2.2.1.** Whitening a correlated sample. (a) Draws with \( \Cov(\Y)=\bSigma \) and the
Mahalanobis contours \( \Delta=1,2,3 \). (b) The symmetric whitening
\( \bSigma^{-1/2}(\Y-\bmu) \). (c) The Cholesky whitening \( \bL^{-1}(\Y-\bmu) \). The contours
become circles in both, and the three coloured points keep their distance from the
origin. Panels (b) and (c) differ by a rotation (@prp-rv-whitening).](whitening.svg){#fig-rv-whitening width=100%}
:::

::: {when-format="pdf"}
![Whitening a correlated sample. (a) Draws with \( \Cov(\Y)=\bSigma \) and the
Mahalanobis contours \( \Delta=1,2,3 \). (b) The symmetric whitening
\( \bSigma^{-1/2}(\Y-\bmu) \). (c) The Cholesky whitening \( \bL^{-1}(\Y-\bmu) \). The contours
become circles in both, and the three coloured points keep their distance from the
origin. Panels (b) and (c) differ by a rotation (@prp-rv-whitening).](whitening.pdf){width=100%}
:::

```{.python .run #cell-whitening-whiten}
import numpy as np
rng = np.random.default_rng(2)

mu = np.array([3.0, 1.0])
Sigma = np.array([[4.0, 2.4],
                  [2.4, 2.25]])                   # correlation 2.4 / (2 * 1.5) = 0.8

lam, U = np.linalg.eigh(Sigma)                    # spectral decomposition
W_sym = U @ np.diag(lam ** -0.5) @ U.T            # Sigma^{-1/2}
L = np.linalg.cholesky(Sigma)                     # Sigma = L L', L lower triangular
W_chol = np.linalg.inv(L)                         # L^{-1}

n = 400
Y = mu + rng.standard_normal((n, 2)) @ L.T        # any draws with Cov = Sigma will do
Z_sym = (Y - mu) @ W_sym.T
Z_chol = (Y - mu) @ W_chol.T

Q = W_chol @ np.linalg.inv(W_sym)                 # L^{-1} Sigma^{1/2}
print("Q'Q =", np.round(Q.T @ Q, 12))
```

Whitening runs in both directions. If \( \Cov(\Z)=\I \), then \( \bmu+\bSigma^{1/2}\Z \), or
\( \bmu+\bL\Z \), has mean \( \bmu \) and covariance \( \bSigma \). This is how correlated vectors are
simulated, and with \( \Z \) standard normal it is how [Chapter 3](../ch03-multivariate-normal/index.html) constructs
\( \Normal_p(\bmu,\bSigma) \). The Mahalanobis distance then appears in the exponent of the
density (@thm-mvn-density). When \( \bSigma \) is singular, no matrix whitens all of \( \Y \),
since a degenerate combination cannot be scaled to variance one. The best one can do is
whiten within the support \( \bmu+\C(\bSigma) \) using the Moore–Penrose inverse square root
(@exr-rv-singular-whitening).

## Second moments in the linear model {#sec-rv-linear-model}

Here is how the rules of this section enter regression. Suppose
\( \Y=\X\bbeta+\be \) with \( \X \) an \( n\times p \) constant matrix of rank \( p \), \( \E(\be)=\bzero \)
and \( \Cov(\be)=\sigma^2\I_n \). Adding the constant \( \X\bbeta \) changes the mean but not the
covariance, so \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\I_n \). The least squares estimator
\( \hbeta=(\X\T\X)^{-1}\X\T\Y \) is \( \A\Y \) with \( \A=(\X\T\X)^{-1}\X\T \), so
\[
\E(\hbeta)=(\X\T\X)^{-1}\X\T\X\bbeta=\bbeta,
\qquad
\Cov(\hbeta)=\sigma^2(\X\T\X)^{-1}\X\T\X(\X\T\X)^{-1}=\sigma^2(\X\T\X)^{-1}.
\]
No distributional assumption was used, and none is needed for the Gauss–Markov theorem
of Chapter 7. Fitted values and residuals are also linear in \( \Y \), and their
covariance matrices are computed the same way in Chapter 5 and
[Chapter 6](../ch06-projections/index.html). If instead \( \Cov(\be)=\sigma^2\V \) with \( \V \) positive definite, then
\( \V^{-1/2}\Y=\V^{-1/2}\X\bbeta+\V^{-1/2}\be \) is a model with \( \Cov(\V^{-1/2}\be)=\sigma^2\I \).
Whitening the data turns correlated errors back into the standard case. This is
generalized least squares, which [Chapter 6](../ch06-projections/index.html) and Chapter 31 develop.

## Exercises

### A. Check your understanding

::: {#exr-rv-differences}
[A1]

Let \( \Cov(\Y)=\sigma^2\I_3 \). Find the covariance and correlation matrices of the successive
differences \( (Y_2-Y_1,\,Y_3-Y_2)\T \), and of \( (Y_2-Y_1,\,Y_3-Y_2,\,Y_3-Y_1)\T \). Which is singular,
and which linear combination is degenerate?
:::

::: {#exr-rv-cholesky-residual}
[A2]

For \( p=2 \) and \( \bSigma=\bL\bL\T \), show that the second component of \( \bL^{-1}(\Y-\bmu) \) is a multiple of
\( Y_2-\mu_2-(\sigma_{12}/\sigma_{11})(Y_1-\mu_1) \). Explain the connection with @thm-rv-blp.
:::

### B. Practice

::: {#exr-rv-moving-average}
[B1]

Let \( \varepsilon_0,\dots,\varepsilon_n \) be uncorrelated with variance \( \sigma^2 \), and let
\( Y_t=\varepsilon_t+\theta\varepsilon_{t-1} \) for \( t=1,\dots,n \). Find \( \Cov(\Y) \), and show that it is
positive definite for every real \( \theta \), including \( \theta=\pm1 \).
:::

::: {.solution}
\( \Y=\A\be \), where \( \A \) is \( n\times(n+1) \) with \( \theta \) in
column \( t \) and \( 1 \) in column \( t+1 \) of row \( t \) (columns indexed by \( \varepsilon_0,\dots,\varepsilon_n \)).
So \( \Cov(\Y)=\sigma^2\A\A\T \) is tridiagonal, with \( \sigma^2(1+\theta^2) \) on the diagonal and \( \sigma^2\theta \)
next to it. The last \( n \) columns of \( \A \) form a lower bidiagonal matrix with unit diagonal, which
is nonsingular. So \( \A \) has full row rank \( n \), and \( \A\A\T \) is positive definite
([Section 2.2](02-linear.html)) whatever the value of \( \theta \).
:::

::: {#exr-rv-cov-rank}
[B2]

Show that \( \C(\A\bSigma\A\T)=\C(\A\bSigma) \) for every nonnegative definite \( \bSigma \) and conformable
\( \A \). Deduce that \( \rank\Cov(\A\Y)\le\rank\Cov(\Y) \).
:::

::: {#exr-rv-whitening-closest}
[B3]

Let \( \bSigma \) be positive definite. Among all \( \W \) with \( \W\bSigma\W\T=\I \), show that
\( \E\norm{\W(\Y-\bmu)-(\Y-\bmu)}^2 \) is minimized uniquely by \( \W=\bSigma^{-1/2} \). *Hint:* use
@prp-rv-whitening(b) and show that \( \tr(\Q\bSigma^{1/2})\le\tr(\bSigma^{1/2}) \) for orthogonal \( \Q \).
:::

::: {.solution}
Write \( \W=\Q\bSigma^{-1/2} \). By
@cor-rv-quadform-special(a) with \( \A=(\W-\I)\T(\W-\I) \),
\[
\E\norm{(\W-\I)(\Y-\bmu)}^2=\tr\bigl((\W-\I)\bSigma(\W-\I)\T\bigr)
=\tr(\W\bSigma\W\T)-2\tr(\W\bSigma)+\tr(\bSigma)=p+\tr(\bSigma)-2\tr(\Q\bSigma^{1/2}).
\]
Let \( \bSigma^{1/2}=\bm\Gamma\bLambda^{1/2}\bm\Gamma\T \) and \( \G=\bm\Gamma\T\Q\bm\Gamma \), which is orthogonal, so each
\( \lvert g_{ii}\rvert\le1 \). Then
\( \tr(\Q\bSigma^{1/2})=\tr(\G\bLambda^{1/2})=\sum_ig_{ii}\lambda_i^{1/2}\le\sum_i\lambda_i^{1/2} \). Since every
\( \lambda_i>0 \), equality forces every \( g_{ii}=1 \). An orthogonal matrix with unit diagonal is \( \I \),
because each column has norm one. So \( \G=\I \), which means \( \Q=\I \).
:::

::: {#exr-rv-singular-whitening}
[B4]

Let \( \bSigma=\Cov(\Y) \) have rank \( r \), with spectral decomposition
\( \bSigma=\bm\Gamma_r\bLambda_r\bm\Gamma_r\T \), where \( \bm\Gamma_r\T\bm\Gamma_r=\I_r \) and \( \bLambda_r \) holds the positive
eigenvalues. Show that \( \Z=\bLambda_r^{-1/2}\bm\Gamma_r\T(\Y-\bmu) \) has \( \Cov(\Z)=\I_r \) and that
\( \Y=\bmu+\bm\Gamma_r\bLambda_r^{1/2}\Z \) with probability one. Show that
\( (\Y-\bmu)\T\bSigma^{+}(\Y-\bmu)=\norm{\Z}^2 \), where \( \bSigma^+ \) is the Moore–Penrose inverse, and that its
expectation is \( r \).
:::

::: {.solution}
By @thm-rv-linear,
\( \Cov(\Z)=\bLambda_r^{-1/2}\bm\Gamma_r\T\bm\Gamma_r\bLambda_r\bm\Gamma_r\T\bm\Gamma_r\bLambda_r^{-1/2}=\I_r \). By
@thm-rv-cov-nnd(d), \( \Y-\bmu\in\C(\bSigma)=\C(\bm\Gamma_r) \) with probability one, and on that event
\( \Y-\bmu=\bm\Gamma_r\bm\Gamma_r\T(\Y-\bmu)=\bm\Gamma_r\bLambda_r^{1/2}\Z \). The Moore–Penrose inverse is
\( \bSigma^+=\bm\Gamma_r\bLambda_r^{-1}\bm\Gamma_r\T \), so
\( (\Y-\bmu)\T\bSigma^+(\Y-\bmu)=\norm{\bLambda_r^{-1/2}\bm\Gamma_r\T(\Y-\bmu)}^2=\norm{\Z}^2 \), with expectation
\( \tr\Cov(\Z)=r \).
:::

### C. Going deeper

::: {#exr-rv-gls-cov}
[C1]

In \( \Y=\X\bbeta+\be \) with \( \X \) of full column rank, let \( \Cov(\be)=\sigma^2\V \) with \( \V \) positive definite.
Find \( \Cov(\hbeta) \) for \( \hbeta=(\X\T\X)^{-1}\X\T\Y \) and \( \Cov(\tilde{\bbeta}) \) for
\( \tilde{\bbeta}=(\X\T\V^{-1}\X)^{-1}\X\T\V^{-1}\Y \). Show that \( \Cov(\hbeta)-\Cov(\tilde{\bbeta}) \) is nonnegative
definite. *Hint:* show that \( \hbeta-\tilde{\bbeta} \) and \( \tilde{\bbeta} \) are uncorrelated.
:::

::: {.solution}
By @thm-rv-linear,
\[
\Cov(\hbeta)=\sigma^2(\X\T\X)^{-1}\X\T\V\X(\X\T\X)^{-1},
\qquad
\Cov(\tilde{\bbeta})=\sigma^2(\X\T\V^{-1}\X)^{-1}.
\]
Write \( \hbeta-\tilde{\bbeta}=\bD\Y \) with \( \bD=(\X\T\X)^{-1}\X\T-(\X\T\V^{-1}\X)^{-1}\X\T\V^{-1} \), and note
\( \bD\X=\I-\I=\bzero \). Then
\( \Cov(\hbeta-\tilde{\bbeta},\tilde{\bbeta})=\sigma^2\bD\V\V^{-1}\X(\X\T\V^{-1}\X)^{-1}=\bzero \), so by
@prp-rv-sums(b),
\( \Cov(\hbeta)=\Cov(\tilde{\bbeta})+\Cov(\hbeta-\tilde{\bbeta}) \), and the last term is nonnegative definite.
:::
