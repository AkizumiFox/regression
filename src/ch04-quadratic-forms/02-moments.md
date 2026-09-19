# Mean and variance of quadratic forms

## The mean needs no normality

@thm-rv-quadform-mean computed \( \E(\Y\T\A\Y) \) for any random vector
with mean \( \bmu \) and covariance \( \bSigma \):
\[
\E(\Y\T\A\Y)=\tr(\A\bSigma)+\bmu\T\A\bmu .
\]{#eq-qf-mean}

The formula has two parts with distinct roles. The
trace term is what the noise contributes. The term \( \bmu\T\A\bmu \) is what the
signal contributes. For nonnegative definite \( \A \) it vanishes exactly when
\( \A\bmu=\bzero \), which describes the mean vectors a well-designed estimator is meant to
ignore.

::: {#exm-qf-drift}
[Two estimators of \( \sigma^2 \) under a drifting mean]

Let \( Y_1,\dots,Y_n \) be uncorrelated with common variance \( \sigma^2 \), and consider
\[
S^2=\frac{1}{n-1}\sum_{i=1}^n(Y_i-\bar{Y})^2,\qquad
Q=\frac{1}{2(n-1)}\sum_{i=1}^{n-1}(Y_{i+1}-Y_i)^2 ,
\]
the sample variance and the **mean square successive difference**. In matrix
form \( S^2=\Y\T\mathbf{C}\Y/(n-1) \) with \( \mathbf{C}=\I-n^{-1}\bone\bone\T \), and
\( Q=\Y\T\bD\T\bD\Y/\{2(n-1)\} \), where \( \bD \) is the \( (n-1)\times n \) first-difference
matrix with rows \( (\dots,-1,1,\dots) \). Both matrices annihilate \( \bone \), so both
estimators ignore a constant mean. Since \( \tr(\mathbf{C})=n-1 \) and
\( \tr(\bD\T\bD)=\tr(\bD\bD\T)=2(n-1) \), @eq-qf-mean with \( \bSigma=\sigma^2\I \) shows
that both are unbiased when the mean is constant.

Now let the mean drift linearly, \( \mu_i=\mu_0+\theta i \). Then \( \mathbf{C}\bmu \) is the
centred trend and \( \bD\bmu=\theta\bone_{n-1} \), so
\[
\E(S^2)=\sigma^2+\theta^2\frac{n(n+1)}{12},\qquad
\E(Q)=\sigma^2+\frac{\theta^2}{2}.
\]
The bias of \( S^2 \) grows quadratically with the length of the series, while that of
\( Q \) does not grow at all. With \( n=20 \), \( \sigma=1 \) and
\( \theta=0.1 \), the means are \( 1.350 \) and
\( 1.005 \). Differencing removes smooth trends, and this is why \( Q \) was
proposed for series whose level may wander (von Neumann 1941).
The price of robustness to drift is paid in variance, which we compute below.
:::

## A canonical representation

Higher moments and distributions need normality, and one representation
answers all such questions at once.

::: {#thm-qf-canonical}
[Canonical form of a normal quadratic form]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) and let \( \A \) be symmetric. Write
\( \bSigma=\bL\bL\T \) with \( \bL \) of size \( n\times k \), and take a spectral decomposition
\( \bL\T\A\bL=\Q\bLambda\Q\T \) with \( \bLambda=\diag(\lambda_1,\dots,\lambda_k) \). Put
\( \bb=\Q\T\bL\T\A\bmu \). Then
\[
\Y\T\A\Y\;\overset{d}{=}\;\bmu\T\A\bmu+\sum_{j=1}^k\bigl(\lambda_jW_j^2+2b_jW_j\bigr),
\qquad \W\sim\Normal_k(\bzero,\I_k).
\]{#eq-qf-canonical}

For \( \lambda_j\ne0 \) the \( j \)th term equals \( \lambda_j(W_j+b_j/\lambda_j)^2-b_j^2/\lambda_j \).
The nonzero \( \lambda_j \) are the nonzero eigenvalues of \( \A\bSigma \).
:::

::: {.proof}
Such an \( \bL \) exists. By @thm-mat-spectral, \( \bSigma=\bU\bD\bU\T \) with \( \bU \)
orthogonal and \( \bD \) diagonal with entries \( d_i\ge0 \), and \( \bL=\bU\bD^{1/2} \) works with
\( k=n \). Let \( \Z\sim\Normal_k(\bzero,\I) \). The vector \( \bmu+\bL\Z \) is normal with mean
\( \bmu \) and covariance \( \bL\bL\T=\bSigma \) (@thm-mvn-linear). A normal distribution
is determined by its mean and covariance (@thm-mvn-mgf), so \( \Y \) and \( \bmu+\bL\Z \)
have the same law, and so do \( \Y\T\A\Y \) and
\[
(\bmu+\bL\Z)\T\A(\bmu+\bL\Z)=\bmu\T\A\bmu+2\bmu\T\A\bL\Z+\Z\T\bL\T\A\bL\Z .
\]
Put \( \W=\Q\T\Z \), again \( \Normal_k(\bzero,\I) \). Then
\( \Z\T\bL\T\A\bL\Z=\W\T\bLambda\W=\sum\lambda_jW_j^2 \) and
\( \bmu\T\A\bL\Z=\bb\T\W \), which is @eq-qf-canonical. The completed square is
algebra. Finally, \( \bL\T\A\bL=\bH\mathbf{F} \) with \( \mathbf{F}=\A\bL \) and \( \bH=\bL\T \), and
\( \mathbf{F}\bH=\A\bSigma \). By @prp-mat-eigen-basic(d) these two products have the same nonzero
eigenvalues, with the same multiplicities.
:::

So a normal quadratic form is a constant plus a sum of independent pieces. Each piece
is either a scaled noncentral \( \chi^2(1) \) variable, with scale \( \lambda_j \), or a
normal variable with mean zero, when \( \lambda_j=0 \) but \( b_j\ne0 \). Because the
\( \lambda_j \) may be negative, \( \Y\T\A\Y \) need not be positive. The linear pieces
are special to singular covariances. If \( \bSigma \) is positive definite, \( \bL \) can be
taken square and nonsingular, and then \( \lambda_j=0 \) gives \( \A\bL\mathbf{q}_j=\bzero \) for
the \( j \)th column \( \mathbf{q}_j \) of \( \Q \), so \( b_j=(\A\bL\mathbf{q}_j)\T\bmu=0 \). A linear piece
needs a direction in which \( \A \) responds to the mean while the noise is absent.

::: {#cor-qf-mgf}
[Moment generating function]

In the setting of @thm-qf-canonical, for \( \lvert t\rvert<1/(2\max_j\lvert\lambda_j\rvert) \)
(any \( t \) if all \( \lambda_j=0 \)),
\[
\E\,e^{t\Y\T\A\Y}=e^{t\bmu\T\A\bmu}\prod_{j=1}^k(1-2\lambda_jt)^{-1/2}
\exp\Bigl(\frac{2b_j^2t^2}{1-2\lambda_jt}\Bigr).
\]
:::

::: {.proof}
The terms of @eq-qf-canonical are independent. If \( \lambda_j=0 \), the term
\( 2b_jW_j \) has mgf \( e^{2b_j^2t^2} \). If \( \lambda_j\ne0 \), apply
@lem-qf-shifted-square to \( W_j+\nu_j \) with \( \nu_j=b_j/\lambda_j \) at the
argument \( \lambda_jt \), and multiply by \( e^{-tb_j^2/\lambda_j} \). The exponent is
\( b_j^2t/\{\lambda_j(1-2\lambda_jt)\}-b_j^2t/\lambda_j=2b_j^2t^2/(1-2\lambda_jt) \).
:::

## Variance and covariances

::: {#thm-qf-mean-var}
[Moments of normal quadratic forms]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \), let \( \A \) and \( \B \) be symmetric \( n\times n \)
matrices, and let \( \mathbf{K} \) be an \( m\times n \) matrix. Then

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\Y\T\A\Y)=\tr(\A\bSigma)+\bmu\T\A\bmu \);

2. \( \Var(\Y\T\A\Y)=2\tr\{(\A\bSigma)^2\}+4\bmu\T\A\bSigma\A\bmu \);

3. \( \Cov(\mathbf{K}\Y,\Y\T\A\Y)=2\mathbf{K}\bSigma\A\bmu \);

4. \( \Cov(\Y\T\A\Y,\Y\T\B\Y)=2\tr(\A\bSigma\B\bSigma)+4\bmu\T\A\bSigma\B\bmu \).
:::

:::

::: {.proof}
(a) is @eq-qf-mean. For (b), use @eq-qf-canonical. Independent terms
have additive variances, and for \( W\sim\Normal(0,1) \),
\( \Var(\lambda W^2+2bW)=2\lambda^2+4b^2 \), because \( \Cov(W^2,W)=\E W^3=0 \). Hence
\[
\Var(\Y\T\A\Y)=2\sum_j\lambda_j^2+4\norm{\bb}^2
=2\tr\{(\bL\T\A\bL)^2\}+4\bmu\T\A\bL\Q\Q\T\bL\T\A\bmu .
\]
Now \( \tr\{(\bL\T\A\bL)^2\}=\tr(\bL\T\A\bSigma\A\bL)=\tr(\A\bSigma\A\bSigma) \) by
@thm-mat-trace-cyclic, and \( \Q\Q\T=\I \), so the second term is
\( 4\bmu\T\A\bSigma\A\bmu \).
For (c), write \( \Y=\bmu+\bL\Z \) as in the proof of @thm-qf-canonical. Then
\( \Y\T\A\Y-\E(\Y\T\A\Y)=\Z\T\bL\T\A\bL\Z-\tr(\A\bSigma)+2\bmu\T\A\bL\Z \) and
\( \mathbf{K}\Y-\E(\mathbf{K}\Y)=\mathbf{K}\bL\Z \). Every entry of \( \E[\Z(\Z\T\mathbf{G}\Z)] \) is a sum of
third moments of a \( \Normal_k(\bzero,\I) \) vector, which vanish. So only the linear
term contributes, and \( \E[\mathbf{K}\bL\Z\cdot2\Z\T\bL\T\A\bmu]=2\mathbf{K}\bL\bL\T\A\bmu \).
For (d), apply (b) to \( \A \), \( \B \) and \( \A+\B \), and use
\( \Var(X+Y)=\Var X+\Var Y+2\Cov(X,Y) \).
:::

For \( \bmu=\bzero \), part (b) can also be read off the fourth moments of a normal vector
(@exr-mvn-fourth-moments). The variance has the same two-part structure as the mean: a pure-noise term, plus a
term that is present only when \( \A \) responds to the mean. Part (c) shows that a
quadratic form and a linear form are uncorrelated when \( \mathbf{K}\bSigma\A=\mathbf{0} \). Under
normality this will turn out to give independence (@thm-qf-indep-linear).

::: {#exm-qf-drift-variance}
[@exm-qf-drift, continued]

With \( \bSigma=\sigma^2\I \) and a constant mean, @thm-qf-mean-var(b) gives
\( \Var(S^2)=2\sigma^4\tr(\mathbf{C}^2)/(n-1)^2=2\sigma^4/(n-1) \). For \( Q \), the matrix
\( \bD\bD\T \) is tridiagonal with \( 2 \) on the diagonal and \( -1 \) beside it, so
\( \tr\{(\bD\T\bD)^2\}=\tr\{(\bD\bD\T)^2\}=4(n-1)+2(n-2)=6n-8 \) and
\[
\Var(Q)=\frac{2\sigma^4(6n-8)}{4(n-1)^2}=\frac{\sigma^4(3n-4)}{(n-1)^2}.
\]
The ratio \( \Var(Q)/\Var(S^2)=(3n-4)/\{2(n-1)\} \) tends to \( 3/2 \). With
\( n=20 \) the two variances are \( 0.1551 \) and
\( 0.1053 \), a ratio of \( 1.474 \). So with no drift,
\( Q \) needs about half as many observations again to match \( S^2 \). Under the drift
\( \theta=0.1 \), the signal term \( 4\bmu\T\A\bSigma\A\bmu \) raises the
variances to \( 0.1552 \) and \( 0.1789 \). Combining variance and bias,
the root mean squared errors are \( 0.394 \) for \( Q \) and
\( 0.549 \) for \( S^2 \). The listing evaluates the formulas, and the script
confirms all eight moments with \( 400{,}000 \) simulated series.
:::

```{.python .run #cell-successive-differences-moments}
import numpy as np
rng = np.random.default_rng(4402)
n, sigma, theta = 20, 1.0, 0.1
t = np.arange(1, n + 1)
mu = 5.0 + theta * t                                   # mean with a linear trend
Sigma = sigma**2 * np.eye(n)

def qf_moments(A, mu, Sigma):
    """Mean and variance of Y^T A Y for Y ~ N(mu, Sigma), A symmetric."""
    AS = A @ Sigma
    mean = np.trace(AS) + mu @ A @ mu
    var = 2 * np.trace(AS @ AS) + 4 * mu @ AS @ A @ mu
    return mean, var

C = np.eye(n) - np.ones((n, n)) / n                    # centring matrix
D = np.diff(np.eye(n), axis=0)                         # (n-1) x n first differences
A_s2 = C / (n - 1)
A_q = D.T @ D / (2 * (n - 1))
for name, A in [("S2", A_s2), ("Q", A_q)]:
    print(name, qf_moments(A, mu, Sigma), qf_moments(A, 0 * mu, Sigma))
```

Normality enters (b)–(d) through the third and fourth moments of \( \Y \). For
independent coordinates with excess kurtosis the variance acquires an extra term
that involves only the diagonal of \( \A \) (@thm-rv-quadform-variance). Variance formulas for
sums of squares derived under normality can be badly wrong for heavy-tailed errors,
a point taken up in [Chapter 19](../ch19-theory-of-departures/index.html) (@prp-dep-sigma-interval, @prp-dep-balance).

## Exercises

### A. Check your understanding

::: {#exr-qf-symmetrize}
[A1]

Let \( \Y\sim\Normal_2(\bzero,\I) \) and \( \A=\begin{psmallmatrix}0&1\\0&0\end{psmallmatrix} \). Compute
\( \Var(\Y\T\A\Y) \) directly, and compare with \( 2\tr(\A^2) \) and with \( 2\tr(\A_s^2) \), where
\( \A_s=\tfrac12(\A+\A\T) \).
:::

### B. Practice

::: {#exr-qf-sample-cov-var}
[B1]

Let \( (U_i,W_i) \), \( i=1,\dots,n \), be independent copies of a bivariate normal pair with
variances \( \sigma_U^2 \), \( \sigma_W^2 \) and covariance \( \sigma_{UW} \), and let
\( s_{UW}=(n-1)^{-1}\sum_i(U_i-\bar{U})(W_i-\bar{W}) \). Write \( s_{UW}=\Y\T\A\Y \) for the stacked vector
\( \Y=(\mathbf{U}\T,\mathbf{W}\T)\T \) and a symmetric \( 2n\times2n \) matrix \( \A \), and use @thm-qf-mean-var
to show that \( \Var(s_{UW})=(\sigma_U^2\sigma_W^2+\sigma_{UW}^2)/(n-1) \). Check the answer against
\( \Var(S^2) \) when \( \mathbf{U}=\mathbf{W} \).
:::

::: {#exr-qf-kurtosis}
[B2]

Let \( Y_1,\dots,Y_n \) be independent with mean \( 0 \), variance \( \sigma^2 \) and \( \E Y_i^4=\mu_4 \). By
@exr-rv-cov-quadforms, for symmetric \( \A \) and \( \B \),
\[
\Cov(\Y\T\A\Y,\Y\T\B\Y)=(\mu_4-3\sigma^4)\sum_ia_{ii}b_{ii}+2\sigma^4\tr(\A\B).
\]
Show that in a one-way layout with \( g\ge2 \) groups, not all of size one, the between-group and
within-group sums of squares are correlated unless \( \mu_4=3\sigma^4 \). Which sign does the
correlation have for heavy-tailed errors?
:::

::: {.solution}
In the one-way layout the between-group matrix is
\( \A=\M-\bP_1 \) and the within-group matrix is \( \B=\I-\M \). They satisfy \( \A\B=\mathbf{0} \), so
\( \tr(\A\B)=0 \) and the covariance reduces to \( (\mu_4-3\sigma^4)\sum_ia_{ii}b_{ii} \). An observation in
group \( \ell \) has \( a_{ii}=1/n_\ell-1/n \) and \( b_{ii}=1-1/n_\ell \), so
\[
\sum_ia_{ii}b_{ii}=\sum_\ell n_\ell\Bigl(\frac1{n_\ell}-\frac1n\Bigr)\Bigl(1-\frac1{n_\ell}\Bigr)
=\sum_\ell\Bigl(1-\frac{n_\ell}{n}\Bigr)\Bigl(1-\frac1{n_\ell}\Bigr).
\]
With \( g\ge2 \) every first factor is positive, and the second is positive for any group with two or
more observations. So the sum is positive, and the covariance vanishes iff \( \mu_4=3\sigma^4 \). For
heavy-tailed errors, \( \mu_4>3\sigma^4 \), the two sums of squares are positively correlated: a
large error inflates both.
:::

::: {#exr-qf-sse-bias}
[B3]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \), and let \( \M \) be symmetric idempotent of rank \( r \) with
\( \M\bmu\neq\bmu \), as when a regression model omits a needed regressor. Show that
\( \norm{(\I-\M)\Y}^2/\sigma^2\sim\chi^2(n-r,\norm{(\I-\M)\bmu}^2/\sigma^2) \), and find the bias of
\( \norm{(\I-\M)\Y}^2/(n-r) \) as an estimator of \( \sigma^2 \). In the setting of @thm-qf-nested-f
without the assumption \( \bP\bmu=\bmu \), which parts of the conclusion survive?
:::

::: {#exr-qf-mgf-det}
[B4]

Let \( \Y\sim\Normal_n(\bzero,\bSigma) \). Use @cor-qf-mgf and @prp-mat-eigen-basic(d) to show that
\( \E\,e^{t\Y\T\A\Y}=\det(\I_n-2t\A\bSigma)^{-1/2} \) for \( \lvert t\rvert \) small.
:::

### C. Going deeper

::: {#exr-qf-best-quadratic}
[C1]

Let \( \Y\sim\Normal_n(\mu\bone,\sigma^2\I) \) with \( \mu \) and \( \sigma^2 \) unknown. Among symmetric \( \A \) with
\( \E(\Y\T\A\Y)=\sigma^2 \) for all \( \mu \) and \( \sigma^2 \), show that \( \A=\mathbf{C}/(n-1) \), which gives \( S^2 \),
uniquely minimizes \( \Var(\Y\T\A\Y) \) at every parameter value. (Compare
@exr-rv-best-quadratic, where the mean is known to be zero.)
:::

::: {.solution}
Since \( \E(\Y\T\A\Y)=\sigma^2\tr\A+\mu^2\bone\T\A\bone \),
unbiasedness means \( \tr\A=1 \) and \( \bone\T\A\bone=0 \). By @thm-qf-mean-var,
\( \Var(\Y\T\A\Y)=2\sigma^4\tr(\A^2)+4\mu^2\sigma^2\norm{\A\bone}^2 \). The map \( \A\mapsto\mathbf{C}\A\mathbf{C} \) is
an orthogonal projection for the inner product \( \tr(\A\B) \) on symmetric matrices, so
\( \tr(\A^2)\ge\tr\{(\mathbf{C}\A\mathbf{C})^2\} \), with equality iff \( \A=\mathbf{C}\A\mathbf{C} \). Also
\( \tr(\mathbf{C}\A\mathbf{C})=\tr\A-\bone\T\A\bone/n=1 \). The symmetric matrix \( \mathbf{C}\A\mathbf{C} \) has rank at most \( n-1 \).
If \( \lambda_1,\dots,\lambda_{n-1} \) are its eigenvalues on \( \bone\perpc \), the Cauchy–Schwarz inequality
gives \( 1=(\sum\lambda_i)^2\le(n-1)\sum\lambda_i^2 \). Hence \( \tr(\A^2)\ge1/(n-1) \), with equality iff
\( \A=\mathbf{C}\A\mathbf{C} \) and all \( \lambda_i=1/(n-1) \), that is, iff \( \A=\mathbf{C}/(n-1) \). That matrix also has
\( \A\bone=\bzero \), so it minimizes both variance terms at once.
:::
