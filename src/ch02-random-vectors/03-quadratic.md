# Expectation of quadratic forms

Sums of squares are quadratic forms in the response. The sample variance is
\( (n-1)^{-1}\Y\T(\I-n^{-1}\mathbf{J})\Y \). The residual sum of squares of a regression is
\( \Y\T(\I-\M)\Y \), where \( \M \) is the projection onto the column space of \( \X \)
([Chapter 6](../ch06-projections/index.html)). Test statistics are ratios of such forms. Their expectations
follow from a short trace argument that needs no distributional assumption.

## The basic formula

::: {#thm-rv-quadform-mean}
[Mean of a quadratic form]

Let \( \Y \) be a \( p\times1 \) random vector with finite second moments, \( \E(\Y)=\bmu \) and
\( \Cov(\Y)=\bSigma \), and let \( \A \) be a constant \( p\times p \) matrix. Then
\[
\E(\Y\T\A\Y)=\tr(\A\bSigma)+\bmu\T\A\bmu .
\]
:::

::: {.proof}
A \( 1\times1 \) matrix equals its trace, and \( \tr(\mathbf{a}\mathbf{b}\T)=\mathbf{b}\T\mathbf{a} \)
(@thm-mat-trace-cyclic), so \( \Y\T\A\Y=\tr(\A\Y\Y\T) \). By
@prp-rv-expectation-linear and @prp-rv-cov-basic(b),
\[
\E(\Y\T\A\Y)=\tr\bigl(\A\,\E(\Y\Y\T)\bigr)=\tr\bigl(\A(\bSigma+\bmu\bmu\T)\bigr)
=\tr(\A\bSigma)+\bmu\T\A\bmu .
\qedhere
\]
:::

The formula has two parts. The trace term is what the form picks up from random
variation. The second term is the value of the form at the mean, which is what it would be
without any randomness. The theorem does not require independence, normality or symmetry
of \( \A \). Since \( \Y\T\A\Y=\Y\T\A_s\Y \) with \( \A_s=\tfrac12(\A+\A\T) \), and
\( \tr(\A\bSigma)=\tr(\A_s\bSigma) \) when \( \bSigma \) is symmetric, we may assume \( \A \) symmetric
whenever that is convenient.

::: {#cor-rv-quadform-special}
Under the conditions of @thm-rv-quadform-mean:

::: {.enumerate options="label=(\alph*)"}
1. \( \E\bigl[(\Y-\bb)\T\A(\Y-\bb)\bigr]=\tr(\A\bSigma)+(\bmu-\bb)\T\A(\bmu-\bb) \) for any constant \( \bb \);

2. if \( \bSigma=\sigma^2\I \), then \( \E(\Y\T\A\Y)=\sigma^2\tr(\A)+\bmu\T\A\bmu \);

3. if \( \bSigma \) is positive definite, the squared Mahalanobis distance
           @eq-rv-mahalanobis has \( \E\,\Delta^2(\Y)=p \);

4. for random vectors \( \bU \) (\( p\times1 \)) and \( \V \) (\( q\times1 \)) with finite second
           moments and a constant \( p\times q \) matrix \( \A \),
           \( \E(\bU\T\A\V)=\tr\bigl(\A\,\Cov(\V,\bU)\bigr)+\E(\bU)\T\A\,\E(\V) \).
:::

:::

::: {.proof}
(a) applies the theorem to \( \Y-\bb \), which has covariance matrix \( \bSigma \). (b) is
immediate. (c) is (a) with \( \bb=\bmu \) and \( \A=\bSigma^{-1} \), so the trace is
\( \tr(\I_p)=p \). For (d), \( \bU\T\A\V=\tr(\A\V\bU\T) \) and
\( \E(\V\bU\T)=\Cov(\V,\bU)+\E(\V)\E(\bU)\T \).
:::

::: {#exm-rv-sample-variance}
[Sample variance with unequal means and variances]

Let \( Y_1,\dots,Y_n \) be uncorrelated with means \( \mu_i \) and variances \( \sigma_i^2 \), and let
\( S^2=(n-1)^{-1}\Y\T\mathbf{C}\Y \) with \( \mathbf{C}=\I-n^{-1}\mathbf{J} \). Here
\( \bSigma=\diag(\sigma_1^2,\dots,\sigma_n^2) \), so \( \tr(\mathbf{C}\bSigma)=(1-1/n)\sum_i\sigma_i^2 \), and
\( \bmu\T\mathbf{C}\bmu=\sum_i(\mu_i-\bar{\mu})^2 \). Therefore
\[
\E(S^2)=\overline{\sigma^2}+\frac{1}{n-1}\sum_{i=1}^n(\mu_i-\bar{\mu})^2,
\qquad \overline{\sigma^2}=\frac1n\sum_{i=1}^n\sigma_i^2 .
\]
With a common mean, \( S^2 \) is unbiased for the *average* variance, even when the
variances differ. Differences in the means inflate it. That is exactly why a residual mean
square estimates \( \sigma^2 \) only when the mean structure is right.
:::

::: {#exm-rv-rss-bias}
[Residual sum of squares under a wrong mean]

Let \( \Cov(\Y)=\sigma^2\I_n \) and \( \E(\Y)=\boldsymbol{\uptheta} \), where \( \boldsymbol{\uptheta} \) need not lie in
\( \C(\X) \). The residual sum of squares is \( \Y\T(\I-\M)\Y \), where \( \M \) is symmetric and
idempotent of rank \( r=\rank(\X) \), so \( \tr(\I-\M)=n-r \) (@thm-mat-idempotent). By
@cor-rv-quadform-special(b),
\[
\E(\text{RSS})=\sigma^2(n-r)+\norm{(\I-\M)\boldsymbol{\uptheta}}^2 .
\]{#eq-rv-rss-mean}

If the model is right, \( \boldsymbol{\uptheta}\in\C(\X) \), the second term vanishes and
\( \text{RSS}/(n-r) \) is unbiased for \( \sigma^2 \). If the model is wrong, the estimate is
biased upwards by the squared distance from the true mean to the model space, divided by
\( n-r \). For instance, fit a straight line at \( x=1,\dots,6 \) when the true mean is the parabola
\( \theta_i=\tfrac14(x_i-3.5)^2 \) and \( \sigma^2=1 \). Then \( \sigma^2(n-r)=4 \),
\( \norm{(\I-\M)\boldsymbol{\uptheta}}^2=2.333 \), and \( \E(\text{RSS})=6.333 \).
[Chapter 19](../ch19-theory-of-departures/index.html) studies such departures systematically (@thm-dep-omitted).
:::

The theorem uses only two moments, so it holds equally for skewed and heavy-tailed data.
The listing checks it for a correlated vector \( \Y=\bmu+\B\bU \), where \( \bU \) has independent,
centred exponential components, which are strongly skewed. With the matrices in the
listing, \( \tr(\A\bSigma)=6.9 \) and \( \bmu\T\A\bmu=17.0 \), so
\( \E(\Y\T\A\Y)=23.9 \). The average over \( 10^6 \) draws is
\( 23.940 \), with Monte Carlo standard error \( 0.018 \). The script also
verifies the formula exactly, by enumerating a discrete distribution with eight atoms.

```{.python .run #cell-quadform-moments-mean}
import itertools
import numpy as np
rng = np.random.default_rng(20260917)

mu = np.array([1.0, -1.0, 2.0])
B = np.array([[1.0, 0.0, 0.0],
              [0.5, 1.0, 0.0],
              [-0.5, 0.8, 0.6]])
Sigma = B @ B.T                                   # Cov(mu + B U) when Cov(U) = I
A = np.array([[2.0, 1.0, 0.0],
              [1.0, 1.0, -1.0],
              [0.0, -1.0, 3.0]])

theory = np.trace(A @ Sigma) + mu @ A @ mu        # Theorem: no normality needed

N = 1_000_000
U = rng.exponential(size=(N, 3)) - 1.0            # skewed, mean 0, variance 1
Y = mu + U @ B.T
Q = np.einsum("ij,jk,ik->i", Y, A, Y)             # Y_i' A Y_i for every draw
print(f"theory {theory:.4f}   simulation {Q.mean():.4f} ± {Q.std() / np.sqrt(N):.4f}")
```

## The variance needs fourth moments

The mean of \( \Y\T\A\Y \) involves only \( \E(Y_iY_j) \). Its second moment involves
\( \E(Y_iY_jY_kY_l) \), so no formula for \( \Var(\Y\T\A\Y) \) can be written in terms of \( \bmu \) and
\( \bSigma \) alone. Two random vectors with the same mean and covariance matrix can give the
same quadratic form very different variances. The next theorem handles independent
components, the case that matters for residual sums of squares under the classical error
assumptions. [Chapter 4](../ch04-quadratic-forms/index.html) gives the normal case with a general \( \bSigma \)
(@thm-qf-mean-var).

::: {#thm-rv-quadform-variance}
[Variance of a quadratic form]

Let \( Y_1,\dots,Y_n \) be independent with means \( \theta_1,\dots,\theta_n \) and common central
moments \( \sigma^2=\E(Y_i-\theta_i)^2 \), \( \mu_3=\E(Y_i-\theta_i)^3 \) and
\( \mu_4=\E(Y_i-\theta_i)^4<\infty \). Let \( \A \) be a constant symmetric \( n\times n \) matrix, and
let \( \mathbf{a}=(a_{11},\dots,a_{nn})\T \) be its diagonal. Then
\[
\Var(\Y\T\A\Y)=(\mu_4-3\sigma^4)\,\mathbf{a}\T\mathbf{a}+2\sigma^4\tr(\A^2)
+4\sigma^2\,\boldsymbol{\uptheta}\T\A^2\boldsymbol{\uptheta}+4\mu_3\,\boldsymbol{\uptheta}\T\A\mathbf{a} .
\]
:::

::: {.proof}
Write \( \Y=\boldsymbol{\uptheta}+\be \), where the \( \varepsilon_i \) are independent with mean zero. Then
\[
\Y\T\A\Y=q+\ell+\boldsymbol{\uptheta}\T\A\boldsymbol{\uptheta},\qquad q=\be\T\A\be,\quad \ell=2\bb\T\be,\quad \bb=\A\boldsymbol{\uptheta},
\]
so \( \Var(\Y\T\A\Y)=\Var(q)+\Var(\ell)+2\Cov(q,\ell) \). We compute the three terms.

*The linear term.* \( \Cov(\be)=\sigma^2\I \), so \( \Var(\ell)=4\sigma^2\bb\T\bb=4\sigma^2\boldsymbol{\uptheta}\T\A^2\boldsymbol{\uptheta} \).

*The cross term.* Since \( \E(\ell)=0 \), \( \Cov(q,\ell)=\E(q\ell)=2\sum_{i,j,k}a_{ij}b_k\E(\varepsilon_i\varepsilon_j\varepsilon_k) \).
By independence and zero means, \( \E(\varepsilon_i\varepsilon_j\varepsilon_k) \) vanishes unless
\( i=j=k \), when it equals \( \mu_3 \). So \( \Cov(q,\ell)=2\mu_3\sum_ia_{ii}b_i=2\mu_3\,\boldsymbol{\uptheta}\T\A\mathbf{a} \).

*The quadratic term.* \( \E(q^2)=\sum_{i,j,k,l}a_{ij}a_{kl}\E(\varepsilon_i\varepsilon_j\varepsilon_k\varepsilon_l) \).
The fourth-order moment is \( \mu_4 \) if all four indices agree. It is \( \sigma^4 \) if the indices
form two distinct pairs, which happens in three patterns: \( i=j\ne k=l \), \( i=k\ne j=l \) and
\( i=l\ne j=k \). Otherwise some index occurs exactly once, and the moment is zero. Using
the symmetry of \( \A \),
\[
\E(q^2)=\mu_4\sum_ia_{ii}^2+\sigma^4\Bigl[\sum_{i\ne k}a_{ii}a_{kk}+2\sum_{i\ne j}a_{ij}^2\Bigr]
=\mu_4\,\mathbf{a}\T\mathbf{a}+\sigma^4\bigl[(\tr\A)^2-\mathbf{a}\T\mathbf{a}+2\tr(\A^2)-2\mathbf{a}\T\mathbf{a}\bigr],
\]
since \( \sum_{i,j}a_{ij}^2=\tr(\A^2) \). By @thm-rv-quadform-mean, \( \E(q)=\sigma^2\tr\A \), and
subtracting its square gives \( \Var(q)=(\mu_4-3\sigma^4)\mathbf{a}\T\mathbf{a}+2\sigma^4\tr(\A^2) \).
Adding the three terms proves the formula.
:::

For normal errors, \( \mu_3=0 \) and \( \mu_4=3\sigma^4 \), and the formula collapses to
\( 2\sigma^4\tr(\A^2)+4\sigma^2\boldsymbol{\uptheta}\T\A^2\boldsymbol{\uptheta} \). Both extra terms describe how
nonnormality leaks into second-order behaviour. The kurtosis term depends on \( \A \) only
through its diagonal. The skewness term needs both a skewed error distribution and a
mean that the form does not annihilate. For the residual sum of squares of a correct
model, \( \A\boldsymbol{\uptheta}=(\I-\M)\boldsymbol{\uptheta}=\bzero \), so the skewness term disappears but the kurtosis
term remains.

::: {#exm-rv-rss-variance}
[Skewed errors and a misspecified line]

Continue @exm-rv-rss-bias, now with independent centred exponential errors, for which
\( \sigma^2=1 \), \( \mu_3=2 \) and \( \mu_4=9 \). @thm-rv-quadform-variance gives
\( \Var(\text{RSS})=29.80 \). Of this, the kurtosis term contributes
\( 16.73 \) and the skewness term \( -4.27 \). Normal errors
with the same mean and variance would give \( 17.33 \), so the exponential
errors multiply the variance by \( 1.72 \) even though the mean of RSS is
unchanged. Over \( 10^6 \) simulated data sets, the variance of RSS is
\( 30.04 \). The script also confirms the formula exactly for a two-point
error distribution.
:::

```{.python .run #cell-quadform-moments-variance}
se = Q.std() / np.sqrt(N)
vals, probs = np.array([-1.0, 2.0]), np.array([2 / 3, 1 / 3])
exact = 0.0
for idx in itertools.product(range(2), repeat=3):
    u = vals[list(idx)]
    y = mu + B @ u
    exact += np.prod(probs[list(idx)]) * (y @ A @ y)
x = np.arange(1.0, 7.0)                           # n = 6 design points
X = np.column_stack([np.ones_like(x), x])
M = X @ np.linalg.solve(X.T @ X, X.T)
A2 = np.eye(6) - M                                # RSS = Y'(I - M)Y for a straight line
theta = 0.25 * (x - 3.5) ** 2                     # true mean is curved: the line is wrong
a = np.diag(A2)

def var_quadform(A, theta, m2, m3, m4):
    """Var(Y'AY) for independent Y_i with means theta_i and common central moments."""
    a = np.diag(A)
    return ((m4 - 3 * m2**2) * a @ a + 2 * m2**2 * np.trace(A @ A)
            + 4 * m2 * theta @ A @ A @ theta + 4 * m3 * theta @ A @ a)

# A2 = I - M for a straight line at x = 1..6; theta = true (curved) mean
mean_rss = np.trace(A2) + theta @ A2 @ theta      # sigma^2 = 1
var_expo = var_quadform(A2, theta, 1.0, 2.0, 9.0) # centred exponential errors
var_norm = var_quadform(A2, theta, 1.0, 0.0, 3.0) # what normal errors would give

E = rng.exponential(size=(N, 6)) - 1.0
R = np.einsum("ij,jk,ik->i", theta + E, A2, theta + E)
print(f"mean  {mean_rss:.3f} vs {R.mean():.3f}")
print(f"var   {var_expo:.3f} vs {R.var():.3f}   (normal-theory value {var_norm:.3f})")
```

## Exercises

### A. Check your understanding

::: {#exr-rv-distance-two-vectors}
[A1]

For random vectors \( \bU \) and \( \V \) of the same size, show that
\( \E\norm{\bU-\V}^2=\tr\Cov(\bU)+\tr\Cov(\V)-2\tr\Cov(\bU,\V)+\norm{\E\bU-\E\V}^2 \).
:::

### B. Practice

::: {#exr-rv-lag-quadform}
[B1]

Let \( \E(Y_t)=\mu \) and \( \Cov(Y_t,Y_s)=\sigma^2\phi^{\lvert t-s\rvert} \) for \( t,s=1,\dots,n \), with
\( \lvert\phi\rvert<1 \). Show that the sample variance satisfies
\[
\E(S^2)=\sigma^2\Bigl[1-\frac{2}{n(n-1)}\sum_{k=1}^{n-1}(n-k)\phi^k\Bigr].
\]
Is \( S^2 \) biased up or down when \( \phi>0 \)? What happens as \( n\to\infty \) with \( \phi \) fixed?
:::

::: {.solution}
With \( \mathbf{C}=\I-n^{-1}\bone\bone\T \) and \( \E(\Y)=\mu\bone \),
\( \mathbf{C}\E(\Y)=\bzero \), so \( \E(\Y\T\mathbf{C}\Y)=\tr(\mathbf{C}\bSigma)=\tr(\bSigma)-n^{-1}\bone\T\bSigma\bone \). Here
\( \tr(\bSigma)=n\sigma^2 \), and \( \bone\T\bSigma\bone \) sums all entries. There are \( n \) diagonal entries and
\( 2(n-k) \) entries at lag \( k \), so
\( \bone\T\bSigma\bone=\sigma^2\bigl[n+2\sum_{k=1}^{n-1}(n-k)\phi^k\bigr] \). Dividing
\( n\sigma^2-n^{-1}\bone\T\bSigma\bone \) by \( n-1 \) gives the formula. For \( \phi>0 \) the sum is positive, so \( S^2 \)
is biased downwards: positively correlated observations cluster, and the sample looks less
variable than the process is. As \( n\to\infty \) the sum is \( O(n) \), so the bias is \( O(1/n) \) and vanishes.
:::

::: {#exr-rv-rss-correlated}
[B2]

Suppose \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\V \). Show that \( \E(\text{RSS})=\sigma^2\tr\bigl((\I-\M)\V\bigr) \).
For the equicorrelated \( \V=(1-\rho)\I+\rho\bone\bone\T \) and a model with \( \bone\in\C(\X) \), show that
\( \E(\text{RSS})/(n-r)=\sigma^2(1-\rho) \).
:::

::: {#exr-rv-var-sample-variance}
[B3]

Let \( Y_1,\dots,Y_n \) be independent and identically distributed with variance \( \sigma^2 \) and finite
fourth central moment \( \mu_4 \). Use @thm-rv-quadform-variance to show that
\( \Var(S^2)=(\mu_4-3\sigma^4)/n+2\sigma^4/(n-1) \). Compare with the normal case.
:::

::: {.solution}
Take \( \A=\mathbf{C}/(n-1) \) and \( \boldsymbol{\uptheta}=\mu\bone \). Then
\( \A\boldsymbol{\uptheta}=\bzero \), so the last two terms of @thm-rv-quadform-variance vanish. Every diagonal
entry of \( \A \) is \( (1-1/n)/(n-1)=1/n \), so \( \mathbf{a}\T\mathbf{a}=1/n \), and
\( \tr(\A^2)=\tr(\mathbf{C})/(n-1)^2=1/(n-1) \). This gives the formula. Normal data have
\( \mu_4=3\sigma^4 \) and \( \Var(S^2)=2\sigma^4/(n-1) \). Heavy tails (\( \mu_4>3\sigma^4 \)) add a term of the same
order \( 1/n \), so the normal-theory variance is wrong even in large samples.
:::

::: {#exr-rv-cov-quadforms}
[B4]

Under the conditions of @thm-rv-quadform-variance with \( \boldsymbol{\uptheta}=\bzero \), show that
\( \Cov(\Y\T\A\Y,\Y\T\B\Y)=(\mu_4-3\sigma^4)\mathbf{a}\T\mathbf{b}+2\sigma^4\tr(\A\B) \), where \( \mathbf{a} \) and \( \mathbf{b} \) are the
diagonals of the symmetric matrices \( \A \) and \( \B \). *Hint:* apply the theorem to \( \A+\B \).
:::

::: {#exr-rv-best-quadratic}
[B5]

Let \( Y_1,\dots,Y_n \) be independent \( \Normal(0,\sigma^2) \). Among symmetric \( \A \) with \( \E(\Y\T\A\Y)=\sigma^2 \),
show that \( \Var(\Y\T\A\Y) \) is minimized by \( \A=n^{-1}\I \). *Hint:* \( \tr(\A^2)\ge(\tr\A)^2/n \).
:::

### C. Going deeper

::: {#exr-rv-fourth-moment-bound}
[C1]

::: {.enumerate options="label=(\alph*)"}
1. Exhibit two random vectors in \( \Real^n \), both with mean \( \bzero \) and covariance \( \I_n \), for
           which \( \Var(\Y\T\Y) \) equals \( 0 \) and \( 2n \). So the variance of a quadratic form, unlike
           its mean in @thm-rv-quadform-mean, is not determined by \( \bmu \) and \( \bSigma \). Check
           both answers against @thm-rv-quadform-variance.

2. Let \( Y_1,\dots,Y_n \) be independent and identically distributed with variance \( \sigma^2>0 \)
           and fourth central moment \( \mu_4 \). Using @exr-rv-var-sample-variance, show that
           \( \Var(S^2)\ge0 \) for every \( n\ge2 \) forces \( \mu_4\ge\sigma^4 \), and prove that inequality
           directly.

3. Show that \( \mu_4=\sigma^4 \) holds only for the two-point law \( \Pr(Y_i=\theta\pm\sigma)=\tfrac12 \),
           that this law gives \( \Var(S^2)=2\sigma^4/\{n(n-1)\} \), and that this is the smallest
           variance the sample variance can have for a given \( \sigma^2 \).
:::

:::

::: {.solution}
(a) Let the \( Y_i \) be independent random signs. Then \( \Y\T\Y=n \) is
constant, so its variance is \( 0 \); here \( \sigma^2=1 \), \( \mu_4=1 \), \( \A=\I \) and
\( \mathbf{a}\T\mathbf{a}=\tr(\A^2)=n \), so @thm-rv-quadform-variance gives \( (1-3)n+2n=0 \). Let the
\( Y_i \) be independent \( \Normal(0,1) \) instead. Then \( \mu_4=3 \) and the formula gives
\( 0\cdot n+2n=2n \), which is the variance of a \( \chi^2(n) \) variable. Both vectors have mean
\( \bzero \) and covariance \( \I_n \).

(b) Writing the formula of @exr-rv-var-sample-variance over a common denominator,
\[
\Var(S^2)=\frac{\mu_4-3\sigma^4}{n}+\frac{2\sigma^4}{n-1}
=\frac{n(\mu_4-\sigma^4)-\mu_4+3\sigma^4}{n(n-1)} .
\]
If \( \mu_4<\sigma^4 \) the numerator is negative for all large \( n \), which is impossible. Directly,
\( \mu_4-\sigma^4=\Var\{(Y_i-\theta)^2\}\ge0 \) with \( \theta=\E Y_i \).

(c) Equality means \( (Y_i-\theta)^2 \) is constant, necessarily equal to its mean \( \sigma^2 \), so
\( Y_i-\theta=\pm\sigma \); the two signs have probability \( \tfrac12 \) each because \( \E(Y_i-\theta)=0 \).
Substituting \( \mu_4=\sigma^4 \) in the display gives \( 2\sigma^4/\{n(n-1)\} \). The display is increasing
in \( \mu_4 \) for \( n\ge2 \), so no law with the same \( \sigma^2 \) does better.
:::

::: {#exr-rv-residual-bias-range}
[C2]

Let \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\V \), where \( \V \) is nonnegative definite with unit
diagonal, let \( \M \) be the orthogonal projection onto \( \C(\X) \) with \( \rank(\X)=r<n \), and let
\( s^2=\text{RSS}/(n-r) \). By @exr-rv-rss-correlated, \( \E(s^2)=\sigma^2\tr\{(\I-\M)\V\}/(n-r) \).

::: {.enumerate options="label=(\alph*)"}
1. Show that \( 0\le\tr\{(\I-\M)\V\}\le n \), so that \( 0\le\E(s^2)\le\sigma^2n/(n-r) \).

2. Show that both bounds are attained. For the lower one take \( \bone\in\C(\X) \) and
           \( \V=\bone\bone\T \); for the upper one take \( \X=\bone \) with \( n \) even and \( \V=\bv\bv\T \),
           where \( \bv=(1,-1,1,-1,\dots)\T \). Describe the data in each case.

3. Show that \( \E(s^2)=\sigma^2 \) for *every* model matrix \( \X \) if and only if \( \V=\I \).
:::

:::

::: {.solution}
(a) Put \( \mathbf{N}=\I-\M \). Then \( \tr(\mathbf{N}\V)=\tr(\mathbf{N}\V\mathbf{N})\ge0 \),
since \( \mathbf{N}\V\mathbf{N} \) is nonnegative definite. For the upper bound write the spectral
decomposition \( \V=\sum_k\lambda_k\bu_k\bu_k\T \) with \( \lambda_k\ge0 \) (@thm-mat-spectral). Then
\( \tr(\mathbf{N}\V)=\sum_k\lambda_k\norm{\mathbf{N}\bu_k}^2\le\sum_k\lambda_k=\tr(\V)=n \), because
\( \norm{\mathbf{N}\bu_k}\le\norm{\bu_k}=1 \).

(b) With \( \V=\bone\bone\T \) and \( \bone\in\C(\X) \) we have \( \mathbf{N}\bone=\bzero \), so
\( \tr(\mathbf{N}\V)=\bone\T\mathbf{N}\bone=0 \) and \( \E(s^2)=0 \). Here \( \Y=\X\bbeta+\sigma Z\bone \) for a
single scalar \( Z \): every observation is displaced by the same amount, the fit absorbs it and
the residuals vanish. With \( \X=\bone \) and \( \V=\bv\bv\T \), \( \bone\T\bv=0 \) for even \( n \), so
\( \mathbf{N}\bv=\bv \) and \( \tr(\mathbf{N}\V)=\norm{\bv}^2=n \), giving \( \E(s^2)=\sigma^2n/(n-1) \). Here
\( Y_i=\mu+\sigma v_iZ \): the errors alternate in sign, none of the variation is along \( \bone \), and
the residual sum of squares carries all \( n \) units of variance while being divided by \( n-1 \).
Both \( \V \) are singular; replacing them by \( (1-\varepsilon)\V+\varepsilon\I \) keeps them correlation
matrices, makes them positive definite, and changes \( \E(s^2) \) by \( O(\varepsilon) \).

(c) If \( \V=\I \) then \( \tr(\mathbf{N}\V)=n-r \) for every \( \X \). Conversely, take \( \X=\mathbf{u} \), a
single column with \( \norm{\mathbf{u}}=1 \), so that \( r=1 \) and \( \M=\mathbf{u}\mathbf{u}\T \). Unbiasedness
requires \( \tr(\V)-\mathbf{u}\T\V\mathbf{u}=n-1 \), that is \( \mathbf{u}\T\V\mathbf{u}=1 \) for every unit \( \mathbf{u} \).
As in @exr-rv-isotropic(a), a symmetric matrix with \( \mathbf{u}\T(\V-\I)\mathbf{u}=0 \) for all \( \mathbf{u} \) is
\( \V=\I \).
:::
