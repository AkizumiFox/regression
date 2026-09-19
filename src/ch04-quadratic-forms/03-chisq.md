# When a quadratic form is chi-squared

A quadratic form is a weighted sum of independent noncentral \( \chi^2(1) \) pieces, and
the weights are eigenvalues. It is itself chi-squared when every nonzero weight
equals one. This section turns that observation into checkable conditions on \( \A \)
and \( \bSigma \), and shows that for the spherical case the condition is necessary as
well as sufficient.

## Spherical normal vectors

::: {#thm-qf-chisq-identity}
[Chi-squared forms in a spherical normal vector]

Let \( \Y\sim\Normal_n(\bmu,\I_n) \), let \( \A \) be symmetric, and let \( r\ge1 \) be an integer
and \( \gamma\ge0 \). Then \( \Y\T\A\Y\sim\chi^2(r,\gamma) \) if and only if \( \A \) is idempotent
with \( \rank(\A)=r \) and \( \gamma=\bmu\T\A\bmu \).
:::

::: {.proof}
*Sufficiency.* By @thm-mat-idempotent, \( \A=\Q_1\Q_1\T \), where the \( n\times r \)
matrix \( \Q_1 \) has orthonormal columns and \( r=\rank(\A)=\tr(\A) \).
Then \( \Y\T\A\Y=\norm{\Q_1\T\Y}^2 \) and \( \Q_1\T\Y\sim\Normal_r(\Q_1\T\bmu,\Q_1\T\Q_1)=\Normal_r(\Q_1\T\bmu,\I_r) \),
so by definition \( \Y\T\A\Y\sim\chi^2(r,\norm{\Q_1\T\bmu}^2) \), and
\( \norm{\Q_1\T\bmu}^2=\bmu\T\A\bmu \).

*Necessity.* Apply @thm-qf-canonical with \( \bL=\I \), so that \( \A=\Q\bLambda\Q\T \).
Let \( \boldsymbol{\upnu}=\Q\T\bmu \). Then \( \bb=\Q\T\A\bmu=\bLambda\boldsymbol{\upnu} \), so \( b_j=\lambda_j\nu_j \). The
constant left after completing squares is
\( \bmu\T\A\bmu-\sum_{\lambda_j\ne0}\lambda_j\nu_j^2=\boldsymbol{\upnu}\T\bLambda\boldsymbol{\upnu}-\sum_j\lambda_j\nu_j^2=0 \).
Hence \( \Y\T\A\Y \) has the law of \( \sum_{j\in J}\lambda_j(W_j+\nu_j)^2 \), where \( J \) indexes
the nonzero eigenvalues. \( J \) is not empty, because a \( \chi^2(r,\gamma) \) variable with
\( r\ge1 \) is not identically zero.

By @lem-qf-shifted-square, near \( t=0 \) the logarithm of the mgf of this sum is
\( \sum_{j\in J}\{-\tfrac12\log(1-2\lambda_jt)+\lambda_j\nu_j^2t/(1-2\lambda_jt)\} \).
Expanding \( -\tfrac12\log(1-x)=\sum_{k\ge1}x^k/(2k) \) and \( x/(1-x)=\sum_{k\ge1}x^k \), the
coefficient of \( t^k \) is \( (2^{k-1}/k)\sum_{j\in J}\lambda_j^k(1+k\nu_j^2) \). The same
expansion of @eq-qf-ncchisq-mgf gives \( (2^{k-1}/k)(r+k\gamma) \). Equal
distributions have equal mgfs, and two power series that agree on an interval have
equal coefficients, so
\[
\sum_{j\in J}\lambda_j^k\,(1+k\nu_j^2)=r+k\gamma\qquad\text{for every }k\ge1 .
\]{#eq-qf-cumulant-match}

Let \( \rho=\max_{j\in J}\lvert\lambda_j\rvert \). For even \( k \) every term on the left is
nonnegative and one of them is at least \( \rho^k \). The right side grows only linearly
in \( k \), so \( \rho\le1 \). Split \( J \) into \( J_+=\{\lambda_j=1\} \), \( J_-=\{\lambda_j=-1\} \) and
\( J_0=\{0<\lvert\lambda_j\rvert<1\} \), and write \( N_\pm=\lvert J_\pm\rvert \) and
\( S_\pm=\sum_{J_\pm}\nu_j^2 \). The terms over \( J_0 \) tend to zero as \( k\to\infty \), so
\[
N_++(-1)^kN_-+k\bigl(S_++(-1)^kS_-\bigr)-r-k\gamma\;\longrightarrow\;0 .
\]
Dividing by \( k \) and letting \( k\to\infty \) through even and through odd values gives
\( S_++S_-=\gamma=S_+-S_- \). Then the same limit without the factor \( k \) gives
\( N_++N_-=r=N_+-N_- \). So \( N_-=S_-=0 \), \( N_+=r \) and \( S_+=\gamma \). The terms over \( J_+ \) in
@eq-qf-cumulant-match now account for all of \( r+k\gamma \), and with \( k=2 \) the
remaining sum \( \sum_{J_0}\lambda_j^2(1+2\nu_j^2) \) is zero. Every term of that sum is
positive, so \( J_0 \) is empty. Thus every eigenvalue of \( \A \) is \( 0 \) or \( 1 \), \( \A \) is
idempotent, \( \rank(\A)=N_+=r \), and \( \bmu\T\A\bmu=\sum_{J_+}\nu_j^2=\gamma \).
:::

The necessity half is the one that makes the theorem a tool for *recognizing*
distributions. If a quadratic form in a spherical normal vector is chi-squared at all,
then it is the squared length of a projection. There is no other way.

## General covariance matrices

::: {#thm-qf-chisq}
[When a quadratic form is chi-squared]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) and let \( \A \) be symmetric.

::: {.enumerate options="label=(\alph*)"}
1. *(Projections.)* If \( \bSigma=\sigma^2\I \) with \( \sigma>0 \) and \( \bP \) is symmetric
           and idempotent, then
           \[
\frac{\norm{\bP\Y}^2}{\sigma^2}=\frac{\Y\T\bP\Y}{\sigma^2}
          \sim\chi^2\Bigl(\rank(\bP),\,\frac{\norm{\bP\bmu}^2}{\sigma^2}\Bigr).
\]

2. *(Positive definite \( \bSigma \).)* If \( \bSigma \) is positive definite, then
           \( \Y\T\A\Y \) has a \( \chi^2(r,\gamma) \) distribution for some \( r\ge1 \) and \( \gamma\ge0 \)
           iff \( \A\bSigma \) is idempotent and \( \A\ne\mathbf{0} \). In that case
           \( r=\rank(\A)=\tr(\A\bSigma) \) and \( \gamma=\bmu\T\A\bmu \).

3. *(Any \( \bSigma \).)* If
           \[
\text{(i) } \bSigma\A\bSigma\A\bSigma=\bSigma\A\bSigma,\qquad
          \text{(ii) } \bSigma\A\bSigma\A\bmu=\bSigma\A\bmu,\qquad
          \text{(iii) } \bmu\T\A\bSigma\A\bmu=\bmu\T\A\bmu,
\]
           then \( \Y\T\A\Y\sim\chi^2\bigl(\tr(\A\bSigma),\bmu\T\A\bmu\bigr) \). The three
           conditions hold in particular when \( \A\bSigma\A=\A \), and when \( \A\bSigma \) is
           idempotent and \( \bmu\in\C(\bSigma) \).
:::

:::

::: {.proof}
(a) \( \sigma^{-1}\Y\sim\Normal_n(\sigma^{-1}\bmu,\I) \). Apply @thm-qf-chisq-identity
with \( \A=\bP \), using \( \bP\T\bP=\bP \) to write \( \Y\T\bP\Y=\norm{\bP\Y}^2 \) and
\( \bmu\T\bP\bmu=\norm{\bP\bmu}^2 \).

(b) Let \( \bSigma^{1/2} \) be the positive definite square root
(@thm-mat-square-root). Then \( \Z=\bSigma^{-1/2}\Y\sim\Normal_n(\bSigma^{-1/2}\bmu,\I) \) and
\( \Y\T\A\Y=\Z\T\B\Z \) with \( \B=\bSigma^{1/2}\A\bSigma^{1/2} \) symmetric. By
@thm-qf-chisq-identity, \( \Y\T\A\Y \) is noncentral chi-squared iff \( \B \) is
idempotent and nonzero. Now \( \B^2=\bSigma^{1/2}\A\bSigma\A\bSigma^{1/2} \), so \( \B^2=\B \)
iff \( \A\bSigma\A=\A \) (multiply by \( \bSigma^{-1/2} \) on both sides) iff
\( \A\bSigma\A\bSigma=\A\bSigma \) (multiply by \( \bSigma \) or \( \bSigma^{-1} \) on the right). Also \( \B\ne\mathbf{0} \) iff \( \A\ne\mathbf{0} \).
The degrees of freedom are \( \rank(\B)=\rank(\A) \), which equals
\( \tr(\B)=\tr(\A\bSigma) \), and the noncentrality is
\( \bmu\T\bSigma^{-1/2}\B\bSigma^{-1/2}\bmu=\bmu\T\A\bmu \).

(c) Take \( \bL \) of size \( n\times k \) with full column rank \( k=\rank(\bSigma) \) and
\( \bL\bL\T=\bSigma \), for instance \( \bL=\bU_1\bD_1^{1/2} \) from the eigenvectors \( \bU_1 \)
for the positive eigenvalues \( \bD_1 \) of \( \bSigma \). (If \( \bSigma=\mathbf{0} \) the claim is
trivial.) As in @thm-qf-canonical, \( \Y\overset{d}{=}\bmu+\bL\Z \) with
\( \Z\sim\Normal_k(\bzero,\I) \). Since \( \bL \) has full column rank, \( \bL\T \) has full row rank, so
\( \bL\mathbf{G}\bL\T=\bL\bH\bL\T \) implies \( \mathbf{G}=\bH \), and \( \bL\mathbf{g}=\bL\mathbf{h} \) implies
\( \mathbf{g}=\mathbf{h} \). Put \( \B=\bL\T\A\bL \) and \( \mathbf{c}=\bL\T\A\bmu \). Condition (i) reads
\( \bL\B^2\bL\T=\bL\B\bL\T \), so \( \B^2=\B \). Condition (ii) reads \( \bL\B\mathbf{c}=\bL\mathbf{c} \), so
\( \B\mathbf{c}=\mathbf{c} \). Condition (iii) reads \( \mathbf{c}\T\mathbf{c}=\bmu\T\A\bmu \). Hence
\[
(\bmu+\bL\Z)\T\A(\bmu+\bL\Z)=\bmu\T\A\bmu+2\mathbf{c}\T\Z+\Z\T\B\Z
=(\Z+\mathbf{c})\T\B(\Z+\mathbf{c}),
\]
because \( (\Z+\mathbf{c})\T\B(\Z+\mathbf{c})=\Z\T\B\Z+2\mathbf{c}\T\B\Z+\mathbf{c}\T\B\mathbf{c} \) and
\( \B\mathbf{c}=\mathbf{c} \). Now \( \Z+\mathbf{c}\sim\Normal_k(\mathbf{c},\I) \), and @thm-qf-chisq-identity
gives the law \( \chi^2(\rank\B,\mathbf{c}\T\B\mathbf{c}) \), where \( \rank\B=\tr\B=\tr(\A\bSigma) \) and
\( \mathbf{c}\T\B\mathbf{c}=\mathbf{c}\T\mathbf{c}=\bmu\T\A\bmu \). If \( \B=\mathbf{0} \), then \( \mathbf{c}=\B\mathbf{c}=\bzero \) and
the form is identically \( 0 \), which is \( \chi^2(0) \). If \( \A\bSigma\A=\A \), all three
conditions follow by substitution. If \( \A\bSigma \) is idempotent and \( \bmu=\bSigma\bv \),
then (i) follows from multiplying \( \A\bSigma\A\bSigma=\A\bSigma \) by \( \bSigma \) on the left,
and (ii) and (iii) from the same identity applied to \( \bv \).
:::

Part (a) is the case used most often in this book. For the linear model
\( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) and the projection \( \M \) onto \( \C(\X) \),
[Chapter 6](../ch06-projections/index.html) shows that the residual sum of squares is
\( \text{SSE}=\norm{(\I-\M)\Y}^2 \). Since \( (\I-\M)\X\bbeta=\bzero \), part (a) gives
\[
\text{SSE}/\sigma^2\sim\chi^2\bigl(n-\rank(\X)\bigr),
\]
a central distribution whatever the value of \( \bbeta \). The projections that define
test statistics, by contrast, do not annihilate the mean under the alternative, and
their noncentrality \( \norm{\bP\X\bbeta}^2/\sigma^2 \) measures how far the truth is from
the hypothesis in the metric that the test can see.

In part (c) idempotence of \( \A\bSigma \) alone is *not* enough when \( \bSigma \) is
singular. A mean with a component outside \( \C(\bSigma) \) is a deterministic shift that
the form may pick up as a constant or as a linear term (@exr-qf-singular). Part (c) also covers the
natural generalization of the Mahalanobis distance.

::: {#cor-qf-mahalanobis}
[Mahalanobis forms]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \).

::: {.enumerate options="label=(\alph*)"}
1. If \( \bSigma \) is positive definite, \( (\Y-\bmu)\T\bSigma^{-1}(\Y-\bmu)\sim\chi^2(n) \) and
           \( \Y\T\bSigma^{-1}\Y\sim\chi^2(n,\bmu\T\bSigma^{-1}\bmu) \).

2. If \( \G \) is a symmetric generalized inverse of \( \bSigma \) and \( \bmu\in\C(\bSigma) \), then
           \( \Y\T\G\Y\sim\chi^2\bigl(\rank(\bSigma),\bmu\T\G\bmu\bigr) \).
:::

:::

::: {.proof}
(a) \( \A=\bSigma^{-1} \) satisfies \( \A\bSigma\A=\A \), and \( \Y-\bmu\sim\Normal_n(\bzero,\bSigma) \). The
central statement was already obtained by whitening in @prp-mvn-mahalanobis.
(b) \( \G\bSigma \) is idempotent because \( \bSigma\G\bSigma=\bSigma \), so part (c) of the theorem
applies and the degrees of freedom are \( \tr(\G\bSigma)=\tr(\bL\T\G\bL) \). With \( \bL \) as in
that proof, \( \bL\bL\T\G\bL\bL\T=\bL\bL\T \) gives \( \bL\T\G\bL=\I_k \), whose trace is
\( k=\rank(\bSigma) \).
:::

A symmetric generalized inverse always exists (@prp-mat-ginverse-props(d)). The value of \( \Y\T\G\Y \) does not
depend on which one is chosen, with probability one (@exr-qf-ginverse-form).

::: {#exm-qf-gls-centring}
[Correlated observations about a common level]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) with \( \bSigma \) positive definite, and suppose we want to
measure how far \( \Y \) is from a common level \( m\bone \). The generalized least squares level is
\( \hat{m}=\bone\T\bSigma^{-1}\Y/\bone\T\bSigma^{-1}\bone \), and the natural discrepancy is
\( (\Y-\hat{m}\bone)\T\bSigma^{-1}(\Y-\hat{m}\bone)=\Y\T\A\Y \) with
\[
\A=\bSigma^{-1}-\frac{\bSigma^{-1}\bone\bone\T\bSigma^{-1}}{\bone\T\bSigma^{-1}\bone}.
\]
A direct multiplication gives \( \A\bSigma\A=\A \), and \( \tr(\A\bSigma)=n-1 \). By
@thm-qf-chisq(b) the form is noncentral chi-squared with \( n-1 \) degrees of freedom and
noncentrality \( \bmu\T\A\bmu \), which is zero exactly when \( \bmu \) is constant. The listing takes \( n=8 \), an
autoregressive correlation \( \bSigma=(\rho^{\lvert i-j\rvert}) \) with
\( \rho=0.6 \), and a linear trend in the mean, for which
\( \gamma=4.200 \). [Figure 4.3.1](03-chisq.html#fig-qf-chisq)(a) compares \( 200{,}000 \)
simulated values with the \( \chi^2(7,4.200) \) density.
:::

```{.python .run #cell-quadform-chisq-gls}
import numpy as np
from scipy import integrate, stats
rng = np.random.default_rng(31415)
n, rho = 8, 0.6
idx = np.arange(n)
Sigma = rho ** np.abs(idx[:, None] - idx[None, :])     # AR(1) correlation matrix
mu = 0.4 * (idx - idx.mean())                          # linear trend, mean zero
one = np.ones(n)

W = np.linalg.inv(Sigma)
A = W - np.outer(W @ one, W @ one) / (one @ W @ one)   # GLS-centred form
AS = A @ Sigma
print("A Sigma idempotent:", np.allclose(AS @ AS, AS))
df, gamma = np.trace(AS), mu @ A @ mu                  # chi^2(df, gamma)

L = np.linalg.cholesky(Sigma)                          # Y = mu + L Z
Z = rng.standard_normal((200_000, n))
Y = mu + Z @ L.T
q_gls = np.einsum("ij,jk,ik->i", Y, A, Y)
print(df, gamma, stats.kstest(q_gls, stats.ncx2(df, gamma).cdf).pvalue)
```

::: {when-format="html"}
![**Figure 4.3.1.** Two quadratic forms measuring spread about a common level, for \( n=8 \)
autoregressive observations with a trending mean. (a) With generalized least squares
centring, \( \A\bSigma \) is idempotent and the form is exactly noncentral chi-squared.
(b) With ordinary centring it is a weighted sum of noncentral \( \chi^2(1) \) terms.
A two-moment scaled chi-squared fits closely. The chi-squared law that would hold
for independent observations does not.](quadform_chisq.svg){#fig-qf-chisq width=100%}
:::

::: {when-format="pdf"}
![Two quadratic forms measuring spread about a common level, for \( n=8 \)
autoregressive observations with a trending mean. (a) With generalized least squares
centring, \( \A\bSigma \) is idempotent and the form is exactly noncentral chi-squared.
(b) With ordinary centring it is a weighted sum of noncentral \( \chi^2(1) \) terms.
A two-moment scaled chi-squared fits closely. The chi-squared law that would hold
for independent observations does not.](quadform_chisq.pdf){width=100%}
:::

## Forms that are not chi-squared

Most quadratic forms fail the conditions of @thm-qf-chisq. The ordinary sum of
squares about the mean, \( \Y\T\mathbf{C}\Y \) with \( \mathbf{C}=\I-n^{-1}\bone\bone\T \), is an example as
soon as the observations are correlated. Then \( \mathbf{C}\bSigma \) is not idempotent, and by
@thm-qf-canonical the law is that of \( \sum_j\lambda_j(W_j+\nu_j)^2 \), with weights
\( \lambda_j \) equal to the nonzero eigenvalues of \( \mathbf{C}\bSigma \). In
@exm-qf-gls-centring these range from \( 0.259 \) to
\( 1.832 \). Two practical ways to handle such a law are the following.

**Exact computation.**  The characteristic function of the weighted sum is known
in closed form (put \( t=iu \) in @cor-qf-mgf), and a distribution function can be
recovered from a characteristic function \( \varphi \) by the inversion formula of
Gil-Pelaez (1951),
\[
\Pr(X>x)=\frac12+\frac1\pi\int_0^\infty\frac{\operatorname{Im}\{e^{-iux}\varphi(u)\}}{u}\,du .
\]
For a weighted sum of \( k \) chi-squared terms the integrand decays like \( u^{-1-k/2} \), so
ordinary quadrature works well. Imhof (1961) gave the real-variable form of
this integral that most software uses. The listing implements it directly and checks it
against simulation.

**Two-moment approximation.**  Match \( a\chi^2(\nu) \) to the mean \( m \) and variance
\( v \) from @thm-qf-mean-var. Since \( a\chi^2(\nu) \) has mean \( a\nu \) and variance
\( 2a^2\nu \), this gives \( a=v/(2m) \) and \( \nu=2m^2/v \), with \( \nu \) usually not an integer
(Satterthwaite 1946; Patnaik 1949).

In @exm-qf-gls-centring, the simulated \( 95 \)th percentile of \( \Y\T\mathbf{C}\Y \) is
\( 26.51 \). The exact tail probability beyond it is
\( 0.0490 \), and the approximation with \( a=2.488 \) and
\( \nu=4.68 \) gives \( 0.0484 \). Treating the observations as
independent, and so using \( \chi^2(7,6.720) \), gives
\( 0.040 \). [Figure 4.3.1](03-chisq.html#fig-qf-chisq)(b) shows the three laws. Here the error from
ignoring correlation is modest at the \( 5\% \) point but visible in the shape. It is not a
monotone function of \( \rho \), but strong correlation makes it severe: with
\( \rho=0.9 \) and the same mean, the naive tail probability at the true
\( 95 \)th percentile is \( 0.197 \).

```{.python .run #cell-quadform-chisq-canonical}
def canonical(A, mu, L):
    """Y^T A Y = sum_j lam_j (W_j + nu_j)^2 + c (all lam_j nonzero here), W ~ N(0, I)."""
    lam, Qm = np.linalg.eigh(L.T @ A @ L)
    b = Qm.T @ L.T @ A @ mu
    keep = np.abs(lam) > 1e-10
    assert np.allclose(b[~keep], 0)                    # no linear terms in this example
    c = mu @ A @ mu - np.sum(b[keep] ** 2 / lam[keep])
    return lam[keep], b[keep] / lam[keep], c

def tail_prob(x, lam, nu, c=0.0):
    """P(sum lam_j (W_j + nu_j)^2 + c > x) by Gil-Pelaez inversion."""
    def integrand(u):
        z = 1 - 2j * lam * u
        phi = np.prod(z ** -0.5 * np.exp(1j * u * lam * nu**2 / z))
        return (np.exp(-1j * u * (x - c)) * phi).imag / u
    edges = np.linspace(0, 400, 401)                   # integrand decays like u^(-1-k/2)
    val = sum(integrate.quad(integrand, a, b)[0] for a, b in zip(edges[:-1], edges[1:]))
    return 0.5 + val / np.pi
```

## Exercises

### A. Check your understanding

::: {#exr-qf-which-idempotent}
[A1]

Let \( \Y\sim\Normal_n(\bzero,\I) \) with \( n\ge2 \), and \( \A=\alpha\I+\beta\bone\bone\T \). Find all
\( (\alpha,\beta) \) for which \( \Y\T\A\Y \) has a chi-squared distribution, and identify each form.
:::

::: {.solution}
\( \A \) has eigenvalue \( \alpha+n\beta \) on \( \bone \) and \( \alpha \)
on \( \bone\perpc \). By @thm-qf-chisq-identity, the form is chi-squared iff \( \A\ne\mathbf{0} \) and both
eigenvalues lie in \( \{0,1\} \). The three solutions are \( (1,0) \), giving \( \norm{\Y}^2\sim\chi^2(n) \);
\( (0,1/n) \), giving \( n\bar{Y}^2\sim\chi^2(1) \); and \( (1,-1/n) \), giving \( \sum(Y_i-\bar{Y})^2\sim\chi^2(n-1) \).
:::

::: {#exr-qf-two-sample}
[A2]

Two independent samples of sizes \( m \) and \( k \) come from \( \Normal(\mu_1,\sigma^2) \) and
\( \Normal(\mu_2,\sigma^2) \). Write the pooled within-sample sum of squares as \( \norm{(\I-\M)\Y}^2 \) for a
suitable projection, and show that the two-sample \( t \) statistic has the
\( t\bigl(m+k-2,(\mu_1-\mu_2)/(\sigma\sqrt{1/m+1/k})\bigr) \) distribution.
:::

### B. Practice

::: {#exr-qf-equicorrelated}
[B1]

Let \( \Y\sim\Normal_n\bigl(\mu\bone,\sigma^2\{(1-\rho)\I+\rho\bone\bone\T\}\bigr) \) with
\( -1/(n-1)<\rho<1 \), and let \( \bar{Y} \) and \( S^2 \) be the sample mean and variance. Show that
\( (n-1)S^2/\{\sigma^2(1-\rho)\}\sim\chi^2(n-1) \) and that \( S^2 \) is independent of \( \bar{Y} \). Deduce
that \( \sqrt n(\bar{Y}-\mu)/S \) has the law of \( \{(1+(n-1)\rho)/(1-\rho)\}^{1/2} \) times a \( t(n-1) \)
variable, and describe what a positive \( \rho \) does to the coverage of the usual \( 95\% \) \( t \) interval
for \( \mu \).
:::

::: {#exr-qf-singular}
[B2]

Let \( \bSigma=\diag(1,0) \), \( \bmu=(0,1)\T \), \( \A=\I_2 \) and \( \Y\sim\Normal_2(\bmu,\bSigma) \). Show that
\( \A\bSigma \) is idempotent but \( \Y\T\A\Y \) is not chi-squared. Which of the conditions of
@thm-qf-chisq(c) fail? What changes if \( \bmu=(m,0)\T \)? Repeat the first part with
\( \A=\begin{psmallmatrix}1&1\\1&0\end{psmallmatrix} \), where the shift enters through a linear term.
:::

::: {.solution}
\( \A\bSigma=\bSigma \) is idempotent. But \( Y_2=1 \) with probability
one, so \( \Y\T\A\Y=Y_1^2+1\ge1 \), while every \( \chi^2(r,\gamma) \) with \( r\ge1 \) gives positive probability
to \( (0,1) \). Conditions (i) and (ii) hold, since \( \bSigma^3=\bSigma \) and \( \bSigma^2\bmu=\bzero=\bSigma\bmu \).
Condition (iii) fails, since \( \bmu\T\bSigma\bmu=0\ne1=\bmu\T\bmu \). The part of the mean outside
\( \C(\bSigma) \) adds a constant. If \( \bmu=(m,0)\T\in\C(\bSigma) \), then \( \Y\T\A\Y=Y_1^2\sim\chi^2(1,m^2) \), as
@thm-qf-chisq(c) predicts. For \( \A=\begin{psmallmatrix}1&1\\1&0\end{psmallmatrix} \) and \( \bmu=(0,1)\T \),
\( \A\bSigma=\begin{psmallmatrix}1&0\\1&0\end{psmallmatrix} \) is again idempotent, and
\( \Y\T\A\Y=Y_1^2+2Y_1Y_2=Y_1^2+2Y_1=(Y_1+1)^2-1 \), which is negative with positive probability. Conditions (i)
and (ii) hold, since \( \bSigma\A\bSigma=\bSigma \) and \( \bSigma\A\bSigma\A\bmu=\bSigma\A\bmu=(1,0)\T \), but
(iii) fails because \( \bmu\T\A\bSigma\A\bmu=1\ne0=\bmu\T\A\bmu \). The shift now enters through the
linear term \( 2Y_1 \).
:::

::: {#exr-qf-ginverse-form}
[B3]

In @cor-qf-mahalanobis(b), show that \( \Y\T\G\Y \) takes the same value for every generalized
inverse \( \G \) of \( \bSigma \), with probability one. *Hint:* @thm-rv-cov-nnd.
:::

::: {.solution}
By @thm-rv-cov-nnd, \( \Y-\bmu\in\C(\bSigma) \) with
probability one. Since \( \bmu\in\C(\bSigma) \), also \( \Y\in\C(\bSigma) \), so \( \Y=\bSigma\bv \) for some
(random) \( \bv \). For any generalized inverse \( \G \),
\( \Y\T\G\Y=\bv\T\bSigma\G\bSigma\bv=\bv\T\bSigma\bv \), which does not involve \( \G \).
:::

::: {#exr-qf-satterthwaite}
[B4]

Let \( \Y\sim\Normal_n(\bzero,\bSigma) \) with \( \bSigma \) positive definite, and let \( \A\ne\mathbf{0} \) be
nonnegative definite. Show that the two-moment approximation \( a\chi^2(\nu) \) to \( \Y\T\A\Y \) has
\( \nu=\{\tr(\A\bSigma)\}^2/\tr\{(\A\bSigma)^2\} \), and that \( \nu\le\rank(\A) \), with equality iff all
nonzero eigenvalues of \( \A\bSigma \) are equal.
:::

### C. Going deeper

::: {#exr-qf-necessity-singular}
[C1]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) with \( \bSigma \) possibly singular, and suppose
\( \Y\T\A\Y\sim\chi^2(r,\gamma) \) with \( r\ge1 \). Show that \( \bSigma\A\bSigma\A\bSigma=\bSigma\A\bSigma \).
*Hint:* in @cor-qf-mgf, compare the coefficients of \( t^k \) for \( k\ge3 \) with those of
@eq-qf-ncchisq-mgf, and adapt the argument of @thm-qf-chisq-identity.
:::
