# Bayesian generalized linear models and MCMC

@prp-prc-quadratic closed off the closed form. What is left is a posterior known only up
to a constant, and two ways to work with it: approximate it by a normal density at its
mode, accurate to order \( n^{-1} \) and the origin of the BIC; or sample from it,
accurate to whatever the Monte Carlo error allows. This section does both, with the
samplers written out rather than delegated.

## The posterior

::: {#def-prc-posterior}
[The Bayesian generalized linear model]

Let \( \y \) follow the generalized linear model of @def-glm-model with link \( g \),
dispersion \( \phi \) and log-likelihood \( \ell(\bbeta) \), and let \( \pi \) be a prior
density on \( \Real^p \). The **posterior density** of \( \bbeta \) is
\[
\pi(\bbeta\mid\y)=\frac{\exp\{\ell(\bbeta)\}\,\pi(\bbeta)}{m(\y)},\qquad
m(\y)=\int_{\Real^p}\exp\{\ell(\bb)\}\,\pi(\bb)\,d\bb ,
\]{#eq-prc-posterior}

whenever \( 0<m(\y)<\infty \); the posterior is then **proper**. The **maximum a
posteriori** estimate \( \hbeta_{\text{MAP}} \) maximizes @eq-prc-posterior, equivalently
maximizes the penalized log-likelihood \( \ell(\bbeta)+\log\pi(\bbeta) \). The default
choice throughout this chapter is the independent normal prior
\( \bbeta\sim\Normal_p(\mathbf{m},\bSigma) \), for which
\[
\log\pi(\bbeta\mid\y)=\ell(\bbeta)
-\tfrac12(\bbeta-\mathbf{m})\T\bSigma^{-1}(\bbeta-\mathbf{m})+\text{const}.
\]{#eq-prc-normal-prior-post}
:::

Because @eq-prc-normal-prior-post is concave whenever \( \ell \) is (@thm-glm-concave),
\( \hbeta_{\text{MAP}} \) is found by the Fisher scoring of @thm-glm-irls with the score
replaced by \( \X\T\W\bD^{-1}(\y-\bmu)/\phi-\bSigma^{-1}(\bbeta-\mathbf{m}) \) and the
information by \( \X\T\W\X/\phi+\bSigma^{-1} \): one extra term in each. What the prior
buys is not only shrinkage but existence.

::: {#prp-prc-propriety}
[When a flat prior gives a posterior]

Let \( \pi \) be a prior on \( \Real^p \).

::: {.enumerate options="label=(\alph*)"}
1. If \( \pi \) is proper and the likelihood is bounded, the posterior is proper. Both the
   logistic and the Poisson likelihoods are bounded.

2. Take the improper prior \( \pi\propto1 \) and the ungrouped logistic model with
   \( \rank(\X)=p \). Then the posterior is proper if and only if the data are in overlap
   in the sense of @def-bin-separation, that is if and only if the maximum likelihood
   estimate exists (@thm-bin-separation).
:::
:::

::: {.proof}
(a) If \( L\le L_{\max}<\infty \) then \( m(\y)\le L_{\max}\int\pi=L_{\max}<\infty \),
and \( m(\y)>0 \) because the integrand is positive on a set of positive measure. The
logistic likelihood is a product of probabilities, so \( L\le1 \); the Poisson likelihood
\( \prod_i\mu_i^{y_i}e^{-\mu_i}/y_i! \), maximized freely over
\( \bmu\in(0,\infty)^n \) at \( \mu_i=y_i \), is bounded by
\( \prod_iy_i^{y_i}e^{-y_i}/y_i! \) (with \( 0^0=1 \)).

(b) Write \( s_i=2y_i-1 \) and use @eq-bin-signed-loglik:
\( -\ell(\bbeta)=\sum_i\log\{1+\exp(-s_i\x_{(i)}\T\bbeta)\} \).

*Overlap implies propriety.* Since \( \log(1+e^{-t})\ge\max(0,-t)=t_- \),
\[
-\ell(\bbeta)\ \ge\ \sum_{i=1}^n\bigl(s_i\x_{(i)}\T\bbeta\bigr)_-\ =:\ c(\bbeta),
\]
with \( c \) continuous and homogeneous of degree one. Under overlap \( c(\bb)>0 \) for
every \( \bb\ne\bzero \), since \( c(\bb)=0 \) forces \( s_i\x_{(i)}\T\bb\ge0 \) for
all \( i \), which overlap permits only at the origin. The unit sphere is compact, so
\( \delta=\min_{\norm{\bb}=1}c(\bb)>0 \) and, by homogeneity,
\( c(\bbeta)\ge\delta\norm{\bbeta} \). Hence
\( \exp\{\ell(\bbeta)\}\le e^{-\delta\norm{\bbeta}} \), integrable over
\( \Real^p \).

*Separation implies impropriety.* Let \( \bb\ne\bzero \) satisfy
\( s_i\x_{(i)}\T\bb\ge0 \) for all \( i \), put \( S=\{i:\x_{(i)}\T\bb=0\} \),
\( \varrho=\min_{i\notin S}s_i\x_{(i)}\T\bb>0 \), \( \omega=\max_i\norm{\x_{(i)}} \),
and fix \( \rho>0 \). For \( \bbeta=t\bb+\bu \) with
\( t\ge t_0:=2\omega\rho/\varrho \) and \( \norm{\bu}\le\rho \), every
\( i\notin S \) has \( s_i\x_{(i)}\T\bbeta\ge t\varrho-\omega\rho\ge0 \) and every
\( i\in S \) has \( s_i\x_{(i)}\T\bbeta=s_i\x_{(i)}\T\bu\ge-\omega\rho \). So
\( -\ell(\bbeta)\le n\log\{1+e^{\omega\rho}\} \) on that whole tube, giving
\( \exp\{\ell(\bbeta)\}\ge\{1+e^{\omega\rho}\}^{-n} \) there; and the tube has
infinite Lebesgue measure, containing a half-infinite cylinder about the ray. Hence
\( m(\y)=\infty \).
:::

::: {.warning}
The second part is not a curiosity. Separation is common in small samples and with rare
outcomes, and a "noninformative" Bayesian analysis of separated binary data summarizes a
measure of infinite mass: whatever the sampler prints is an artefact of where it wandered.
@exm-prc-separated-bayes makes the divergence visible.
:::

## The Laplace approximation

For large \( n \) the integral in @eq-prc-posterior is dominated by a shrinking
neighbourhood of the mode, where the integrand looks normal. Write
\( q_n(\bbeta)=\ell(\bbeta)+\log\pi(\bbeta) \) for the unnormalized log posterior,
\( \hbeta_n \) for its maximizer and \( \bH_n=-\nabla^2q_n(\hbeta_n) \).

::: {#thm-prc-laplace}
[Laplace approximation]

Let \( \rho_n=n^{-3/8} \) and \( B_n=\{\bbeta:\norm{\bbeta-\hbeta_n}\le\rho_n\} \), and
suppose there are constants \( 0<c\le C<\infty \) and \( a,K>0 \), free of \( n \), with:

::: {.enumerate options="label=(L\arabic*)"}
1. \( q_n \) is four times continuously differentiable on \( B_n \), with
   \( \nabla q_n(\hbeta_n)=\bzero \) and \( \bH_n \) positive definite;

2. \( cn\le\lambda_{\min}(\bH_n)\le\lambda_{\max}(\bH_n)\le Cn \), and every third and
   fourth partial derivative of \( q_n \) is bounded in absolute value by \( Cn \) on
   \( B_n \);

3. \( \int_{\Real^p\setminus B_n}\exp\{q_n(\bbeta)-q_n(\hbeta_n)\}\,d\bbeta
   \le K\exp(-an^{1/4}) \).
:::

Then
\[
\int_{\Real^p}e^{q_n(\bbeta)}\,d\bbeta
=(2\pi)^{p/2}\det(\bH_n)^{-1/2}e^{q_n(\hbeta_n)}\bigl\{1+O(n^{-1})\bigr\}.
\]{#eq-prc-laplace}
:::

::: {.proof}
Write \( I_n \) for the integral and \( J_n=I_n\det(\bH_n)^{1/2}e^{-q_n(\hbeta_n)} \), so
the claim is \( J_n=(2\pi)^{p/2}\{1+O(n^{-1})\} \). Substitute
\( \bbeta=\hbeta_n+\bH_n^{-1/2}\bu \), Jacobian \( \det(\bH_n)^{-1/2} \):
\[
J_n=\int_{\Real^p}\exp\bigl\{q_n(\hbeta_n+\bH_n^{-1/2}\bu)-q_n(\hbeta_n)\bigr\}\,d\bu .
\]
Because \( \nabla q_n(\hbeta_n)=\bzero \), Taylor's theorem on \( B_n \) gives
\[
q_n(\hbeta_n+\bH_n^{-1/2}\bu)-q_n(\hbeta_n)=-\tfrac12\norm{\bu}^2+A_n(\bu)+r_n(\bu),
\]
where \( A_n \) is the cubic form built from the third derivatives at \( \hbeta_n \) and
\( r_n \) is the fourth-order remainder. By (L2) the derivatives are bounded by \( Cn \)
and \( \norm{\bH_n^{-1/2}\bu}\le(cn)^{-1/2}\norm{\bu} \), so
\[
|A_n(\bu)|\le c_3n^{-1/2}\norm{\bu}^3,\qquad
|r_n(\bu)|\le c_4n^{-1}\norm{\bu}^4,
\]{#eq-prc-taylor-bounds}

with \( c_3=p^3C/(6c^{3/2}) \) and \( c_4=p^4C/(24c^2) \), valid whenever
\( \hbeta_n+\bH_n^{-1/2}\bu\in B_n \).

Let \( S_n=\{\bu:\norm{\bu}\le\sqrt c\,n^{1/8}\} \) and \( U_n \) the image of
\( B_n \). Every \( \bu\in S_n \) has
\( \norm{\bH_n^{-1/2}\bu}\le(cn)^{-1/2}\norm{\bu}\le n^{-3/8}=\rho_n \), so
\( S_n\subseteq U_n \); and every \( \bu\in U_n \) has
\( \norm{\bu}\le\lambda_{\max}(\bH_n)^{1/2}\rho_n\le\sqrt C\,n^{1/8} \). Split
\( J_n \) over \( S_n \), over \( U_n\setminus S_n \) and outside \( U_n \).

*Outside \( U_n \).* Changing back to \( \bbeta \) and using (L3) with
\( \det(\bH_n)^{1/2}\le(Cn)^{p/2} \) bounds that piece by \( (Cn)^{p/2}Ke^{-an^{1/4}} \),
which is \( o(n^{-1}) \).

*Between \( S_n \) and \( U_n \).* There \( \norm{\bu}\le\sqrt Cn^{1/8} \), so
@eq-prc-taylor-bounds gives \( |A_n|+|r_n|\le c_3C^{3/2}n^{-1/8}+c_4C^2n^{-1/2}\le1 \) for
large \( n \); the integrand is then at most \( e^{1-\norm{\bu}^2/2} \), whose integral
over \( \norm{\bu}>\sqrt c\,n^{1/8} \) is \( O(e^{-cn^{1/4}/4})=o(n^{-1}) \).

*Inside \( S_n \).* Put \( R_n=A_n+r_n \), so \( |R_n|\le1 \) there for large \( n \).
From \( |e^x-1-x-\tfrac12x^2|\le\tfrac16|x|^3e^{|x|} \) for \( |x|\le1 \),
\[
\begin{aligned}
\int_{S_n}e^{-\norm{\bu}^2/2}e^{R_n}d\bu
&=\int_{S_n}e^{-\norm{\bu}^2/2}\Bigl\{1+R_n+\tfrac12R_n^2\Bigr\}d\bu+\mathcal E_n,\\
|\mathcal E_n|&\le\tfrac e6\int_{S_n}e^{-\norm{\bu}^2/2}|R_n|^3d\bu .
\end{aligned}
\]
Four pieces. First,
\( \int_{S_n}e^{-\norm{\bu}^2/2}d\bu=(2\pi)^{p/2}+O(e^{-cn^{1/4}/2}) \). Second,
\( A_n \) is homogeneous of degree three, hence odd, and \( S_n \) is a ball about the
origin, so \( \int_{S_n}e^{-\norm{\bu}^2/2}A_n(\bu)d\bu=0 \) exactly. Third,
\( \int_{S_n}e^{-\norm{\bu}^2/2}|r_n|d\bu\le c_4n^{-1}\int e^{-\norm{\bu}^2/2}
\norm{\bu}^4d\bu=O(n^{-1}) \). Fourth, by @eq-prc-taylor-bounds
\( R_n^2\le2c_3^2n^{-1}\norm{\bu}^6+2c_4^2n^{-2}\norm{\bu}^8 \), of Gaussian integral
\( O(n^{-1}) \), and \( |R_n|^3=O(n^{-3/2}\norm{\bu}^9)+O(n^{-3}\norm{\bu}^{12}) \), of
Gaussian integral \( O(n^{-3/2}) \). Adding the three regions gives
\( J_n=(2\pi)^{p/2}+O(n^{-1}) \), which is @eq-prc-laplace.
:::

::: {.remark}
[When the conditions hold]

For a generalized linear model with a proper, positive, twice continuously differentiable
prior and the conditions (G1)–(G4) of @thm-glm-asymptotics, (L1) and (L2) hold with
probability tending to one, because \( \bH_n=\X\T\hat{\W}\X/\phi+O(1) \) and (G4) makes
its smallest eigenvalue grow linearly. Condition (L3) says the posterior concentrates;
verifying it in general is a separate argument this book does not give. Kass and Raftery
(1995, section 4) discuss the regularity conditions and the accuracy of the approximation.
:::

::: {#cor-prc-bic}
[Where the BIC comes from]

Under the conditions of @thm-prc-laplace, suppose in addition that \( \bH_n/n\to\bar{\bH} \)
positive definite, that the prior is continuous and positive at the limit of
\( \hbeta_n \), and that \( \hbeta \) maximizes \( \ell \) with
\( \ell(\hbeta)-\ell(\hbeta_n)=O(n^{-1}) \). Then
\[
-2\log m(\y)=-2\ell(\hbeta)+p\log n+O(1),
\]{#eq-prc-bic}

where the \( O(1) \) term is
\( \log\det(\bH_n/n)-p\log(2\pi)-2\log\pi(\hbeta_n)+o(1) \). The leading two terms are the
BIC of @def-sel-aic-bic.
:::

::: {.proof}
Take logarithms in @eq-prc-laplace:
\[
2\log m(\y)=2\ell(\hbeta_n)+2\log\pi(\hbeta_n)+p\log(2\pi)-\log\det\bH_n+O(n^{-1}).
\]
Write \( \log\det\bH_n=p\log n+\log\det(\bH_n/n) \) and replace \( \ell(\hbeta_n) \) by
\( \ell(\hbeta) \) at a cost of \( O(n^{-1}) \); rearranging gives @eq-prc-bic with the
stated remainder, which converges because \( \bH_n/n\to\bar{\bH} \).
:::

The remainder in @eq-prc-bic is \( O(1) \), not \( o(1) \): the BIC stays a bounded
distance from \( -2\log m(\y) \), and that distance carries the whole influence of the
prior, so a BIC difference is only a crude approximation to a log Bayes
factor (@prp-sel-aic-bic). What makes it useful is that a difference in dimension
contributes \( \log n \), which grows, while the neglected term does not. A third
criterion belongs to the same family: the deviance information criterion of Spiegelhalter,
Best, Carlin and van der Linde (2002) replaces \( p \) by an effective number of
parameters estimated from posterior draws, which suits a model whose parameters are
themselves shrunk by a prior.

::: {#exm-prc-laplace}
[How good is the approximation]

Take a logistic regression with an intercept and one covariate, a
\( \Normal(0,5^2) \) prior on each coefficient, and data simulated at
\( \bbeta=(-0.4,1.2)\T \). Adaptive Gauss–Hermite quadrature on a tensor grid centred at
the mode computes @eq-prc-posterior to machine accuracy, so the Laplace error can be
measured exactly.

At \( n=25 \) the Laplace marginal likelihood is too small by
\( 0.10700 \), or eleven per cent; the relative error falls to
\( 0.001114 \) at \( n=1600 \). Regressing the log error on
\( \log n \) over seven sample sizes gives a slope of \( -1.06 \),
and over those seven sizes the product \( n\times\text{error} \) stays between
\( 1.530 \) and \( 2.675 \): the
\( O(n^{-1}) \) of @eq-prc-laplace, over a decade and a half of sample sizes
([Figure 39.4.1](04-bayesian-glms-mcmc.html#fig-prc-laplace), panel (b)).

Panel (a) shows what that error looks like at \( n=25 \): the marginal posterior of the
slope is visibly right-skewed, mean \( 2.3291 \) against mode
\( 1.9271 \), with total variation distance
\( 0.1827 \) to the normal approximation. Reporting the mode with the
curvature as a standard error would put the estimate in the wrong place and the interval
on the wrong side. The corollary is visible too: over the same seven sizes the gap between
the BIC and \( -2\log m(\y) \) stays between \( -3.41 \) and
\( -1.61 \), an \( O(1) \) that does not grow.
:::

::: {when-format="html"}
![**Figure 39.4.1.** The Laplace approximation in a two-parameter logistic regression.
(a) The marginal posterior of the slope at 25 observations against its normal
approximation; the vertical line is the mode. (b) Relative error of the Laplace marginal
likelihood against sample size, on log scales, with a reference line of slope
\( -1 \).](laplace_error.svg){#fig-prc-laplace width=100%}
:::

::: {when-format="pdf"}
![The Laplace approximation in a two-parameter logistic regression.
(a) The marginal posterior of the slope at 25 observations against its normal
approximation; the vertical line is the mode. (b) Relative error of the Laplace marginal
likelihood against sample size, on log scales, with a reference line of slope
\( -1 \).](laplace_error.pdf){width=100%}
:::

```{.python .run #cell-laplace-laplace}
import numpy as np
from scipy import special, stats
TAU = 5.0                                       # prior standard deviation, N(0, TAU^2 I)
BETA0 = np.array([-0.4, 1.2])                   # the truth used to simulate


def simulate(n, seed):
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.ones(n), rng.uniform(-2, 2, n)])
    y = rng.binomial(1, special.expit(X @ BETA0)).astype(float)
    return X, y


def log_post(beta, X, y):
    """Log likelihood plus log prior, up to a constant free of beta."""
    eta = X @ beta
    return np.sum(y * eta - np.logaddexp(0.0, eta)) - np.sum(beta**2) / (2 * TAU**2)


def mode_and_hessian(X, y):
    """Posterior mode and the negative Hessian there, by Newton's method."""
    beta = np.zeros(X.shape[1])
    for _ in range(100):
        pi = special.expit(X @ beta)
        H = X.T @ ((pi * (1 - pi))[:, None] * X) + np.eye(len(beta)) / TAU**2
        g = X.T @ (y - pi) - beta / TAU**2
        beta = beta + np.linalg.solve(H, g)
    return beta, H


def log_marginal_laplace(X, y):
    """log of the Laplace value of the integral of exp(log_post)."""
    beta, H = mode_and_hessian(X, y)
    p = len(beta)
    return (log_post(beta, X, y) + p / 2 * np.log(2 * np.pi)
            - 0.5 * np.linalg.slogdet(H)[1] - p / 2 * np.log(2 * np.pi * TAU**2))


def log_marginal_quadrature(X, y, nodes=60):
    """The same integral by adaptive Gauss-Hermite quadrature: the reference value."""
    beta, H = mode_and_hessian(X, y)
    p = len(beta)
    L = np.linalg.cholesky(np.linalg.inv(H))
    z, w = np.polynomial.hermite_e.hermegauss(nodes)            # weight exp(-z^2/2)
    Z = np.stack(np.meshgrid(z, z, indexing="ij"), -1).reshape(-1, p)
    logw = np.log(w)[:, None] + np.log(w)[None, :]
    pts = beta + Z @ L.T
    vals = np.array([log_post(b, X, y) for b in pts]) + (Z**2).sum(1) / 2 + logw.ravel()
    return (special.logsumexp(vals) + np.linalg.slogdet(L)[1]
            - p / 2 * np.log(2 * np.pi * TAU**2))


sizes = [25, 50, 100, 200, 400, 800, 1600]
rel_err = []
for k, n in enumerate(sizes):
    X, y = simulate(n, 3910 + k)
    a, b = log_marginal_laplace(X, y), log_marginal_quadrature(X, y)
    rel_err.append(abs(np.expm1(a - b)))
rel_err = np.array(rel_err)
slope = np.polyfit(np.log(sizes), np.log(rel_err), 1)[0]
print("n, relative error of the Laplace marginal likelihood:")
for n, e in zip(sizes, rel_err):
    print(f"  {n:5d}  {e:.3e}   n * error = {n * e:.3f}")
print(f"log-log slope {slope:.3f}")
```

Two refinements are quoted without proof. For a positive smooth \( f \), applying
@thm-prc-laplace to \( \exp\{q_n+\log f\} \) in the numerator of
\( \E\{f(\bbeta)\mid\y\}=\int fe^{q_n}/\int e^{q_n} \) raises the accuracy of the
ratio to \( O(n^{-2}) \), the leading errors cancelling (Tierney and Kadane, 1986); and
the Bernstein–von Mises theorem says the posterior of
\( \bH_n^{1/2}(\bbeta-\hbeta_n) \) converges in total variation to
\( \Normal_p(\bzero,\I) \), so a credible set is asymptotically a confidence set
(van der Vaart, 1998, chapter 10).

## Sampling the posterior

When \( n \) is small, or the posterior skewed as in @exm-prc-laplace, an approximation
at the mode is not enough. The alternative is a Markov chain whose stationary distribution
is the posterior.

::: {#def-prc-mcmc}
[Metropolis–Hastings, Gibbs, and the two diagnostics]

Let \( \pi(\cdot\mid\y) \) be a density known up to a constant.

::: {.enumerate options="label=(\alph*)"}
1. The **Metropolis–Hastings** chain with proposal density \( q(\cdot\mid\cdot) \) moves
   from \( \bbeta^{(t)} \) by drawing \( \bbeta^{*}\sim q(\cdot\mid\bbeta^{(t)}) \) and
   setting \( \bbeta^{(t+1)}=\bbeta^{*} \) with probability
   \[
\alpha(\bbeta^{*}\mid\bbeta^{(t)})
=\min\Bigl\{1,\ \frac{\pi(\bbeta^{*}\mid\y)\,q(\bbeta^{(t)}\mid\bbeta^{*})}
{\pi(\bbeta^{(t)}\mid\y)\,q(\bbeta^{*}\mid\bbeta^{(t)})}\Bigr\},
\]{#eq-prc-mh}

   and \( \bbeta^{(t+1)}=\bbeta^{(t)} \) otherwise. The normalizing constant cancels. A
   symmetric proposal, \( q(\bb\mid\bb')=q(\bb'\mid\bb) \), gives the **random-walk
   Metropolis** algorithm, in which \( \alpha \) is the ratio of posterior densities
   capped at one.

2. The **Gibbs sampler** partitions \( \boldsymbol{\uppsi}=(\boldsymbol{\uppsi}_1,\dots,
   \boldsymbol{\uppsi}_K) \) and replaces each block in turn by a draw from its **full
   conditional** \( \pi(\boldsymbol{\uppsi}_k\mid\boldsymbol{\uppsi}_{-k},\y) \). It is the
   special case of (a) in which every proposal is accepted.

3. For \( c \) chains of length \( N \), split each in half and let \( W \) be the mean of
   the \( 2c \) within-sequence variances and \( B \) the variance of the \( 2c \)
   sequence means. The **potential scale reduction** is
   \[
\hat R=\Bigl\{\frac{(N/2-1)W/(N/2)+B}{W}\Bigr\}^{1/2}.
\]{#eq-prc-rhat}

4. With \( \rho_t \) the lag-\( t \) autocorrelation of the stationary chain, the
   **effective sample size** of \( cN \) draws is
   \( \text{ESS}=cN/(1+2\sum_{t\ge1}\rho_t) \), the sum truncated by Geyer's rule in the
   form used here: add the autocorrelations in consecutive pairs
   \( \rho_{2k-1}+\rho_{2k} \), \( k=1,2,\dots \), and stop at the first pair that is
   negative.
:::
:::

::: {#prp-prc-mh-invariance}
[The chain has the right stationary distribution]

Let \( \pi \) be a density on \( \Real^p \) and \( q(\cdot\mid\cdot) \) a transition
density with \( q(\bb'\mid\bb)>0 \) whenever \( q(\bb\mid\bb')>0 \). The
Metropolis–Hastings kernel of @eq-prc-mh satisfies detailed balance with respect to
\( \pi \):
\[
\pi(\bb)\,q(\bb'\mid\bb)\,\alpha(\bb'\mid\bb)
=\pi(\bb')\,q(\bb\mid\bb')\,\alpha(\bb\mid\bb')
\qquad\text{for all }\bb\ne\bb' ,
\]
and consequently \( \pi \) is a stationary distribution of the chain.
:::

::: {.proof}
Fix \( \bb\ne\bb' \) and write \( u=\pi(\bb)q(\bb'\mid\bb) \),
\( v=\pi(\bb')q(\bb\mid\bb') \). By @eq-prc-mh, \( \alpha(\bb'\mid\bb)=\min(1,v/u) \)
and \( \alpha(\bb\mid\bb')=\min(1,u/v) \), so both sides of the identity equal
\( \min(u,v) \). For stationarity let \( P(\bb,\cdot) \) be the kernel; for any
measurable \( A \),
\[
\int P(\bb,A)\pi(\bb)\,d\bb
=\int_A\!\!\int\pi(\bb)q(\bb'\mid\bb)\alpha(\bb'\mid\bb)\,d\bb\,d\bb'
+\int_A\pi(\bb')\,\varrho(\bb')\,d\bb' ,
\]
where \( \varrho(\bb')=1-\int q(\bb''\mid\bb')\alpha(\bb''\mid\bb')d\bb'' \) is the
probability of staying put. By detailed balance the inner integral of the first term is
\( \pi(\bb')\int q(\bb\mid\bb')\alpha(\bb\mid\bb')d\bb=\pi(\bb')\{1-\varrho(\bb')\} \), so
the two terms add to \( \int_A\pi(\bb')d\bb' \).
:::

Detailed balance is all that is checked here; convergence *to* \( \pi \) from an
arbitrary start needs irreducibility and aperiodicity, and a central limit theorem for the
ergodic averages more still (Tierney, 1994). The one design decision is the proposal. A
covariance \( \varsigma^2\bSigma_{\text{prop}} \) mixes well only when
\( \bSigma_{\text{prop}} \) resembles the posterior covariance, and the Laplace matrix
\( \bH_n^{-1} \) of @thm-prc-laplace is the obvious candidate, at the cost of one Newton
solve. Roberts, Gelman and Gilks (1997) showed that in a high-dimensional limit the
optimal scale is \( \varsigma=2.38/\sqrt p \), with acceptance rate
near \( 0.234 \) — derived for a product target, used everywhere as a rule of thumb.

## Data augmentation for the probit model

A Gibbs sampler needs full conditionals that can be drawn from, and @prp-prc-quadratic
says the conditional of \( \bbeta \) is not normal. Albert and Chib (1993) observed that
enlarging the parameter space restores normality, through the latent variable
representation of the probit model (@def-bin-latent).

::: {#prp-prc-albert-chib}
[The Albert–Chib sampler]

Consider the ungrouped probit model \( \Pr(Y_i=1)=\Phi(\x_{(i)}\T\bbeta) \) with prior
\( \bbeta\sim\Normal_p(\mathbf{m},\bSigma) \). Introduce latent variables
\( Z_i\mid\bbeta\sim\Normal(\x_{(i)}\T\bbeta,1) \), independent, and set
\( Y_i=\mathbf{1}\{Z_i>0\} \). Then:

::: {.enumerate options="label=(\alph*)"}
1. the marginal distribution of \( \Y \) given \( \bbeta \) is the probit model, so the
   posterior of \( \bbeta \) is unchanged;

2. given \( \bz \) and \( \y \), \( \bbeta\sim\Normal_p(\V(\X\T\bz+\bSigma^{-1}\mathbf{m}),
   \V) \) with \( \V=(\X\T\X+\bSigma^{-1})^{-1} \);

3. given \( \bbeta \) and \( \y \), the \( Z_i \) are independent, \( Z_i \) being
   \( \Normal(\x_{(i)}\T\bbeta,1) \) truncated to \( (0,\infty) \) if \( y_i=1 \) and to
   \( (-\infty,0) \) if \( y_i=0 \).
:::
:::

::: {.proof}
(a) \( \Pr(Y_i=1\mid\bbeta)=\Pr(Z_i>0\mid\bbeta)=\Phi(\x_{(i)}\T\bbeta) \), and the
\( Z_i \) are independent given \( \bbeta \), hence so are the \( Y_i \).

(b) The joint density of \( (\bbeta,\bz,\y) \) is proportional to
\[
\exp\Bigl\{-\tfrac12(\bbeta-\mathbf{m})\T\bSigma^{-1}(\bbeta-\mathbf{m})
-\tfrac12\norm{\bz-\X\bbeta}^2\Bigr\}
\prod_i\mathbf{1}\{\text{sign}(z_i)\text{ matches }y_i\}.
\]
With \( \bz \) fixed the indicator is constant in \( \bbeta \), and the exponent is that
of the conjugate normal update of @exr-opt-bayes-known-variance with \( \sigma^2=1 \),
response \( \bz \) and prior precision \( \bSigma^{-1} \); completing the square as in
the proof of @thm-opt-bayes-conjugate gives the stated normal law.

(c) With \( \bbeta \) fixed, the joint density factors over \( i \) into
\( \phi(z_i-\x_{(i)}\T\bbeta) \) times the indicator, which is the stated truncated normal
density up to its normalizing constant.
:::

The price is \( n \) extra draws per sweep; the gain is that every step is accepted and
nothing needs tuning. The same device works for the logit link by a different route — the
logistic distribution is a scale mixture of normals (Holmes and Held, 2006), the binomial
likelihood a mixture over Pólya–gamma variables (Polson, Scott and Windle, 2013) — and,
more intricately, for multinomial and ordinal models.

::: {#exm-prc-anes-mcmc}
[The 1996 survey, sampled]

Return to the logistic regression of @exm-bin-anes: the intended vote of
\( 944 \) respondents on party identification, age in decades,
education and income. Put an independent \( \Normal(0,10^2) \) prior on each of the
\( 5 \) coefficients — nearly flat on this scale — and run four
random-walk Metropolis chains of \( 8000 \) draws from starts
\( 8 \) posterior standard deviations away in every coordinate, with the Laplace
covariance as the proposal shape and scale \( 2.38/\sqrt5 \). The acceptance rate is
\( 0.298 \), and discarding the first half of each chain leaves
\( 16000 \) draws.

Panel (a) of [Figure 39.4.2](04-bayesian-glms-mcmc.html#fig-prc-mcmc) shows the four
chains for the party coefficient: two of them climb far above the mode before turning
back, and all four are indistinguishable once the first few hundred iterations are past.
Panel (b) applies the whole recipe — discard the first half, then split what is left — at
each run length, and plots \( \hat R-1 \) on a log scale. It starts at
\( 5.63 \) with
\( 200 \) draws per chain, is below \( 1.05 \)
by \( 1000 \) and below \( 1.02 \) by
\( 1400 \), and its own Monte Carlo variability keeps it brushing the
\( 1.01 \) line until about \( 6000 \); at the full eight thousand
the largest value over the five coefficients is \( 1.0066 \).
Effective sample size is the binding constraint: \( 1079 \), or
\( 0.067 \) effective draws per iteration.

The posterior mean of that coefficient is \( 1.2300 \) with posterior
standard deviation \( 0.0705 \), against the maximum likelihood
\( 1.2188 \) and standard error \( 0.0709 \):
with \( n=944 \) and a prior this weak the two analyses agree to the second decimal, as
@thm-prc-laplace and the Bernstein–von Mises theorem say they must. What the simulation
buys is the ability to ask for any posterior summary at all.

The probit version, sampled by @prp-prc-albert-chib with \( 4000 \)
draws per chain, gives posterior mean \( 0.6780 \) against the maximum
likelihood \( 0.6733 \); the ratio of the logistic coefficient to it
is \( 1.814 \), against the \( \pi/\sqrt3=1.8138 \)
of @prp-bin-latent-scale. Its effective sample size per iteration is
\( 0.092 \), half again as good and with no tuning, though each
iteration costs \( n \) truncated normal draws. Panel (c) compares the two
autocorrelation functions.
:::

::: {when-format="html"}
![**Figure 39.4.2.** Diagnostics for the samplers of @exm-prc-anes-mcmc. (a) The first
1200 draws of four Metropolis chains for the party-identification coefficient, from starts
eight posterior standard deviations away in every coordinate. (b) \( \hat R-1 \), on a log
scale, against draws generated per chain, recomputed from scratch at each length after
discarding the first half; the reference line is \( \hat R=1.01 \). (c) Autocorrelation
functions of the Metropolis and Albert–Chib
chains.](mcmc_diagnostics.svg){#fig-prc-mcmc width=100%}
:::

::: {when-format="pdf"}
![Diagnostics for the samplers of @exm-prc-anes-mcmc. (a) The first
1200 draws of four Metropolis chains for the party-identification coefficient, from starts
eight posterior standard deviations away in every coordinate. (b) \( \hat R-1 \), on a log
scale, against draws generated per chain, recomputed from scratch at each length after
discarding the first half; the reference line is \( \hat R=1.01 \). (c) Autocorrelation
functions of the Metropolis and Albert–Chib
chains.](mcmc_diagnostics.pdf){width=100%}
:::

```{.python .run #cell-mcmc-data}
import numpy as np
import statsmodels.api as sm
from scipy import special, stats
anes = sm.datasets.anes96.load_pandas().data
y = anes["vote"].to_numpy(float)                          # 1 if the respondent expected to vote Dole
X = np.column_stack([np.ones(len(y)), anes["PID"], anes["age"] / 10.0,
                     anes["educ"], anes["income"]])
n, p = X.shape
TAU = 10.0                                                # prior: beta ~ N(0, TAU^2 I)

mle = sm.GLM(y, X, family=sm.families.Binomial()).fit()
w = mle.fittedvalues * (1 - mle.fittedvalues)             # the working weights at the mode
V_prop = np.linalg.inv(X.T @ (w[:, None] * X) + np.eye(p) / TAU**2)   # Laplace covariance
CHOL = np.linalg.cholesky(V_prop)
print("maximum likelihood:", mle.params.round(4))
```

```{.python .run #cell-mcmc-metropolis}
def log_post(beta):
    eta = X @ beta
    return np.sum(y * eta - np.logaddexp(0.0, eta)) - np.sum(beta**2) / (2 * TAU**2)


def metropolis(n_draws, seed, start, step=2.38 / np.sqrt(5)):
    """Random-walk Metropolis with the Laplace covariance as the proposal shape."""
    rng = np.random.default_rng(seed)
    beta, lp = start.copy(), log_post(start)
    out, accepted = np.empty((n_draws, p)), 0
    for t in range(n_draws):
        prop = beta + step * (CHOL @ rng.standard_normal(p))
        lp_prop = log_post(prop)
        if np.log(rng.random()) < lp_prop - lp:           # the acceptance probability
            beta, lp, accepted = prop, lp_prop, accepted + 1
        out[t] = beta
    return out, accepted / n_draws


def split_rhat(draws):
    """Gelman-Rubin statistic on the two halves of every chain."""
    c, N = draws.shape
    s = draws.reshape(2 * c, N // 2)
    means, variances = s.mean(axis=1), s.var(axis=1, ddof=1)
    B, W = means.var(ddof=1), variances.mean()
    return np.sqrt(((N // 2 - 1) / (N // 2) * W + B) / W)


def ess(draws):
    """Effective sample size, with Geyer's initial positive sequence."""
    c, N = draws.shape
    centred = draws - draws.mean(axis=1, keepdims=True)
    var = draws.var(ddof=1)
    rho = np.array([np.mean([np.dot(x[:N - t], x[t:]) / N for x in centred]) / var
                    for t in range(1, min(N, 400))])
    pairs = rho[: 2 * (len(rho) // 2)].reshape(-1, 2).sum(axis=1)
    k = np.argmax(pairs < 0) if np.any(pairs < 0) else len(pairs)
    return c * N / (1 + 2 * pairs[:k].sum())


short, short_rates = zip(*[metropolis(1500, 3960 + c, mle.params) for c in range(2)])
short = np.array(short)[:, 750:, :]
print(f"two short chains: acceptance {np.mean(short_rates):.3f}, "
      f"party id mean {short[:, :, 1].mean():.4f}, R-hat {split_rhat(short[:, :, 1]):.4f}, "
      f"ESS {ess(short[:, :, 1]):.0f}")
```

```{.python .run #cell-mcmc-gibbs}
def albert_chib(n_draws, seed, start):
    """Gibbs for the probit model: latent normals, then a normal draw for beta."""
    rng = np.random.default_rng(seed)
    V = np.linalg.inv(X.T @ X + np.eye(p) / TAU**2)
    L = np.linalg.cholesky(V)
    beta, out = start.copy(), np.empty((n_draws, p))
    for t in range(n_draws):
        eta = X @ beta
        lower = stats.norm.cdf(-eta)                      # P(latent < 0)
        u = np.where(y > 0, lower + (1 - lower) * rng.random(n), lower * rng.random(n))
        z = eta + stats.norm.ppf(np.clip(u, 1e-12, 1 - 1e-12))
        mean = V @ (X.T @ z)
        beta = mean + L @ rng.standard_normal(p)
        out[t] = beta
    return out


gibbs_short = np.array([albert_chib(600, 3970 + c, np.zeros(p)) for c in range(2)])[:, 300:, :]
print(f"short Gibbs run: probit party id mean {gibbs_short[:, :, 1].mean():.4f}, "
      f"R-hat {split_rhat(gibbs_short[:, :, 1]):.4f}")
```

## Checking the fit from the posterior

::: {#def-prc-ppc}
[Posterior predictive distribution and check]

The **posterior predictive** distribution of a replicate data set \( \Y^{\text{rep}} \) is
\[
p(\y^{\text{rep}}\mid\y)=\int f(\y^{\text{rep}}\mid\bbeta)\,\pi(\bbeta\mid\y)\,d\bbeta ,
\]{#eq-prc-ppd}

which is sampled by drawing \( \bbeta \) from the posterior and then
\( \y^{\text{rep}} \) from the model. For a statistic \( T \) the **posterior predictive
\( p \)-value** is \( \Pr\{T(\Y^{\text{rep}})\ge T(\y)\mid\y\} \), or the
lower-tail version \( \Pr\{T(\Y^{\text{rep}})\le T(\y)\mid\y\} \) when small values
of \( T \) are the departure of interest, estimated by the proportion of replicates on
the relevant side of the observed value.
:::

The check is only as good as the statistic. One the model was fitted to match — the total
number of successes in a model with an intercept — has a \( p \)-value near \( 1/2 \) by
construction; a useful statistic is one the model does *not* control, such as a subgroup
total, a count of zeros or a maximum.

::: {.warning}
The data are used twice: once to form the posterior and again to compute \( T(\y) \). The
resulting \( p \)-value is conservative — its distribution under the model is concentrated
towards \( 1/2 \) rather than uniform — so a middling value is weak evidence of nothing
while an extreme one is strong. Bayarri and Berger (2000) quantify the conservatism; Meng
(1994) and Gelman, Meng and Stern (1996) give the general framework, including
discrepancies \( T(\y,\bbeta) \) that depend on the parameter too.
:::

::: {#exm-prc-ppc}
[Where the survey model fails]

The model of @exm-prc-anes-mcmc treats party identification as a single linear term. Take
\( T \) to be the number of Dole votes among the \( 108 \)
respondents who called themselves independent Democrats — a subgroup the model has no
parameter for. Observed: \( 7 \); across four thousand posterior
predictive replicates the mean is \( 14.4 \), and only
\( 0.034 \) of them are as low as seven
([Figure 39.4.3](04-bayesian-glms-mcmc.html#fig-prc-ppc), panel (a)). Given the
conservatism just described, that is a real signal.

It is not one a quadratic term would catch: adding \( \text{PID}^2 \) drops the deviance
by \( 0.019 \), while treating party identification as a seven-level
factor drops it by \( 12.68 \) on five degrees of freedom,
\( p=0.027 \) — the same finding from the likelihood side. Panel (b)
shows why: at average covariates the fitted curve passes well above the observed
proportion at party identification two and below it at four. The monotone logistic curve
is not flexible enough near the middle of this scale.
:::

::: {when-format="html"}
![**Figure 39.4.3.** (a) Posterior predictive distribution of the number of Dole votes
among the 108 independent Democrats, with the observed value marked. (b) Fitted
probability of a Dole vote against party identification at the average of the other
covariates, with a 95 per cent pointwise credible band and the raw category proportions
(unadjusted).](posterior_predictive.svg){#fig-prc-ppc width=100%}
:::

::: {when-format="pdf"}
![(a) Posterior predictive distribution of the number of Dole votes
among the 108 independent Democrats, with the observed value marked. (b) Fitted
probability of a Dole vote against party identification at the average of the other
covariates, with a 95 per cent pointwise credible band and the raw category proportions
(unadjusted).](posterior_predictive.pdf){width=100%}
:::

```{.python .run #cell-mcmc-ppc}
group = anes["PID"].to_numpy() == 2                       # the independent Democrats
observed = y[group].sum()
short4 = np.array([metropolis(2000, 3980 + c, mle.params)[0] for c in range(4)])[:, 1000:, :]
draws_short = short4.reshape(-1, p)
rng2 = np.random.default_rng(3990)
idx2 = rng2.choice(len(draws_short), 2000, replace=False)
rep_short = np.array([rng2.binomial(1, special.expit(X @ draws_short[i]))[group].sum()
                      for i in idx2])
print(f"short run: {int(observed)} of {int(group.sum())} against {rep_short.mean():.1f} "
      f"replicated; posterior predictive p-value {np.mean(rep_short <= observed):.3f}")
```

## Exercises

### A. Check your understanding

::: {#exr-prc-mh-symmetric}
[A1]

Show that for a symmetric proposal the acceptance probability @eq-prc-mh reduces to
\( \min\{1,\exp(q_n(\bbeta^{*})-q_n(\bbeta^{(t)}))\} \), and that an uphill move is always
accepted.
:::

::: {.solution}
Symmetry gives \( q(\bbeta^{(t)}\mid\bbeta^{*})=q(\bbeta^{*}\mid\bbeta^{(t)}) \), so those
factors cancel and \( \alpha=\min\{1,\pi(\bbeta^{*}\mid\y)/\pi(\bbeta^{(t)}\mid\y)\} \).
The normalizing constant cancels too, leaving the ratio of the unnormalized densities,
which is \( \exp\{q_n(\bbeta^{*})-q_n(\bbeta^{(t)})\} \). If
\( q_n(\bbeta^{*})\ge q_n(\bbeta^{(t)}) \) the ratio is at least one and \( \alpha=1 \).
:::

::: {#exr-prc-rhat-reading}
[A2]

A chain is run once, from a start far from the posterior mode, and \( \hat R \) computed
from its two halves is \( 1.00 \). Is that evidence of convergence? What does the
definition in @eq-prc-rhat require?
:::

::: {.solution}
No. \( \hat R \) near one says only that the pieces compared agree; a single chain stuck
in one region has two halves that agree perfectly. The definition compares \( 2c \)
sequences, and the value of the diagnostic comes from starting the \( c \) chains
*over-dispersed relative to the posterior*, so that disagreement is detectable. A single
chain also cannot reveal a second mode it never visited.
:::

### B. Practice

::: {#exr-prc-laplace-onedim}
[B1]

Take \( q_n(\beta)=n\,g(\beta) \) for a fixed smooth \( g \) with a unique interior
maximum at \( \beta^{*} \) and \( g''(\beta^{*})<0 \), and suppose \( g \) is bounded away
from its maximum outside any neighbourhood of \( \beta^{*} \). Verify (L1)–(L3) of
@thm-prc-laplace and write out the approximation explicitly.
:::

::: {#exr-prc-bic-prior}
[B2]

@cor-prc-bic says the prior enters only the \( O(1) \) remainder. Show that if the prior
is itself allowed to depend on \( n \) — say
\( \bbeta\sim\Normal_p(\bzero,n\bSigma_0) \) with \( \bSigma_0 \) fixed, so that the
prior is made vaguer as data accumulate — then the remainder grows like
\( p\log n \), and say what that does to the BIC.
:::

::: {.solution}
With \( \pi \) the \( \Normal_p(\bzero,n\bSigma_0) \) density,
\( -2\log\pi(\hbeta_n)=p\log(2\pi)+p\log n+\log\det\bSigma_0
+\hbeta_n\T\bSigma_0^{-1}\hbeta_n/n \), so the term \( -2\log\pi(\hbeta_n) \) in the
remainder grows like \( p\log n \) instead of staying bounded. The expansion then reads
\( -2\log m(\y)=-2\ell(\hbeta)+2p\log n+O(1) \): a prior whose spread grows with \( n \)
doubles the dimension penalty. This is Lindley's paradox in miniature — a prior that is
made vaguer as data accumulate penalizes extra parameters more and more heavily — and it is
why the BIC's derivation assumes a fixed prior.
:::

::: {#exr-prc-ess-thinning}
[B3]

A chain of \( 16\,000 \) draws has effective sample size \( 1000 \). A colleague thins it
by keeping every sixteenth draw, obtaining \( 1000 \) nearly independent draws, and says
nothing has been lost. Is that right? What is thinning good for?
:::

### C. Going deeper

::: {#exr-prc-gibbs-as-mh}
[C1]

Show that the Gibbs update of a block \( \boldsymbol{\uppsi}_k \) is a Metropolis–Hastings
step with proposal
\( q(\boldsymbol{\uppsi}^{*}\mid\boldsymbol{\uppsi})
=\pi(\boldsymbol{\uppsi}^{*}_k\mid\boldsymbol{\uppsi}_{-k},\y) \) and that its acceptance
probability is one. Deduce that @prp-prc-mh-invariance covers the Gibbs sampler, and
explain why a Gibbs sampler can nevertheless fail to be irreducible.
:::
