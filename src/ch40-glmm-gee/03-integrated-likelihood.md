# Integrated likelihood approximations

Everything the last section promised rests on the integral @eq-gmm-integrated, which has no
closed form except in the normal linear case. Fitting a GLMM is therefore a numerical
problem before it is a statistical one, and the choice of method is not innocent: three
methods in common use can differ in the second significant figure of \( \hat\bbeta \) and
in the first of \( \hat\tau \).

## Why it is hard, and what makes it easier

For the normal linear mixed model \( \Y\mid\bu \) and \( \bu \) are both normal, so
@eq-gmm-integrated is a Gaussian integral and gives @eq-mix-marginal in closed form. Nothing
of the sort happens otherwise: a Poisson likelihood in \( e^{u} \) against a normal density
in \( u \) is not a recognizable integral in any dimension. What *does* help is
independence. If the random effects are attached to clusters that are independent of each
other, the \( q \)-dimensional integral splits into \( m \) integrals of dimension
\( q_0 \) each, and with a random intercept \( q_0=1 \). A one-dimensional integral
repeated \( m \) times is tractable; a \( q \)-dimensional one with \( q \) in the
hundreds is not.

::: {#thm-gmm-likelihood}
[The integrated likelihood and its approximations]

In @def-gmm-glmm write
\( \psi_i(\bu_i)=\log f(\y_i\mid\bu_i;\bbeta,\phi)+\log p(\bu_i;\boldsymbol{\uptheta}) \).

::: {.enumerate options="label=(\alph*)"}
1. **Factorization.** If the random effects split into independent blocks
   \( \bu_1,\dots,\bu_m \) with \( \bu_i \) entering only the observations of cluster
   \( i \), then
   \[
   L(\bbeta,\boldsymbol{\uptheta},\phi)
   =\prod_{i=1}^{m}\int_{\Real^{q_0}}\exp\{\psi_i(\bu_i)\}\,d\bu_i ,
   \]{#eq-gmm-factorized}

   a product of \( m \) integrals of dimension \( q_0 \).

2. **Laplace approximation.** Fix \( i \) and take \( q_0=1 \). Let \( \psi_i \) be
   four times continuously differentiable with a unique interior maximizer
   \( \hat u_i \), and \( c_i=-\psi_i''(\hat u_i)>0 \). Then, formally,
   \[
   \begin{aligned}
   \int e^{\psi_i(u)}\,du
   &=\sqrt{\frac{2\pi}{c_i}}\;e^{\psi_i(\hat u_i)}\,
   \bigl\{1+\delta_i+\dots\bigr\},\\
   \delta_i&=\frac{\psi_i^{(4)}(\hat u_i)}{8c_i^{2}}
   +\frac{5\{\psi_i'''(\hat u_i)\}^{2}}{24c_i^{3}} .
   \end{aligned}
   \]{#eq-gmm-laplace}

   Suppose \( \psi_i=n_ik_i+s_i \), where \( k_i \) does not depend on \( n_i \) and has
   a nondegenerate maximum and \( s_i \) has derivatives bounded uniformly in \( n_i \).
   (In a GLMM \( n_ik_i \) is the conditional log-likelihood and \( s_i \) the log prior
   density of \( \bu_i \), which does not grow with \( n_i \).) Then both terms of
   \( \delta_i \) are \( O(n_i^{-1}) \): the Laplace approximation has *relative* error
   \( O(n_i^{-1}) \) and the error in the log-likelihood is
   \( O\bigl(\sum_i n_i^{-1}\bigr) \).

3. **Adaptive Gauss–Hermite quadrature.** Let \( (z_k,w_k)_{k=1}^{K} \) be the nodes and
   weights of \( K \)-point Gauss–Hermite quadrature, so that
   \( \int e^{-z^{2}}P(z)\,dz=\sum_k w_kP(z_k) \) for every polynomial \( P \) of degree at
   most \( 2K-1 \). With \( \hat\sigma_i=c_i^{-1/2} \), the substitution
   \( u=\hat u_i+\sqrt2\,\hat\sigma_i z \) gives the approximation
   \[
   \int e^{\psi_i(u)}\,du
   \approx\sqrt2\,\hat\sigma_i\sum_{k=1}^{K}w_k\,e^{z_k^{2}}\,
   e^{\psi_i(\hat u_i+\sqrt2\,\hat\sigma_i z_k)} .
   \]{#eq-gmm-aghq}

   It is exact whenever \( e^{\psi_i(\hat u_i+\sqrt2\hat\sigma_i z)+z^{2}} \) is a
   polynomial of degree at most \( 2K-1 \) in \( z \), and \( K=1 \) reproduces @eq-gmm-laplace
   without its correction terms.

4. **Penalized quasi-likelihood** avoids the integral altogether, by a penalized
   maximization that is cheaper than either of the above and carries a bias that does not
   vanish. It is stated and proved separately as @prp-gmm-pql.
:::

:::

::: {.proof}
(a) Under the stated independence the integrand factorizes and Fubini's theorem applies,
every factor being nonnegative.

(b) Substitute \( u=\hat u_i+s/\sqrt{c_i} \) and expand \( \psi_i \) about \( \hat u_i \),
using \( \psi_i'(\hat u_i)=0 \):
\[
\begin{aligned}
\psi_i(u)&=\psi_i(\hat u_i)-\tfrac12 s^{2}+a_3s^{3}+a_4s^{4}+\dots,\\
a_3&=\frac{\psi_i'''(\hat u_i)}{6c_i^{3/2}},\qquad
a_4=\frac{\psi_i^{(4)}(\hat u_i)}{24c_i^{2}} .
\end{aligned}
\]
Hence
\( \int e^{\psi_i}\,du=c_i^{-1/2}e^{\psi_i(\hat u_i)}\int
e^{-s^{2}/2}e^{a_3s^{3}+a_4s^{4}+\dots}\,ds \).
Expanding the second exponential as \( 1+(a_3s^{3}+a_4s^{4})+\tfrac12 a_3^{2}s^{6}+\dots \)
and using the standard normal moments \( \E S^{3}=0 \), \( \E S^{4}=3 \), \( \E S^{6}=15 \)
gives the factor \( \sqrt{2\pi}\,\{1+3a_4+\tfrac{15}{2}a_3^{2}+\dots\} \), which is @eq-gmm-laplace.
If \( \psi_i=n_ik_i+s_i \) then \( c_i=n_i\lvert k_i''\rvert+O(1) \),
\( \psi_i'''=n_ik_i'''+O(1) \) and \( \psi_i^{(4)}=n_ik_i^{(4)}+O(1) \), so
\( a_4=O(n_i^{-1}) \) and \( a_3^{2}=O(n_i^{-1}) \).

(c) The substitution is a change of variable, and
\( e^{\psi_i(u)}=e^{-z^{2}}\bigl\{e^{z^{2}}e^{\psi_i(\hat u_i+\sqrt2\hat\sigma_iz)}\bigr\} \),
so the quadrature rule applies to the term in braces. For \( K=1 \) the node is
\( z_1=0 \) and the weight is \( w_1=\sqrt\pi \), giving
\( \sqrt2\,\hat\sigma_i\sqrt\pi\,e^{\psi_i(\hat u_i)}=\sqrt{2\pi/c_i}\,e^{\psi_i(\hat u_i)} \).
:::

::: {.remark}
[What "formally" means in (b)]

The expansion in the proof interchanges an infinite series with an integral over the whole
line, which the argument does not justify, and the remainder can also be affected by the
behaviour of \( \psi_i \) far from its mode. The rigorous statement — that the expansion is
asymptotic in \( n_i \), with an explicit remainder, under conditions on the derivatives of
\( k_i \) and on the tails — is classical, and is given in the form used here by
Tierney and Kadane (1986) and Barndorff-Nielsen and Cox (1989, chapter 3).
[Chapter 39](../ch39-glms-in-practice-bayes/index.html) proves the companion
statement for the posterior and the marginal likelihood (@thm-prc-laplace). We take the
error order as given; @exm-gmm-compare checks it numerically.
:::

Two readings of (b). First, the accuracy of the Laplace approximation is governed by the
observations *per cluster*, not by the sample size: a study with
\( 10\,000 \) patients seen twice each is the hard case, and no number of
extra patients repairs it. Second, the relevant \( n_i \) is an *effective* one — a cluster
of twenty Bernoulli responses says much less about its own \( u_i \) than a cluster of
twenty gamma responses — which is why the same approximation behaves so differently in
@exm-gmm-compare and @exm-gmm-grunfeld.

Part (c) matters because plain Gauss–Hermite quadrature places its nodes where the *prior*
is large. When the data are informative the posterior of \( u_i \) sits elsewhere and is
much narrower, so most nodes fall where the integrand is negligible. The adaptive version
puts them at each cluster's own mode and curvature — the same \( \hat u_i \) and
\( \hat\sigma_i \) that @prp-gmm-prediction predicts the random effect with — at the cost
of one Newton iteration per cluster per likelihood evaluation.

```{.python .run #cell-approximations-quadrature}
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize, minimize_scalar

BETA = np.array([-0.5, 1.0])
TAU = 1.0

def expit(z):
    return 1.0 / (1.0 + np.exp(-z))

def cluster_mode(y, X, beta, tau):
    """Mode and curvature of log f(y_i | u) + log p(u), cluster by cluster."""
    eta0, u = X @ beta, np.zeros(y.shape[0])
    for _ in range(100):
        mu = expit(eta0 + u[:, None])
        step = ((y - mu).sum(1) - u / tau**2) / (-(mu * (1 - mu)).sum(1) - 1 / tau**2)
        u = u - step
        if np.max(np.abs(step)) < 1e-12:
            break
    mu = expit(eta0 + u[:, None])
    return u, 1.0 / np.sqrt((mu * (1 - mu)).sum(1) + 1 / tau**2)

def loglik(y, X, beta, tau, K, adapt=True):
    """Integrated log-likelihood by K-node Gauss-Hermite quadrature.

    Adaptive: centre and scale the nodes at each cluster's own mode and curvature.
    K = 1 with adapt=True is exactly the Laplace approximation."""
    node, w = hermgauss(K)
    if adapt:
        centre, scale = cluster_mode(y, X, beta, tau)
    else:
        centre, scale = np.zeros(y.shape[0]), np.full(y.shape[0], tau)
    pts = centre[:, None] + np.sqrt(2) * scale[:, None] * node            # (m, K)
    eta = (X @ beta)[:, :, None] + pts[:, None, :]                        # (m, n, K)
    cond = (y[:, :, None] * eta - np.logaddexp(0.0, eta)).sum(1)          # (m, K)
    prior = -0.5 * (pts / tau) ** 2 - np.log(tau) - 0.5 * np.log(2 * np.pi)
    terms = np.log(w) + node**2 + cond + prior
    return float(np.sum(np.log(np.sqrt(2) * scale) + np.logaddexp.reduce(terms, axis=1)))

def fit_ml(y, X, K, start=None):
    p = X.shape[2]
    obj = lambda par: -loglik(y, X, par[:p], np.exp(par[p]), K)
    x0 = np.r_[np.zeros(p), 0.0] if start is None else start
    r = minimize(obj, x0, method="Nelder-Mead",
                 options=dict(maxiter=8000, maxfev=8000, xatol=1e-8, fatol=1e-10))
    return r.x[:p], float(np.exp(r.x[p])), -r.fun

M_demo, N_demo = 60, 5
rng_demo = np.random.default_rng(40_007)
x_demo = rng_demo.normal(size=(M_demo, N_demo))
X_demo = np.stack([np.ones((M_demo, N_demo)), x_demo], axis=2)
u_demo = rng_demo.normal(0, TAU, M_demo)
y_demo = (rng_demo.uniform(size=(M_demo, N_demo))
          < expit(X_demo @ BETA + u_demo[:, None])).astype(float)

for K in [1, 5, 15]:
    b, t_hat, ll = fit_ml(y_demo, X_demo, K)
    print(f"AGHQ K = {K:2d}:  beta = {b.round(4)}   tau = {t_hat:.4f}   loglik = {ll:.4f}")
```

## Penalized quasi-likelihood

The oldest practical method avoids the integral by treating the random effects as
parameters with a penalty.

::: {#prp-gmm-pql}
[Penalized quasi-likelihood]

In @def-gmm-glmm let \( \ell(\bbeta,\bu) \) be the conditional log-likelihood of the data
given \( \bu \), and consider the **penalized quasi-likelihood**
\[
\ell^{\text{pen}}(\bbeta,\bu)=\ell(\bbeta,\bu)-\tfrac12\bu\T\G^{-1}\bu .
\]{#eq-gmm-pql-objective}

::: {.enumerate options="label=(\alph*)"}
1. For fixed \( \boldsymbol{\uptheta} \), the stationary equations of @eq-gmm-pql-objective
   are
   \[
   \X\T\W\bD^{-1}(\y-\bmu)=\bzero,\qquad
   \Z\T\W\bD^{-1}(\y-\bmu)=\phi\,\G^{-1}\bu ,
   \]
   with \( \W \) and \( \bD \) the working weights and derivatives of @thm-glm-score
   evaluated at \( \boldsymbol{\upeta}=\X\bbeta+\Z\bu \). Equivalently, with the working
   response \( \bz=\boldsymbol{\upeta}+\bD^{-1}(\y-\bmu) \) of @eq-glm-working-response, they
   are Henderson's mixed-model equations (@thm-mix-henderson) for the linear mixed model
   \( \bz\mid\bu\sim\Normal_n(\X\bbeta+\Z\bu,\ \phi\W^{-1}) \),
   \( \bu\sim\Normal_q(\bzero,\G) \).

2. Iterating — update \( (\bbeta,\bu) \) by those equations, then update
   \( \boldsymbol{\uptheta} \) by maximum likelihood or REML in the working linear mixed
   model, and repeat — is the Breslow–Clayton PQL algorithm. It is therefore the Laplace
   approximation @eq-gmm-laplace with the further simplification that \( \W \) is treated
   as fixed when \( \boldsymbol{\uptheta} \) is updated, so that \( c_i \) contributes no
   derivative terms.

3. The resulting \( \hat\bbeta \) and \( \hat{\boldsymbol{\uptheta}} \) are biased towards
   zero, and the bias does **not** vanish as \( m\to\infty \) with cluster sizes fixed: PQL
   is inconsistent. It is worst for binary responses with small \( n_i \) and large
   \( \tau \), and negligible for responses that carry much information per observation.
:::

:::

::: {.proof}
(a) Differentiating \( \ell \) with respect to \( \bbeta \) and \( \bu \) reproduces @eq-glm-score
with model matrix \( [\X,\Z] \), because \( \bu \) enters the linear predictor
exactly as \( \bbeta \) does; the penalty contributes \( -\G^{-1}\bu \) to the second block
and nothing to the first. Multiplying by \( \phi \) gives the displayed pair. For the second
form, \( \W\bD^{-1}(\y-\bmu)=\W(\bz-\boldsymbol{\upeta}) \) by the definition of \( \bz \),
so the equations read
\( \X\T\W(\bz-\X\bbeta-\Z\bu)=\bzero \) and
\( \Z\T\W(\bz-\X\bbeta-\Z\bu)=\phi\G^{-1}\bu \), which are @eq-mix-henderson with
\( \R=\phi\W^{-1} \).

(b) is the observation that maximizing @eq-gmm-pql-objective over \( \bu \) locates
\( \hat\bu \) at the mode of the integrand of @eq-gmm-integrated, which is the first
ingredient of @eq-gmm-laplace; the remaining ingredient, \( \log\lvert c_i\rvert \), is what
the \( \boldsymbol{\uptheta} \) update drops.

(c) We do not prove this. The asymptotic bias expansions, and corrections built from them,
are given by Breslow and Lin (1995) and Lin and Breslow (1996); @exm-gmm-compare measures the
bias by simulation.
:::

::: {.remark}
[PQL is empirical Bayes]

Since \( \exp\{\ell^{\text{pen}}\} \) is proportional to the joint density of \( \y \) and
\( \bu \), the PQL estimate of \( \bu \) is the posterior mode of @prp-gmm-prediction with
\( \bbeta \) and \( \boldsymbol{\uptheta} \) replaced by estimates: an empirical Bayes
predictor, and nothing about it is wrong. What is wrong is using the same objective to
estimate \( \boldsymbol{\uptheta} \), because it has been maximized over \( \bu \)
rather than integrated over it, and the two differ by exactly the \( \log c_i \) term that
carries the information about \( \tau \).
:::

## Monte Carlo alternatives

When \( q_0 \) exceeds two or three, a product grid needs \( K^{q_0} \) nodes and stops
being practical. **Monte Carlo EM** treats \( \bu \) as missing data: the E step's
expectation over the posterior of \( \bu \) is replaced by an average over Markov chain
Monte Carlo draws, and the M step is an ordinary weighted GLM fit. EM's monotonicity
survives only up to Monte Carlo error, so the number of draws must grow with the iterations;
[Chapter 41](../ch41-missing-data/index.html) develops the algorithm and proves that the
exact version never decreases the likelihood. **Simulated maximum likelihood** estimates each
factor of @eq-gmm-factorized by importance sampling with the Laplace approximation as
proposal — @eq-gmm-aghq with random nodes. **Fully Bayesian inference** samples
\( \bbeta \), \( \bu \) and \( \boldsymbol{\uptheta} \) together, as
[Section 32.6](../ch32-linear-mixed-models/06-bayes.html) did for the linear mixed
model (@prp-mix-bayes); without conjugacy the full conditionals need Metropolis–Hastings
steps or data augmentation, which
[Chapter 39](../ch39-glms-in-practice-bayes/index.html) builds (@def-prc-mcmc).

## How different are the answers?

::: {#exm-gmm-compare}
[Three methods on clustered binary data]

Simulate \( 150 \) clusters of \( 5 \) binary
observations from a random-intercept logistic model with
\( \bbeta=(-0.5,1.0)\T \), a covariate varying within the cluster, and
\( \tau=1 \). Over \( 200 \) replications the average estimates are

| | \( \hat\beta_1 \) (truth \( 1.0 \)) | \( \hat\tau \) (truth \( 1.0 \)) |
|:---|---:|---:|
| PQL | 0.8970 | 0.8168 |
| Laplace | 1.0022 | 0.9353 |
| adaptive Gauss–Hermite, \( K=9 \) | 1.0087 | 0.9909 |

PQL attenuates the slope by \( 10.3 \) per cent and the
standard deviation by \( 18.3 \) per cent. The Laplace
approximation is nearly unbiased for the slope but still shrinks \( \hat\tau \) by about
six per cent, which is consistent with the \( O(n_i^{-1}) \) error of
@thm-gmm-likelihood(b) at \( n_i=5 \); the direction and size of the induced bias are
measured here, not derived, since part (b) bounds the error in the likelihood and says
nothing about where its maximizer moves. Nine adaptive nodes leave a bias of about one per
cent, which is the ordinary finite-sample bias of maximum likelihood rather than a numerical
error.

Panel (c) of [Figure 40.3.1](#fig-gmm-approximations) separates the numerical question from
the statistical one. At a fixed parameter value, the error in the integrated log-likelihood
with \( K=5 \) is \( 3.65136\times 10^{-3} \) for the adaptive rule and
\( 1.01404\times 10^{-1} \) for the plain one — a factor of about
\( 28 \) — and the adaptive rule is down to
\( 3.759\times 10^{-10} \) by \( K=13 \), while the plain one is still
above \( 10^{-6} \), at \( 2.460\times 10^{-6} \), when
\( K=19 \).
:::

::: {when-format="html"}
![**Figure 40.3.1.** Approximations to the integrated likelihood, on clustered binary data
with \( m=150 \), \( n_i=5 \) and \( \tau=1 \). (a), (b) Sampling distributions of
\( \hat\beta_1 \) and \( \hat\tau \) over \( 200 \) replications under three methods;
the horizontal line is the truth. (c) Error in the integrated log-likelihood against the
number of quadrature nodes, adaptive and
plain.](approximations.svg){#fig-gmm-approximations width=100%}
:::

::: {when-format="pdf"}
![Approximations to the integrated likelihood, on clustered binary data
with \( m=150 \), \( n_i=5 \) and \( \tau=1 \). (a), (b) Sampling distributions of
\( \hat\beta_1 \) and \( \hat\tau \) over \( 200 \) replications under three methods;
the horizontal line is the truth. (c) Error in the integrated log-likelihood against the
number of quadrature nodes, adaptive and plain.](approximations.pdf){width=100%}
:::

```{.python .run #cell-approximations-pql}
def weighted_lmm(z, X, W):
    """REML fit of a weighted random-intercept linear mixed model: z_ij = x_ij'b + u_i + e_ij
    with Var(e_ij) = 1/W_ij and Var(u_i) = t2. Returns b, t2 and the predicted u."""
    def solve(t2):
        s, WX = W.sum(1), W[..., None] * X
        den = 1 + t2 * s
        XtVX = (np.einsum("mnp,mnq->pq", WX, X)
                - t2 * np.einsum("mp,mq->pq", WX.sum(1), WX.sum(1) / den[:, None]))
        XtVz = (np.einsum("mnp,mn->p", WX, z)
                - t2 * np.einsum("mp,m->p", WX.sum(1), (W * z).sum(1) / den))
        return np.linalg.solve(XtVX, XtVz), XtVX, den

    def neg_reml(log_t2):
        t2 = np.exp(log_t2)
        b, XtVX, den = solve(t2)
        r = z - X @ b
        q = np.sum(W * r * r) - t2 * np.sum((W * r).sum(1) ** 2 / den)
        ld = np.sum(np.log(den)) - np.sum(np.log(W))
        return 0.5 * (ld + np.linalg.slogdet(XtVX)[1] + q)

    opt = minimize_scalar(neg_reml, bounds=(-12.0, 6.0), method="bounded",
                          options=dict(xatol=1e-10))
    t2 = float(np.exp(opt.x))
    b, _, den = solve(t2)
    r = z - X @ b
    return b, t2, t2 * (W * r).sum(1) / den

def fit_pql(y, X, maxit=300):
    """Breslow-Clayton penalized quasi-likelihood: a linear mixed model on the working
    response, refitted until it stops changing."""
    beta, u = np.zeros(X.shape[2]), np.zeros(y.shape[0])
    for _ in range(maxit):
        eta = X @ beta + u[:, None]
        mu = expit(eta)
        v = np.clip(mu * (1 - mu), 1e-8, None)          # working weights
        z = eta + (y - mu) / v                          # working response
        b_new, t2, u_new = weighted_lmm(z, X, v)
        done = max(np.max(np.abs(b_new - beta)), np.max(np.abs(u_new - u))) < 1e-10
        beta, u = b_new, u_new
        if done:
            break
    return beta, float(np.sqrt(t2)), u

b_pql, t_pql, _ = fit_pql(y_demo, X_demo)
print(f"PQL        :  beta = {b_pql.round(4)}   tau = {t_pql:.4f}")
```

The contrast with @exm-gmm-grunfeld is instructive. That fit has \( n_i=20 \) continuous
responses per firm, and refitting it with one adaptive node instead of nine changes no
parameter by more than \( 5\times10^{-4} \) and the maximized log-likelihood, which is
\( -937.100 \), by less than \( 0.01 \). What separates the two is not
the software but the information each cluster carries about its own random effect, exactly
as @thm-gmm-likelihood(b) says.

::: {.idea}
Report which approximation was used. "A random-intercept logistic model was fitted" is not a
complete description when PQL and adaptive quadrature can differ by a tenth of the
coefficient. With binary responses in small clusters, use quadrature with enough nodes that
the answer stops moving, and say how many.
:::

## Exercises

### A. Check your understanding

::: {#exr-gmm-nodes}
[A1]

A GLMM with two correlated random effects per cluster is fitted by non-adaptive
Gauss–Hermite quadrature on a product grid. How many likelihood evaluations per cluster does
\( K=20 \) cost? What if there are four random effects per cluster? What does this say about
crossed random effects, where \( q_0 \) is effectively the whole vector?
:::

::: {#exr-gmm-which-method}
[A2]

For each of the following, say whether PQL is likely to be adequate and why: (a) counts with
means in the hundreds, ten per cluster; (b) binary responses, three per cluster, with a
large estimated \( \tau \); (c) grouped binomial responses with \( N_i=200 \) trials, one per
cluster.
:::

### B. Practice

::: {#exr-gmm-one-node}
[B1]

Verify @thm-gmm-likelihood(c) for \( K=1 \) directly: show that the one-point Gauss–Hermite
rule has node \( z_1=0 \) and weight \( w_1=\sqrt\pi \), and that @eq-gmm-aghq
then reduces to @eq-gmm-laplace without its correction terms. Then explain why
the *non*-adaptive rule with \( K=1 \) gives something quite different, and say what.
:::

::: {.solution}
The one-point rule must integrate \( e^{-z^{2}}P(z) \) exactly for \( P \) of degree
\( \le1 \). Taking \( P=1 \) gives \( w_1=\int e^{-z^{2}}dz=\sqrt\pi \); taking \( P=z \)
gives \( w_1z_1=\int ze^{-z^{2}}dz=0 \), so \( z_1=0 \). Substituting in @eq-gmm-aghq
gives \( \sqrt2\hat\sigma_i\sqrt\pi e^{\psi_i(\hat u_i)} \), which is
\( \sqrt{2\pi/c_i}\,e^{\psi_i(\hat u_i)} \). The non-adaptive rule puts its single node at
\( u=0 \) with scale \( \tau \), giving \( \sqrt{2\pi}\,\tau\,e^{\psi_i(0)} \): the
integrand evaluated at the *prior* mode with the prior's width, which is the right answer
only when the data say nothing.
:::

::: {#exr-gmm-effective-n}
[B2]

For a cluster of \( n \) Bernoulli responses with common success probability \( \pi \) and
a random intercept, show that \( c_i=n\pi(1-\pi)+\tau^{-2} \). Compare with the gamma case
of @exm-gmm-grunfeld, where \( c_i=n/\phi+\tau^{-2} \), and use @thm-gmm-likelihood(b) to
explain the difference between the two fits in terms of an effective cluster size.
:::

::: {.solution}
For the logit, \( \psi_i''(u)=-\sum_j\mu_{ij}(1-\mu_{ij})-\tau^{-2} \), which at a common
\( \pi \) is \( -n\pi(1-\pi)-\tau^{-2} \). For the gamma with a log link,
\( \psi_i''(u)=-\nu\sum_jy_{ij}e^{-\eta_{ij}}-\tau^{-2} \), which at the mode is about
\( -n/\phi-\tau^{-2} \) since \( \nu=1/\phi \) and \( y_{ij}e^{-\eta_{ij}}\approx1 \). The
Laplace error is governed by \( c_i\tau^{2} \), so the effective cluster size is
\( n\pi(1-\pi)\tau^{2} \) against \( n\tau^{2}/\phi \). With \( n=5 \), \( \pi=1/2 \),
\( \tau=1 \) the first is \( 1.25 \); with \( n=20 \), \( \phi=0.0930 \), \( \tau=0.52 \) the
second is about \( 58 \). The two fits are in completely different regimes.
:::
