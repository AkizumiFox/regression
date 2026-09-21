# Working correlation

@thm-gmm-gee says the working correlation may be wrong without costing consistency. It does
not say the choice is free. A wrong one costs *efficiency*, changes the small-sample
behaviour of the sandwich, and — because it changes the estimate — can be turned into a
diagnostic for the mean model.

## The catalogue

The structures in use are those of @def-cls-covariance, written as correlation matrices and
estimated by moments rather than by likelihood. For a cluster with \( n_i \) observations at
occasions \( j=1,\dots,n_i \), and Pearson residuals
\( r_{ij}=(y_{ij}-\hat\mu_{ij})/\sqrt{V(\hat\mu_{ij})} \) with
\( \hat\phi=\sum_{ij}r_{ij}^{2}/(n-p) \):

| Structure | \( (\R_i)_{jk} \), \( j\ne k \) | Parameters |
|:---|:---|---:|
| independence | \( 0 \) | \( 0 \) |
| exchangeable | \( \alpha \) | \( 1 \) |
| first-order autoregressive | \( \alpha^{\lvert j-k\rvert} \) | \( 1 \) |
| \( q \)-dependent | \( \alpha_{\lvert j-k\rvert} \) if \( \lvert j-k\rvert\le q \), else \( 0 \) | \( q \) |
| unstructured | \( \alpha_{jk} \) | \( n(n-1)/2 \) |

The parameters are estimated by matching averages of residual products. Writing
\( N_{\text{p}}=\sum_i\binom{n_i}{2} \) for the within-cluster pairs and
\( N_{\text{a}}=\sum_i(n_i-1) \) for the adjacent ones,
\[
\hat\alpha_{\text{exch}}
=\frac{1}{\hat\phi\,(N_{\text{p}}-p)}\sum_i\sum_{j<k}r_{ij}r_{ik},
\qquad
\hat\alpha_{\text{ar}}
=\frac{1}{\hat\phi\,(N_{\text{a}}-p)}\sum_i\sum_{j<n_i}r_{ij}r_{i,j+1},
\]
and \( \hat\alpha_{jk}=\sum_ir_{ij}r_{ik}/(\hat\phi\,m) \), rescaled to have unit
diagonal, for the unstructured case; the \( q \)-dependent estimates are lag-by-lag
averages of the same kind.

The design narrows the list at once. The unstructured matrix needs a common set of occasions
across clusters and costs \( n(n-1)/2 \) parameters; the autoregressive structure needs an
ordering, so it is for longitudinal data and not for pupils in a school; the exchangeable
structure is the one a random intercept induces (@prp-gmm-moments), and is the default for
clustered data with no natural order. And every estimated \( \R_i \) must be a valid
correlation matrix, which the moment estimates above do not guarantee — @exm-gmm-modechoice
showed the exchangeable estimate landing almost exactly on the boundary.

```{.python .run #cell-working-correlation-working}
import numpy as np

KINDS = ["independence", "exchangeable", "ar1", "unstructured"]

def correlation(kind, r, p, phi):
    """Moment estimates of the working correlation from Pearson residuals r (m by n)."""
    m, n = r.shape
    if kind == "independence":
        return np.eye(n)
    if kind == "exchangeable":
        off = (np.sum(r.sum(1) ** 2) - np.sum(r**2)) / 2
        a = off / (phi * (m * n * (n - 1) / 2 - p))
        return np.eye(n) + a * (1 - np.eye(n))
    if kind == "ar1":
        a = np.sum(r[:, :-1] * r[:, 1:]) / (phi * (m * (n - 1) - p))
        return a ** np.abs(np.arange(n)[:, None] - np.arange(n))
    if kind == "unstructured":
        R = (r.T @ r) / (phi * m)
        return R / np.sqrt(np.outer(np.diag(R), np.diag(R)))
    raise ValueError(kind)

def gee(y, X, kind, family="poisson", maxit=100):
    m, n, p = X.shape
    beta = np.zeros(p)
    beta[0] = np.log(max(y.mean(), 1e-3)) if family != "binomial" else 0.0
    for _ in range(maxit):
        eta = X @ beta
        mu = np.exp(eta)
        dmu, V = (mu, mu) if family == "poisson" else (mu, mu**2)
        r = (y - mu) / np.sqrt(V)
        phi = np.sum(r**2) / (m * n - p)
        R = correlation(kind, r, p, phi)
        Rinv = np.linalg.pinv(R)
        SX = (dmu / np.sqrt(V))[..., None] * X
        A = np.einsum("mjp,jk,mkq->pq", SX, Rinv, SX)
        step = np.linalg.solve(A, np.einsum("mjp,jk,mk->p", SX, Rinv, r))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-10:
            break
    u = np.einsum("mjp,jk,mk->mp", SX, Rinv, r)
    Ainv = np.linalg.inv(A)
    return beta, Ainv @ (u.T @ u) @ Ainv, phi, R
```

## Efficiency, and where it comes from

::: {#prp-gmm-working}
[Choosing the working correlation]

In the setting of @thm-gmm-gee:

::: {.enumerate options="label=(\alph*)"}
1. **Efficiency.** The asymptotic covariance \( \A^{-1}\B\A^{-1} \) of @eq-gmm-gee-asymptotics
   is minimized, in the nonnegative-definite ordering and over all
   unbiased estimating functions linear in \( \Y_i-\bmu_i \), by
   \( \V_i=\bSigma_i^{0} \) (@exr-gmm-efficiency-bound). Any other working covariance costs
   efficiency, but the cost depends on the covariate: it is substantial for a covariate that
   varies *within* clusters and negligible for one that is constant within a cluster, for
   which the independence choice is in fact hard to beat.

2. **Small samples.** The sandwich @eq-gmm-sandwich is biased downwards by a term of order
   \( p/m \), so Wald tests built on it reject too often when \( m \) is small; the
   bias-corrected forms of [Section 40.4](04-gee.html) and a \( t_{m-p} \) reference repair
   most of it. The bias grows with the number of parameters in \( \R_i \), because
   \( \hat{\boldsymbol{\upalpha}} \) is then noisier and condition (E5) of @thm-gmm-gee
   is stretched.

3. **QIC.** Since there is no likelihood, AIC is unavailable. Pan's **quasi-likelihood
   information criterion** replaces the log-likelihood by the quasi-likelihood
   \( Q(\bmu;\y)=\sum_{ij}\int_{y_{ij}}^{\mu_{ij}}(y_{ij}-t)/\{\phi V(t)\}\,dt \) of @def-ql-quasi
   computed under the *independence* working correlation at the fitted \( \hat{\bmu} \), with
   \( \phi \) the Pearson estimate of @prp-ql-dispersion held fixed across the models
   compared, and replaces the parameter count by a sandwich-based trace:
   \[
   \text{QIC}=-2Q(\hat{\bmu};\y)
   +2\operatorname{tr}\bigl(\hat{\boldsymbol{\Omega}}_{\text{I}}\hat{\V}_{\text{R}}\bigr),
   \]{#eq-gmm-qic}

   where \( \hat{\boldsymbol{\Omega}}_{\text{I}}=\sum_i\hat\bD_i\T(\hat{\A}^{V}_i)^{-1}\hat\bD_i/\phi \)
   is the independence-model information and \( \hat{\V}_{\text{R}} \) is @eq-gmm-sandwich.
   If the working model is correct the trace is \( p \) and @eq-gmm-qic
   is the Akaike form @def-sel-aic-bic. It is a sound criterion for choosing the
   *mean model*, and a weakly justified one for choosing \( \R_i \).
:::

:::

::: {.proof}
(a) is @exr-gmm-efficiency-bound; the covariate-dependence is illustrated in @exm-gmm-efficiency
and explained after it. (b) is the cluster version of the bias
calculation in @exr-cls-crve-bias, applied to \( \hat\bD_i\T\hat\V_i^{-1}\br_i \) in place of
\( \X_g\T\he_g \). (c) is a definition; the claim that the trace equals \( p \) when the
working model is correct is @thm-gmm-gee(d), which gives
\( \hat{\V}_{\text{R}}\to\A^{-1} \), together with
\( \boldsymbol{\Omega}_{\text{I}}=\A \) when \( \R_i=\I \) is correct.
:::

::: {.remark}
[Why QIC is weak as a chooser of correlation]

@eq-gmm-qic measures fit by a quasi-likelihood that ignores the correlation altogether, and
penalizes by a trace that is smallest when the sandwich is closest to the model-based
covariance. Changing \( \R_i \) therefore moves QIC through the penalty almost alone, and
rewards whichever structure makes the sandwich small — not whichever is right. Pan (2001)
proposed it primarily for the mean model, and Hin and Wang (2009) give a variant (CIC) that
drops the fit term for correlation selection. Neither has the standing AIC has for
likelihoods.
:::

::: {#exm-gmm-efficiency}
[How much efficiency is at stake]

Simulate \( 60 \) clusters of \( 6 \) Poisson counts with a log
link, one covariate varying within the cluster and one constant within it, and generate the
dependence in two ways: a shared random intercept, which makes the true correlation
exchangeable, and a latent Gaussian AR(1) process with \( \rho=0.8 \). Fit by GEE with each
of the four working correlations and record the sampling variance of each coefficient over
\( 400 \) replications.

For the **within-cluster** covariate under exchangeable truth, the efficiency of the
exchangeable working correlation relative to independence is
\( 1.789 \); AR(1) and unstructured reach
\( 1.583 \) and \( 1.547 \). Under AR(1)
truth the figures are \( 1.282 \) for AR(1),
\( 1.183 \) for exchangeable and
\( 1.160 \) for unstructured. Matching the truth helps, mismatching it
still helps, and the gain is worth having without being dramatic.

For the **cluster-level** covariate the ranking reverses: under exchangeable truth the
relative efficiencies are \( 0.889 \),
\( 0.938 \) and \( 0.834 \), all below
one, so independence is best of the four. The reason is @exr-gmm-two-clusters: a covariate
constant within a cluster is estimated from cluster means, on which an exchangeable
reweighting has no effect at all, so all a nontrivial \( \R_i \) contributes is the noise
in \( \hat{\boldsymbol{\upalpha}} \).
:::

::: {when-format="html"}
![**Figure 40.5.1.** Efficiency of four working correlations relative to the independence
choice, for a within-cluster and a cluster-level covariate, over \( 400 \) replications with
\( m=60 \) clusters of \( n_i=6 \) Poisson counts. (a) The true dependence is exchangeable
(a shared random intercept). (b) The true dependence is a latent AR(1) process with
\( \rho=0.8 \). Bars above one mean more efficient than
independence.](efficiency.svg){#fig-gmm-efficiency width=100%}
:::

::: {when-format="pdf"}
![Efficiency of four working correlations relative to the independence
choice, for a within-cluster and a cluster-level covariate, over \( 400 \) replications with
\( m=60 \) clusters of \( n_i=6 \) Poisson counts. (a) The true dependence is exchangeable
(a shared random intercept). (b) The true dependence is a latent AR(1) process with
\( \rho=0.8 \). Bars above one mean more efficient than
independence.](efficiency.pdf){width=100%}
:::

## A diagnostic hiding in the efficiency

By @thm-gmm-gee(a), every working correlation estimates the *same* \( \bbeta^{0} \),
provided the mean model @eq-gmm-gee-mean is right. Two estimates of the same quantity
differing by more than their sampling error therefore accuse the mean model, not the
correlation. This is the Hausman idea, and it has real diagnostic value because the
independence and exchangeable estimates weight within- and between-cluster variation
differently: independence is dominated by between-cluster comparisons when clusters are
large, and a nontrivial \( \R_i \) shifts weight towards within-cluster ones.

::: {#exm-gmm-qic}
[Grunfeld's panel by GEE]

Fit the marginal version of @exm-gmm-grunfeld — gamma variance function
\( V(\mu)=\mu^{2} \), log link, investment on log value and log capital — to the
\( 11 \) firms by GEE. The AR(1) working correlation estimates a lag-one
correlation of \( 0.894 \), as one expects of annual investment figures. Yet
QIC, with the dispersion held at the independence estimate of the larger model, ranks the
four structures

| Structure | QIC, value and capital | QIC, value only |
|:---|---:|---:|
| independence | 8450.1 | 8657.6 |
| exchangeable | 8524.2 | 8647.4 |
| AR(1) | 8574.1 | 8657.8 |
| unstructured | 8754.7 | 8704.1 |

and prefers independence. As a comparison of *mean* models the table is informative: at three
of the four structures the model with capital wins, by between about \( 80 \) and
\( 210 \) QIC units, and that is the comparison QIC was designed for. Under the
unstructured correlation the ranking reverses, which with \( 11 \) clusters and a
\( 20\times20 \) matrix to estimate is a comment on the structure rather than on the model.
As a comparison of structures the column mostly reports that the sandwich is smallest under
independence — the remark above in numerical form.

The estimates themselves are the more interesting output. The independence fit gives
\( \hat\beta_{\text{value}}=0.8421 \) and
\( \hat\beta_{\text{capital}}=0.3386 \); the exchangeable fit gives
\( 0.7214 \) and \( 0.2113 \). The gap in
the first coefficient is more than its sandwich standard error of
\( 0.1030 \), and the two gaps point the same way. By @thm-gmm-gee(a)
that accuses the mean model, and two standard explanations fit equally well: firm effects
correlated with the covariates, as @exm-gmm-grunfeld found, or feedback from investment to
the following years' market value, which violates the full-covariate conditional mean
condition and makes the exchangeable fit — not the independence one — the inconsistent
member of the pair. Either way the remedy is to separate the two sources of variation.
Adding the firm means of log value and log capital, as @exr-gmm-within-between asks, makes
the independence and exchangeable estimates agree to machine precision, at within-firm
slopes \( 0.6142 \) and
\( 0.2298 \), both well below the independence figures above.
:::

```{.python .run #cell-working-correlation-qic}
import statsmodels.api as sm

d = sm.datasets.grunfeld.load_pandas().data.sort_values(["firm", "year"])
m, n = d["firm"].nunique(), 20
y = d["invest"].to_numpy().reshape(m, n)
lv = np.log(d["value"].to_numpy()).reshape(m, n)
lc = np.log(d["capital"].to_numpy()).reshape(m, n)

def qic(y, X, kind, phi):
    """Pan's QIC: the quasi-likelihood of the independence model at the fitted mean, plus a
    penalty built from the sandwich. For V(mu) = mu^2 the quasi-likelihood is -y/mu - log mu.
    The dispersion phi is held fixed across the models compared, as in Mallows's Cp."""
    beta, robust, _, _ = gee(y, X, kind, family="gamma")
    mu = np.exp(X @ beta)
    Q = np.sum(-y / mu - np.log(mu)) / phi
    omega = np.einsum("mjp,mjq->pq", X, X) / phi     # independence information, log link
    return -2 * Q + 2 * np.trace(omega @ robust)


X_full = np.stack([np.ones((m, n)), lv, lc], axis=2)
X_small = np.stack([np.ones((m, n)), lv], axis=2)
phi_ref = gee(y, X_full, "independence", family="gamma")[2]     # from the largest model
for name, XX in [("value + capital", X_full), ("value only", X_small)]:
    for kind in KINDS:
        print(f"{name:16s} {kind:14s} QIC {qic(y, XX, kind, phi_ref):10.2f}")
```

## Practical advice

- **Count the clusters first.** Everything here that uses a sandwich needs \( m \) large.
  Below about forty clusters use a bias correction and a \( t \) reference; below about
  fifteen, prefer a likelihood-based model with a covariance structure chosen by
  REML (@prp-cls-selection) and accept the assumptions.
- **Match the structure to the design, not to the data.** Exchangeable for unordered
  clusters, autoregressive for equally spaced longitudinal data, unstructured only when the
  occasions are common to all clusters and \( m \) is large.
- **Use independence when the covariates are time-varying, or the covariate of interest is
  cluster-level.** In the first case a non-diagonal \( \R_i \) needs the full-covariate
  conditional mean condition of [Section 40.4](04-gee.html), which feedback destroys; in the
  second it loses nothing (@exm-gmm-efficiency). Independence needs neither, estimates
  nothing extra, and cannot produce an invalid \( \R_i \).
- **Compare structures to check the mean model.** A large difference between the independence
  and exchangeable estimates is a warning about \( \X \), not about \( \R_i \).
- **Say which covariance you report**, and never put a GEE coefficient and a GLMM
  coefficient in the same table as if they estimated the same thing. The first costs
  coverage (@exm-gmm-coverage); the second is the commonest error in this literature.

::: {.remark}
[Beyond correlations]

For binary responses a correlation is an awkward measure of dependence, because the Fréchet
bounds on the joint cell probability confine it to an interval that depends on the two
means, while a pairwise odds ratio is free in \( (0,\infty) \) whatever the means are.
**Alternating logistic regressions**, due to Carey, Zeger and Diggle (1993), exploit this:
they replace \( \R_i \) by a model for the pairwise odds ratios and alternate between the
mean equation @eq-gmm-gee and a logistic regression for the association, whose parameters
are unconstrained and interpretable in their own right.
:::

## Exercises

### A. Check your understanding

::: {#exr-gmm-which-structure}
[A1]

For each design, name the working correlation you would use and say why: (a) pupils within
\( 300 \) schools, class size varying; (b) five annual measurements on each of
\( 500 \) patients; (c) left and right eyes of \( 80 \) patients; (d) irregularly timed
clinic visits, between two and eleven per patient, \( 1200 \) patients.
:::

::: {#exr-gmm-qic-reading}
[A2]

In @exm-gmm-qic, QIC prefers independence among the working correlations but the AR(1)
parameter is estimated at \( 0.894 \). Are these two statements in conflict?
Explain what each of them is measuring.
:::

### B. Practice

::: {#exr-gmm-exchangeable-bound}
[B1]

Show that the exchangeable correlation matrix
\( \R=(1-\alpha)\I+\alpha\bone\bone\T \) of size \( n \) has eigenvalues
\( 1+(n-1)\alpha \) (once) and \( 1-\alpha \) (\( n-1 \) times), so that it is positive
definite exactly when \( -1/(n-1)<\alpha<1 \). Verify the value
\( 1+3\hat\alpha=0.0956 \) quoted in @exm-gmm-modechoice from
\( \hat\alpha=-0.3015 \), and say what a practitioner should do when a
moment estimate falls outside the range.
:::

::: {.solution}
\( \R=(1-\alpha)\I+\alpha\bone\bone\T \); \( \bone \) is an eigenvector with eigenvalue
\( (1-\alpha)+n\alpha=1+(n-1)\alpha \), and any vector orthogonal to \( \bone \) has
eigenvalue \( 1-\alpha \). Positivity of both gives the range. With
\( \hat\alpha=-0.3015 \) and \( n=4 \),
\( 1+3(-0.3015)=1-0.9045=0.0955 \), which matches the quoted
\( 0.0956 \) to rounding. Outside the range the working covariance is not a
covariance at all and the update @eq-gmm-gee-update can diverge; the usual remedies are to
shrink \( \hat\alpha \) towards zero, to project \( \R \) onto the nearest correlation
matrix, or — better — to use the independence structure, which is never invalid.
:::

::: {#exr-gmm-hausman}
[B2]

Formalize the diagnostic of @exm-gmm-qic. Let \( \hbeta_{\text{I}} \) and
\( \hbeta_{\text{E}} \) be the independence and exchangeable GEE estimates. Explain why
\( \hbeta_{\text{I}}-\hbeta_{\text{E}} \) has mean zero to first order when the mean model is
correct, and why its covariance is *not* the difference of the two covariances in general.
What resampling scheme would give a usable reference distribution?
:::

::: {.solution}
Both estimates satisfy @eq-gmm-gee-asymptotics with the same \( \bbeta^{0} \), so their
difference is \( o_p(1) \) and, after scaling, asymptotically normal with mean zero. The
classical Hausman simplification \( \Cov(\hat\bbeta_1-\hat\bbeta_2)
=\Cov(\hat\bbeta_2)-\Cov(\hat\bbeta_1) \) needs one of the two to be *efficient*, which
neither is here, so the covariance must be obtained from the joint asymptotic distribution of
the two estimating functions or by resampling. The cluster bootstrap of @exr-bs-wild-cluster
— resample whole clusters with replacement and refit both — gives the
reference distribution directly.
:::
