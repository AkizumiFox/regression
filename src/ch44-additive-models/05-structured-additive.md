# Structured additive regression

Four kinds of term have appeared: a smooth function of a covariate, a coefficient varying
with a modifier, an effect attached to a region of a map, an effect attached to a cluster.
Each was a pair \( (\Z_j,\bP_j) \) — a design block and a nonnegative definite
penalty — estimated by the same penalized criterion. This section says that out loud, shows
that the framework contains the mixed models of
[Chapter 32](../ch32-linear-mixed-models/index.html) exactly, and takes from them the thing
still missing: a principled way to choose the smoothing parameters.

## One predictor

::: {#def-add-star}
[Structured additive regression]

A **structured additive predictor** is
\[
\boldsymbol{\upeta}=\X\bbeta+\sum_{j=1}^{J}\Z_j\bgamma_j ,
\]{#eq-add-star}

where \( \X \) collects the unpenalized columns and each term \( j \) is a pair
\( (\Z_j,\bP_j) \) with \( \Z_j \) an \( n\times d_j \) design block and
\( \bP_j \) a symmetric nonnegative definite \( d_j\times d_j \) penalty matrix,
identified by the constraint \( \bone\T\Z_j\bgamma_j=0 \) of
@prp-add-identifiability(c). A **structured additive regression model** (STAR) fits
@eq-add-star by minimizing the penalized criterion @eq-add-criterion, with one smoothing
parameter \( \lambda_j\ge0 \) per term.
:::

The content of the definition is the catalogue of pairs. Everything this chapter and
[Chapter 43](../ch43-smoothing/index.html) have built is in
[Table 44.5.1](05-structured-additive.html#tab-add-terms), and the list is open: a new kind
of effect enters as soon as somebody writes down its design and its penalty.

[**Table 44.5.1.** The catalogue of terms. Each row is a design block and a penalty; the
last column is the dimension of the part of the term that the penalty never touches, which
is what survives as \( \lambda_j\to\infty \) and which the centring constraint reduces by
one. In the last row \( m_1 \) and \( m_2 \) are the null-space dimensions of the two
marginal penalties, so \( 1 \) for first differences in both coordinates and \( 4 \) for
second differences.]{#tab-add-terms}

| term | design block \( \Z_j \) | penalty \( \bP_j \) | \( \dim\Null(\bP_j) \) |
|---|---|---|---|
| polynomial of degree \( d \) | powers of \( z \) | \( \bzero \) | \( d+1 \) |
| P-spline, difference order \( r \) | B-splines at \( z \) | \( \mathbf{D}_r\T\mathbf{D}_r \) | \( r \) |
| smoothing spline | natural spline basis | \( \int B_l''B_m'' \) | \( 2 \) |
| varying coefficient | B-splines times \( u \) | as for the P-spline | as for the P-spline |
| Markov random field | region indicators | @eq-add-mrf | components of the map |
| kriging term | radial basis at knots | \( \{r(\norm{\kappa_l-\kappa_m})\} \) | \( 0 \) |
| random intercept | cluster indicators | \( \I \) | \( 0 \) |
| random slope | cluster indicators times \( u \) | \( \I \) | \( 0 \) |
| tensor-product surface | products of B-splines | @eq-add-tensor | \( m_1m_2 \) |

::: {#prp-add-star-identified}
[When a STAR model is identified]

In @def-add-star, assume every \( \Z_j \) has full column rank. The penalized criterion has
a unique minimizer for every \( \lambda_j>0 \) if and only if the only vector
\( \bgamma=(\bbeta\T,\bgamma_1\T,\dots)\T \) with
\( \X\bbeta+\sum_j\Z_j\bgamma_j=\bzero \) and \( \bP_j\bgamma_j=\bzero \) for every
\( j \) is \( \bgamma=\bzero \). Equivalently: there is no exact concurvity among
\( \C(\X) \) and the *unpenalized parts* \( \Z_j\Null(\bP_j) \) of the terms.
:::

::: {.proof}
The criterion is the convex quadratic \( \norm{\y-\Z\bgamma}^2+\bgamma\T\bP\bgamma \)
with \( \bP \) the block penalty @eq-add-blockpenalty, and its stationarity condition
is @eq-add-normal with \( \A=\Z\T\Z+\bP \). If \( \A \) is positive definite the
criterion is strictly convex and @thm-reg-existence gives a unique minimizer. Conversely, if
\( \A \) is singular, pick \( \bgamma_0\ne\bzero \) with \( \A\bgamma_0=\bzero \); expanding
the criterion at a minimizer \( \hat{\bgamma} \) gives
\( Q(\hat{\bgamma}+t\bgamma_0)=Q(\hat{\bgamma})+t^2\bgamma_0\T\A\bgamma_0=Q(\hat{\bgamma}) \)
for every \( t \), so the minimizer is not unique. Uniqueness for every \( \lambda_j>0 \) is
therefore positive definiteness of \( \A \), that is
\( \bgamma\T\A\bgamma=\norm{\Z\bgamma}^2+\sum_j\lambda_j\bgamma_j\T\bP_j\bgamma_j=0 \)
forces \( \bgamma=\bzero \). Both summands are nonnegative, so the condition is exactly the
one stated, and \( \bgamma_j\T\bP_j\bgamma_j=0 \) is equivalent to
\( \bP_j\bgamma_j=\bzero \) by @thm-mat-pd-characterizations.
:::

Two P-spline terms can therefore be identified even when their spaces nearly intersect,
because only their *linear* parts — the null spaces of the second-difference penalties —
need be independent. That is far weaker than the requirement for the unpenalized model in
@prp-add-concurvity(a), and it explains why penalized additive fits rarely fail outright:
they fail softly, through inflated variances.

## The mixed-model representation

::: {#prp-add-mixed}
[A STAR model is a linear mixed model]

Let \( \bP_j \) have rank \( r_j \) and spectral decomposition
\( \bP_j=\mathbf{P}_j^{+}\bLambda_j\mathbf{P}_j^{+\top} \) with
\( \bLambda_j=\diag(\nu_{j1},\dots,\nu_{jr_j}) \) positive, and let
\( \mathbf{P}_j^{0} \) be an orthonormal basis of \( \Null(\bP_j) \). Put
\[
\X_j^{0}=\Z_j\mathbf{P}_j^{0},\qquad
\Z_j^{+}=\Z_j\mathbf{P}_j^{+}\bLambda_j^{-1/2},\qquad
\bu_j=\bLambda_j^{1/2}\mathbf{P}_j^{+\top}\bgamma_j .
\]{#eq-add-mixed-transform}

::: {.enumerate options="label=(\alph*)"}
1. \( \Z_j\bgamma_j=\X_j^{0}\mathbf{b}_j+\Z_j^{+}\bu_j \) with
   \( \mathbf{b}_j=\mathbf{P}_j^{0\top}\bgamma_j \), and
   \( \bgamma_j\T\bP_j\bgamma_j=\norm{\bu_j}^2 \).

2. With \( \tau_j^2=\sigma^2/\lambda_j \), the penalized criterion @eq-add-criterion
   divided by \( \sigma^2 \) is the criterion @eq-mix-penalized of @thm-mix-henderson for
   the linear mixed model @eq-mix-model with fixed design \( [\X,\X_1^{0},\dots,\X_J^{0}] \),
   random design \( [\Z_1^{+},\dots,\Z_J^{+}] \),
   \( \G=\diag(\tau_1^2\I_{r_1},\dots,\tau_J^2\I_{r_J}) \) and \( \R=\sigma^2\I \).

3. Consequently the penalized estimate is the BLUP of @thm-mix-blup; the smoothing
   parameters are variance ratios, \( \lambda_j=\sigma^2/\tau_j^2 \); and the restricted
   likelihood @eq-mix-reml-loglik of the transformed model is a likelihood for
   \( (\sigma^2,\tau_1^2,\dots,\tau_J^2) \), hence for
   \( (\lambda_1,\dots,\lambda_J) \).
:::

:::

::: {.proof}
(a) The columns of \( \mathbf{P}_j^{+} \) and \( \mathbf{P}_j^{0} \) form an orthonormal
basis of \( \Real^{d_j} \), so
\( \bgamma_j=\mathbf{P}_j^{0}\mathbf{b}_j+\mathbf{P}_j^{+}\mathbf{P}_j^{+\top}\bgamma_j \)
and \( \mathbf{P}_j^{+}\mathbf{P}_j^{+\top}\bgamma_j
=\mathbf{P}_j^{+}\bLambda_j^{-1/2}\bu_j \). Multiplying by \( \Z_j \) gives the first
claim. For the second,
\( \bgamma_j\T\bP_j\bgamma_j
=\bgamma_j\T\mathbf{P}_j^{+}\bLambda_j\mathbf{P}_j^{+\top}\bgamma_j=\norm{\bu_j}^2 \).
(b) Substituting (a) into @eq-add-criterion turns it into
\[
\Bigl\lVert\y-\X\bbeta-\sum_j\X_j^{0}\mathbf{b}_j-\sum_j\Z_j^{+}\bu_j\Bigr\rVert^2
+\sum_j\lambda_j\norm{\bu_j}^2 ,
\]
and dividing by \( \sigma^2 \) and writing \( \lambda_j/\sigma^2=1/\tau_j^2 \) gives
\( \R^{-1}=\sigma^{-2}\I \) in the first term and \( \G^{-1} \) in the second,
which is @eq-mix-penalized.
(c) Immediate from @thm-mix-henderson(a)–(b), which identifies the minimizer of
@eq-mix-penalized with \( (\hbeta,\hat{\bu}) \), and from @def-mix-reml applied to the
transformed model.
:::

::: {.idea}
A penalty is a prior, and a smoothing parameter is a variance ratio. That one
identification does three jobs: the null space of a penalty becomes *fixed* effects, which
is why [Table 44.5.1](05-structured-additive.html#tab-add-terms) records
\( \dim\Null(\bP_j) \); it supplies a likelihood for
\( \boldsymbol{\uplambda} \), which cross-validation does not; and it puts every result of
[Chapter 32](../ch32-linear-mixed-models/index.html) — BLUP, REML, Henderson's equations,
the posterior of @prp-mix-bayes — at the service of smoothing.
:::

## Choosing the smoothing parameters

Two routes are now open. The first is prediction error: @def-sel-gcv and the information
criteria @prp-sel-aic-bic apply verbatim with \( \tr\bS \) in place of the number of
parameters, minimized over \( \boldsymbol{\uplambda} \) numerically. The second is the
restricted likelihood of @prp-add-mixed(c), usually preferred: it is less prone to the
occasional catastrophic undersmoothing cross-validation shows when its criterion is nearly
flat, and it comes with the variance interpretation. Both choose from the data the fit uses,
so @prp-sel-selection-bias applies to both.

Maximizing the restricted likelihood directly is awkward, since each evaluation needs a
determinant, and the updates below avoid it. What is proved is that they have the stated
fixed points; that those are stationary points of the restricted likelihood is due to Schall
(1991), with the general version in Wood and Fasiolo (2017).

::: {#prp-add-reml-update}
[Restricted-likelihood updates]

In the mixed-model representation of @prp-add-mixed, the expectation-maximization algorithm
for the variance components, with \( (\mathbf{b},\bu) \) as the missing data, updates
\[
\tau_j^{2,\text{new}}
=\frac{1}{r_j}\Bigl[\hat{\bgamma}_j\T\bP_j\hat{\bgamma}_j
+\sigma^2\tr\bigl\{(\A^{-1})_{jj}\bP_j\bigr\}\Bigr] ,
\]{#eq-add-em-update}

where \( \hat{\bgamma}_j \) and \( \A \) are computed at the current
\( \boldsymbol{\uplambda} \). A pair \( (\boldsymbol{\uptau}^2,\sigma^2) \) is a fixed point
of @eq-add-em-update together with the update
\( \sigma^{2,\text{new}}=\norm{\y-\hat{\y}}^2/(n-\mathrm{df}) \) if and only if
\[
\hat{\tau}_j^{2}=\frac{\hat{\bgamma}_j\T\bP_j\hat{\bgamma}_j}{\mathrm{df}_j-m_j},
\qquad
\hat{\sigma}^2=\frac{\norm{\y-\hat{\y}}^2}{n-\mathrm{df}},
\qquad
\hat\lambda_j=\frac{\hat\sigma^2}{\hat\tau_j^2} ,
\]{#eq-add-reml-update}

with \( m_j=d_j-r_j=\dim\Null(\bP_j) \) and \( \mathrm{df}_j \), \( \mathrm{df} \)
the effective degrees of freedom of @prp-add-df.
:::

::: {.proof}
@eq-add-em-update is the expectation-maximization step for a variance component in a normal
linear mixed model: the complete-data estimate of \( \tau_j^2 \) is
\( \norm{\bu_j}^2/r_j \), whose conditional expectation given \( \y \) is
\( \norm{\hat{\bu}_j}^2/r_j+\tr\{\Cov(\bu_j\mid\y)\}/r_j \). By @prp-add-mixed(a),
\( \norm{\hat{\bu}_j}^2=\hat{\bgamma}_j\T\bP_j\hat{\bgamma}_j \), and by
@prp-mix-bayes(a) the conditional covariance of \( \bgamma_j \) is
\( \sigma^2(\A^{-1})_{jj} \), so @eq-add-mixed-transform gives

\[
\tr\Cov(\bu_j\mid\y)
=\sigma^2\tr\bigl\{\bLambda_j^{1/2}\mathbf{P}_j^{+\top}(\A^{-1})_{jj}
\mathbf{P}_j^{+}\bLambda_j^{1/2}\bigr\}
=\sigma^2\tr\bigl\{(\A^{-1})_{jj}\bP_j\bigr\},
\]

the trace term as written. (The EM step itself is quoted from the mixed-model literature;
what is proved here is its fixed point.)

For the fixed point, @prp-add-df(b) gives
\( \mathrm{df}_j=d_j-\lambda_j\tr\{(\A^{-1})_{jj}\bP_j\} \), that is
\[
\tr\bigl\{(\A^{-1})_{jj}\bP_j\bigr\}
=\frac{d_j-\mathrm{df}_j}{\lambda_j}
=\frac{(d_j-\mathrm{df}_j)\,\tau_j^2}{\sigma^2}.
\]
Substituting into @eq-add-em-update and demanding
\( \tau_j^{2,\text{new}}=\tau_j^2 \),
\[
r_j\tau_j^2=\hat{\bgamma}_j\T\bP_j\hat{\bgamma}_j+(d_j-\mathrm{df}_j)\tau_j^2 ,
\]
so \( \tau_j^2(r_j-d_j+\mathrm{df}_j)=\hat{\bgamma}_j\T\bP_j\hat{\bgamma}_j \), and
\( r_j-d_j=-m_j \) gives the first equation of @eq-add-reml-update. The statement for
\( \sigma^2 \) is its own update, and \( \lambda_j=\sigma^2/\tau_j^2 \) is
@prp-add-mixed(c).
:::

Read @eq-add-reml-update as an estimate of a variance: the numerator is the roughness the
term shows, the denominator \( \mathrm{df}_j-m_j \) the number of dimensions it used in
showing it, discounting the unpenalized part the penalty cannot see. The solver is one
assembly and one solve,

```{.python .run #cell-psplines-fit}
import numpy as np

def assemble(blocks, X, lambdas):
    """Stack the parametric columns and the terms into one design and one block penalty."""
    Z = np.column_stack([X] + [Zj for Zj, _ in blocks])
    K = np.zeros((Z.shape[1], Z.shape[1]))
    start = X.shape[1]
    for (_, Kj), lam in zip(blocks, lambdas):
        d = Kj.shape[0]
        K[start:start + d, start:start + d] = lam * Kj
        start += d
    return Z, K

class Fit:
    """One penalized solve: coefficients, components, operators and degrees of freedom."""

    def __init__(self, y, blocks, X, lambdas, W=None, z=None, operators=True):
        Z, K = assemble(blocks, X, lambdas)
        w = np.ones(len(y)) if W is None else W
        Zw = Z * w[:, None]
        B = Z.T @ Zw                                           # the weighted cross-product
        A = B + K
        Ainv = np.linalg.inv(A)
        self.Z, self.K, self.A, self.w = Z, K, A, w
        self.gamma = Ainv @ (Zw.T @ (y if z is None else z))
        self.fitted = Z @ self.gamma
        self.edf_total = np.trace(Ainv @ B)                    # tr(S) without forming S
        self.parts, self.edf, self.comp = [], [], []
        start = X.shape[1]
        for Zj, _ in blocks:
            sl = slice(start, start + Zj.shape[1])
            self.parts.append(Zj @ self.gamma[sl])
            self.edf.append(np.trace(Ainv[sl] @ B[:, sl]))
            start += Zj.shape[1]
        if operators:                                          # the n x n smoothers
            C = Ainv @ Zw.T
            self.S = Z @ C
            start = X.shape[1]
            for Zj, _ in blocks:
                self.comp.append(Zj @ C[start:start + Zj.shape[1]])
                start += Zj.shape[1]
        self.beta = self.gamma[:X.shape[1]]
```

and the updates themselves are a handful more:

```{.python .run #cell-psplines-reml}
def null_dim(K, tol=1e-9):
    """The dimension of the penalty null space: the part of a term that is never penalized."""
    e = np.linalg.eigvalsh(K)
    return int(np.sum(e < tol * max(e.max(), 1.0)))

def reml_lambdas(y, blocks, X, lambdas=None, iters=80, tol=1e-8):
    """Restricted-likelihood updates for a normal response: refit, then set each variance
    component to tau_j^2 = gamma_j' K_j gamma_j / (df_j - m_j) and sigma^2 to the penalized
    residual sum of squares over n - df, and put lambda_j = sigma^2 / tau_j^2."""
    lam = np.ones(len(blocks)) if lambdas is None else np.asarray(lambdas, float)
    for _ in range(iters):
        fit = Fit(y, blocks, X, lam, operators=False)
        sigma2 = np.sum((y - fit.fitted) ** 2) / (len(y) - fit.edf_total)
        new = lam.copy()
        start = X.shape[1]
        for j, (_, Kj) in enumerate(blocks):
            d = Kj.shape[0]
            g = fit.gamma[start:start + d]
            tau2 = max(g @ Kj @ g, 1e-12) / max(fit.edf[j] - null_dim(Kj), 1e-6)
            new[j] = min(max(sigma2 / tau2, 1e-8), 1e10)
            start += d
        done = np.max(np.abs(np.log(new) - np.log(lam))) < tol
        lam = new
        if done:
            break
    return lam, Fit(y, blocks, X, lam), sigma2
```

::: {#exm-add-star}
[Four kinds of term at once]

A simulated data set of \( 1200 \) observations has a binary covariate with
a linear effect, a smooth effect of \( z_1 \), a coefficient of \( u \) varying with
\( z_2 \), an areal effect on an \( 8\times8 \) map of \( 64 \) districts and a
random intercept for \( 40 \) clusters: four terms,
\( 132 \) coefficients, smoothing parameters from @eq-add-reml-update.

The fit recovers the residual standard deviation as \( 0.713 \) (true
\( 0.7 \)) and the linear effect as \( 0.770 \) (true \( 0.8 \)). The four
terms take \( 7.86 \), \( 5.13 \),
\( 42.92 \) and \( 36.48 \) effective degrees of
freedom, \( 94.40 \) out of \( 132 \), with implied variance
components \( 0.1609 \), \( 0.0170 \),
\( 0.2002 \) and \( 0.2646 \) — the last the
variance of the random intercept, truly \( 0.25 \). Root mean squared errors against the
truth are \( 0.031 \), \( 0.045 \),
\( 0.129 \) and \( 0.119 \), against a noise
standard deviation of \( 0.7 \).
[Figure 44.5.1](05-structured-additive.html#fig-add-star) shows all four.
:::

::: {when-format="html"}
![**Figure 44.5.1.** One structured additive fit with four kinds of term. (a) the smooth term
against the truth. (b) the estimated coefficient function. (c) and (d) fitted district and
cluster effects against the truth, with the line of equality; both are shrunk, as a best
linear unbiased predictor must be.](star_terms.svg){#fig-add-star width=100%}
:::

::: {when-format="pdf"}
![One structured additive fit with four kinds of term. (a) the smooth term
against the truth. (b) the estimated coefficient function. (c) and (d) fitted district and
cluster effects against the truth, with the line of equality; both are shrunk, as a best
linear unbiased predictor must be.](star_terms.pdf){width=100%}
:::

The shrinkage in panels (c) and (d) is not a defect but @thm-mix-blup: a predictor of a
random effect trades bias for variance, and the slope of predicted against true falls below
one by exactly the amount the estimated variance ratio dictates.

## The Bayesian reading

@prp-add-mixed already contains it. Give each term the prior
\[
\pi(\bgamma_j\mid\tau_j^2)\ \propto\
(\tau_j^2)^{-r_j/2}\exp\Bigl\{-\frac{1}{2\tau_j^2}\bgamma_j\T\bP_j\bgamma_j\Bigr\},
\]{#eq-add-prior}

normal on \( \Null(\bP_j)\perpc \) and flat on the null space, give \( \bbeta \)
a flat prior, and the penalized estimate is the posterior mode (@def-prc-posterior); with
\( \boldsymbol{\uplambda} \) fixed the posterior of \( \bgamma \) is normal with
covariance \( \sigma^2\A^{-1} \) by @prp-mix-bayes(a). That covariance produces the bands
in this chapter's figures: the Wahba–Nychka bands of @thm-smo-mixed, with their caveat that
they are honest about bias averaged over the curve, not pointwise, and condition on
\( \boldsymbol{\uplambda} \) as though it were known. Treating
\( \boldsymbol{\uptau}^2 \) and \( \sigma^2 \) as unknown too gives a fully Bayesian
model, sampled as in
[Section 39.4](../ch39-glms-in-practice-bayes/04-bayesian-glms-mcmc.html) (@def-prc-mcmc),
whose advantage is precisely that it propagates the uncertainty in
\( \boldsymbol{\uplambda} \) into the bands instead of ignoring it.

::: {.warning}
[Improper priors on variance components]

The prior @eq-add-prior is improper on \( \Null(\bP_j) \), which is harmless: that
part is a fixed effect with a flat prior. The prior on \( \tau_j^2 \) is not harmless. An
inverse-gamma prior with both parameters tending to zero, or a flat prior on
\( \tau_j^2 \), can produce an improper posterior, and the symptom in a sampler is not a
crash but a chain that drifts towards zero variance and never settles (@prp-prc-priors). Use
a proper prior whose scale is set on the scale of the data, and check the sampler, never the
fit alone.
:::

## Exercises

### A. Check your understanding

::: {#exr-add-catalogue}
[A1]

For each of the following, give the design block and the penalty: (i) a random slope of a
covariate \( u \) within \( K \) clusters; (ii) a linear effect that is not to be penalized
at all; (iii) a smooth term that is to be penalized towards a straight line rather than
towards a constant.
:::

### B. Practice

::: {#exr-add-random-is-ridge}
[B1]

Take one term with \( \Z=\mathbf{D} \) the \( n\times K \) matrix of cluster indicators and
\( \bP=\I \). Show that the penalized estimate of the cluster effects is
\( \hat\gamma_k=m_k\bar r_k/(m_k+\lambda) \), where \( m_k \) is the cluster size and
\( \bar r_k \) the mean residual in cluster \( k \), and identify this with @thm-mix-blup
for the one-way random effects model of @def-mix-oneway.
:::

::: {.solution}
With the other effects fixed the criterion is
\( \sum_k\sum_{i\in k}(r_i-\gamma_k)^2+\lambda\sum_k\gamma_k^2 \), which separates over
\( k \) and gives \( \hat\gamma_k=m_k\bar r_k/(m_k+\lambda) \). With
\( \lambda=\sigma^2/\tau^2 \) that is
\( \bar r_k\,m_k\tau^2/(m_k\tau^2+\sigma^2) \), the shrinkage factor of @thm-mix-blup for
the one-way model.
:::

::: {#exr-add-gcv-vs-reml}
[B2]

Write the generalized cross-validation criterion @def-sel-gcv for a STAR fit as
\( \mathrm{GCV}(\boldsymbol{\uplambda})=n\norm{\y-\bS\y}^2/(n-\tr\bS)^2 \). Show that it is
minimized where \( \sigma^2 \) of @eq-add-reml-update is stationary in a particular sense,
and explain why the two criteria nevertheless disagree in general.
:::

### C. Going deeper

::: {#exr-add-reml-not-monotone}
[C1]

The updates @eq-add-reml-update are a fixed-point iteration, not a descent method, and they
are not guaranteed to converge from every start. Construct a two-term model in which the
restricted likelihood has a stationary point at \( \tau_1^2=0 \) (the term is smoothed
away), verify that this is a fixed point of @eq-add-reml-update, and explain how a run that
lands there should be reported.
:::

::: {.solution}
Let the first term be a P-spline in a covariate with no effect and the second anything. If
the data show no roughness in that term,
\( \hat{\bgamma}_1\T\bP_1\hat{\bgamma}_1\to0 \) while \( \mathrm{df}_1\to m_1 \),
so @eq-add-reml-update sends \( \tau_1^2\to0 \) and \( \lambda_1\to\infty \), and there
the update reproduces itself because the numerator stays zero. This is the boundary of the
parameter space, where @thm-mix-boundary says the usual asymptotics for a test of
\( \tau_1^2=0 \) fail; the honest report is "the term was smoothed to its null space", not
"significant at \( 0 \) degrees of freedom".
:::
