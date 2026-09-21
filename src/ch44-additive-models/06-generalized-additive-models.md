# Generalized additive models

Everything so far assumed a normal response and an identity link. Lifting that costs one
paragraph, because [Part VIII](../ch34-exponential-families-glm/index.html) did the work: a
generalized linear model is fitted by repeated weighted least squares (@thm-glm-irls), and
repeated *penalized* weighted least squares is what this chapter has been doing throughout.
Put the two together and the model class that carries the name appears.

## The model and its fit

::: {#thm-add-gam}
[Generalized additive models]

Let \( y_1,\dots,y_n \) be independent from an exponential dispersion family in the sense
of @def-glm-model, with dispersion \( \phi \), variance function \( V \), link \( g \) and
mean \( \mu_i=g^{-1}(\eta_i) \), and let \( \boldsymbol{\upeta}=\Z\bgamma \) be the
structured additive predictor @eq-add-star with block penalty \( \bP \)
as in @eq-add-blockpenalty. Define the **penalized deviance**
\[
Q(\bgamma)=D(\bgamma)+\bgamma\T\bP\bgamma ,
\]{#eq-add-penlik}

with \( D \) the deviance of @def-glm-deviance, and let \( \hat{\bgamma} \) minimize it.

::: {.enumerate options="label=(\alph*)"}
1. Minimizing @eq-add-penlik is the same as maximizing the **penalized log-likelihood**
   \( \ell(\bgamma)-(2\phi)^{-1}\bgamma\T\bP\bgamma \). Neither the minimizer nor
   the iteration below depends on \( \phi \).

2. *(Penalized IRLS.)* With \( \W \) the working weights and \( \bz \) the working
   response of @eq-glm-working-response, both evaluated at \( \bgamma^{(t)} \), the Fisher
   scoring step for @eq-add-penlik is
   \[
   \bgamma^{(t+1)}=\bigl(\Z\T\W\Z+\bP\bigr)^{-1}\Z\T\W\bz ,
   \]{#eq-add-pirls}

   that is, the penalized weighted least squares fit of the working response on the same
   design. At a fixed point the unpenalized columns satisfy
   \( \X\T\W\bD^{-1}(\y-\bmu)=\bzero \), where \( \bD=\diag\{h'(\eta_i)\} \); for the
   canonical link with prior weights \( w_i\equiv1 \) this is \( \X\T(\y-\bmu)=\bzero \).

3. *(Degrees of freedom.)* At convergence, with \( \A=\Z\T\W\Z+\bP \), the
   effective degrees of freedom are \( \mathrm{df}=\tr(\A^{-1}\Z\T\W\Z) \) and split over
   terms exactly as in @prp-add-df(b). The approximate covariance of \( \hat{\bgamma} \)
   is \( \phi\A^{-1}\Z\T\W\Z\A^{-1} \) in its frequentist form and \( \phi\A^{-1} \) in the
   Bayesian form of @prp-add-mixed.

4. *(Inference for a term.)* Let \( \V_j=\phi(\A^{-1})_{jj} \) and let \( \V_j^{r-} \) be
   the rank-\( r \) pseudo-inverse of \( \V_j \) with \( r \) the integer nearest to
   \( \mathrm{df}_j \). Under \( f_j\equiv0 \), and treating
   \( \boldsymbol{\uplambda} \) as fixed and \( \phi \) as known,
   \( T_j=\hat{\bgamma}_j\T\V_j^{r-}\hat{\bgamma}_j \) has approximately the
   \( \chi^2(r) \) distribution; with \( \phi \) estimated, \( T_j/r \) is referred to
   \( F(r,\,n-\mathrm{df}) \). This is an approximation, stated here without proof; see the
   notes.
:::

:::

::: {.proof}
(a) By @def-glm-deviance, \( D(\bgamma)=2\phi\{\ell_{\mathrm{sat}}-\ell(\bgamma)\} \) with
\( \ell_{\mathrm{sat}} \) free of \( \bgamma \); dividing @eq-add-penlik by \( 2\phi \)
and dropping the constant gives the stated criterion, and rescaling an objective by a
positive constant does not move its minimizer.

(b) Let \( \bU=\partial\ell/\partial\bgamma \) and
\( \boldsymbol{\mathcal I}=\Z\T\W\Z/\phi \) be the score and expected information of the
unpenalized model (@thm-glm-score), so \( \partial D/\partial\bgamma=-2\phi\bU \) and
\( \E(\partial^2D/\partial\bgamma\partial\bgamma\T)=2\Z\T\W\Z \). Fisher scoring for
@eq-add-penlik therefore takes
\[
\begin{aligned}
\bgamma^{(t+1)}
&=\bgamma^{(t)}-\bigl(2\Z\T\W\Z+2\bP\bigr)^{-1}
\bigl(-2\phi\bU+2\bP\bgamma^{(t)}\bigr)\\
&=\A^{-1}\bigl\{\A\bgamma^{(t)}+\phi\bU-\bP\bgamma^{(t)}\bigr\},
\end{aligned}
\]
and \( \A\bgamma^{(t)}-\bP\bgamma^{(t)}=\Z\T\W\Z\bgamma^{(t)}=\Z\T\W\boldsymbol{\upeta}^{(t)} \),
while \( \phi\bU=\Z\T\W\bD^{-1}(\y-\bmu) \) by @thm-glm-score. Adding,
\[
\bgamma^{(t+1)}=\A^{-1}\Z\T\W\bigl\{\boldsymbol{\upeta}^{(t)}+\bD^{-1}(\y-\bmu)\bigr\}
=\A^{-1}\Z\T\W\bz^{(t)} ,
\]
which is @eq-add-pirls; \( \phi \) has cancelled, proving the last claim of (a). At a fixed
point the gradient of @eq-add-penlik vanishes, that is
\( \Z\T\W\bD^{-1}(\y-\bmu)=\bP\bgamma \); the rows of \( \bP \) belonging to
\( \X \) are zero, so those components read \( \X\T\W\bD^{-1}(\y-\bmu)=\bzero \). For the
canonical link \( \W\bD^{-1}=\diag\{w_ih'(\eta_i)/V(\mu_i)\} \) reduces to
\( \diag(w_i) \), since \( h'(\eta)=V(\mu) \) there (@thm-glm-score(d)), and with
\( w_i\equiv1 \) the identity becomes \( \X\T(\y-\bmu)=\bzero \).

(c) At convergence the fitted linear predictor is \( \Z\A^{-1}\Z\T\W\bz \), of the same
form as in @prp-add-df; its trace is
\( \tr(\Z\A^{-1}\Z\T\W)=\tr(\A^{-1}\Z\T\W\Z) \) (@thm-mat-trace-cyclic), and the block
decomposition of @prp-add-df(b) goes through with \( \Z\T\W\Z \) in place of
\( \Z\T\Z \). Treating \( \bz \) as a response with covariance \( \phi\W^{-1} \), as in
@thm-glm-irls, gives \( \Cov(\hat{\bgamma})\approx\phi\A^{-1}\Z\T\W\Z\A^{-1} \); the
Bayesian form is @prp-add-mixed with \( \R=\phi\W^{-1} \).

(d) Not proved here; see the notes.
:::

::: {.idea}
Penalized IRLS is the whole algorithm: a generalized additive model is a generalized linear
model whose weighted least squares step has acquired a penalty. Deviance, residuals,
dispersion, offsets, the canonical-link orthogonality all survive, because the working
problem is the same working problem.
:::

What Hastie and Tibshirani (1990) called **local scoring** differs only in solving the
inner penalized weighted least squares problem by backfitting (@thm-add-backfitting) rather
than by one direct solve — essential when a \( d\times d \) solve was expensive, a matter
of taste now. @exr-add-local-scoring works out the correspondence.

The inner loop needs one addition to the toolkit: the working weights and the working
response. The first cell below repeats, unchanged, the term constructor of
[Section 44.2](02-backfitting.html) and the solver of
[Section 44.5](05-structured-additive.html), so that this page runs on its own; the second
is the new part.

```{.python .run #cell-psplines-toolkit}
import numpy as np
from scipy.interpolate import BSpline

def bspline_basis(x, n_basis=15, degree=3, lo=None, hi=None):
    """n_basis B-splines of the given degree on equally spaced knots covering the data."""
    lo = x.min() if lo is None else lo
    hi = x.max() if hi is None else hi
    inner = np.linspace(lo, hi, n_basis - degree + 1)
    h = inner[1] - inner[0]
    knots = np.r_[lo - h * np.arange(degree, 0, -1), inner, hi + h * np.arange(1, degree + 1)]
    return BSpline.design_matrix(np.clip(x, lo, hi), knots, degree, extrapolate=False).toarray()

def difference_penalty(d, order=2):
    """K = D'D for the order-th difference matrix D: the P-spline penalty."""
    D = np.diff(np.eye(d), n=order, axis=0)
    return D.T @ D

def spline_term(x, n_basis=15, order=2, weight=None):
    """A P-spline term: its design, its penalty, and the map from new covariate values to
    design rows. The constraint sum_i f(x_i) = 0 is absorbed into the basis, and
    weight = u turns the term into the varying-coefficient term f(x) u."""
    lo, hi = x.min(), x.max()
    B = bspline_basis(x, n_basis, lo=lo, hi=hi)
    Q, _ = np.linalg.qr(B.sum(axis=0)[:, None], mode="complete")
    U = Q[:, 1:]                                  # a basis of {gamma : sum_i f(x_i) = 0}
    Z = B if weight is None else B * np.asarray(weight)[:, None]

    def design(xnew, wnew=None):
        Bn = bspline_basis(np.asarray(xnew, float), n_basis, lo=lo, hi=hi)
        return (Bn if wnew is None else Bn * np.asarray(wnew)[:, None]) @ U

    return Z @ U, U.T @ difference_penalty(n_basis, order) @ U, design

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

def null_dim(K, tol=1e-9):
    """The dimension of the penalty null space: the part of a term that is never penalized."""
    e = np.linalg.eigvalsh(K)
    return int(np.sum(e < tol * max(e.max(), 1.0)))
```

```{.python .run #cell-psplines-pirls}
def penalized_irls(y, blocks, X, lambdas, family="poisson", offset=0.0, iters=60, tol=1e-10):
    """Penalized IRLS: at each step a penalized weighted least squares fit of the working
    response on the same design -- Chapter 34's IRLS with the penalty added."""
    eta = np.log(np.maximum(y, 0.5)) if family == "poisson" else np.zeros(len(y))
    for _ in range(iters):
        if family == "poisson":
            mu = np.exp(eta + offset)
            w = dmu = mu                            # w = (dmu/deta)^2 / V(mu) = mu
        else:                                       # binomial with the logit link
            mu = 1 / (1 + np.exp(-(eta + offset)))
            w = dmu = mu * (1 - mu)
        w = np.maximum(w, 1e-10)
        z = eta + (y - mu) / dmu                    # the working response
        fit = Fit(z, blocks, X, lambdas, W=w, z=z, operators=False)
        done = np.max(np.abs(fit.fitted - eta)) < tol
        eta = fit.fitted
        if done:
            break
    fit.eta, fit.mu = eta, mu
    return fit

def pearson_dispersion(y, fit, family="poisson"):
    """The Pearson estimate of the dispersion, with the effective degrees of freedom."""
    v = fit.mu if family == "poisson" else fit.mu * (1 - fit.mu)
    return np.sum((y - fit.mu) ** 2 / v) / (len(y) - fit.edf_total)

def poisson_deviance(y, mu):
    """The Poisson deviance of Chapter 34, with the 0 log 0 convention."""
    t = np.where(y > 0, y * np.log(np.maximum(y, 1e-300) / mu), 0.0)
    return 2 * np.sum(t - (y - mu))

def irls_reml(y, blocks, X, family="poisson", offset=0.0, outer=40, tol=1e-7,
              dispersion=False):
    """Penalized IRLS with the smoothing parameters updated at every step: each working
    problem gets the update of reml_lambdas, with the dispersion of Chapter 38 in place of
    the nominal phi = 1 when dispersion is True."""
    lam = np.ones(len(blocks))
    for _ in range(outer):
        fit = penalized_irls(y, blocks, X, lam, family=family, offset=offset)
        phi = pearson_dispersion(y, fit, family) if dispersion else 1.0
        new = lam.copy()
        start = X.shape[1]
        for j, (_, Kj) in enumerate(blocks):
            d = Kj.shape[0]
            g = fit.gamma[start:start + d]
            tau2 = max(g @ Kj @ g, 1e-12) / max(fit.edf[j] - null_dim(Kj), 1e-6)
            new[j] = min(max(phi / tau2, 1e-8), 1e10)
            start += d
        done = np.max(np.abs(np.log(new) - np.log(lam))) < tol
        lam = new
        if done:
            break
    fit = penalized_irls(y, blocks, X, lam, family=family, offset=offset)
    fit.phi = pearson_dispersion(y, fit, family)
    return lam, fit
```

## Choosing the smoothing parameters, and the dispersion

@eq-add-reml-update was derived for a normal response. The standard device, the
*performance iteration* of Gu (1992), applies it to the working problem: at each penalized
IRLS step treat \( \bz \) as a normal response with covariance \( \phi\W^{-1} \), update
\( \boldsymbol{\uplambda} \) by @eq-add-reml-update with the weights in place, and
continue. That is what `irls_reml` does above. Recomputing \( \boldsymbol{\uplambda} \) by
maximizing a properly defined restricted likelihood for the generalized model at every step
is more stable and is what modern software does; Wood (2011) gives it.

One choice there matters more than it looks. For a Poisson or binomial response the nominal
dispersion is \( \phi=1 \), and taking that at face value when the data are overdispersed
makes \( \lambda_j=\hat\phi/\hat\tau_j^2 \) too small by the factor \( \hat\phi \), so
every curve is undersmoothed. [Chapter 38](../ch38-quasi-likelihood/index.html) built the repair: estimate \( \phi \) by the
Pearson statistic over the residual degrees of freedom (@prp-ql-dispersion) and use *that*
in the update, which is the quasi-likelihood position of @def-ql-quasi.

::: {#exm-add-gam}
[Physician visits, with two smooth terms]

@exm-cnt-visits fitted a Poisson log-linear model to the RAND Health Insurance Experiment,
\( n=20190 \) person-years, with every regressor entering linearly — among
them a general disease index, which has no reason to be linear on the log scale. Take that
index and one further variable of the same kind, the log of the participation incentive
payment, and give each a P-spline term of \( 12 \) basis functions, keeping the
coinsurance, deductible-plan, physical-limitation and health-status regressors linear: that
is a generalized additive model, fitted by @eq-add-pirls. Its baseline throughout is the
**linear fit**, meaning the same list of regressors with the disease index and the incentive
payment entering linearly too.

The disease term takes \( 9.11 \) effective degrees of freedom and the
incentive term \( 4.48 \), out of \( 20.59 \) for the whole
model. The estimated dispersion is \( 6.22 \), so these data are as
overdispersed as [Chapter 37](../ch37-counts/index.html) found them; the smoothing
parameters, \( 1.16 \) and \( 103.94 \), were chosen with
that dispersion in place. The coefficient of log coinsurance is
\( -0.1146 \), a rate ratio of \( 0.8917 \), against
\( -0.0874 \) in the linear fit and
\( -0.0712 \) in the smaller model of @exm-cnt-visits, which carried no
incentive-payment term at all; self-reported physical limitation still multiplies the rate
by \( e^{0.2777} \).

The deviance falls from \( 84378 \) for the linear fit to
\( 83619 \), at a cost of \( 11.59 \) degrees of freedom. Scaling
by the dispersion gives a quasi-\( F \) statistic of \( 10.54 \) on
\( 11.59 \) and \( 20169 \) degrees of freedom, whose nominal
\( p \)-value is \( 1.6\times 10^{-20} \). The deviance difference is the
statistic of @thm-glm-deviance(a), but its reference distribution there assumed a fixed
nested pair, and here the larger model was smoothed. Read as a rough diagnostic rather
than as a test, it says what panel (a) of
[Figure 44.6.1](06-generalized-additive-models.html#fig-add-gam) shows: the effect of the
disease index is far from linear — it rises steeply over the first half of the range and
then flattens — and no straight line would do.

The same panel is a warning. Above a disease index of \( 40 \) lie
\( 55 \) person-years, \( 0.3\% \) of the sample, and the curve there
swings by more than a unit on the log scale with a standard error of
\( 0.347 \), against \( 0.011 \) at the median. That wiggle is
an artefact of a handful of observations and a globally chosen smoothing parameter; report
it as "no information above 40", not as a finding.
:::

::: {when-format="html"}
![**Figure 44.6.1.** A Poisson generalized additive model for physician visits. (a) and (b) the two
smooth terms, with two-standard-error bands and a rug of covariate quantiles. (c) and (d)
observed mean, and observed variance, against fitted mean over twenty bins, the last with
the Poisson and quasi-Poisson lines.](gam_visits.svg){#fig-add-gam width=100%}
:::

::: {when-format="pdf"}
![A Poisson generalized additive model for physician visits. (a) and (b) the two
smooth terms, with two-standard-error bands and a rug of covariate quantiles. (c) and (d)
observed mean, and observed variance, against fitted mean over twenty bins, the last with
the Poisson and quasi-Poisson lines.](gam_visits.pdf){width=100%}
:::

```{.python .run #cell-gam-visits-fit}
import numpy as np
import statsmodels.api as sm

data = sm.datasets.randhie.load_pandas().data          # 20190 person-years, public domain
y = data["mdvis"].to_numpy(float)                      # outpatient physician visits in the year
disea = data["disea"].to_numpy(float)                  # a general disease index
lpi = data["lpi"].to_numpy(float)                      # log participation incentive payment
linear = ["lncoins", "idp", "physlm", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in linear])

Z1, K1, design1 = spline_term(disea, n_basis=12)       # the two smooth terms
Z2, K2, design2 = spline_term(lpi, n_basis=12)
lambdas, fit = irls_reml(y, [(Z1, K1), (Z2, K2)], X, family="poisson", dispersion=True)

print("smoothing parameters", lambdas.round(2))
print(f"edf: disease index {fit.edf[0]:.2f}, incentive payment {fit.edf[1]:.2f}, "
      f"model {fit.edf_total:.2f}")
print(f"dispersion {fit.phi:.2f}")
for name, b in zip(linear, fit.beta[1:]):
    print(f"  {name:8s} {b:8.4f}   rate ratio {np.exp(b):.4f}")
```

```{.python .run #cell-gam-visits-compare}
X_linear = np.column_stack([X, disea, lpi])            # the same regressors, both smooths linear
fit_linear = penalized_irls(y, [], X_linear, [], family="poisson")
D_gam = poisson_deviance(y, fit.mu)
D_linear = poisson_deviance(y, fit_linear.mu)
ddf = fit.edf_total - fit_linear.edf_total
F = (D_linear - D_gam) / ddf / fit.phi                 # a quasi-F, not an exact null statistic
print(f"deviance {D_linear:.0f} (linear) against {D_gam:.0f} (smooth)")
print(f"quasi-F {F:.2f} on {ddf:.2f} and {len(y) - fit.edf_total:.0f} degrees of freedom")
```

## What the p-values are worth

::: {.warning}
[Inference after choosing the smoothing parameters]

The statistic in @thm-add-gam(d) and the quasi-\( F \) of @exm-add-gam both treat
\( \boldsymbol{\uplambda} \) as a constant. It is not: it was chosen from the same data,
by minimizing a criterion that rewards fitting those data well — @prp-sel-selection-bias in
a new costume. Three consequences, in order of importance:

1. The reference distribution is approximate even when everything else is right. How good
   the approximation of @thm-add-gam(d) is depends on the term; Wood (2013) reports both the
   cases where it works and the cases where it fails, and simulation is the only way to know
   which one is at hand.
2. A term that the smoothing-parameter selection has pushed to the boundary
   \( \lambda_j=\infty \) — smoothed to a straight line, or to zero — sits on the edge of
   the parameter space, where @thm-mix-boundary says the usual asymptotics fail.
3. A \( p \)-value for a term that was *selected* into the model by the same fit is not a
   \( p \)-value at all.

Use these numbers to rank terms and to screen, nothing more. For a claim that matters, hold
out data, or put the whole fitting procedure — smoothing-parameter selection included —
inside the resampling loop of
[Chapter 23](../ch23-resampling-inference/index.html).
:::

## Checking the fit

The model is still a generalized linear model at heart, so the diagnostics of
[Chapter 39](../ch39-glms-in-practice-bayes/index.html) transfer unchanged; panels (c) and
(d) of [Figure 44.6.1](06-generalized-additive-models.html#fig-add-gam) show two of them,
binned by fitted mean. Panel (c) is the binned analogue of the exact recalibration identity
@prp-bin-calibration, which a penalized fit no longer satisfies exactly, since the score
identity of @thm-add-gam(b) holds only for the unpenalized columns. Panel (d) puts the
observed variance in the same bins against the fitted mean: far above the Poisson line
\( \sigma^2=\mu \) and close to \( \hat\phi\mu \), so the mean model is adequate, the
variance model is not, and every standard error must carry the dispersion (@def-ql-quasi).
One check is specific to these models: a term whose effective degrees of freedom sit close to
\( d_j-1 \), the most its basis allows, is asking for more basis functions; refit with more
and see whether the curve moves.

A smoothing parameter shrinks a term towards its penalty null space but not through it:
however large \( \lambda_j \) grows, a P-spline term with a second-difference penalty keeps
a straight line and never leaves the model. Both standard repairs are penalties, so both live
inside @def-add-star: a second penalty on the null space, letting the restricted likelihood
shrink the term to zero (Marra and Wood, 2011; @exr-add-select); or a penalty that is not
differentiable at zero, the group lasso of
[Section 30.3](../ch30-regularization-boosting/03-group-penalties.html) with each term a
group, which removes terms exactly.

## Where this leaves the book

Part IX set out to drop the linear predictor, and it has. A modern additive fit is a linear
model in a basis, with a quadratic penalty per term, a variance ratio per penalty and a
mixed model underneath; the response may come from any exponential dispersion family, and
the terms may be curves, surfaces, maps or clusters.

One assumption has survived it all. Every model in this book, from
[Chapter 5](../ch05-model-and-least-squares/index.html) to this page, has modelled the
*mean* of the response and treated the rest of its distribution as a nuisance —
\( \sigma^2 \), or \( \phi \), or a variance function fixed by the family. When the
spread or the shape is itself the object of interest, that is not enough.
[Chapter 45](../ch45-quantile-gamlss/index.html), the last of the book, gives up the mean: it
estimates quantiles directly, and then lets *every* parameter of a response distribution
carry a structured additive predictor of its own (@def-qnt-gamlss) — this chapter's
machinery, run once per parameter.

## Exercises

### A. Check your understanding

::: {#exr-add-gam-phi}
[A1]

In @thm-add-gam(a) the dispersion cancels from the iteration but not from the
smoothing-parameter update. Explain the difference, and say what happens to the fitted
curves if \( \phi=1 \) is assumed when in truth \( \phi=6 \).
:::

::: {.solution}
Given \( \boldsymbol{\uplambda} \), \( \phi \) only rescales @eq-add-penlik, so the
minimizer is unchanged. The update \( \lambda_j=\phi/\tau_j^2 \) compares the penalty term
with the noise level, so getting the noise level wrong by a factor of six makes every
\( \lambda_j \) six times too small: undersmoothed curves, too many effective degrees of
freedom, bands too narrow.
:::

### B. Practice

::: {#exr-add-local-scoring}
[B1]

*Local scoring.* Show that replacing the direct solve in @eq-add-pirls by one or more
backfitting cycles on the weighted penalized normal equations gives the local scoring
algorithm of Hastie and Tibshirani, and that @thm-add-backfitting applies to the inner loop
with \( \Z_j\T\Z_j \) replaced by \( \Z_j\T\W\Z_j \). What must be true of \( \W \) for the
convergence proof to go through?
:::

::: {.solution}
@eq-add-pirls is \( \A\bgamma=\Z\T\W\bz \) with \( \A=\Z\T\W\Z+\bP \): the
structure of @eq-add-normal with the inner product
\( \inner{\bu}{\bv}_{\W}=\bu\T\W\bv \) in place of the Euclidean one, so backfitting on
it is block Gauss–Seidel, which is local scoring. The proof of @thm-add-backfitting needs
\( \A \) symmetric nonnegative definite with positive definite diagonal blocks, hence
\( \W \) nonnegative definite — working weights are, being
\( h'(\eta_i)^2/V(\mu_i)\ge0 \) — and positive on enough observations to keep each
\( \Z_j\T\W\Z_j+\lambda_j\bP_j \) nonsingular.
:::

::: {#exr-add-select}
[B2]

*Shrinking a whole term away.* Let \( \bP_j \) have null space of dimension
\( m_j \) with orthonormal basis \( \mathbf{P}_j^{0} \), and replace the penalty by
\( \bP_j+\epsilon\mathbf{P}_j^{0}\mathbf{P}_j^{0\top} \) for a small \( \epsilon>0 \).
Show that the new penalty is positive definite, that
\( \lambda_j\to\infty \) now shrinks the whole term to zero, and that
\( \mathrm{df}_j\to0 \).
:::

::: {.solution}
Split \( \bgamma=\bgamma^{+}+\bgamma^{0} \) over the range and null space of
\( \bP_j \); the new quadratic form is
\( \bgamma^{+\top}\bP_j\bgamma^{+}+\epsilon\norm{\bgamma^{0}}^2 \), zero only if
both pieces vanish, so it is positive definite with null space \( \{\bzero\} \). By
@prp-add-df(b) the term's degrees of freedom then fall to \( 0 \) as
\( \lambda_j\to\infty \), with \( \hat{\bgamma}_j\to\bzero \): the restricted likelihood
can select the term out, which is what @eq-add-reml-update does when it explains nothing.
:::

### C. Going deeper

::: {#exr-add-gam-boundary-count}
[C1]

In @exm-add-gam the curve above a disease index of \( 40 \) rests on
\( 55 \) observations. Suppose all of them were deleted. Using the fact that
the P-spline penalty of @def-smo-pspline has the straight lines as its null space, describe
exactly what the fitted curve does over the empty range, and compute the limit of its
variance as the range of emptiness grows.
:::

