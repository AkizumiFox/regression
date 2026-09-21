# Mixed-model and Bayesian representations

The criteria of [Section 43.5](05-choosing-lambda.html) treat \( \lambda \) as a
tuning constant to be optimized. There is another way to think about it. A
quadratic penalty is, up to a constant, the logarithm of a normal density, so
penalizing \( \bgamma\T\bP\bgamma \) is the same as believing in advance that
\( \bgamma \) is a draw from a normal distribution with precision proportional to
\( \bP \). Then \( \lambda \) is a ratio of two variances, and a variance is not a
tuning constant but a parameter the data can estimate. The machinery of
[Chapter 32](../ch32-linear-mixed-models/index.html) becomes available, with a
second dividend: an interval whose width knows about the bias.

## Splitting the penalty

The obstacle is that \( \bP \) is singular. A normal prior with precision matrix
\( \bP/\tau^2 \) puts no constraint at all on the \( k \)-dimensional null space
\( \Null(\bP) \) — the straight lines, for a second-order penalty — so it is not a
proper distribution. The fix is to split \( \bgamma \) into the part the penalty
does not see, which becomes a *fixed* effect, and the rest, which becomes a
*random* effect with a proper distribution.

::: {#thm-smo-mixed}
[A penalized spline is a linear mixed model]

Let \( \bP \) be nonnegative definite of order \( m \) and rank \( r=m-k \), write
\( \bP=\bL\bL\T \) with \( \bL \) of size \( m\times r \) and full column rank, let
\( \X_\gamma \) be \( m\times k \) with columns spanning \( \Null(\bP) \), and set
\( \Z_\gamma=\bL(\bL\T\bL)^{-1} \).

::: {.enumerate options="label=(\alph*)"}
1. *(The transformation.)* \( [\X_\gamma,\Z_\gamma] \) is nonsingular, so every
   \( \bgamma \) is uniquely \( \bgamma=\X_\gamma\bbeta+\Z_\gamma\bu \); moreover
   \( \bu=\bL\T\bgamma \) and \( \bgamma\T\bP\bgamma=\bu\T\bu \).

2. *(The mixed model.)* With \( \X=\B\X_\gamma \) and \( \Z=\B\Z_\gamma \), the
   P-spline criterion @eq-smo-pspline is
   \[
   \norm{\y-\X\bbeta-\Z\bu}^2+\lambda\,\bu\T\bu ,
   \]{#eq-smo-mixed}

   which is the penalized criterion @eq-mix-penalized of the linear mixed model
   \( \Y=\X\bbeta+\Z\bu+\be \), \( \bu\sim\Normal_r(\bzero,\sigma_u^2\I) \),
   \( \be\sim\Normal_n(\bzero,\sigma^2\I) \) independent, with
   \[
   \lambda=\sigma^2/\sigma_u^2 .
   \]{#eq-smo-lambda-ratio}

   Its solution is Henderson's (@thm-mix-henderson): \( \hat{\bbeta} \) is the
   generalized least squares estimate and \( \hat{\bu} \) the BLUP
   (@thm-mix-blup).

3. *(The Bayesian reading.)* Suppose \( \sigma^2 \) is known and \( \bgamma \) has
   the prior with density proportional to
   \( \exp\{-\bgamma\T\bP\bgamma/(2\tau^2)\} \), \( \tau^2=\sigma^2/\lambda \) —
   equivalently \( \bbeta \) flat and \( \bu\sim\Normal_r(\bzero,\tau^2\I) \). Then
   \[
   \bgamma\mid\Y=\y\ \sim\
   \Normal_m\bigl\{\hat{\bgamma},\ \sigma^2(\B\T\B+\lambda\bP)^{-1}\bigr\} ,
   \]{#eq-smo-posterior}

   with \( \hat{\bgamma} \) the penalized estimate of @prp-smo-pspline-ridge(a).

4. *(REML.)* The marginal distribution is
   \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\V_\lambda) \) with
   \( \V_\lambda=\I+\lambda^{-1}\Z\Z\T \). Profiling \( \sigma^2 \) out of the
   restricted likelihood (@def-mix-reml, @thm-mix-reml) leaves
   \[
   -2\ell_R(\lambda)=(n-k)\log\hat\sigma^2_\lambda
   +\log\lvert\V_\lambda\rvert+\log\lvert\X\T\V_\lambda^{-1}\X\rvert+c,
   \]{#eq-smo-reml}

   where
   \( \hat\sigma^2_\lambda=\y\T\mathbf{A}_\lambda\y/(n-k) \) with
   \( \mathbf{A}_\lambda=\V_\lambda^{-1}-\V_\lambda^{-1}\X(\X\T\V_\lambda^{-1}\X)^{-1}\X\T\V_\lambda^{-1} \),
   and \( c \) does not depend on \( \lambda \). Minimizing @eq-smo-reml over
   \( \lambda>0 \) estimates the smoothing parameter.
:::

:::

::: {.proof}
(a) From \( \bP=\bL\bL\T \) we get \( \bP\bgamma=\bzero \) iff
\( \norm{\bL\T\bgamma}^2=\bgamma\T\bP\bgamma=0 \) iff \( \bL\T\bgamma=\bzero \), so
\( \Null(\bP)=\Null(\bL\T)=\C(\bL)^{\perp} \) and
\( \C(\Z_\gamma)=\C(\bL)=\Null(\bP)^{\perp} \). Hence
\( \C(\X_\gamma) \) and \( \C(\Z_\gamma) \) are orthogonal subspaces of dimensions
\( k \) and \( r \) that together span \( \Real^m \), so \( [\X_\gamma,\Z_\gamma] \)
has rank \( m \) and the decomposition is unique. Multiplying
\( \bgamma=\X_\gamma\bbeta+\Z_\gamma\bu \) by \( \bL\T \) gives
\( \bL\T\bgamma=\bzero+\bL\T\bL(\bL\T\bL)^{-1}\bu=\bu \), and then
\( \bgamma\T\bP\bgamma=\norm{\bL\T\bgamma}^2=\bu\T\bu \).

(b) Substituting \( \B\bgamma=\X\bbeta+\Z\bu \) and (a) into @eq-smo-pspline
gives @eq-smo-mixed. In the mixed model with \( \G=\sigma_u^2\I \) and
\( \R=\sigma^2\I \), Henderson's criterion is
\( \sigma^{-2}\norm{\y-\X\bbeta-\Z\bu}^2+\sigma_u^{-2}\bu\T\bu \), which is
\( \sigma^{-2} \) times @eq-smo-mixed when \( \lambda=\sigma^2/\sigma_u^2 \); the
two have the same minimizers.

(c) The prior is flat in \( \bbeta \) and normal in \( \bu \), so the joint density
of \( (\y,\bgamma) \) is proportional to
\( \exp[-\{\norm{\y-\B\bgamma}^2+\lambda\bgamma\T\bP\bgamma\}/(2\sigma^2)] \). By
the algebra in the proof of @prp-smo-pspline-ridge(a), the exponent is
\( -\{(\bgamma-\hat{\bgamma})\T(\B\T\B+\lambda\bP)(\bgamma-\hat{\bgamma})+C\}/(2\sigma^2) \)
with \( C \) free of \( \bgamma \). As a function of \( \bgamma \) this is the
normal density @eq-smo-posterior (@def-prc-posterior).

(d) \( \Cov(\Y)=\sigma_u^2\Z\Z\T+\sigma^2\I=\sigma^2\V_\lambda \)
by @eq-mix-marginal, and @thm-mix-reml gives the restricted log likelihood
\( -\tfrac12\{(n-k)\log\sigma^2+\log\lvert\V_\lambda\rvert
+\log\lvert\X\T\V_\lambda^{-1}\X\rvert+\y\T\mathbf{A}_\lambda\y/\sigma^2\} \) up to a
constant. Differentiating in \( \sigma^2 \) gives
\( \hat\sigma^2_\lambda=\y\T\mathbf{A}_\lambda\y/(n-k) \), and substituting it back
gives @eq-smo-reml.
:::

::: {.idea}
The same algebra has now been read three ways. A quadratic penalty is a
shrinkage rule (@thm-shr-ridge), a normal prior (@prp-prc-priors) and a variance
component. Only the third makes the tuning constant estimable, because only the
third supplies a probability model for the coefficients and with it a likelihood
in \( \lambda \). Nothing is obtained for nothing: that model says the true
\( f \) is a draw from a particular distribution of curves, a stronger assumption
than "\( f \) is smooth".
:::

::: {.remark}
[What the transformation looks like]

For a second-order difference penalty, \( \Null(\bD_2) \) is spanned by
\( \bone \) and \( (1,\dots,m)\T \), so the two fixed-effect columns of
\( \X=\B\X_\gamma \) span the constant and the linear
function (@prp-smo-penalty-limit(c)): the mixed model has a straight line as its
fixed part and the deviations from it as its random part — the decomposition a
data analyst would choose by hand, forced by the penalty rather than assumed.
Other choices of \( \bL \) and \( \X_\gamma \) give the same fitted values, since
they differ by a nonsingular reparameterization (@thm-proj-reparam).
:::

## Bands that know about the bias

@prp-smo-linear-smoother(e) left a problem: the obvious interval is centred at
\( \E\{\hat f(x_i)\} \), which is not \( f(x_i) \). The posterior
@eq-smo-posterior offers a different interval, and it is wider.

::: {#prp-smo-bayes-bands}
[Wahba–Nychka bands]

Let \( \bS_\lambda=\B(\B\T\B+\lambda\bP)^{-1}\B\T \) with diagonal entries
\( s_{ii} \).

::: {.enumerate options="label=(\alph*)"}
1. Under @thm-smo-mixed(c), \( f(x_i)=\B_{i\cdot}\bgamma \) has posterior
   distribution \( \Normal(\hat f(x_i),\sigma^2s_{ii}) \), giving the **Bayesian
   band**
   \[
   \hat f(x_i)\ \pm\ z_{1-\alpha/2}\,\hat\sigma\sqrt{s_{ii}} .
   \]{#eq-smo-band}

2. \( s_{ii}\ge(\bS_\lambda\bS_\lambda\T)_{ii} \), so @eq-smo-band is never
   narrower than the interval of @prp-smo-linear-smoother(e), and is strictly
   wider unless \( \bS_\lambda \) is a projection.

3. *(Stated, not proved.)* Treating \( f \) as fixed, the extra width in (b)
   accounts approximately for the squared bias, in the sense that
   \( n^{-1}\sum_i\sigma^2s_{ii} \) approximates
   \( n^{-1}\sum_i[\,\Var\{\hat f(x_i)\}+\text{bias}\{\hat f(x_i)\}^2\,] \) when
   \( \lambda \) is near its risk-optimal value. Consequently the *average* over
   the design points of the coverage of @eq-smo-band is close to \( 1-\alpha \),
   while the coverage at any single \( x_i \) is not: it is below nominal where
   \( f \) curves sharply and above it where \( f \) is nearly straight. This is
   **across-the-function** coverage (Wahba, 1983; Nychka, 1988).
:::

:::

::: {.proof}
(a) By @eq-smo-posterior, \( \B_{i\cdot}\bgamma \) is normal with mean
\( \B_{i\cdot}\hat{\bgamma}=\hat f(x_i) \) and variance
\( \sigma^2\B_{i\cdot}(\B\T\B+\lambda\bP)^{-1}\B_{i\cdot}\T=\sigma^2s_{ii} \).

(b) \( \bS_\lambda \) is symmetric with eigenvalues \( s_j\in(0,1] \)
by @prp-smo-pspline-ridge(b), so \( \bS_\lambda-\bS_\lambda^2 \) is symmetric with
eigenvalues \( s_j(1-s_j)\ge0 \) and hence nonnegative definite; its \( i \)th
diagonal entry \( s_{ii}-(\bS_\lambda\bS_\lambda\T)_{ii} \) is therefore
nonnegative, and zero for all \( i \) only if the matrix is zero, that is only if
every \( s_j\in\{0,1\} \).

Part (c) is not proved here; see the notes.
:::

::: {.warning}
"Across-the-function" is a weaker guarantee than it sounds. It says that if you
draw the band and count how often it contains the true curve across the design
points, the fraction is about \( 1-\alpha \). It does not say that the band
contains the curve *everywhere* with probability \( 1-\alpha \) — that is a
simultaneous band, which is wider (Sun and Loader, 1994) — nor that the interval
at any one chosen point has the nominal coverage. If a particular \( x \) is the
question, @eq-smo-band is not the answer to it.
:::

::: {#exm-smo-coverage}
[Measuring the coverage]

Take \( n=120 \) equally spaced points, the mean function
of @exm-smo-bandwidth and \( \sigma=0.25 \), with a
cubic P-spline on \( K=20 \) knots, so \( m=24 \)
coefficients split into \( k=2 \) fixed and
\( r=22 \) random ones. On one data set the
transformation of @thm-smo-mixed reproduces the penalized fit to machine
precision, REML chooses \( \lambda=5.01\times 10^{0} \),
which is \( 8.12 \) effective degrees of freedom, against
\( 7.10 \) for GCV; the REML estimate of
\( \sigma \) is \( 0.254 \).

Over \( 400 \) replicates, with \( \lambda \) re-estimated by
REML each time, the nominal \( 95 \) per cent band @eq-smo-band covers
\( f(x_i) \) at a fraction \( 0.954 \) of the design
points on average, the pointwise coverages running from
\( 0.900 \) at the narrow dip near
\( x=0.8 \), where the fit is most biased, to
\( 0.985 \) on the flat
stretches ([Figure 43.6.1](06-mixed-model-and-bayes.html#fig-smo-bands)). The
Bayesian standard errors are on average \( 1.12 \)
times the naive ones of @prp-smo-linear-smoother(e), whose bands cover only
\( 0.927 \) of the design points: a twelve per cent
widening buys back nearly three percentage points of average coverage, and what
is left over is the bias that no single width can equalize.
:::

::: {when-format="html"}
![**Figure 43.6.1.** (a) One data set, the REML fit, and the nominal \( 95 \) per
cent Bayesian band @eq-smo-band. (b) Coverage of that band at each design point
over 400 replicates, \( \lambda \) re-estimated each time; dashed line nominal,
dotted line the average over the design points.](bayes_bands.svg){#fig-smo-bands width=100%}
:::

::: {when-format="pdf"}
![(a) One data set, the REML fit, and the nominal \( 95 \) per
cent Bayesian band @eq-smo-band. (b) Coverage of that band at each design point
over 400 replicates, \( \lambda \) re-estimated each time; dashed line nominal,
dotted line the average over the design points.](bayes_bands.pdf){width=100%}
:::

```{.python .run #cell-mixed-bayes-transform}
import numpy as np


def f_true(x):
    """The mean function of Section 43.1."""
    return 2 * x + np.exp(-25 * (x - 0.35) ** 2) - 0.75 * np.exp(-50 * (x - 0.80) ** 2)


def uniform_knots(a, b, K, d):
    """K interior knots equally spaced in (a, b), continued d further steps at each end."""
    delta = (b - a) / (K + 1)
    return a + delta * np.arange(-d, K + 2 + d)


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for j in range(C.shape[1]):
            if kn[j + deg] > kn[j]:
                C[:, j] += (x - kn[j]) / (kn[j + deg] - kn[j]) * B[:, j]
            if kn[j + deg + 1] > kn[j + 1]:
                C[:, j] += (kn[j + deg + 1] - x) / (kn[j + deg + 1] - kn[j + 1]) * B[:, j + 1]
        B = C
    return B


def split_penalty(P, tol=1e-9):
    """Bases X_g of N(P) and Z_g of its complement with Z_g' P Z_g = I and X_g' P = 0."""
    e, U = np.linalg.eigh(P)
    pos = e > tol * e.max()
    L = U[:, pos] * np.sqrt(e[pos])                    # P = L L' with L of full column rank
    Xg = U[:, ~pos]                                    # the null space of the penalty
    Zg = L @ np.linalg.inv(L.T @ L)
    return Xg, Zg


n, sigma, K, d, k = 120, 0.25, 20, 3, 2
x = (np.arange(1, n + 1) - 0.5) / n
kn = uniform_knots(0.0, 1.0, K, d)
B = bspline_basis(x, kn, d)
m = B.shape[1]
Dk = np.diff(np.eye(m), n=k, axis=0)
P = Dk.T @ Dk
Xg, Zg = split_penalty(P)
X, Z = B @ Xg, B @ Zg
r = Z.shape[1]

rng = np.random.default_rng(431)
g = rng.normal(size=m)                                           # any coefficient vector
beta_u = np.linalg.solve(np.column_stack([Xg, Zg]), g)
print("the penalty, in the old and the new coordinates:", g @ P @ g, beta_u[k:] @ beta_u[k:])

y = f_true(x) + sigma * rng.normal(size=n)
lam = 1.0
fit_pen = B @ np.linalg.solve(B.T @ B + lam * P, B.T @ y)
C = np.block([[X.T @ X, X.T @ Z], [Z.T @ X, Z.T @ Z + lam * np.eye(r)]])   # Henderson's equations
sol = np.linalg.solve(C, np.concatenate([X.T @ y, Z.T @ y]))
fit_mix = X @ sol[:k] + Z @ sol[k:]
print("penalized fit versus mixed-model fit:", np.max(np.abs(fit_pen - fit_mix)))
```

```{.python .run #cell-mixed-bayes-reml}
Uz, sz, _ = np.linalg.svd(Z, full_matrices=False)     # the design is fixed, so do this once


def reml(y, lam):
    """Profiled REML criterion (to be minimized) for V = I + Z Z' / lam."""
    w = sz**2 / (lam + sz**2)
    def Vi(M):                                         # multiplication by V^{-1}
        M2 = M if M.ndim == 2 else M[:, None]
        out = M2 - Uz @ (w[:, None] * (Uz.T @ M2))
        return out if M.ndim == 2 else out[:, 0]

    logdetV = np.sum(np.log1p(sz**2 / lam))
    ViX = Vi(X)
    A = X.T @ ViX
    resid = Vi(y) - ViX @ np.linalg.solve(A, ViX.T @ y)
    s2 = y @ resid / (n - k)
    return 0.5 * ((n - k) * np.log(s2) + logdetV + np.linalg.slogdet(A)[1]), s2


def reml_lambda(y, grid):
    return grid[int(np.argmin([reml(y, l)[0] for l in grid]))]


grid = np.geomspace(1e-6, 1e6, 121)
lam_reml = reml_lambda(y, grid)
S_reml = B @ np.linalg.solve(B.T @ B + lam_reml * P, B.T)
df_reml = np.trace(S_reml)
gcvs = []
for l in grid:
    S = B @ np.linalg.solve(B.T @ B + l * P, B.T)
    gcvs.append(np.mean((y - S @ y) ** 2) / (1 - np.trace(S) / n) ** 2)
lam_gcv = grid[int(np.argmin(gcvs))]
df_gcv = np.trace(B @ np.linalg.solve(B.T @ B + lam_gcv * P, B.T))
print(f"REML: lambda = {lam_reml:.3g}, df = {df_reml:.2f}; GCV: df = {df_gcv:.2f}")
```

## REML against the criteria, and full Bayes

REML and GCV estimate different things. GCV aims at prediction error; REML
maximizes a likelihood under the assumption that \( f \) really is a draw from
the prior — an assumption nobody believes, since nobody believes the deviations
from a straight line are exchangeable normal increments. In practice the two
rarely differ by much: in @exm-smo-coverage REML chose
\( 8.12 \) effective degrees of freedom against
GCV's \( 7.10 \). But REML is the more stable of
the two, its criterion better determined near the optimum, and it almost never
produces the interpolating fits that GCV occasionally does (Reiss and Ogden,
2009). It has no answer when the model for the errors is wrong: with the AR(1)
errors of @exm-smo-correlated it overfits just as GCV does, unless the
correlation is put into \( \Cov(\be) \), which the mixed-model form at least
makes possible (@def-cls-covariance).

::: {.remark}
[Full Bayes]

@thm-smo-mixed(c) fixes \( \lambda \). A fully Bayesian treatment puts a prior on
\( \tau^2 \) and \( \sigma^2 \) too — conventionally inverse-gamma, which is
conditionally conjugate — and samples \( (\bgamma,\tau^2,\sigma^2) \) by the Gibbs
sampler of @def-prc-mcmc: each sweep draws \( \bgamma \) from @eq-smo-posterior,
then \( \tau^2 \) from an inverse-gamma with shape \( a+r/2 \) and scale
\( b+\bu\T\bu/2 \), then \( \sigma^2 \) likewise. The posterior for \( f \) then
averages over the smoothing parameter instead of conditioning on an estimate of
it, which is the cleanest answer to the selection problem of
[Section 43.5](05-choosing-lambda.html). It is also where the prior shows its
weakness: the posterior for \( \tau^2 \) is sensitive to the hyperparameters, as
@prp-prc-priors warns, because a variance component with \( r \) terms is not well
identified. [Chapter 44](../ch44-additive-models/index.html) uses this machinery
for several smooth terms at once.
:::

## What this book does not do with smoothing

Three threads end here rather than continue, and it is only fair to name them.

**Adaptive smoothing.** A single \( \lambda \) is a compromise when \( f \) is
rough in one place and flat in another; letting it vary with \( x \) — a second
smooth for \( \log\lambda \), free-knot splines (Denison, Mallick and Smith,
1998), wavelet thresholding (Donoho and Johnstone, 1994) — is a subject this book
does not enter.

**Simultaneous bands.** @eq-smo-band is pointwise, or at best
across-the-function; a band containing the whole curve with a stated probability
needs the tube formula of Sun and Loader (1994) or a bootstrap over the whole
fit.

**Shape constraints.** Monotone, convex or unimodal fits add linear inequality
constraints to @eq-smo-pspline, turning it into a quadratic program; Ramsay
(1988) treats monotone regression splines.

What does continue is the penalty. Everything after
[Section 43.2](02-b-splines.html) rests on one idea — a rich basis, a quadratic
penalty, a trace that counts what was spent, and a variance ratio that can be
estimated — and that idea is indifferent to how many covariates there are, to
whether the response is normal, and to which parameter of the response
distribution is modelled.
[Chapter 44](../ch44-additive-models/index.html) adds covariates and
[Chapter 45](../ch45-quantile-gamlss/index.html) parameters, both by writing down
more terms of the form @eq-smo-pspline.

## Exercises

### A. Check your understanding

::: {#exr-smo-variance-ratio}
[A1]

Using @eq-smo-lambda-ratio, say what \( \sigma_u^2\to0 \) and
\( \sigma_u^2\to\infty \) correspond to
in @prp-smo-pspline-ridge, and what a REML estimate \( \hat\sigma_u^2=0 \) would mean
about the fitted curve. Is that estimate on the boundary of the parameter space,
and does @thm-mix-boundary apply?
:::

::: {.solution}
\( \sigma_u^2\to0 \) gives \( \lambda\to\infty \): the random part is annihilated
and only the fixed part, the straight line, survives.
\( \sigma_u^2\to\infty \) gives \( \lambda\to0 \), the unpenalized spline fit. A
REML estimate of zero means the data give no evidence of departure from a
straight line, and the reported fit *is* the least squares line. It is on the
boundary, so the null distribution of a likelihood ratio test of
\( \sigma_u^2=0 \) is not the usual \( \chi^2 \); @thm-mix-boundary gives the
mixture, and the same caution applies to testing whether a smooth term is needed
at all.
:::

### B. Practice

::: {#exr-smo-reml-vs-gcv}
[B1]

Run the simulation of @exm-smo-criteria again with \( \lambda \) chosen by REML
instead, and compare the median selected degrees of freedom, the average squared
error, and the upper tail of the distribution of df, with the values for CV, GCV
and AIC\( _c \). Which criterion is safest, and at what cost?
:::

::: {#exr-smo-coverage-sigma}
[B2]

@exm-smo-coverage uses \( \hat\sigma \) from @eq-smo-sigma-hat in @eq-smo-band.
Repeat the coverage simulation with \( \sigma \) known and with \( \hat\sigma \)
estimated, and with \( z_{1-\alpha/2} \) replaced by a \( t \) quantile on
\( n-2\tr(\bS)+\tr(\bS\bS\T) \) degrees of freedom. How much of the departure from
nominal coverage is due to estimating \( \sigma \), and how much to the bias?
:::

### C. Going deeper

::: {#exr-smo-blup-shrinkage}
[C1]

Show that in the mixed model of @thm-smo-mixed(b), the BLUP satisfies
\( \hat{\bu}=(\Z\T\Z+\lambda\I)^{-1}\Z\T(\y-\X\hat{\bbeta}) \), and hence that
\( \hat{\bu} \) is a ridge estimate (@thm-shr-ridge) on the residuals from the
fixed part. Explain, using @thm-mix-blup, why the shrinkage factor is
\( \sigma_u^2/(\sigma_u^2+\sigma^2/d_j^2) \) in the \( j \)th singular direction of
\( \Z \), and what happens to it as \( n \) grows with \( r \) fixed.
:::

::: {.solution}
The second block of Henderson's equations is
\( \Z\T\X\hat{\bbeta}+(\Z\T\Z+\lambda\I)\hat{\bu}=\Z\T\y \), which rearranges to
the displayed formula; it is @eq-shr-ridge-criterion with response
\( \y-\X\hat{\bbeta} \) and ridge parameter \( \lambda=\sigma^2/\sigma_u^2 \). With
\( \Z=\sum_jd_j\bu_j\bv_j\T \), @thm-shr-ridge(a) gives the factor
\( d_j^2/(d_j^2+\lambda)=\sigma_u^2/\{\sigma_u^2+\sigma^2/d_j^2\} \), which is
@eq-mix-shrinkage: the random effect is shrunk according to how much information
the data carry about it. As \( n \) grows with the basis fixed, every
\( d_j^2 \) grows like \( n \) and the factors tend to one, which is why a
low-rank smoother with a fixed basis is asymptotically a parametric fit unless
\( m \) grows with \( n \).
:::

::: {#exr-smo-improper-prior}
[C2]

The prior of @thm-smo-mixed(c) is improper on \( \Null(\bP) \). Show that the
posterior @eq-smo-posterior is nevertheless proper, provided
\( \B\X_\gamma \) has full column rank. By @exr-smo-rank-deficient that holds for
\( k\le2 \) whenever the design has at least \( k \) distinct points, however the
knots fall; show that it fails when it has fewer — for \( k=2 \), a design
concentrated at a single \( x \) — and say what the posterior then does in the
unidentified direction. Compare with the discussion of improper priors in
[Section 39.5](../ch39-glms-in-practice-bayes/05-priors-as-regularizers.html).
:::
