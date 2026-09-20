# Iteratively reweighted least squares

The likelihood equations @eq-glm-score have no closed-form solution. Any general-purpose
optimizer will solve them, but there is a better way, and it is why generalized linear
models became practical: the natural Newton-type iteration is a *weighted least squares fit*,
repeated. The apparatus of [Chapter 6](../ch06-projections/index.html) and
[Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html)
is therefore reusable, including the stable routines of
[Chapter 10](../ch10-computation/index.html).

## Fisher scoring is weighted least squares

Newton's method for maximizing \( \ell \) replaces the observed information by its inverse:
\( \bbeta^{(t+1)}=\bbeta^{(t)}+\{-\nabla^2\ell\}^{-1}\bU \). **Fisher scoring** replaces the
observed information by the expected information \( \boldsymbol{\mathcal I} \) instead.
By @thm-glm-score(d) the two are the same under the canonical link; in general they differ, and
the expected version is both simpler and guaranteed nonnegative definite.

::: {#thm-glm-irls}
[Fisher scoring is iteratively reweighted least squares]

Let \( \bbeta^{(t)} \) be in the interior of \( \mathcal B \) with \( \X\T\W^{(t)}\X \) nonsingular,
where \( \W^{(t)} \) is the working-weight matrix of @thm-glm-score evaluated at
\( \bbeta^{(t)} \). Define the **working response**
\[
z_i^{(t)}=\eta_i^{(t)}+\frac{y_i-\mu_i^{(t)}}{h'(\eta_i^{(t)})}
=\eta_i^{(t)}+\bigl(y_i-\mu_i^{(t)}\bigr)\,g'\bigl(\mu_i^{(t)}\bigr).
\]{#eq-glm-working-response}

Then the Fisher scoring step
\( \bbeta^{(t+1)}=\bbeta^{(t)}+\boldsymbol{\mathcal I}(\bbeta^{(t)})^{-1}\bU(\bbeta^{(t)}) \)
is exactly the weighted least squares estimate
\[
\bbeta^{(t+1)}=\bigl(\X\T\W^{(t)}\X\bigr)^{-1}\X\T\W^{(t)}\bz^{(t)} ,
\]{#eq-glm-irls}

that is, the regression of the working response \( \bz^{(t)} \) on \( \X \) with weights
\( \W^{(t)} \). The dispersion \( \phi \) cancels and never enters the iteration.
:::

::: {.proof}
By @eq-glm-information and @eq-glm-score, with everything evaluated at \( \bbeta^{(t)} \),
\[
\begin{aligned}
\boldsymbol{\mathcal I}^{-1}\bU
&=\phi\bigl(\X\T\W\X\bigr)^{-1}\cdot\frac1\phi\X\T\W\bD^{-1}(\y-\bmu)\\
&=\bigl(\X\T\W\X\bigr)^{-1}\X\T\W\bD^{-1}(\y-\bmu),
\end{aligned}
\]
so \( \phi \) has already cancelled. Adding \( \bbeta^{(t)}=(\X\T\W\X)^{-1}\X\T\W\X\bbeta^{(t)} \)
and using \( \X\bbeta^{(t)}=\boldsymbol{\upeta}^{(t)} \),
\[
\bbeta^{(t+1)}=\bigl(\X\T\W\X\bigr)^{-1}\X\T\W\bigl\{\boldsymbol{\upeta}^{(t)}+\bD^{-1}(\y-\bmu)\bigr\},
\]
and the braced vector is \( \bz^{(t)} \) by @eq-glm-working-response, since
\( \bD^{-1}=\diag\{1/h'(\eta_i)\} \) and \( g'(\mu)=1/h'(\eta) \).
:::

::: {.idea}
The working response is a first-order Taylor expansion of the link applied to the data:
\( g(y_i)\approx g(\mu_i)+(y_i-\mu_i)g'(\mu_i)=z_i \). So each step regresses "the data,
mapped to the linear-predictor scale" on \( \X \). The working weight is the reciprocal of
that quantity's approximate variance, because
\[
\Var(z_i)\approx g'(\mu_i)^2\Var(Y_i)=\frac{\phi}{w_i}\frac{V(\mu_i)}{h'(\eta_i)^2}
=\frac{\phi}{\W_{ii}} .
\]
Both the response and the weights depend on the current fit, which is why the regression
has to be iterated.
:::

Written out, the algorithm is short.

::: {.algorithm}
**Iteratively reweighted least squares.**

1. Choose starting fitted means \( \mu_i^{(0)} \) and set \( \eta_i^{(0)}=g(\mu_i^{(0)}) \).
2. For \( t=0,1,2,\dots \):
   compute \( \W_{ii}=w_ih'(\eta_i^{(t)})^2/V(\mu_i^{(t)}) \) and \( z_i^{(t)} \) from
   @eq-glm-working-response; solve the weighted least squares problem @eq-glm-irls for
   \( \bbeta^{(t+1)} \); set \( \boldsymbol{\upeta}^{(t+1)}=\X\bbeta^{(t+1)} \) and
   \( \mu_i^{(t+1)}=h(\eta_i^{(t+1)}) \).
3. Stop when \( \bbeta \) (or the deviance of
   [Section 34.5](05-deviance-and-inference.html)) stops changing; report
   \( \hbeta=\bbeta^{(t+1)} \) and \( (\X\T\hat{\W}\X)^{-1} \), with \( \hat{\W} \) the weights at
   the solution.
:::

Because the map is a fixed-point iteration whose fixed points solve @eq-glm-score, any limit
is a solution of the likelihood equations; under the canonical link with
\( \rank(\X)=p \), @thm-glm-concave then makes it *the* maximum likelihood estimate.

**Starting values.** The usual device is to start from the data themselves: set
\( \mu_i^{(0)}=y_i \), adjusted if necessary to keep it inside \( \mathcal M \). For a Poisson
response take \( \mu_i^{(0)}=y_i+0.1 \); for a binomial proportion,
\( \mu_i^{(0)}=(w_iy_i+\tfrac12)/(w_i+1) \), the empirical-logit correction
of @exr-tr-empirical-logit; for a gamma response, \( \mu_i^{(0)}=y_i \). The first step is then
the one-step weighted least squares fit that preceded generalized linear models
historically, and the iteration refines it.

**Relation to the least squares chapters.** Each step of @eq-glm-irls is the weighted least
squares problem of @thm-het-wls with *known* weights, so its convergence theory is numerical,
not statistical: the final \( (\X\T\hat{\W}\X)^{-1} \) is the inverse expected information
of the likelihood, not the covariance of a fixed-weight estimate, and the two coincide
asymptotically for a different reason (@thm-glm-asymptotics). Compute the step by the QR or
Cholesky routines of [Section 6.10](../ch06-projections/10-computation.html) applied to
\( \W^{1/2}\X \) and \( \W^{1/2}\bz \), never by inverting \( \X\T\W\X \).

::: {.warning}
The name "iteratively reweighted least squares" is also used for the robust M-estimation
algorithm of @prp-res-irls, and the two are different. There the weights come from a
\( \rho \) function chosen to resist outliers, and the iteration provably decreases the
objective at every step. Here the weights come from the variance function, the iteration is
Fisher scoring, and there is no monotonicity guarantee unless the link is canonical.
:::

## The algorithm in practice

The listing below implements @eq-glm-irls in a dozen lines. A family is nothing but
\( (V,g,h) \) together with \( h' \) and a starting rule; the solver never needs to know which
distribution it is fitting, which is the computational content of the whole chapter.

```{.python .run #cell-irls-families}
import numpy as np
from scipy import stats

# A family is the variance function V, the link g, its inverse h, the derivative
# h' = d(mu)/d(eta), and a rule for the starting fitted means.
def poisson_log():
    return dict(name="Poisson/log", V=lambda mu: mu, g=np.log, h=np.exp,
                dmu=lambda eta: np.exp(eta), start=lambda y, w: y + 0.1)

def binomial_logit():
    return dict(name="binomial/logit", V=lambda mu: mu * (1 - mu),
                g=lambda mu: np.log(mu / (1 - mu)), h=lambda eta: 1 / (1 + np.exp(-eta)),
                dmu=lambda eta: 1 / (2 + np.exp(eta) + np.exp(-eta)),
                start=lambda y, w: (w * y + 0.5) / (w + 1.0))

def binomial_probit():
    return dict(name="binomial/probit", V=lambda mu: mu * (1 - mu),
                g=stats.norm.ppf, h=stats.norm.cdf, dmu=stats.norm.pdf,
                start=lambda y, w: (w * y + 0.5) / (w + 1.0))

def gamma_log():
    return dict(name="gamma/log", V=lambda mu: mu ** 2, g=np.log, h=np.exp,
                dmu=lambda eta: np.exp(eta), start=lambda y, w: y)

def irls(X, y, family, w=None, tol=1e-12, max_iter=60, trace=False):
    """Fisher scoring for a generalized linear model, as weighted least squares.

    At each step the working response is z = eta + (y - mu)/h'(eta) and the working
    weight is W = w h'(eta)^2 / V(mu); beta is the weighted least squares coefficient
    of z on X. The dispersion phi never enters, because it cancels.
    """
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    mu = family["start"](y, w)
    eta = family["g"](mu)
    beta, path = np.zeros(X.shape[1]), []
    for _ in range(max_iter):
        d = family["dmu"](eta)                             # h'(eta)
        W = w * d ** 2 / family["V"](mu)                   # working weights
        z = eta + (y - mu) / d                             # working response
        XtW = X.T * W
        beta_new = np.linalg.solve(XtW @ X, XtW @ z)       # one weighted least squares fit
        path.append(beta_new)
        eta = X @ beta_new
        mu = family["h"](eta)
        if np.max(np.abs(beta_new - beta)) < tol * (1 + np.max(np.abs(beta_new))):
            beta = beta_new
            break
        beta = beta_new
    W = w * family["dmu"](eta) ** 2 / family["V"](mu)      # weights at the solution
    cov_unscaled = np.linalg.inv((X.T * W) @ X)            # (X^T W X)^{-1}
    return (beta, cov_unscaled, np.array(path)) if trace else (beta, cov_unscaled)
```

```{.python .run #cell-irls-check}
import statsmodels.api as sm

rng = np.random.default_rng(3401)
n = 300
Xs = np.column_stack([np.ones(n), rng.normal(size=n), rng.binomial(1, 0.4, n)])
beta_true = np.array([0.7, 0.5, -0.4])
counts = rng.poisson(np.exp(Xs @ beta_true))

beta_hat, cov = irls(Xs, counts.astype(float), poisson_log())
sm_fit = sm.GLM(counts, Xs, family=sm.families.Poisson()).fit()
print("from scratch:", np.round(beta_hat, 6))
print("statsmodels :", np.round(sm_fit.params, 6))
print("standard errors:", np.round(np.sqrt(np.diag(cov)), 6), np.round(sm_fit.bse, 6))
```

The script behind this section runs the same solver on four models — Poisson/log,
binomial/logit, binomial/probit and gamma/log — and checks the coefficients and standard
errors against `statsmodels` to at least eight decimal places in every case. Two details
matter. For the binomial fits the response passed to the solver is the *proportion*
\( s_i/m_i \) with prior weight \( w_i=m_i \), which is @exm-glm-binomial in code. For the
gamma fit the solver returns \( (\X\T\hat{\W}\X)^{-1} \), so the standard errors need the
extra factor \( \hat\phi \), estimated separately; on Engel's data that
gives \( \hbeta=(0.5068,\,0.8629) \) with
\( \hat\phi=0.0177 \) and a standard error of
\( 0.0198 \) for the elasticity.

## Convergence and its failures

::: {#exm-glm-convergence}
[Two convergence rates]

[Figure 34.4.1](04-irls.html#fig-glm-convergence) records
\( \max_j\lvert\beta_j^{(t)}-\hat\beta_j\rvert \) against the iteration count for the four
models. The canonical-link fits reach full double precision in
\( 5 \) and \( 4 \) steps, with error curves that bend
downwards on the logarithmic scale: the correct digits roughly double at each step. The
non-canonical fits take \( 5 \) and \( 7 \)
steps, along straight lines, so their correct digits grow by a fixed amount per step.

This is what the theory predicts. Under the canonical link Fisher scoring *is* Newton's
method (@thm-glm-score(d)), quadratically convergent near a nondegenerate maximum. Under any
other link the two information matrices differ by \( \phi^{-1}\X\T\bD_r\X \), which vanishes
only in expectation, so the iteration is Newton's method with a perturbed Hessian and
converges linearly — in practice, one or two extra iterations.
:::

::: {when-format="html"}
![**Figure 34.4.1.** Convergence of iteratively reweighted least squares for four models.
Solid curves are canonical links, dashed curves are not. On a logarithmic scale a quadratic
rate bends downwards and a linear rate is
straight.](irls_convergence.svg){#fig-glm-convergence width=72%}
:::

::: {when-format="pdf"}
![Convergence of iteratively reweighted least squares for four models.
Solid curves are canonical links, dashed curves are not. On a logarithmic scale a quadratic
rate bends downwards and a linear rate is
straight.](irls_convergence.pdf){width=72%}
:::

Four things go wrong often enough to be worth recognizing.

**The estimate does not exist.** When the maximum of @exm-glm-no-maximum is not attained
the iteration cannot converge, and its behaviour is characteristic: the coefficients grow
without limit, the deviance falls towards its infimum, the working weights collapse to zero.
Software that stops on a small change in the deviance then reports a large but finite
estimate with a large standard error and no warning.

::: {#exm-glm-separation-irls}
[What separation looks like from inside the algorithm]

Ten binary responses were generated with \( y_i=1 \) exactly when \( x_i>0 \), so the data are
completely separated. Running the solver on them gives

| step | \( \hat\beta_1 \) | deviance | largest working weight |
|---:|---:|---:|---:|
| 1 | \( 1.36 \) | \( 4.589 \) | \( 0.245 \) |
| 4 | \( 5.74 \) | \( 0.931 \) | \( 0.159 \) |
| 7 | \( 16.95 \) | \( 0.0574 \) | \( 0.0142 \) |
| 8 | — | — | — |

At the eighth step every fitted probability has reached \( 0 \) or \( 1 \) in double precision,
so \( V(\hat\mu_i)=0 \), the working weights are \( 0/0 \), and the iteration returns nothing.
The slope was never going to settle: the supremum of the likelihood is approached only as
\( \beta_1\to\infty \). [Chapter 35](../ch35-binary-responses/index.html) characterizes
exactly when this happens and what to do about it (@thm-bin-separation).
:::

```{.python .run #cell-irls-separation}
x_sep = np.array([-2.4, -1.7, -1.1, -0.6, -0.2, 0.3, 0.8, 1.4, 1.9, 2.6])
y_sep = (x_sep > 0).astype(float)                        # the two groups do not overlap
X_sep = np.column_stack([np.ones(10), x_sep])

for iters in range(1, 9):
    with np.errstate(divide="ignore", invalid="ignore"):
        b, _ = irls(X_sep, y_sep, binomial_logit(), tol=0.0, max_iter=iters)
        mu = 1 / (1 + np.exp(-(X_sep @ b)))
        dev = -2 * np.sum(np.where(y_sep == 1, np.log(mu), np.log(1 - mu)))
    print(f"{iters} steps: slope {b[1]:9.4f}   deviance {dev:.3e}"
          f"   largest working weight {np.max(mu * (1 - mu)):.3e}")
```

**The step leaves the admissible set.** With a link whose range is not all of \( \Real \), a
full step can produce \( \eta_i^{(t+1)}\notin g(\mathcal M) \) — a negative fitted mean for an
identity-link Poisson model, where the variance function is then negative and the "weights"
meaningless. The repair is **step halving**: if the proposed \( \bbeta^{(t+1)} \) is
inadmissible, or the deviance increases, halve the step and retry. A few halvings almost
always suffice, and software applies them automatically.

**The information matrix is nearly singular.** The matrix inverted at each step is
\( \X\T\W\X \), not \( \X\T\X \), so full column rank is no longer enough: if the weights are
tiny wherever a regressor varies, the weighted matrix can be numerically singular even for a
well-conditioned \( \X \). The diagnostics of
[Chapter 26](../ch26-collinearity/index.html) apply to \( \W^{1/2}\X \), not \( \X \).

**The iteration converges to the wrong stationary point.** For non-canonical links the
log-likelihood need not be concave, and Fisher scoring from a bad start can settle at a
local maximum; the defence is several starting values compared by maximized
log-likelihood.

## Exercises

### A. Check your understanding

::: {#exr-glm-working-variance}
[A1]

Verify the approximation \( \Var(z_i)\approx\phi/\W_{ii} \) stated in the key-idea box, and
explain which two approximations are being made.
:::

### B. Practice

::: {#exr-glm-one-step}
[B1]

Take grouped binomial data with all \( 0<s_i<m_i \) and start the iteration at
\( \mu_i^{(0)}=s_i/m_i \). Show that the first working response is the empirical logit
\( \log\{s_i/(m_i-s_i)\} \) and the first working weight is \( s_i(m_i-s_i)/m_i \). Compare
with @exr-tr-empirical-logit.
:::

::: {.solution}
With \( \mu^{(0)}_i=s_i/m_i \) the starting linear predictor is already
\( \eta^{(0)}_i=\log\{\mu^{(0)}_i/(1-\mu^{(0)}_i)\}=\log\{s_i/(m_i-s_i)\} \), and
\( y_i-\mu^{(0)}_i=0 \), so \( z^{(0)}_i=\eta^{(0)}_i \) is the empirical logit. The working
weight is \( \W_{ii}=m_i\mu^{(0)}_i(1-\mu^{(0)}_i)=s_i(m_i-s_i)/m_i \), which is the weight
of @exr-tr-empirical-logit with \( \hat p_i=s_i/m_i \). The first iterate is therefore exactly
the empirical-logit weighted least squares estimate; later iterates are not, which is why
the two methods give different answers in small samples.
:::

::: {#exr-glm-irls-offset}
[B2]

Modify @eq-glm-irls for a model with an offset
\( \eta_i=o_i+\x_{(i)}\T\bbeta \), as in @exr-glm-offset. Show that the working response becomes
\( z_i=\eta_i-o_i+(y_i-\mu_i)g'(\mu_i) \) and nothing else changes.
:::

### C. Going deeper

::: {#exr-glm-irls-monotone}
[C1]

Show that under the canonical link, Fisher scoring is Newton's method applied to a strictly
concave function, and deduce that when the maximum exists, the iteration converges to it
from any starting point at which the step is well defined *provided the step is damped*: for
each \( t \) there is \( s_t\in(0,1] \) with \( \ell(\bbeta^{(t)}+s_t\mathbf{d}^{(t)})>\ell(\bbeta^{(t)}) \),
where \( \mathbf{d}^{(t)} \) is the Newton direction. Why does the undamped step not carry the
same guarantee? (Compare @prp-res-irls, where the undamped step *is* monotone, because
the argument there is majorization rather than Newton.)
:::

::: {.solution}
By @thm-glm-concave the Hessian \( -\boldsymbol{\mathcal I} \) is negative definite, so the
Newton direction \( \mathbf{d}=\boldsymbol{\mathcal I}^{-1}\bU \) satisfies
\( \bU\T\mathbf{d}=\bU\T\boldsymbol{\mathcal I}^{-1}\bU>0 \) whenever \( \bU\ne\bzero \): it is an
ascent direction, so the directional derivative of \( \ell \) at \( \bbeta^{(t)} \) along
\( \mathbf{d} \) is positive and some small enough step increases \( \ell \). The undamped step
has length fixed by the local quadratic approximation, which can overshoot badly when \( \ell \)
is far from quadratic, so it may decrease \( \ell \) or leave \( \mathcal B \) altogether.
:::

