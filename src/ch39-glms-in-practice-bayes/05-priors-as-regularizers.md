# Priors as regularizers

@prp-reg-map already said that every penalty is a prior and every prior a penalty: the
maximizer of a penalized log-likelihood is the mode of a posterior whose log prior is
minus the penalty. This section takes that identity in the other direction and asks what a
prior *does* in a generalized linear model. The answer is more than shrinkage: a prior
carrying almost no information can still be the difference between an estimate and no
estimate at all.

## Two correspondences and a warning about scale

::: {#prp-prc-priors}
[What a prior does]

Let \( \ell \) be the log-likelihood of a generalized linear model with
\( \rank(\X)=p \) and a concave log-likelihood (@thm-glm-concave), and write
\( \bbeta=(\beta_0,\bgamma\T)\T \) with \( \beta_0 \) an intercept.

::: {.enumerate options="label=(\alph*)"}
1. *(Normal prior, ridge penalty.)* Under
   \( \gamma_j\sim\Normal(0,\tau^2) \) independently and a flat prior on \( \beta_0 \), the
   posterior mode maximizes \( \ell(\bbeta)-\tfrac12\lambda\norm{\bgamma}^2 \) with
   \( \lambda=1/\tau^2 \). If \( \ell \) is bounded above, as it is for the binomial
   and the Poisson, and the response is not constant — for the binomial, both outcomes
   occur — the mode exists and is unique whether or not the data are separated. The
   hypothesis on the response cannot be dropped: with a flat prior on \( \beta_0 \) the
   objective is not coercive in \( \beta_0 \) when every \( y_i \) is \( 1 \).
   Penalizing the intercept as well makes existence unconditional.

2. *(Laplace prior, lasso penalty.)* Under \( \gamma_j \) independent with density
   \( (2b)^{-1}e^{-|\gamma_j|/b} \) and a flat prior on \( \beta_0 \), the posterior mode
   maximizes \( \ell(\bbeta)-\lambda\sum_{j}|\gamma_j| \) with \( \lambda=1/b \), and is
   therefore the lasso-penalized estimate of @def-shr-lasso, with exact zeros. The
   *posterior*, by contrast, is never sparse: whenever it is proper,
   \( \Pr(\gamma_j=0\mid\y)=0 \) for every \( j \).

3. *(What a normal prior says about a probability.)* In a binary model with logit link,
   suppose the linear predictor at some covariate value has the prior
   \( \eta\sim\Normal(0,\sigma^2) \). The induced prior density of
   \( \pi=e^{\eta}/(1+e^{\eta}) \) on \( (0,1) \) is unimodal with mode \( 1/2 \) if
   \( \sigma^2\le2 \), and has exactly two modes, symmetric about \( 1/2 \), with a local
   minimum at \( 1/2 \), if \( \sigma^2>2 \).

4. *(Existence.)* Under the prior of (a), or under independent Cauchy priors on the
   \( \gamma_j \) with a flat prior on \( \beta_0 \), the logistic posterior mode exists
   for every data set in which both outcomes occur, separated or not; under (a) it is also
   unique. Under any proper prior the posterior
   itself is proper (@prp-prc-propriety(a)), so posterior means and credible intervals
   exist whether or not the maximum likelihood estimate does.
:::
:::

::: {.proof}
(a) The log posterior is \( \ell(\bbeta)-\norm{\bgamma}^2/(2\tau^2) \) up to a constant,
which is the stated penalized log-likelihood with \( \lambda=1/\tau^2 \); this is
@prp-reg-map for the normal prior. Existence and uniqueness for the logistic model with
the intercept left out of the penalty are @exr-bin-ridge-exists, under exactly the stated
hypothesis; the same argument — a concave \( \ell \) bounded above, coercive in
\( \beta_0 \) at fixed \( \bgamma \), plus a strictly concave coercive quadratic —
applies to any family whose log-likelihood is concave and bounded above. When
\( \bgamma \) is penalized and \( \beta_0 \) is not, and every \( y_i \) equals
\( 1 \), the supremum \( 0 \) is approached only along \( \bgamma=\bzero \),
\( \beta_0\to+\infty \) and is not attained.

(b) The log prior is \( -\sum_j|\gamma_j|/b \) up to a constant, giving the penalty with
\( \lambda=1/b \); this is @prp-reg-map for the Laplace prior. For the second claim, the
unnormalized posterior density is positive and continuous on \( \Real^p \), so when it
integrates the marginal posterior of each \( \gamma_j \) is absolutely continuous and
\( \Pr(\gamma_j=0\mid\y)=0 \). The sparsity of the lasso is a property of the *mode* of a
density with a kink, not of the distribution; @exm-reg-laplace-mean makes the same point in
the orthonormal linear model. Note what is *not* claimed: the posterior mean can be
exactly zero, as it is whenever the posterior happens to be symmetric in \( \gamma_j \)
about the origin.

(c) Write \( t=\log\{\pi/(1-\pi)\} \), so \( dt/d\pi=1/\{\pi(1-\pi)\}>0 \). The induced
density satisfies
\[
f(\pi)\ \propto\ \exp\Bigl(-\frac{t^2}{2\sigma^2}\Bigr)\frac{1}{\pi(1-\pi)},
\]
with a constant of proportionality free of \( \pi \). Since \( t\mapsto\pi \) is
strictly increasing, the stationary points of \( f \) in \( \pi \) are the stationary
points of \( \psi(t)=\log f(\pi(t)) \) in \( t \). Using
\( \pi(1-\pi)=e^{t}/(1+e^{t})^2 \),
\[
\begin{aligned}
\psi(t)&=-\frac{t^2}{2\sigma^2}+2\log(1+e^{t})-t+\text{const},\\
\psi'(t)&=-\frac{t}{\sigma^2}+\frac{2e^t}{1+e^t}-1
=\tanh(t/2)-\frac{t}{\sigma^2}.
\end{aligned}
\]
Call the last expression \( G(t) \); it is odd and \( G(0)=0 \). Now
\( G'(t)=\tfrac12\operatorname{sech}^2(t/2)-1/\sigma^2 \), which is even, strictly
decreasing in \( |t| \), and has maximum \( \tfrac12-1/\sigma^2 \) at \( t=0 \).

If \( \sigma^2\le2 \) then \( G'\le0 \) everywhere, with equality only at \( t=0 \) in the
boundary case, so \( G \) is strictly decreasing and \( t=0 \) is its only root: \( f \) has
one stationary point. It is a maximum because \( f(\pi)\to0 \) as \( \pi\to0 \) or
\( \pi\to1 \), which follows from \( \exp\{-t^2/(2\sigma^2)+|t|\}\to0 \).

If \( \sigma^2>2 \) then \( G'(0)>0 \) and \( G' \) has exactly two zeros, so \( G \) has
at most three roots; it has at least three, because \( G(0)=0 \), \( G>0 \) just to the
right of zero while \( G(t)\to-\infty \) as \( t\to\infty \), and by oddness on the left.
So there are roots \( t_-<0=t_0<t_+ \) with \( t_-=-t_+ \). Since \( f\to0 \) at both ends
and \( f \) is continuous and positive, the outer two are maxima and the middle one, at
\( \pi=1/2 \), is a minimum.

(d) With the normal prior the objective is coercive and strictly concave, as in (a). With
independent Cauchy priors of scales \( s_j \), the log posterior is
\( \ell(\bbeta)-\sum_j\log\{1+(\gamma_j/s_j)^2\} \) up to a constant. For the logistic
model \( \ell\le0 \), and the penalty tends to \( -\infty \) as \( \norm{\bgamma}\to\infty \);
if \( \bgamma \) stays bounded while \( |\beta_0|\to\infty \) then \( \ell\to-\infty \),
provided both outcomes occur, by the argument of @exr-bin-ridge-exists. So the objective is
coercive and continuous and attains its maximum on a compact sublevel set. It need not be
unique: the Cauchy log density is not concave. The last sentence is
@prp-prc-propriety(a).
:::

Part (c) is the warning. A prior that is "noninformative" about a coefficient is a strong
statement about a probability: a \( \Normal(0,10^2) \) prior on the linear predictor puts
half of its mass on probabilities below \( 0.0012 \) or above
\( 1-0.0012 \), saying before any data that the event is nearly
certain or nearly impossible and that the analyst does not know which. Panel (a) of
[Figure 39.5.1](05-priors-as-regularizers.html#fig-prc-priors) shows the three shapes: a
gentle hump at \( \sigma=1 \), flat at the boundary \( \sigma=\sqrt2 \), a deep U at
\( \sigma=4 \).

This is why the default priors recommended for logistic regression are much tighter than
"vague". Gelman, Jakulin, Pittau and Su (2008) rescale every covariate to mean zero and
standard deviation \( 1/2 \), then put an independent Cauchy prior of scale
\( 2.5 \) on each coefficient and \( 10 \) on the intercept. On that scale a coefficient
of \( 5 \) multiplies the odds by \( e^{2.5}\approx12 \) for a one standard deviation
change, beyond anything social or biomedical data produce; the Cauchy's heavy tail
nonetheless leaves a genuinely large coefficient nearly alone, which a normal prior of the
same scale would not.

## Priors as penalties, checked

::: {.remark}
[Scaling is not optional]

Both correspondences put the same \( \lambda \) on every coefficient, which is a
statement about the units of the covariates: multiplying \( x_j \) by ten divides
\( \gamma_j \) by ten and changes its penalty by a factor of a hundred. Every penalized
or weakly-informative fit therefore standardizes first (@prp-shr-centring,
@prp-reg-rescaled). A prior of fixed scale is meaningful only once the covariates have one.
:::

On the \( 944 \) respondents of @exm-bin-anes, with the four covariates
standardized, maximizing \( \ell(\bbeta)-\tfrac12\lambda\norm{\bgamma}^2 \) by Newton's
method and maximizing \( \ell(\bbeta)-\lambda\sum_j|\gamma_j| \) by coordinate descent
gives the two paths. At \( \lambda=64 \), which is a normal prior with
\( \tau=1/8 \), the ridge estimates of the four slopes are
\( 1.465 \), \( 0.105 \),
\( 0.032 \) and \( 0.153 \); the lasso
estimates are \( 1.773 \) and three exact zeros. Both have shrunk the
party-identification coefficient from its unpenalized value of
\( 2.769 \) on this standardized scale, the ridge further than the
lasso, because the lasso has removed the competitors and the party term absorbs what
they explained.

```{.python .run #cell-priors-penalties}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats
anes = sm.datasets.anes96.load_pandas().data
Xf = np.column_stack([anes["PID"], anes["age"], anes["educ"], anes["income"]])
Xf = (Xf - Xf.mean(0)) / Xf.std(0)                      # standardize, as any penalty requires
Xf = np.column_stack([np.ones(len(Xf)), Xf])
yf = anes["vote"].to_numpy(float)


def map_estimate(penalty, lam, steps=30):
    """Posterior mode: maximize the log-likelihood minus lam times a penalty on the slopes.

    A normal prior gives the ridge penalty and is solved by Newton's method; a Laplace
    prior gives the lasso penalty and is solved by coordinate descent on the working
    least squares problem of each Newton step.
    """
    beta = np.zeros(Xf.shape[1])
    for _ in range(steps):
        pi = special.expit(Xf @ beta)
        w = np.maximum(pi * (1 - pi), 1e-8)
        z = Xf @ beta + (yf - pi) / w                   # the working response
        if penalty == "ridge":
            P = np.diag([0.0] + [1.0] * (Xf.shape[1] - 1))
            beta = np.linalg.solve(Xf.T @ (w[:, None] * Xf) + lam * P, Xf.T @ (w * z))
            continue
        for _ in range(40):                             # coordinate descent for the lasso
            for j in range(Xf.shape[1]):
                r = z - Xf @ beta + Xf[:, j] * beta[j]
                num, den = np.sum(w * Xf[:, j] * r), np.sum(w * Xf[:, j] ** 2)
                beta[j] = (num / den if j == 0
                           else np.sign(num) * max(abs(num) - lam, 0.0) / den)
    return beta


unpenalized = map_estimate("ridge", 0.0)                # the maximum likelihood fit
lams = np.array([1.0, 4.0, 16.0, 64.0, 256.0])
ridge_path = np.array([map_estimate("ridge", l) for l in lams])
lasso_path = np.array([map_estimate("lasso", l) for l in lams])
print(f"unpenalized slopes: {np.round(unpenalized[1:], 3)}")
print("normal prior with sd tau = 1/sqrt(lam)  ->  ridge;  Laplace prior with b = 1/lam  ->  lasso")
for k, l in enumerate(lams):
    print(f"  lam {l:6.1f}  ridge {np.round(ridge_path[k, 1:], 3)}   "
          f"lasso {np.round(lasso_path[k, 1:], 3)}")
```

The Bayesian reading does not stop at the mode. Park and Casella (2008) noted that the
Laplace prior is a scale mixture of normals —
\( \gamma_j\mid\tau_j^2\sim\Normal(0,\tau_j^2) \) with \( \tau_j^2 \) exponential —
which turns the lasso into a hierarchical model with normal full conditionals, and so into
a Gibbs sampler whose output is the whole posterior rather than a sparse point. Priors
designed to be *more* selective make the mixing distribution heavy-tailed at both ends: a
spike-and-slab prior mixes a point mass at zero with a diffuse component, and the
horseshoe of Carvalho, Polson and Scott (2010) takes \( \tau_j \) half-Cauchy, so small
coefficients are shrunk almost to zero and large ones almost not at all. The other
penalties of [Chapter 30](../ch30-regularization-boosting/index.html) carry over by the
same substitution, penalized least squares becoming penalized log-likelihood: the elastic
net of @def-reg-elastic-net, the group lasso of @def-reg-group-lasso for a factor's dummy
variables, and the componentwise boosting of @def-reg-l2boost with the negative
log-likelihood as loss and one Newton step per round.

## A prior that rescues an estimate

::: {#exm-prc-separated-bayes}
[The separated subgroup, three priors]

@exm-bin-separated fitted the intended vote on party identification among the
\( 13 \) respondents whose schooling stopped at or before the eighth
grade, and found complete separation: no maximum likelihood estimate exists, and software
reports whatever coefficient its iteration limit produced. With only two parameters the
posterior can be computed by quadrature, under any prior at all — over the whole plane,
by the substitution \( \beta=s\tan u \), since under separation a Cauchy posterior has a
polynomial tail and any fixed box would discard a large part of it.

Under the flat prior the posterior is not a distribution at all, and the arithmetic shows
it: integrating the likelihood over a square box centred at the origin gives
\( 8.7 \) for half-width ten and \( 778.6 \) for
half-width eighty, each doubling multiplying the mass by at least
\( 4.13 \). That is the area of the separating cone growing, as in the
proof of @prp-prc-propriety(b), and no amount of sampling will reveal it — a chain started
near the data's own scale looks well behaved while wandering up a ridge of infinite mass.

An independent \( \Normal(0,2.5^2) \) prior on the slope and
\( \Normal(0,10^2) \) on the intercept gives a proper posterior with mean
\( 3.077 \), median \( 2.931 \) and a
central \( 95 \) per cent credible interval of
\( (1.046,\ 5.922) \). Independent Cauchy priors of the same scales give
mean \( 6.091 \), median \( 4.001 \) and
interval \( (1.127,\ 22.995) \). Both exclude zero decisively. Both
lower endpoints agree to within a tenth, and the upper endpoints differ by a factor of
\( 3.9 \): the data bound the effect from below, and everything above
that is prior. Panel (b) of
[Figure 39.5.1](05-priors-as-regularizers.html#fig-prc-priors) shows the two posteriors.

Set these beside the frequentist answers of
[Section 35.3](../ch35-binary-responses/03-separation.html). The ordinary likelihood ratio
gives the one-sided interval \( (0.902,\infty) \): the most honest
statement available without a prior, and the least useful. Firth's penalty, the mode under
the Jeffreys prior (@prp-bin-firth), gives \( 1.362 \) with the
penalized likelihood ratio interval \( (0.364,\,4.045) \) — smaller
than either posterior median, and a mode rather than a mean. The Jeffreys penalty falls
off exponentially along the separating ray (@prp-bin-firth) while a Cauchy of scale
\( 2.5 \) falls off only polynomially, so it bites hardest on exactly the direction the
data do not constrain. There is no way to choose among these numbers from the data.
Thirteen observations that separate perfectly contain a bound and nothing more, and every
method reports that bound plus whatever it assumed. The discipline is to say which
assumption produced which number, and to notice that the flat prior is not the
assumption-free option but the one whose assumption is inconsistent.
:::

::: {when-format="html"}
![**Figure 39.5.1.** (a) The prior induced on a probability by a
\( \Normal(0,\sigma^2) \) prior on the logit; the middle curve is the boundary case
\( \sigma^2=2 \) of @prp-prc-priors(c). (b) Posterior densities of the
party-identification slope in the separated subgroup of @exm-bin-separated under a normal
and a Cauchy prior, with Firth's estimate
marked.](prior_separation.svg){#fig-prc-priors width=100%}
:::

::: {when-format="pdf"}
![(a) The prior induced on a probability by a
\( \Normal(0,\sigma^2) \) prior on the logit; the middle curve is the boundary case
\( \sigma^2=2 \) of @prp-prc-priors(c). (b) Posterior densities of the
party-identification slope in the separated subgroup of @exm-bin-separated under a normal
and a Cauchy prior, with Firth's estimate
marked.](prior_separation.pdf){width=100%}
:::

```{.python .run #cell-priors-separated}
sub = anes[anes["educ"] == 1]                          # schooling stopped at or before grade eight
Xs = np.column_stack([np.ones(len(sub)), sub["PID"]])
ys = sub["vote"].to_numpy(float)


def loglik(b0, b1):
    """Logistic log-likelihood at one intercept and a whole grid of slopes."""
    eta = b0 + np.outer(b1, Xs[:, 1])
    return np.sum(ys * eta - np.logaddexp(0.0, eta), axis=1)


def log_prior(b0, b1, kind, scale=2.5):
    if kind == "flat":
        return 0.0 * b0
    if kind == "normal":
        return -(b0 / 10.0) ** 2 / 2 - (b1 / scale) ** 2 / 2
    return -np.log1p((b0 / 10.0) ** 2) - np.log1p((b1 / scale) ** 2)   # Cauchy


def posterior_of_slope(kind, k=1200, s0=30.0, s1=6.0):
    """Marginal posterior of the slope, by quadrature over the whole plane."""
    t, wt = np.polynomial.legendre.leggauss(k)
    u, wu = t * np.pi / 2, wt * np.pi / 2
    b0, b1 = s0 * np.tan(u), s1 * np.tan(u)
    j0, j1 = s0 / np.cos(u) ** 2, s1 / np.cos(u) ** 2
    logq = np.array([loglik(a, b1) for a in b0]) + log_prior(b0[:, None], b1[None, :], kind)
    q = np.exp(logq - logq.max())
    marg = ((j0 * wu)[:, None] * q).sum(axis=0)        # integrate out the intercept
    dens = marg / np.sum(marg * j1 * wu)               # density with respect to b1
    return u, b1, dens, j1, wu


def summarize(u, b1, dens, j1, wu, probs=(0.025, 0.5, 0.975)):
    """Posterior mean and quantiles of the slope from the quadrature grid."""
    f = dens * j1                                      # the density in the u coordinate
    cdf = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(u))])
    return np.sum(b1 * f * wu), np.interp(probs, cdf / cdf[-1], b1)


post_n, post_c = posterior_of_slope("normal"), posterior_of_slope("cauchy")
mean_n, (lo_n, med_n, hi_n) = summarize(*post_n)
mean_c, (lo_c, med_c, hi_c) = summarize(*post_c)
print(f"posterior mean of the slope: normal {mean_n:.3f}, Cauchy {mean_c:.3f}")
```

::: {.idea}
A penalty is a prior and a prior a penalty, but the two readings answer different
questions. The penalty reading asks what estimate to report and tunes \( \lambda \) by
cross-validation; the prior reading asks what is believed before the data, and its scale
must be defensible on the scale of the parameter — for a binary response, as a statement
about probabilities.
:::

## Exercises

### A. Check your understanding

::: {#exr-prc-tau-lambda}
[A1]

A colleague reports a ridge-penalized logistic fit with \( \lambda=100 \) on standardized
covariates. What normal prior does that correspond to, and what odds ratio for a one
standard deviation change does the prior regard as, say, two standard deviations from
zero?
:::

::: {.solution}
\( \lambda=1/\tau^2 \) gives \( \tau=0.1 \). Two prior standard deviations is
\( \gamma_j=0.2 \), an odds ratio of \( e^{0.2}=1.22 \). The prior therefore says that a
one standard deviation change in any covariate almost certainly multiplies the odds by
less than about \( 1.22 \): a strong belief, not a weak one, which may be appropriate in
a high-dimensional screening problem and absurd in a four-covariate model.
:::

::: {#exr-prc-sigma-two}
[A2]

Explain in words why the induced prior of @prp-prc-priors(c) becomes U-shaped for large
\( \sigma \), without using the calculus of the proof.
:::

### B. Practice

::: {#exr-prc-cauchy-mixture}
[B1]

Show that a Cauchy prior with scale \( s \) is a scale mixture of normals:
\( \gamma\mid\lambda\sim\Normal(0,\lambda) \) with \( \lambda \) having an inverse gamma
distribution with shape \( 1/2 \) and rate \( s^2/2 \). Explain how this turns a Cauchy
prior into a Gibbs sampler with normal full conditionals in the probit model.
:::

::: {.solution}
If \( \lambda\sim\text{IG}(1/2,s^2/2) \) then the marginal density of \( \gamma \) is
\[
\begin{aligned}
&\int_0^\infty\frac{1}{\sqrt{2\pi\lambda}}e^{-\gamma^2/(2\lambda)}
\frac{(s^2/2)^{1/2}}{\Gamma(1/2)}\lambda^{-3/2}e^{-s^2/(2\lambda)}\,d\lambda\\
&\qquad\propto\int_0^\infty\lambda^{-2}e^{-(\gamma^2+s^2)/(2\lambda)}d\lambda
\propto\frac{1}{\gamma^2+s^2},
\end{aligned}
\]
the Cauchy density of scale \( s \), the last step by the substitution
\( v=1/\lambda \). In the Albert–Chib scheme of @prp-prc-albert-chib, conditioning on the
\( \lambda_j \) as well as on \( \bz \) makes the prior on \( \bgamma \) normal with
covariance \( \diag(\lambda_j) \), so the full conditional of \( \bbeta \) is normal,
and that of each \( \lambda_j \) is inverse gamma with shape \( 1 \) and rate
\( (\gamma_j^2+s^2)/2 \): every step is a standard draw.
:::

::: {#exr-prc-prior-probability}
[B2]

For \( \sigma=10 \), compute the prior probability that
\( \pi=e^{\eta}/(1+e^{\eta})<0.01 \) when \( \eta\sim\Normal(0,\sigma^2) \), and compare
with \( \sigma=1 \). Relate the answers to @prp-prc-priors(c).
:::

::: {#exr-prc-firth-jeffreys}
[B3]

@prp-bin-firth adds \( \tfrac12\log\det(\X\T\W\X) \) to the log-likelihood. Show that
this is the log of the Jeffreys prior density for the logistic model, so Firth's estimate
is a posterior mode. Why does that mode exist even under separation, when the maximum
likelihood estimate does not?
:::

::: {.solution}
The Jeffreys prior is \( \propto\det\{\boldsymbol{\mathcal I}(\bbeta)\}^{1/2} \) for
the expected information \( \boldsymbol{\mathcal I} \), which for the logistic model is
\( \X\T\W\X \) with
\( \W=\diag\{\pi_i(1-\pi_i)\} \) at \( \bbeta \) (@prp-bin-irls); its logarithm is
\( \tfrac12\log\det(\X\T\W\X) \) up to a constant, so adding it to \( \ell \) gives
Firth's penalized log-likelihood and its maximizer is the posterior mode under that prior.
It exists under separation for the reason @prp-bin-firth proves: along a separating ray
\( t\bb \) the weights decay like \( e^{-t\delta} \), so the log prior falls linearly in
\( t \) while the log-likelihood, already at most zero, has nothing left to gain. The
likelihood alone is flat along that ray; the prior supplies the curvature.
:::

### C. Going deeper

::: {#exr-prc-bimodal-probit}
[C1]

Redo @prp-prc-priors(c) for the probit link: with \( \eta\sim\Normal(0,\sigma^2) \) and
\( \pi=\Phi(\eta) \), show that the induced density of \( \pi \) is unimodal when
\( \sigma^2\le1 \) and bimodal when \( \sigma^2>1 \), and say why the threshold differs
from the logit's.
:::
