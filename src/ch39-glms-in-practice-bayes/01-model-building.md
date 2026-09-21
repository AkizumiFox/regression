# Building a model

A generalized linear model has three parts that can be chosen separately: a response
family, which fixes the variance function \( V(\mu) \); a link \( g \); and a linear
predictor \( \eta=\x\T\bbeta \).
[Chapter 34](../ch34-exponential-families-glm/index.html) fits whichever combination is
handed to it; this section is about how the combination is chosen, and how honest one can
be afterwards.

## The family follows the variance

Two things are known about a response before any model is fitted: what values it can take,
and how its spread changes with its level. The first rules out most families at once — a
count is not negative, a proportion is bounded, a duration is positive. The second
separates those that remain, because within the exponential dispersion family the variance
function determines the family (@exr-glm-variance-determines). For the families
of @def-glm-variance the variance is a power of the mean,
\[
\Var(Y_i)=\phi\,\mu_i^{\zeta}/w_i ,
\]{#eq-prc-power}

with \( \zeta=0 \) for the normal, \( 1 \) for the Poisson, \( 2 \) for the gamma and
\( 3 \) for the inverse Gaussian; in between are the Tweedie families
of @exr-glm-tweedie. So "which family?" is largely the question "which \( \zeta \)?",
and with replicates it can be answered before any model is fitted. A second choice arises
at once: the scale on which the systematic part is additive. Transforming the response and
transforming the mean are different operations (@exr-glm-link-vs-transform), and models on
different response scales are not comparable by their likelihoods unless the change of
variable is accounted for. The proposition collects the three facts that make these
choices more than a matter of taste.

::: {#prp-prc-model-building}
[Three facts about building a generalized linear model]

::: {.enumerate options="label=(\alph*)"}
1. *(Reading the variance function off replicates.)* Suppose the observations fall into
   \( G\ge2 \) groups, group \( g \) holding \( m_g\ge2 \) independent replicates with
   mean \( \mu_g>0 \), variance \( \phi\mu_g^{\zeta} \) and finite fourth moment. Let
   \( \bar y_g \) be the mean of the replicates in group \( g \), let \( s_g^2 \) be
   their variance with divisor \( m_g-1 \), and let
   \( \hat\zeta \) be the slope of the weighted least squares regression of
   \( \log s_g^2 \) on \( \log\bar y_g \) with any fixed positive weights \( c_g \). If
   the \( \log\mu_g \) are not all equal, then \( \hat\zeta\to\zeta \) in probability as
   \( \min_gm_g\to\infty \) with \( G \), the \( \mu_g \) and \( \phi \) fixed.

2. *(Additivity is a property of the link, not of the data.)* Let
   \( g(\mu(x_1,x_2))=a(x_1)+b(x_2) \) on \( I_1\times I_2 \), with \( a \) and \( b \)
   continuous and nonconstant, and let \( J \) be the interval of values taken by
   \( a+b \). Let \( \tilde g \) be a second link and \( s=\tilde g\circ g^{-1} \), twice
   differentiable on \( J \). Then \( \tilde g(\mu(x_1,x_2)) \) is again a sum of a
   function of \( x_1 \) and a function of \( x_2 \) if and only if \( s \) is affine on
   \( J \).

3. *(Comparing models across response scales.)* Let \( t \) be strictly increasing and
   continuously differentiable with \( t'>0 \), and let a model for the transformed
   responses \( \tilde y_i=t(y_i) \) have maximized log-likelihood \( \tilde\ell \) and
   \( d \) estimated parameters. As a model for \( \y \) it has maximized log-likelihood
   \( \tilde\ell+\sum_i\log t'(y_i) \), so its AIC on the scale of \( \y \) is
   \[
\text{AIC}=-2\tilde\ell+2d-2\sum_{i=1}^n\log t'(y_i).
\]{#eq-prc-jacobian}
:::
:::

::: {.proof}
(a) With \( m_g \) replicates of finite fourth moment, \( \bar y_g\to\mu_g \) and
\( s_g^2\to\phi\mu_g^{\zeta} \) in probability as \( m_g\to\infty \), the second because
\( \Var(s_g^2)\to0 \) when the fourth moment is finite. Both limits are positive, and
\( \log \) is continuous there, so \( L_g=\log\bar y_g\to\log\mu_g \) and
\( V_g=\log s_g^2\to\log\phi+\zeta\log\mu_g \) in probability. The weighted slope
\[
\hat\zeta=\frac{\sum_gc_g(L_g-\bar L_c)(V_g-\bar V_c)}{\sum_gc_g(L_g-\bar L_c)^2},
\qquad
\bar L_c=\frac{\sum_gc_gL_g}{\sum_gc_g},
\]
is a continuous function of \( (L_1,V_1,\dots,L_G,V_G) \) wherever the denominator is
nonzero. At the limit point the denominator is
\( \sum_gc_g(\log\mu_g-\overline{\log\mu}_c)^2>0 \), because the \( \log\mu_g \) are not
all equal, and the numerator is \( \zeta \) times the same quantity. The continuous
mapping theorem gives \( \hat\zeta\to\zeta \).

(b) If \( s(u)=\alpha+\beta u \) then \( \tilde g(\mu)=\alpha+\beta a(x_1)+\beta b(x_2) \),
a sum. Conversely suppose \( \tilde g(\mu(x_1,x_2))=A(x_1)+B(x_2) \). If
\( a(x_1)=a(x_1') \) then \( \mu(x_1,x_2)=\mu(x_1',x_2) \) for every \( x_2 \), so
\( A(x_1)=A(x_1') \); hence \( A^{*}(u)=A(x_1) \) is well defined on the range \( U \) of
\( a \), and likewise \( B^{*} \) on the range \( V \) of \( b \). Both are
nondegenerate intervals, being continuous images of intervals under nonconstant maps, and
\[
s(u+v)=A^{*}(u)+B^{*}(v),\qquad u\in U,\ v\in V .
\]
Fix \( u_0,v_0 \) interior to \( U,V \). For any other such \( u,v \) the four terms
cancel:
\[
s(u+v)-s(u_0+v)-s(u+v_0)+s(u_0+v_0)=0 .
\]
Differentiating in \( u \) and then in \( v \) gives \( s''(u+v)=0 \). As \( u+v \)
ranges over an open subinterval of \( J \) whose closure is \( J \), and \( s'' \) is
continuous, \( s''\equiv0 \): \( s \) is affine.

(c) If \( \tilde Y=t(Y) \) has density \( f_{\tilde Y} \), then \( Y \) has density
\( f_Y(y)=f_{\tilde Y}(t(y))\,t'(y) \), because \( t \) is a strictly increasing
differentiable bijection onto its image. Summing logs over independent observations,
\( \ell(\cdot;\y)=\tilde\ell(\cdot;\tilde{\y})+\sum_i\log t'(y_i) \) for every parameter
value, so the maximizers coincide and the maximized values differ by the stated constant.
Substituting into \( \text{AIC}=-2\ell+2d \) gives @eq-prc-jacobian.
:::

Part (c) is worth stating because the omitted term is enormous. The AIC that software
prints for a normal linear model fitted to \( \log\y \) compares it with other models
*for \( \log\y \)*; against a gamma model for \( \y \) the comparison is meaningless
until \( 2\sum_i\log y_i \), a sum of order \( n \), is added back.

::: {.warning}
The same applies to cross-validation, to the deviance and to any other criterion built on a
likelihood. A predictive score is a score *for a particular quantity on a particular
scale*; predicting \( \log y \) well and predicting \( y \) well are different
achievements, and the second is usually what was wanted.
:::

## A worked choice

::: {#exm-prc-strikes}
[How long does a strike last]

The durations of \( 62 \) contract strikes in United States manufacturing, from one day
to \( 216 \), are recorded with a measure of unanticipated industrial production at the
time. @exm-ql-strikes-fit read the same durations as counts and asked only for a variance
function; here they are read as what they are, positive and continuous, and the question is
which family to use. The response is strongly right-skewed, so a normal linear model is out;
which of the positive families fits the mean–variance relationship?

The covariate takes only \( 9 \) distinct values, so the data arrive already
in replicate groups, of sizes \( 2 \) to \( 18 \). The weighted
regression of @prp-prc-model-building(a) gives
\( \hat\zeta=1.606 \) with a standard error of
\( 0.434 \), computed from the exact normal-theory variance
\( \Var(\log s_g^2)=\psi'\{(m_g-1)/2\} \) — the familiar \( 2/(m_g-1) \) is only its
large-group limit, and understates it badly for the two groups of size two. Panel (a) of
[Figure 39.1.1](01-model-building.html#fig-prc-variance-power) shows the nine points with
the fitted line and the Poisson and gamma references: \( \zeta=2 \) is
\( 0.9 \) standard errors away and \( \zeta=1 \) is
\( 1.4 \), which is weak evidence for the gamma and against anything
Poisson-like.

Fit three models with a log link and the same linear predictor, maximizing each over its
dispersion parameter as well as over \( \bbeta \), and count \( d=3 \) for all three. AIC
on the scale of the durations is \( 585.4 \) for the gamma,
\( 590.9 \) for the log-normal and
\( 603.3 \) for the inverse Gaussian. The log-normal figure is the AIC
of the least squares fit to \( \log\y \), \( 206.8 \), plus the
Jacobian term \( 2\sum_i\log y_i=384.1 \) of @eq-prc-jacobian. Without
that term the log-normal appears to beat the gamma by
\( 378.7 \): not a close call, and in the wrong direction, because
the uncorrected number answers a different question. Leave-one-out cross-validation with
the log score puts the three in the same order, at \( -4.7192 \),
\( -4.7713 \) and \( -4.8664 \) per observation.

The gamma slope is to the last decimal the estimate that
@exm-ql-strikes-fit obtained from \( V(\mu)=\mu^{2} \) alone: the variance function
fixes the estimating equations, and choosing a family on top of it buys a likelihood, not a
different fit. The maximum likelihood shape of the gamma fit is \( \hat\nu=1.003 \),
so the durations are close to exponential given the covariate, which is what \( \zeta=2 \)
with a coefficient of variation near one means. (Software reports instead the moment
estimate from the Pearson dispersion \( \hat\phi=0.941 \), namely
\( 1.063 \); the two agree here, and for the inverse Gaussian they do
not.) The fits also disagree about what they estimate: at the mean covariate value the
gamma model puts the *mean* duration at \( 39.4 \) days, while the
log-normal model puts the *median* at \( 22.2 \) and, after the
retransformation correction of @prp-tr-retransformation, the mean at
\( 46.7 \).
:::

::: {when-format="html"}
![**Figure 39.1.1.** Choosing a family for the strike durations. (a) Log group variance
against log group mean, point size proportional to group size; the fitted slope estimates
\( \zeta \) in @eq-prc-power, and the reference lines have the Poisson and gamma slopes.
(b) The three fitted densities at the mean covariate value, over a histogram of the
durations.](variance_power.svg){#fig-prc-variance-power width=100%}
:::

::: {when-format="pdf"}
![Choosing a family for the strike durations. (a) Log group variance
against log group mean, point size proportional to group size; the fitted slope estimates
\( \zeta \) in @eq-prc-power, and the reference lines have the Poisson and gamma slopes.
(b) The three fitted densities at the mean covariate value, over a histogram of the
durations.](variance_power.pdf){width=100%}
:::

```{.python .run #cell-family-choice-groups}
import numpy as np
import statsmodels.api as sm
strikes = sm.datasets.strikes.load_pandas().data       # 62 contract strikes, public domain
y = strikes["duration"].to_numpy(float)                # days
x = strikes["iprod"].to_numpy(float)                   # unanticipated industrial production

groups = strikes.groupby("iprod")["duration"].agg(["count", "mean", "var"])
groups = groups[groups["count"] >= 2]                  # nine groups, sizes 2 to 18
w = groups["count"].to_numpy() - 1.0                   # weights: m_g - 1 degrees of freedom
L = np.column_stack([np.ones(len(groups)), np.log(groups["mean"])])
zeta_fit = np.linalg.solve(L.T @ (w[:, None] * L), L.T @ (w * np.log(groups["var"])))

print(groups.round(2))
print(f"variance-function power  zeta = {zeta_fit[1]:.3f}")
```

```{.python .run #cell-family-choice-families}
from scipy import optimize, special, stats
X = np.column_stack([np.ones(len(y)), x])
n, p = X.shape
d = p + 1                                              # slopes plus one dispersion parameter

gam = sm.GLM(y, X, family=sm.families.Gamma(sm.families.links.Log())).fit()
igauss = sm.GLM(y, X, family=sm.families.InverseGaussian(sm.families.links.Log())).fit()
lognorm = sm.OLS(np.log(y), X).fit()                   # a normal linear model for log y


def gamma_loglik(yy, mu, nu):
    return np.sum(stats.gamma.logpdf(yy, nu, scale=mu / nu))


def gamma_shape(yy, mu):
    """Maximum likelihood shape: the root of log nu - digamma(nu) + mean log(y/mu)."""
    c = np.mean(np.log(yy / mu))
    return optimize.brentq(lambda v: np.log(v) - special.digamma(v) + c, 1e-4, 1e4)


def igauss_lambda(yy, mu):
    """Maximum likelihood precision of the inverse Gaussian, in closed form."""
    return 1.0 / np.mean((yy - mu) ** 2 / (mu**2 * yy))


nu_hat = gamma_shape(y, gam.fittedvalues)
lam_hat = igauss_lambda(y, igauss.fittedvalues)
ll_gam = gamma_loglik(y, gam.fittedvalues, nu_hat)
ll_ig = np.sum(stats.invgauss.logpdf(y, igauss.fittedvalues / lam_hat, scale=lam_hat))
sigma_ml = np.sqrt(lognorm.ssr / n)                    # the maximizing sigma, divisor n
ll_ln = np.sum(stats.norm.logpdf(np.log(y), lognorm.fittedvalues, sigma_ml) - np.log(y))

for name, ll in [("gamma, log link", ll_gam), ("inverse Gaussian", ll_ig),
                 ("log-normal (on the y scale)", ll_ln)]:
    print(f"{name:28s} loglik {ll:9.2f}   AIC {-2 * ll + 2 * d:9.2f}")
print(f"gamma:  slope {gam.params[1]:8.3f},  maximum likelihood shape {nu_hat:.3f}")
```

## Links, scales and interactions

With the family fixed, the link decides what "additive" means. The canonical
link (@def-glm-link) keeps the log-likelihood concave (@thm-glm-concave) and makes observed
and expected information agree (@thm-glm-score(d)); neither is a reason to prefer it when
its range constrains \( \bbeta \). The gamma's canonical link is the reciprocal, which
forces \( \x\T\bbeta>0 \) and makes every coefficient a statement about \( 1/\mu \);
the log link has no constraint, keeps the concavity (@exr-glm-concavity-links), and makes
each coefficient a multiplicative effect on the mean. Almost all gamma regressions use the
log link, for interpretation rather than fit.

For the covariates, @thm-tr-box-tidwell transfers verbatim with the linear predictor in
place of the mean: to ask whether \( x_j \) should enter as \( x_j^{\lambda} \), add
\( x_j\log x_j \) and test its coefficient. Centring matters for a different reason.
Under a nonlinear link the intercept is the linear predictor at \( \x=\bzero \), a point
no observation resembles when the covariates are uncentred; centring each at a
representative value makes \( e^{\hat\beta_0} \) the fitted mean of a typical case.

Interactions are where the link stops being cosmetic. A two-factor Poisson model with a
log link and no interaction says the factors multiply; with an identity link it says they
add; by @prp-prc-model-building(b) at most one of these is the additive model
of @def-tw-interaction. "No interaction" is therefore never a property of the response
alone, and removing an interaction by changing the link is a legitimate simplification,
provided the link is reported as part of the finding and was not chosen by looking at the
interaction test (@prp-sel-selection-bias). On the log scale an interaction coefficient is
a ratio of ratios: with binary \( x_1,x_2 \) and
\( \log\mu=\beta_0+\beta_1x_1+\beta_2x_2+\beta_{12}x_1x_2 \), \( e^{\beta_{12}} \)
is the factor by which the effect of \( x_1 \) changes when \( x_2 \) switches on; on
the logit scale it is a ratio of odds ratios, harder to think about because odds ratios do
not collapse over omitted covariates (@exm-bin-noncollapsible). Offsets need no new
theory (@prp-cnt-offsets): a known exposure enters with coefficient fixed at one, and the
model describes a rate.

## Comparing and choosing

Because a generalized linear model is fitted by maximum likelihood, everything in
[Chapter 29](../ch29-model-selection/index.html) applies: AIC is
\( -2\ell(\hbeta)+2d \) with \( d=p \) when \( \phi \) is known and \( p+1 \) when it
is estimated, BIC replaces \( 2 \) by \( \log n \) (@def-sel-aic-bic), and their
different targets (@prp-sel-aic-bic) are unchanged. Two cautions are specific here. First,
by @prp-prc-model-building(c) the likelihoods compared must be for the same random
variable, and each must be maximized over the dispersion as well as over \( \bbeta \).
Second, a quasi-likelihood model has no likelihood, so no AIC; the substitute is the
quasi-likelihood information criterion built from the independence-model quasi-likelihood
and a sandwich penalty, which [Chapter 40](../ch40-glmm-gee/index.html) defines. The
dispersion estimate is not a likelihood quantity either (@prp-ql-dispersion), which is why
an overdispersed fit is compared by quasi-\( F \) statistics.

Cross-validation needs a loss, and the natural one is the log score: hold out an
observation, fit to the rest, record the predictive log density of what was held out. That
is @def-sel-cv with the deviance as loss, so @prp-sel-cv-target applies. For a
\( \phi=1 \) family the log score and the deviance rank models identically; with a
dispersion parameter they need not, because the deviance ignores the term in \( \phi \).

::: {.idea}
Three questions, three tools. *Which family?* The support and the mean–variance
relationship. *Which link?* What should be additive, checked
afterwards ([Section 39.2](02-diagnostics.html)). *Which terms?* A criterion, plus
honesty about what the search cost.
:::

## What the search costs

Every warning of [Section 29.5](../ch29-model-selection/05-selection-bias.html) carries
over unchanged, because none used normality: after a search the standard errors are too
small, the \( p \)-values too small and the intervals too
short (@prp-sel-selection-bias). Two remarks are specific here. A drop in deviance is
approximately \( \chi^2(q) \) under the smaller model (@thm-glm-deviance), so a maximum
over many candidate terms is not, and the largest always looks significant; and with small
\( n \) or sparse cells that approximation is itself poor (@exm-glm-hauck-donner,
@prp-bin-fit), so the nominal levels being abused were never right to begin with. Sample
splitting (@prp-sel-splitting) and simultaneous inference over the whole
search (@prp-sel-simultaneous) are the honest ways out.

## Exercises

### A. Check your understanding

::: {#exr-prc-scale-aic}
[A1]

A colleague fits a normal linear model to \( \sqrt{y} \) and a Poisson log-linear model to
\( y \), and compares the two AIC values printed by the software. What is the correction
term of @prp-prc-model-building(c) here, and in which direction does omitting it push the
comparison?
:::

::: {.solution}
With \( t(y)=\sqrt y \), \( t'(y)=1/(2\sqrt y) \), so the correction term is
\[
-2\sum_i\log t'(y_i)=2\sum_i\log\bigl(2\sqrt{y_i}\bigr)=2n\log2+\sum_i\log y_i ,
\]
and the AIC of the square-root model on the scale of \( \y \) is its printed value plus
that amount. It is positive whenever \( \prod_i4y_i>1 \), as for counts of any
reasonable size it is, so the printed comparison flatters the square-root model. The
Poisson model needs no correction, being already a model for \( \y \) — though it is a
mass function rather than a density, which makes the comparison dubious anyway.
:::

::: {#exr-prc-interaction-scale}
[A2]

A two-by-two table of counts has fitted means \( \mu_{11}=10 \), \( \mu_{12}=20 \),
\( \mu_{21}=15 \) and \( \mu_{22}=30 \). Is there an interaction on the log scale? On the
identity scale? Reconcile the two answers with @prp-prc-model-building(b).
:::

::: {.solution}
On the log scale, \( \log\mu_{11}+\log\mu_{22}-\log\mu_{12}-\log\mu_{21}
=\log(10\cdot30/(20\cdot15))=\log1=0 \): no interaction. On the identity scale,
\( 10+30-20-15=5\ne0 \): there is one. The two links differ by the map \( s=\exp \), which
is not affine, so @prp-prc-model-building(b) says there is no reason for both to be
additive, and here only one is.
:::

### B. Practice

::: {#exr-prc-gamma-canonical}
[B1]

The canonical link of the gamma family is \( g(\mu)=-1/\mu \) (@exm-glm-gamma). Show that
with this link and a model matrix containing an intercept, the parameter space
\( \{\bbeta:\x_{(i)}\T\bbeta<0 \text{ for all } i\} \) is an open convex cone, and that the
maximum likelihood estimate may fail to exist by running to its boundary. Contrast with the
log link.
:::

::: {.solution}
The set is an intersection of \( n \) open half-spaces, hence open and convex, and a cone
because each constraint is homogeneous. The log-likelihood tends to \( -\infty \) as any
\( \x_{(i)}\T\bbeta\to0^- \), since \( \mu_i\to\infty \); but the supremum may be
approached only as \( \norm{\bbeta}\to\infty \) inside the cone, in which case no
maximizer exists, exactly as in @thm-bin-separation. With the log link the parameter space
is all of \( \Real^p \) and the log-likelihood is
concave (@exr-glm-concavity-links), so the only obstacle is the analogue of separation.
:::

::: {#exr-prc-logscore-deviance}
[B2]

Show that for a family with \( \phi=1 \) the leave-one-out log score and the leave-one-out
deviance differ by a constant that does not depend on the model, so the two rank models
identically. Show that this fails when \( \phi \) is estimated.
:::

### C. Going deeper

::: {#exr-prc-additivity-three}
[C1]

Extend @prp-prc-model-building(b) to three covariates: if
\( g(\mu)=a_1(x_1)+a_2(x_2)+a_3(x_3) \) with all three nonconstant and continuous, and
\( \tilde g(\mu) \) is also a sum of three functions of the separate covariates, must
\( s=\tilde g\circ g^{-1} \) be affine? Does the conclusion change if only *two* of the
three are required to separate?
:::

::: {#exr-prc-power-link}
[C2]

Let the true link belong to the power family
\( g_\lambda(\mu)=(\mu^{\lambda}-1)/\lambda \), which has \( g_0=\log \), and let the
model actually fitted use the log link. Following the constructed-variable argument of
@prp-tr-constructed, expand \( g_\lambda \) about \( \lambda=0 \) to show that a one
degree of freedom test of \( \lambda=0 \) is obtained by adding a multiple of
\( \hat\eta^2 \) to the linear predictor, and identify the constructed variable and the
multiple exactly.
:::
