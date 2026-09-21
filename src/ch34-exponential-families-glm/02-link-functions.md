# Link functions

The linear model says \( \E(Y_i)=\x_{(i)}\T\bbeta \). For a count that must be nonnegative, or
a probability that must lie in \( (0,1) \), that equation is a promise the right-hand side
cannot keep: a linear function of unbounded regressors takes every real value. The
resolution is to put a monotone function between the mean and the linear predictor.

## The model

::: {#def-glm-model}
[Generalized linear model]

Let \( \X \) be an \( n\times p \) model matrix with rows \( \x_{(i)}\T \), and let
\( w_1,\dots,w_n>0 \) be known prior weights. A **generalized linear model** consists of

::: {.enumerate options="label=(\alph*)"}
1. a **random component**: \( Y_1,\dots,Y_n \) are independent, and \( Y_i \) has the density
   @eq-glm-edf with natural parameter \( \theta_i \), weight \( w_i \) and a dispersion \( \phi \)
   that is the same for all \( i \);

2. a **linear predictor** \( \eta_i=\x_{(i)}\T\bbeta \), with
   \( \boldsymbol{\upeta}=\X\bbeta\in\Real^n \);

3. a **link function** \( g:\mathcal M\to\Real \), strictly monotone and twice continuously
   differentiable, with
   \[
   g(\mu_i)=\eta_i,\qquad \mu_i=\E(Y_i)=b'(\theta_i).
   \]
:::

The inverse \( h=g^{-1} \) is the **response function**, so \( \mu_i=h(\eta_i) \). The
**admissible set** of the model is
\( \mathcal B=\{\bbeta\in\Real^p:\ \x_{(i)}\T\bbeta\in g(\mathcal M)\ \text{for all }i\} \).
:::

The linear model of [Chapter 5](../ch05-model-and-least-squares/index.html) is the case:
normal random component, identity link, \( w_i=1 \), \( \phi=\sigma^2 \), \( \mathcal B=\Real^p \).
Every classical result of Parts II and III is about that one cell of a large table.

Three remarks. The dispersion is common to all observations: the regressors move the mean,
the variance follows through \( V \), and nothing else may vary (modelling \( \phi \) too is
the subject of [Chapter 45](../ch45-quantile-gamlss/index.html), @def-qnt-gamlss). The weights \( w_i \) are known, not estimated, and usually record
group sizes. And \( g \) transforms the mean, not the data:
[Chapter 22](../ch22-transformations/index.html) transformed \( Y \) itself, whereas here the
response is left alone and the link acts inside the expectation.

::: {#def-glm-link}
[Canonical link]

The **canonical link** of a family is \( g=(b')^{-1} \), so that \( \theta_i=\eta_i \): the
natural parameter *is* the linear predictor. Equivalently, \( g'(\mu)=1/V(\mu) \).
:::

| Family | canonical link | \( g(\mu) \) | range of \( \eta \) | other links in use |
|---|---|---|---|---|
| normal | identity | \( \mu \) | \( \Real \) | log, inverse |
| binomial | logit | \( \log\{\mu/(1-\mu)\} \) | \( \Real \) | probit, complementary log–log, identity |
| Poisson | log | \( \log\mu \) | \( \Real \) | identity, square root |
| gamma | inverse | \( -1/\mu \) | \( (-\infty,0) \) | log, identity |
| inverse Gaussian | inverse square | \( -1/(2\mu^2) \) | \( (-\infty,0) \) | log |
| negative binomial | \( \log\{\mu/(\mu+k)\} \) | as shown | \( (-\infty,0) \) | log |

The sign conventions for the gamma and inverse Gaussian are cosmetic: \( 1/\mu \) and
\( -1/\mu \) give the same fitted values with \( \bbeta \) negated. The column headed "range of
\( \eta \)" is not. For the normal, binomial and Poisson the canonical link maps
\( \mathcal M \) onto all of \( \Real \), so \( \mathcal B=\Real^p \) and no constraint is ever
active. For the gamma and inverse Gaussian it maps onto a half-line, and \( \mathcal B \) is
the polyhedron \( \{\bbeta:\x_{(i)}\T\bbeta<0\ \text{for all }i\} \) — an open convex set
with a boundary a fitting algorithm can walk into.

## What a link choice costs and buys

::: {.idea}
The canonical link is a mathematical convenience, not a scientific claim. It makes the
likelihood concave (@thm-glm-concave), makes observed and expected information
agree (@thm-glm-score), reduces the data to \( p \) sufficient statistics
([Section 34.1](01-exponential-dispersion-families.html)), and forces the fitted values to
reproduce the observed weighted column totals of \( \X \) exactly (@cor-glm-marginals). None
of that is a reason to believe the mean depends on the regressors in that particular way.
:::

Four considerations decide a link in practice.

**Range.** A link whose range is all of \( \Real \) removes constraints from the fitting
problem. The log link for a Poisson response guarantees \( \hat\mu_i>0 \) for every
\( \bbeta\in\Real^p \), and the logit link for a binomial response guarantees
\( \hat\mu_i\in(0,1) \). An identity link for either does not.

::: {#exm-glm-identity-poisson}
[An identity link close to its boundary]

Eight counts were simulated from a Poisson model with mean \( 12-1.2x \) at
\( x=1,\dots,8 \). An identity-link Poisson fit recovers the declining mean, and every fitted
value is positive — but the fitted line reaches zero at \( x=9.89 \), barely
outside the range of the data. Predicting at \( x=10 \) would return a negative mean, and a
slightly different sample would push the boundary inside the observed range, where the
likelihood is \( -\infty \). The log link has no such boundary.
:::

**Interpretation.** A log link makes the regressors act multiplicatively,
\( \mu=e^{\beta_0}\prod_je^{\beta_jx_j} \), so \( \beta_j \) is a log rate ratio and
\( e^{\beta_j} \) the factor multiplying the mean per unit of \( x_j \). A logit link makes
\( \beta_j \) a log odds ratio, the subject
of [Chapter 35](../ch35-binary-responses/index.html) (@def-bin-logistic). An identity link
keeps \( \beta_j \) a difference of means, often what a policy question asks about. These are
different scientific statements, and the choice belongs on the scale where effects are
plausibly additive, not on convenience.

**Fit and stability.** Different links give genuinely different mean curves, and a badly
chosen one shows up in the residuals like any other misspecification. Non-canonical links
can also produce a log-likelihood with several stationary points, and links with a
restricted range need the fitting algorithm to be protected; [Section 34.4](04-irls.html)
returns to this.

::: {#exm-glm-engel}
[Food expenditure: three models, two decisions]

Engel's classic budget survey records annual income and food expenditure for
\( 235 \) working-class households. Food expenditure is positive, right-skewed, and
visibly more variable at high incomes. Write \( u=\log(\text{income}) \) and fit three models
with the same linear predictor \( \beta_0+\beta_1u \): a normal random component with the
identity link (ordinary least squares), a gamma random component with the identity link,
and a gamma random component with the log link. The first two differ only in the variance
function, the last two only in the link.

The identity-link fits give \( \hat\beta_1=571.02 \) (normal) and
\( 474.54 \) (gamma) francs per unit of \( \log \) income. The gamma/log fit
gives \( \hat\beta_1=0.8629 \), an elasticity: a one per cent rise in income is
associated with a \( 0.86 \) per cent rise in expected food expenditure. Its estimated
dispersion is \( \hat\phi=0.0177 \), so by @exm-glm-gamma the fitted
coefficient of variation is \( \sqrt{\hat\phi}=0.133 \), the same at every
income.

[Figure 34.2.1](02-link-functions.html#fig-glm-links) separates the two decisions. In panel
(a) the identity-link mean curve is concave in income and flattens, while the log-link curve
is a power law that keeps climbing; at the largest income in the sample the two fitted means
differ by a factor of \( 1.786 \). Panel (b) holds the link at the identity
and changes the variance function, so its bands come from the identity-link fits, the gamma
one using \( \hat\phi=0.0240 \) and a coefficient of variation
of \( 0.155 \). The normal band has constant width and at low
incomes dips below zero; the gamma band widens in proportion to the mean, as the data do.
:::

```{.python .run #cell-links-fits}
import numpy as np
import statsmodels.api as sm

engel = sm.datasets.engel.load_pandas().data
y = engel["foodexp"].to_numpy()                  # food expenditure
X = np.column_stack([np.ones(len(y)), np.log(engel["income"].to_numpy())])

gauss = sm.GLM(y, X, family=sm.families.Gaussian()).fit()            # identity link
gam_id = sm.GLM(y, X, family=sm.families.Gamma(sm.families.links.Identity())).fit()
gam_log = sm.GLM(y, X, family=sm.families.Gamma(sm.families.links.Log())).fit()

for name, fit in [("normal/identity", gauss), ("gamma/identity", gam_id), ("gamma/log", gam_log)]:
    print(f"{name:16s} coefficients {fit.params[0]:9.3f} {fit.params[1]:8.3f}"
          f"   dispersion {fit.scale:8.4f}")
```

::: {when-format="html"}
![**Figure 34.2.1.** Engel's household budget data, \( n=235 \). (a) Two links
with the same gamma random component and the same regressors: the identity link
(solid) and the log link (dashed). (b) Two variance functions with the same identity link:
the normal band \( \hat\mu\pm2\hat\sigma \) has constant width and dips below zero, while the
gamma band \( \hat\mu(1\pm2\sqrt{\hat\phi}) \), with \( \hat\phi \) from the identity-link gamma
fit, widens with the mean.](links_engel.svg){#fig-glm-links width=100%}
:::

::: {when-format="pdf"}
![Engel's household budget data, \( n=235 \). (a) Two links
with the same gamma random component and the same regressors: the identity link
(solid) and the log link (dashed). (b) Two variance functions with the same identity link:
the normal band \( \hat\mu\pm2\hat\sigma \) has constant width and dips below zero, while the
gamma band \( \hat\mu(1\pm2\sqrt{\hat\phi}) \), with \( \hat\phi \) from the identity-link gamma
fit, widens with the mean.](links_engel.pdf){width=100%}
:::

## Links and transformations are not the same thing

Transforming the response and choosing a link both replace \( \mu \) by \( g(\mu) \), but in
different places. Modelling \( \E\{g(Y)\}=\x\T\bbeta \) changes the quantity modelled: by
Jensen's inequality \( \E\{g(Y)\}\ne g\{\E(Y)\} \) unless \( g \) is affine, so the fitted curve
estimates the mean of a transformed variable, and returning to the mean of \( Y \) needs the
retransformation correction of @prp-tr-retransformation. Modelling
\( g\{\E(Y)\}=\x\T\bbeta \) leaves the response alone, and \( \hat\mu_i=h(\x_{(i)}\T\hbeta) \)
estimates the mean of \( Y \) directly.

The transformation also had a second job the link does not do: stabilizing the
variance (@thm-tr-variance-stabilizing). A generalized linear model does not need that,
because the variance function is already specified, and splitting one overloaded decision
into two independent ones is the framework's main practical advantage. The weighted
regression on empirical logits of @exr-tr-empirical-logit is the historical compromise: the
right link, but one weighted least squares step in place of the likelihood.

## Exercises

### A. Check your understanding

::: {#exr-glm-admissible}
[A1]

Describe the admissible set \( \mathcal B \) of @def-glm-model for (a) a Poisson model with
the identity link, (b) a binomial model with the identity link, (c) a gamma model with the
log link. Which of the three is all of \( \Real^p \)?
:::

### B. Practice

::: {#exr-glm-link-vs-transform}
[B1]

Let \( Y \) be positive with \( \E(Y)=\mu \) and \( \Var(Y)=\phi\mu^2 \), as in the gamma family.
Using a second-order Taylor expansion of \( \log Y \) about \( \mu \), show that
\( \E(\log Y)\approx\log\mu-\phi/2 \). Deduce that a least squares regression of \( \log y_i \) on
\( \X \) and a gamma/log generalized linear model estimate intercepts that differ by
approximately \( \phi/2 \), while their slopes agree to this order. For the fit
of @exm-glm-engel, how large is the discrepancy?
:::

::: {.solution}
With \( g(y)=\log y \), \( g(Y)\approx g(\mu)+(Y-\mu)/\mu-(Y-\mu)^2/(2\mu^2) \); taking
expectations gives \( \E(\log Y)\approx\log\mu-\Var(Y)/(2\mu^2)=\log\mu-\phi/2 \). Since the
correction is the same constant for every observation, it is absorbed by the intercept when
the model has one, and the slopes are unaffected to this order. For Engel's data
\( \hat\phi=0.0177 \), so the two intercepts differ by about \( \hat\phi/2=0.0177/2=0.009 \), less
than one per cent on the expenditure scale. (For the gamma, \( \E(\log Y) \) can be computed
exactly with the digamma function; the approximation is good when \( \phi \) is small.)
:::

::: {#exr-glm-sqrt-link}
[B2]

A Poisson model is sometimes fitted with the square-root link \( g(\mu)=\sqrt\mu \),
because @thm-tr-variance-stabilizing makes \( \sqrt{Y} \) approximately variance-stabilized. Find
\( \mathcal B \) for this link, find \( d\mu/d\eta \), and explain why the link being
variance-stabilizing for \( Y \) is irrelevant to whether it is a good link.
:::

::: {.solution}
\( g(\mathcal M)=(0,\infty) \), so \( \mathcal B=\{\bbeta:\x_{(i)}\T\bbeta>0\ \forall i\} \), and
\( \mu=\eta^2 \) gives \( d\mu/d\eta=2\eta=2\sqrt\mu \). Variance stabilization matters when the
model is fitted by *least squares on the transformed data*, because least squares wants
constant variance. A generalized linear model already knows \( V(\mu)=\mu \) and weights each
observation accordingly (@thm-glm-score), so the only question left about the link is
whether \( \sqrt{\E Y} \) is plausibly linear in the regressors.
:::

::: {#exr-glm-cloglog-range}
[B3]

The complementary log–log link for a binomial response is
\( g(\mu)=\log\{-\log(1-\mu)\} \). Show that it is a bijection from \( (0,1) \) onto \( \Real \),
compute \( d\mu/d\eta \), and show that it is *not* symmetric: \( g(\mu)\ne-g(1-\mu) \), unlike
the logit and probit links. Which tail does it approach more slowly?
:::

::: {.solution}
\( \mu\mapsto-\log(1-\mu) \) maps \( (0,1) \) onto \( (0,\infty) \) increasingly, and \( \log \) maps
\( (0,\infty) \) onto \( \Real \), so \( g \) is an increasing bijection with inverse
\( h(\eta)=1-\exp(-e^{\eta}) \) and \( d\mu/d\eta=e^{\eta}\exp(-e^{\eta})=(1-\mu)\{-\log(1-\mu)\} \).
Symmetry would require \( h(-\eta)=1-h(\eta) \), which fails, for example at \( \eta=0 \), where
\( h(0)=1-e^{-1}=0.632 \ne\tfrac12 \). As \( \eta\to-\infty \), \( \mu\approx e^{\eta}\to0 \)
exponentially, while as \( \eta\to+\infty \), \( 1-\mu=\exp(-e^{\eta})\to0 \) doubly
exponentially: the approach to \( 1 \) is far faster than the approach to \( 0 \).
[Chapter 35](../ch35-binary-responses/index.html) derives the link from a proportional-hazards argument (@prp-bin-links).
:::


### C. Going deeper

::: {#exr-glm-power-links}
[C1]

Every link in the table except the two for a binomial response is a power of \( \mu \) or a logarithm. That
is not an accident. For a positive response call a link \( g \) on \( (0,\infty) \)
**scale-equivariant** if for every \( c>0 \) there are constants \( a(c)\ne0 \) and \( b(c) \) with
\( g(c\mu)=a(c)g(\mu)+b(c) \) for all \( \mu>0 \).

::: {.enumerate options="label=(\alph*)"}
1. Check that \( g(\mu)=\mu^{\lambda} \) with \( \lambda\ne0 \) and \( g(\mu)=\log\mu \) are
   scale-equivariant, and identify \( a(c) \) and \( b(c) \).

2. Conversely, let \( g \) be strictly monotone and twice continuously differentiable on \( (0,\infty) \), as
   @def-glm-model requires, and scale-equivariant. Show that \( g(\mu)=A\mu^{\lambda}+B \) for some
   \( A\ne0 \), \( \lambda\ne0 \) and \( B \), or \( g(\mu)=A\log\mu+B \) with \( A\ne0 \). (Put
   \( G(t)=g(e^{t}) \), differentiate the defining relation in \( t \), and solve the resulting equation for
   \( G' \).)

3. Explain what this means for a model with an intercept: if the response is measured in other units,
   \( Y\mapsto cY \), the class of fitted mean curves \( x\mapsto h(\x\T\bbeta) \) is unchanged exactly when
   the link is scale-equivariant, the coefficients being relabelled as
   \( \bbeta\mapsto a(c)\bbeta+b(c)\mathbf{e}_1 \).

4. Show that \( g(\mu)=\mu-1/\mu \), a perfectly legal increasing bijection from \( (0,\infty) \) onto
   \( \Real \), is not scale-equivariant, so that the model it defines is a different model for expenditure
   in francs and for expenditure in centimes.
:::

:::

::: {.solution}
(a) \( (c\mu)^{\lambda}=c^{\lambda}\mu^{\lambda} \), so \( a(c)=c^{\lambda} \) and \( b(c)=0 \); and
\( \log(c\mu)=\log\mu+\log c \), so \( a(c)=1 \) and \( b(c)=\log c \).

(b) With \( G(t)=g(e^{t}) \) the relation reads \( G(t+s)=a(e^{s})G(t)+b(e^{s}) \) for all real \( t,s \).
Differentiating in \( t \) gives \( G'(t+s)=a(e^{s})G'(t) \); setting \( t=0 \) gives
\( G'(s)=a(e^{s})G'(0) \). If \( G'(0)=0 \) then \( G'\equiv0 \) and \( g \) is constant, contradicting strict
monotonicity, so \( G'(0)\ne0 \) and \( \psi=G'/G'(0) \) satisfies \( \psi(t+s)=\psi(t)\psi(s) \),
\( \psi(0)=1 \). Since \( G \) is twice continuously differentiable, \( \psi \) is differentiable;
differentiating in \( s \) at \( s=0 \) gives \( \psi'(t)=\lambda\psi(t) \) with \( \lambda=\psi'(0) \), so
\( \psi(t)=e^{\lambda t} \) and \( G'(t)=G'(0)e^{\lambda t} \). Integrating,
\( G(t)=G(0)+G'(0)(e^{\lambda t}-1)/\lambda \) if \( \lambda\ne0 \), and \( G(t)=G(0)+G'(0)t \) if
\( \lambda=0 \). Writing \( \mu=e^{t} \) gives \( g(\mu)=A\mu^{\lambda}+B \) with \( A=G'(0)/\lambda \), or
\( g(\mu)=A\log\mu+B \) with \( A=G'(0) \).

(c) If \( g(\mu_i)=\x_{(i)}\T\bbeta \) for all \( i \), then \( g(c\mu_i)=a(c)\x_{(i)}\T\bbeta+b(c) \), which
is \( \x_{(i)}\T\{a(c)\bbeta+b(c)\mathbf{e}_1\} \) when the first column of \( \X \) is \( \bone \). So the
rescaled means satisfy the same model with relabelled coefficients, and the set of mean curves the model can
produce is the same set. If \( g \) is not scale-equivariant, the two sets differ, and which model fits
better depends on the units — a property no scientific model should have.

(d) Suppose \( c\mu-1/(c\mu)=a(\mu-1/\mu)+b \) for all \( \mu>0 \). The three functions \( \mu \),
\( 1/\mu \) and \( 1 \) are linearly independent on \( (0,\infty) \), so \( a=c \), \( a=1/c \) and
\( b=0 \), forcing \( c^2=1 \). Only \( c=1 \) works. By (b) this had to happen: \( \mu-1/\mu \) is a sum of
two powers and not of the permitted form.

The conclusion also explains why the binomial links escape: a probability is a pure number, so there are no
units to change, and the logit and complementary log–log are free to be what the science suggests.
:::

::: {#exr-glm-exact-transform}
[C2]

@exr-glm-link-vs-transform compares "transform the data" with "transform the mean" to first order in
\( \phi \). For some families the comparison is exact.

::: {.enumerate options="label=(\alph*)"}
1. Suppose \( Y=\mu W \) where the distribution of \( W \) does not depend on \( \mu \), \( \E W=1 \) and
   \( \E|\log W|<\infty \) — a **scale family**. Show that \( \E(\log Y)=\log\mu+\kappa \) with
   \( \kappa=\E\log W \), a constant, and that \( \kappa<0 \) unless \( W\equiv1 \).

2. Deduce that if \( \log\mu_i=\x_{(i)}\T\bbeta \) and \( \X \) has an intercept, the population least
   squares regression of \( \log Y \) on \( \X \) has coefficients \( \bbeta \) with the intercept lowered by
   exactly \( |\kappa| \): the slopes agree *exactly*, not merely to first order.

3. Show that the gamma of @exm-glm-gamma with shape \( \nu=1/\phi \) is a scale family, and that
   \( \kappa=\psi(\nu)-\log\nu \), where \( \psi=\Gamma'/\Gamma \) is the digamma function. (Differentiate
   \( \int_0^\infty g^{\nu-1}e^{-g}\,dg=\Gamma(\nu) \) with respect to \( \nu \).) The standard expansion
   \( \psi(\nu)=\log\nu-1/(2\nu)-1/(12\nu^2)+\dots \) then returns \( \kappa\approx-\phi/2 \).

4. Show that no such exact statement is available for a Poisson response. Show that the Poisson is not a
   scale family, that \( \E(\log Y)=-\infty \) whenever \( \Pr(Y=0)>0 \), and that even the patched
   \( \log(1+Y) \) has \( \E\log(1+Y)-\log\mu\to\infty \) as \( \mu\to0 \) and
   \( \E\log(1+Y)-\log\mu\le\log(1+1/\mu)\to0 \) as \( \mu\to\infty \), so the gap is far from constant.
:::

:::

::: {.solution}
(a) \( \log Y=\log\mu+\log W \) and \( \log W \) has a law free of \( \mu \), so
\( \E\log Y=\log\mu+\kappa \) with \( \kappa=\E\log W \). By Jensen's inequality for the strictly concave
\( \log \), \( \kappa=\E\log W<\log\E W=0 \) unless \( W \) is degenerate, in which case \( W\equiv1 \) and
\( \kappa=0 \).

(b) The vector with entries \( \E\log Y_i=\x_{(i)}\T\bbeta+\kappa \) is \( \X\bbeta^{*} \) with
\( \bbeta^{*}=\bbeta+\kappa\mathbf{e}_1 \), because the first column is \( \bone \). The population least
squares coefficients are those of the best linear predictor of \( \log Y \), which here reproduces the mean
exactly; by @prp-lm-misspecified they are \( \bbeta^{*} \). Every coefficient but the intercept is
\( \beta_j \) itself, and \( \kappa<0 \).

(c) A gamma response with mean \( \mu \) and shape \( \nu \) is \( Y=(\mu/\nu)G \) with
\( G\sim\text{Gamma}(\nu,1) \), so \( W=G/\nu \) has \( \E W=1 \) and a law free of \( \mu \).
Differentiating \( \int_0^\infty g^{\nu-1}e^{-g}\,dg=\Gamma(\nu) \) under the integral sign gives
\( \int_0^\infty(\log g)g^{\nu-1}e^{-g}\,dg=\Gamma'(\nu) \), that is \( \E\log G=\psi(\nu) \). Hence
\( \kappa=\E\log G-\log\nu=\psi(\nu)-\log\nu \), which is negative, as (a) requires. With
\( \nu=1/\phi \) the expansion gives \( \kappa=-\phi/2-\phi^2/12-\dots \), the leading term of
@exr-glm-link-vs-transform, and the next term is smaller than it by the factor \( 6/\phi \), which for
Engel's fit is several hundred.

(d) For a scale family \( \Var(Y)/(\E Y)^2=\Var(W) \) does not depend on \( \mu \); for the Poisson it is
\( \mu/\mu^2=1/\mu \), so the Poisson is not one. Since \( \Pr(Y=0)=e^{-\mu}>0 \) and \( \log0=-\infty \),
\( \E\log Y=-\infty \) and the transformed regression is not even defined. For the patch, on the one hand
\( 0\le\E\log(1+Y)\le\E Y=\mu\to0 \) as \( \mu\to0 \) by \( \log(1+y)\le y \), while \( \log\mu\to-\infty \),
so the gap diverges; on the other hand Jensen's inequality gives
\( \E\log(1+Y)\le\log(1+\mu) \), so the gap is at most \( \log(1+1/\mu) \), which tends to \( 0 \). A gap
that ranges over the whole positive half-line cannot be absorbed into an intercept, so the two regressions
have genuinely different slopes, and only the link version answers a question about \( \E(Y) \).
:::
