# Generalized linear mixed models

Replace the scalar \( U_i \) of [Section 40.1](01-marginal-versus-conditional.html) by a
vector of random effects with a covariance matrix, and the single covariate by a model
matrix, and the result is the generalized linear mixed model: the linear mixed model of
[Chapter 32](../ch32-linear-mixed-models/index.html) with a link function and an exponential
dispersion family in place of the normal error.

That one change costs more than it looks. In @def-mix-model the marginal distribution was
still normal, with covariance \( \V=\Z\G\Z\T+\R \), so everything Part VII knew about
generalized least squares applied unchanged. Here it is an integral with no closed form, and
the marginal moments have to be computed one family at a time.

## The model

::: {#def-gmm-glmm}
[Generalized linear mixed model]

Let \( \bu\sim\Normal_q\{\bzero,\G(\boldsymbol{\uptheta})\} \) be a vector of **random
effects** with covariance depending on unknown variance components
\( \boldsymbol{\uptheta} \). Given \( \bu \), let the responses \( Y_1,\dots,Y_n \) be
independent, each from the exponential dispersion family @def-glm-edf with prior weight
\( w_i \), dispersion \( \phi \) and variance function \( V \), and with conditional mean
\[
\mu_i(\bu)=\E(Y_i\mid\bu)=h(\eta_i),\qquad
\eta_i=\x_{(i)}\T\bbeta+\bz_{(i)}\T\bu ,
\]{#eq-gmm-model}

so that \( \Var(Y_i\mid\bu)=\phi V\{\mu_i(\bu)\}/w_i \). In matrix form
\( \boldsymbol{\upeta}=\X\bbeta+\Z\bu \) with \( \X \) of size \( n\times p \) and
\( \Z \) of size \( n\times q \), both known. The parameters are the **fixed effects**
\( \bbeta \), the variance components \( \boldsymbol{\uptheta} \) and, if it is not fixed by
the family, the dispersion \( \phi \).
:::

Setting \( \Z=\mathbf{0} \) recovers the generalized linear model @def-glm-model exactly;
taking the family normal and \( h \) the identity, with \( \R=\phi\I \),
recovers @def-mix-model. Two structures cover most applications.

**Clustered data.** Units \( i=1,\dots,m \) with \( n_i \) observations each, one random
effect vector \( \bu_i\sim\Normal_{q_0}(\bzero,\G_0) \) per unit, independent across units.
Then \( \Z \) is block diagonal, \( q=mq_0 \), and both the likelihood and the marginal
covariance break into \( m \) independent pieces. The **random-intercept** model takes
\( q_0=1 \) and \( \bz_{(ij)}=1 \):
\[
g(\mu_{ij})=\x_{(ij)}\T\bbeta+u_i,\qquad u_i\iid\Normal(0,\tau^{2}).
\]{#eq-gmm-random-intercept-model}

**Crossed or nested effects.** Pupils in classes in schools; raters and items. Here
\( \G \) is block diagonal in the *effects* rather than the units and the likelihood does
not factorize. [Chapter 32](../ch32-linear-mixed-models/index.html) described the linear
versions (@exm-mix-nested-crossed); the difficulty a link function adds is computational,
and severe, because the integral in [Section 40.3](03-integrated-likelihood.html) is then
genuinely high-dimensional.

::: {.remark}
[Two extensions, and why they are rarer]

A cluster-level random effect can be put into the linear predictor of a baseline-category
model @def-mlt-baseline or of a cumulative-logit model @def-mlt-cumulative, giving mixed
models for nominal and ordinal responses; in the ordinal case the effect shifts every
cutpoint together, leaving the proportional-odds structure intact. Conditional independence
given \( \bu_i \) can also be relaxed, by letting the random effects themselves follow an
AR(1) process within the cluster — the mixed-model counterpart of the non-diagonal working
correlation of [Section 40.5](05-working-correlation.html). That is rare, because it raises
each cluster's integral from dimension \( q_0 \) to \( n_i \).
:::

::: {#exm-gmm-random-intercept-logit}
[A random-intercept logistic model]

\( m \) physicians each treat \( n_i \) patients; \( Y_{ij}=1 \) if patient \( j \) of
physician \( i \) recovers. With the logistic link of @def-bin-logistic,
\( \operatorname{logit}\Pr(Y_{ij}=1\mid u_i)=\x_{(ij)}\T\bbeta+u_i \) and
\( u_i\sim\Normal(0,\tau^{2}) \), the coefficient of a patient-level covariate compares two
patients *of the same physician*, and \( \tau \) measures how much physicians differ. A
useful summary is the **median odds ratio** \( \exp(\sqrt2\,\tau\,z_{0.75}) \), the median
odds ratio between two randomly chosen physicians treating identical patients (@exr-gmm-mor).
Replacing \( u_i \) by \( u_{0i}+u_{1i}x_{ij} \) gives a random-slope
model, at the price of three parameters in \( \G_0 \) and, by @exr-gmm-probit-general, no
single marginal coefficient at all.
:::

::: {#exm-gmm-poisson-lognormal}
[The Poisson–lognormal model]

With \( Y_i\mid u_i \) Poisson (@def-cnt-poisson) of mean
\( \exp(\x_{(i)}\T\bbeta+u_i) \) and \( u_i\sim\Normal(0,\tau^2) \), the conditional mean
is multiplied by the lognormal variable \( e^{u_i} \). Unlike the binary case this model is
identified with a *single* observation per unit, because the Poisson variance is determined
by its mean and any excess must come from \( u_i \): it is an overdispersion model in the
sense of [Chapter 38](../ch38-quasi-likelihood/index.html) (@prp-ql-overdispersion), and a
direct competitor to the negative binomial of @def-cnt-negbin, which arises the same way
with a gamma mixing variable (@lem-cnt-gamma-mixture).
:::

## The moments a random effect induces

::: {#prp-gmm-moments}
[Marginal moments of a GLMM]

In @def-gmm-glmm write \( \mu_i(\bu)=h(\eta_i) \) and
\( \mu_i^{M}=\E\{\mu_i(\bu)\} \).

::: {.enumerate options="label=(\alph*)"}
1. \( \E(Y_i)=\mu_i^{M} \),
   \[
   \Var(Y_i)=\frac{\phi}{w_i}\E\bigl[V\{\mu_i(\bu)\}\bigr]+\Var\{\mu_i(\bu)\},
   \]{#eq-gmm-marginal-variance}

   and \( \Cov(Y_i,Y_k)=\Cov\{\mu_i(\bu),\mu_k(\bu)\} \) for \( i\ne k \). The marginal
   variance always exceeds the average conditional variance.

2. **Binomial.** If \( Y_i\mid\bu\sim\text{Bin}\{N_i,\pi_i(\bu)\} \) and
   \( \bar\pi_i=\E\pi_i(\bu) \), then \( \E(Y_i)=N_i\bar\pi_i \) and
   \[
   \Var(Y_i)=N_i\bar\pi_i(1-\bar\pi_i)
   \Bigl\{1+(N_i-1)\,\frac{\Var\{\pi_i(\bu)\}}{\bar\pi_i(1-\bar\pi_i)}\Bigr\}.
   \]{#eq-gmm-binomial-overdispersion}

   The overdispersion factor is the beta-binomial one, and it equals \( 1 \) when
   \( N_i=1 \): a *single* Bernoulli response carries no trace of the random effect.

3. **Poisson with a log link**, \( \phi=1 \), and a random intercept
   \( u_i\sim\Normal(0,\tau^{2}) \) shared by the observations of a cluster. Then
   \( \mu_{ij}^{M}=\exp(\x_{(ij)}\T\bbeta+\tau^{2}/2) \),
   \[
   \begin{aligned}
   \Var(Y_{ij})&=\mu_{ij}^{M}+\bigl(\mu_{ij}^{M}\bigr)^{2}\bigl(e^{\tau^{2}}-1\bigr),\\
   \Cov(Y_{ij},Y_{ik})&=\mu_{ij}^{M}\mu_{ik}^{M}\bigl(e^{\tau^{2}}-1\bigr).
   \end{aligned}
   \]{#eq-gmm-poisson-moments}

   The marginal variance function is *exactly* the NB2 function
   \( V(\mu)=\mu+\mu^{2}/\kappa \) of @def-cnt-negbin, with
   \( \kappa=1/(e^{\tau^{2}}-1) \) — although the marginal distribution is not negative
   binomial and has no closed form.
:::

:::

::: {.proof}
(a) Condition on \( \bu \) and use the variance decomposition
\( \Var(Y)=\E\{\Var(Y\mid\bu)\}+\Var\{\E(Y\mid\bu)\} \) and, for \( i\ne k \),
\( \Cov(Y_i,Y_k)=\E\{\Cov(Y_i,Y_k\mid\bu)\}+\Cov\{\E(Y_i\mid\bu),\E(Y_k\mid\bu)\} \), whose
first term vanishes by conditional independence. The final claim is
\( \Var\{\mu_i(\bu)\}\ge0 \).

(b) With \( \pi=\pi_i(\bu) \), @eq-gmm-marginal-variance gives
\( \Var(Y_i)=N_i\E\{\pi(1-\pi)\}+N_i^{2}\Var(\pi) \). Now
\( \E\{\pi(1-\pi)\}=\bar\pi_i-\E(\pi^{2})=\bar\pi_i(1-\bar\pi_i)-\Var(\pi) \), so
\( \Var(Y_i)=N_i\bar\pi_i(1-\bar\pi_i)+N_i(N_i-1)\Var(\pi) \),
which is @eq-gmm-binomial-overdispersion.

(c) For a lognormal factor, \( \E e^{su}=e^{s^{2}\tau^{2}/2} \). Hence
\( \mu_{ij}^{M}=e^{\eta_{ij}^{0}}\E e^{u_i}=e^{\eta_{ij}^{0}+\tau^{2}/2} \), writing
\( \eta_{ij}^{0}=\x_{(ij)}\T\bbeta \). By @eq-gmm-marginal-variance with
\( V(\mu)=\mu \) and \( \phi=w_{ij}=1 \),
\[
\Var(Y_{ij})=\mu_{ij}^{M}+e^{2\eta_{ij}^{0}}\bigl(\E e^{2u_i}-(\E e^{u_i})^{2}\bigr)
=\mu_{ij}^{M}+e^{2\eta_{ij}^{0}}\bigl(e^{2\tau^{2}}-e^{\tau^{2}}\bigr),
\]
and \( e^{2\eta^{0}}(e^{2\tau^{2}}-e^{\tau^{2}})
=(e^{\eta^{0}+\tau^{2}/2})^{2}(e^{\tau^{2}}-1) \). The covariance is the same computation
with \( \eta_{ij}^{0}+\eta_{ik}^{0} \) in the exponent. Matching
\( \mu+\mu^{2}(e^{\tau^{2}}-1) \) with \( \mu+\mu^{2}/\kappa \) gives \( \kappa \).
:::

Part (b) explains why a random-intercept model for *ungrouped* binary data needs repeated
observations: with \( N_i=1 \) the marginal distribution is Bernoulli whatever \( \tau \)
is, and \( \tau \) is identified only through the within-cluster *correlation*. Part (c)
explains why the Poisson–lognormal and the negative binomial are so hard to tell apart: they
agree in their first two moments exactly and differ only in shape, so a likelihood ratio
between them compares tails, not overdispersion. The listing checks @eq-gmm-poisson-moments
by simulation at \( \tau=0.6 \), where the matching NB2 shape is
\( \kappa=1/(e^{0.36}-1)=2.3077 \).

```{.python .run #cell-glmm-grunfeld-poisson}
import numpy as np

rng = np.random.default_rng(40_002)
tau, eta1, eta2 = 0.6, 0.4, -0.2
u = rng.normal(0, tau, 1_000_000)
y1 = rng.poisson(np.exp(eta1 + u))
y2 = rng.poisson(np.exp(eta2 + u))
mu1, mu2 = np.exp(eta1 + tau**2 / 2), np.exp(eta2 + tau**2 / 2)

print(f"E Y1      : {y1.mean():.4f}  (theory {mu1:.4f})")
print(f"Var Y1    : {y1.var():.4f}  (theory {mu1 + mu1**2 * (np.exp(tau**2) - 1):.4f})")
print(f"Cov(Y1,Y2): {np.cov(y1, y2)[0, 1]:.4f}  (theory {mu1 * mu2 * (np.exp(tau**2) - 1):.4f})")
```

## What a random intercept cannot do

::: {#lem-gmm-association}
[A shared random intercept induces nonnegative correlation]

Let \( U \) be a random variable and \( f,g:\Real\to\Real \) both nondecreasing, with
\( \E f(U)^{2} \) and \( \E g(U)^{2} \) finite. Then
\( \Cov\{f(U),g(U)\}\ge0 \). Consequently, in a random-intercept GLMM @eq-gmm-random-intercept-model
with any inverse link \( h \) that is nondecreasing,
\[
\Cov(Y_{ij},Y_{ik})=\Cov\bigl\{h(\eta_{ij}^{0}+u_i),\,h(\eta_{ik}^{0}+u_i)\bigr\}\ge0
\qquad (j\ne k).
\]
:::

::: {.proof}
Let \( U' \) be an independent copy of \( U \). Because both functions are nondecreasing,
\( \{f(U)-f(U')\}\{g(U)-g(U')\}\ge0 \) pointwise. Taking expectations and expanding,
\[
0\le\E\bigl[\{f(U)-f(U')\}\{g(U)-g(U')\}\bigr]=2\bigl[\E\{f(U)g(U)\}-\E f(U)\,\E g(U)\bigr]
=2\Cov\{f(U),g(U)\}.
\]
The second claim is @prp-gmm-moments(a) with \( f(u)=h(\eta_{ij}^{0}+u) \) and
\( g(u)=h(\eta_{ik}^{0}+u) \), both nondecreasing.
:::

Every link in this book has a nondecreasing inverse, so the consequence is a real
restriction: **a model whose only random effect is a cluster-level intercept cannot describe
negatively correlated responses.** Such data exist. If exactly one of a traveller's four
options is chosen, or exactly one of a litter survives a fixed food supply, the responses
compete and their correlation is negative by construction.
[Section 40.4](04-gee.html) analyses a data set of that kind, and the estimating-equation
approach has no difficulty with it precisely because it never commits to a mechanism for the
dependence.

## Estimation, inference and prediction

Fitting is the subject of [Section 40.3](03-integrated-likelihood.html). The random effects
cannot be estimated by maximizing a joint likelihood, because their number grows with the
sample size — the incidental-parameter trap of @exm-opt-neyman-scott. They are integrated
out, giving the **integrated likelihood**
\[
L(\bbeta,\boldsymbol{\uptheta},\phi)
=\int_{\Real^{q}}\prod_{i=1}^{n}f\bigl(y_i\mid\bu;\bbeta,\phi\bigr)\,
p(\bu;\boldsymbol{\uptheta})\,d\bu ,
\]{#eq-gmm-integrated}

which is maximized numerically. Inference for \( \bbeta \) is then the usual likelihood
inference — Wald intervals from the inverse observed information, likelihood ratio tests
between nested mean models — as \( m\to\infty \) with cluster sizes bounded. That claim is
imported, not proved here: @thm-glm-asymptotics covers a generalized linear model, not a
GLMM, and the conditions under which the integrated likelihood behaves like an ordinary one
are set out by McCulloch, Searle and Neuhaus (2008, chapter 14). Two warnings carry over
from Part VII unchanged. First, a test of \( \tau^{2}=0 \), or of any variance component on
the boundary of its range, is *not* a \( \chi^{2}(1) \) test: the limiting null distribution
is the mixture \( \tfrac12\chi^{2}(0)+\tfrac12\chi^{2}(1) \) under the conditions of
@thm-mix-boundary, so the naive \( p \)-value is twice what it should be. Second, REML has
no clean analogue, because @def-mix-reml was built on the distribution of
\( \mathbf{K}\T\Y \) for \( \C(\mathbf{K})=\C(\X)\perpc \), a construction that needs the
normal linear model; the approximate versions in use inherit the bias of the working linear
model they are built on.

In the linear mixed model @thm-mix-blup produced a linear predictor of the random effects,
and @eq-mix-shrinkage showed it shrinks the cluster's own estimate towards zero. Neither
linearity nor an explicit formula survives the link function, but the shrinkage does.

::: {#prp-gmm-prediction}
[Empirical Bayes prediction and shrinkage]

Consider @eq-gmm-random-intercept-model with \( \bbeta \), \( \tau^{2} \) and \( \phi \)
known, write \( \ell_i(u)=\sum_{j}\log f(y_{ij}\mid u) \) for cluster \( i \)'s conditional
log-likelihood, and assume \( \ell_i \) is concave in \( u \) — as it is for a canonical
link by @thm-glm-concave, and for the gamma with a log link by direct differentiation.

::: {.enumerate options="label=(\alph*)"}
1. The posterior density of \( u_i \) given \( \y_i \), proportional to
   \( \exp\{\ell_i(u)-u^{2}/(2\tau^{2})\} \), is log-concave and its mode \( \hat u_i \) is
   unique.

2. If \( \ell_i \) attains its maximum at a finite \( \hat a_i \ne0 \), then
   \( \hat u_i \) lies strictly between \( 0 \) and \( \hat a_i \): the prediction has the
   same sign as the cluster's own estimate and is strictly smaller in magnitude. If
   \( \ell_i \) has no finite maximizer — as when every response in a binary cluster is
   \( 0 \) — then \( \hat u_i \) is finite all the same.

3. Replacing \( \ell_i \) by its quadratic expansion at \( \hat a_i \) gives
   \( \hat u_i\approx\{\tau^{2}/(\tau^{2}+I_i^{-1})\}\,\hat a_i \) with
   \( I_i=-\ell_i''(\hat a_i) \). In the two cases used in this chapter
   \( I_i^{-1}=\phi/n_i \) exactly, which is @eq-mix-shrinkage: clusters with few
   observations are shrunk harder, and all clusters are shrunk harder when \( \tau^{2} \)
   is small relative to the noise.
:::

:::

::: {.proof}
(a) Adding the strictly concave \( -u^{2}/(2\tau^2) \) to the concave \( \ell_i \) keeps it
concave and makes it strictly so, and the sum tends to \( -\infty \) in both directions, so
a unique maximizer exists.

(b) Write \( \psi(u)=\ell_i(u)-u^{2}/(2\tau^{2}) \) and suppose \( \hat a_i>0 \). Then
\( \ell_i'(\hat a_i)=0 \) gives \( \psi'(\hat a_i)=-\hat a_i/\tau^{2}<0 \), while
\( \psi'(0)=\ell_i'(0)>0 \) because \( \ell_i \) is concave with maximum to the right of
\( 0 \). A continuous derivative positive at \( 0 \) and negative at \( \hat a_i \)
vanishes strictly between them, and by strict concavity that zero is the unique maximizer;
the case \( \hat a_i<0 \) is symmetric. If \( \ell_i \) increases everywhere,
\( \psi' \) is still eventually negative because of the penalty, so \( \hat u_i \) is
finite.

(c) With \( \ell_i(u)\approx\ell_i(\hat a_i)-\tfrac12 I_i(u-\hat a_i)^{2} \), the penalized
maximizer is \( \hat u_i=I_i\hat a_i/(I_i+\tau^{-2})=\{\tau^{2}/(\tau^{2}+I_i^{-1})\}\hat
a_i \). For a canonical link \( I_i=\sum_j w_{ij}V(\hat\mu_{ij})/\phi \), so
\( I_i^{-1}=\phi/(n_i\bar V_i) \) with \( \bar V_i \) the average of the
\( w_{ij}V(\hat\mu_{ij}) \); in the normal case \( V\equiv1 \), \( w\equiv1 \) and
\( \bar V_i=1 \). For the gamma with a log link \( \ell_i''(u)=-\nu\sum_j y_{ij}e^{-\eta_{ij}} \)
with \( \nu=1/\phi \), and \( \ell_i'(\hat a_i)=0 \) forces
\( \sum_j y_{ij}e^{-\eta_{ij}}=n_i \) there, so again \( I_i=n_i/\phi \) exactly.
:::

::: {.remark}
[The predictions are not estimates]

\( \hat u_i \) predicts a random variable, so its error is a prediction error variance and
the interval covering \( u_i \) is a *prediction* interval (@thm-ci-prediction-interval),
exactly as in @thm-mix-blup(b). The usual intervals ignore the estimation of \( \bbeta \)
and \( \boldsymbol{\uptheta} \) and are therefore too short when \( m \) is small.
:::

## An investment panel with a multiplicative mean

::: {#exm-gmm-grunfeld}
[Grunfeld's firms, again]

[Chapter 32](../ch32-linear-mixed-models/index.html) fitted a linear mixed model to
Grunfeld's panel of \( 11 \) US firms observed in each of
\( 20 \) years (@exm-mix-grunfeld). Investment there ranges over three orders of magnitude
and is far more variable for large firms than for small ones, which a constant-variance
model cannot represent; a gamma response with a log link can, because the mean is
multiplicative and the conditional coefficient of variation is constant. With
\( Y_{ij} \) the gross investment of firm \( i \) in year \( j \),
\[
Y_{ij}\mid u_i\sim\text{Gamma}\ \text{with}\quad
\log\E(Y_{ij}\mid u_i)=\beta_0+\beta_1\log(\text{value}_{ij})
+\beta_2\log(\text{capital}_{ij})+u_i ,
\]
\( u_i\sim\Normal(0,\tau^{2}) \), and \( \Var(Y_{ij}\mid u_i)=\phi\,\E(Y_{ij}\mid u_i)^{2} \).
Fitting by adaptive Gauss–Hermite quadrature
([Section 40.3](03-integrated-likelihood.html)) gives

| | \( \beta_1 \) (value) | \( \beta_2 \) (capital) | \( \phi \) | \( \tau \) |
|:---|---:|---:|---:|---:|
| marginal gamma GLM | 0.8421 | 0.3386 | 0.2644 | — |
| gamma GLMM | 0.6825 | 0.2315 | 0.0930 | 0.5188 |

The firm-to-firm standard deviation \( \hat\tau=0.5188 \) is large: firms
a standard deviation apart on the log scale differ by a factor of
\( e^{0.5188}\approx1.68 \) in investment at the same value and capital
stock. The conditional coefficient of variation is
\( \sqrt{\hat\phi}=0.3049 \). For a gamma with shape \( 1/\hat\phi \) the conditional
variance of \( \log Y \) is the trigamma function of the shape,
\( 0.0974 \), so the fraction of the variance of \( \log Y \)
attributable to the firm is
\( \hat\tau^{2}/(\hat\tau^{2}+0.0974)=0.734 \).

The fitted conditional intercept is \( -1.3032 \); by @prp-gmm-marginal(b)
the log link causes no attenuation of the slopes, and the marginal
intercept is \( -1.3032+\hat\tau^{2}/2=-1.1686 \). The difference between the
two rows is therefore *not* the marginal–conditional distinction. It is the older
within-versus-between problem of [Section 6.6](../ch06-projections/06-fwl.html): the
independence fit lets differences between firms contribute to the slopes, while the mixed
model assigns most of them to \( u_i \) and estimates the slopes mainly from movement
within a firm. A gap of this size accuses the mean model, and
[Section 40.5](05-working-correlation.html) returns to it with a diagnostic and two
candidate explanations, @exr-gmm-within-between with a remedy.
:::

::: {when-format="html"}
![**Figure 40.2.1.** Shrinkage of the firm effects in the gamma GLMM of @exm-gmm-grunfeld.
(a) \( \hat u_i \) against the firm's own estimate \( \hat a_i \), for 20 years (circles)
and for the first 3 only (triangles); the diagonal is no shrinkage. (b) The factor
\( \tau^{2}/(\tau^{2}+\phi/n_i) \) of @prp-gmm-prediction(c) at
\( \hat\phi \).](shrinkage.svg){#fig-gmm-shrinkage width=100%}
:::

::: {when-format="pdf"}
![Shrinkage of the firm effects in the gamma GLMM of @exm-gmm-grunfeld.
(a) \( \hat u_i \) against the firm's own estimate \( \hat a_i \), for 20 years (circles)
and for the first 3 only (triangles); the diagonal is no shrinkage. (b) The factor
\( \tau^{2}/(\tau^{2}+\phi/n_i) \) of @prp-gmm-prediction(c) at
\( \hat\phi \).](shrinkage.pdf){width=100%}
:::

With \( 20 \) years per firm and \( \hat\phi=0.0930 \) small there is
almost nothing to shrink: the empirical slope of \( \hat u_i \) on \( \hat a_i \) is
\( 0.9830 \), and @prp-gmm-prediction(c) predicts
\( 0.9830 \). Using only each firm's first three years, the same
parameters give \( 0.8951 \). Shrinkage measures how much a cluster's own data are worth,
and is not a fixed property of the model.

## Exercises

### A. Check your understanding

::: {#exr-gmm-identify}
[A1]

For each of the following, say whether the random-intercept variance \( \tau^{2} \) is
identified, and why: (a) one Bernoulli response per cluster; (b) one binomial response with
\( N_i=40 \) trials per cluster; (c) one Poisson count per cluster; (d) five Bernoulli
responses per cluster.
:::

::: {#exr-gmm-kappa}
[A2]

A Poisson–lognormal model is fitted and gives \( \hat\tau^{2}=0.5 \). What NB2 shape
parameter \( \kappa \) would produce the same marginal variance function? Would you expect
the two models to give noticeably different fitted means?
:::

### B. Practice

::: {#exr-gmm-mor}
[B1]

In @exm-gmm-random-intercept-logit, let two physicians be drawn at random and let
\( \text{OR} \) be the odds ratio between them for identical patients, taken in the
direction that makes it at least one. Show that the median of \( \text{OR} \) is
\( \exp(\sqrt2\,\tau\,z_{0.75}) \), where \( z_{0.75}=0.6745 \) is the upper
quartile of the standard normal. Evaluate it for \( \tau=1 \).
:::

::: {.solution}
The log odds ratio between physicians \( i \) and \( k \) is \( u_i-u_k \), which is
\( \Normal(0,2\tau^{2}) \). Taking it in the positive direction gives
\( \lvert u_i-u_k\rvert \), whose median is the upper quartile of
\( \Normal(0,2\tau^{2}) \), namely \( \sqrt2\,\tau z_{0.75} \). Exponentiating gives the
claim. For \( \tau=1 \), \( \sqrt2\times0.6745=0.9539 \) and
\( e^{0.9539}\approx2.6 \): a typical pair of physicians differs by a factor of about
\( 2.6 \) in the odds of recovery, for identical patients.
:::

::: {#exr-gmm-negative}
[B2]

Give a two-level model with a *vector* random effect whose induced within-cluster
correlation is negative, and explain why @lem-gmm-association does not apply. (Hint: let
\( \bz_{(ij)} \) differ in sign across \( j \).)
:::

::: {.solution}
Take \( n_i=2 \), \( q_0=1 \), \( z_{i1}=+1 \) and \( z_{i2}=-1 \), so
\( \eta_{i1}=\eta^{0}_{i1}+u_i \) and \( \eta_{i2}=\eta^{0}_{i2}-u_i \). Then
\( f(u)=h(\eta^{0}_{i1}+u) \) is nondecreasing and
\( g(u)=h(\eta^{0}_{i2}-u) \) is nonincreasing, so the association inequality applies with
\( -g \) and gives \( \Cov\{f(U),g(U)\}\le0 \). The lemma assumed both functions
nondecreasing, which is automatic only when every \( z_{(ij)} \) has the same sign — as for
an intercept.
:::

### C. Going deeper

::: {#exr-gmm-within-between}
[C1]

In @exm-gmm-grunfeld the mixed-model slopes differ markedly from the independence ones.
Add the firm means \( \overline{\log\text{value}}_{i} \) and
\( \overline{\log\text{capital}}_{i} \) to the model as extra fixed effects, so that the
original covariates carry only within-firm variation. Explain why this makes the random
effect orthogonal to the covariates in the population, refit, and compare. What test does
the difference between the two specifications suggest?
:::
