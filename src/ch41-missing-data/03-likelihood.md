# Likelihood under missingness at random

A complete-case analysis throws away recorded values: a unit whose income is missing
still told us its age. Likelihood throws nothing away. It asks what the *observed*
data were likely to be, and the observed data include the partial records. This
section writes that likelihood down, proves the theorem that makes it usable, and
says what the theorem costs.

## The observed-data likelihood

The full model is @eq-mis-selection: a density for the complete array and a
mechanism. Since \( \mathbf{D}_{\mathrm{mis}} \) was never seen, the likelihood of
what was seen integrates it out,
\[
L_{\mathrm{full}}(\boldsymbol{\uptheta},\boldsymbol{\uppsi})
=\int f(\mathbf{D}_{\mathrm{obs}},\mathbf{D}_{\mathrm{mis}}\mid\boldsymbol{\uptheta})\,
 f(\R\mid\mathbf{D}_{\mathrm{obs}},\mathbf{D}_{\mathrm{mis}},\boldsymbol{\uppsi})\,
 d\mathbf{D}_{\mathrm{mis}} ,
\]{#eq-mis-full-likelihood}

a function of both parameters. Beside it stands the much simpler
\[
L(\boldsymbol{\uptheta})=f(\mathbf{D}_{\mathrm{obs}}\mid\boldsymbol{\uptheta})
=\int f(\mathbf{D}_{\mathrm{obs}},\mathbf{D}_{\mathrm{mis}}\mid\boldsymbol{\uptheta})\,
 d\mathbf{D}_{\mathrm{mis}} ,
\]{#eq-mis-observed-likelihood}

which ignores the mechanism entirely. The **observed-data likelihood**
@eq-mis-observed-likelihood is what software computes when told that some entries are
absent. The question is when using it is the same as using @eq-mis-full-likelihood.

::: {#thm-mis-ignorable}
[Ignorability]

Suppose that

::: {.enumerate options="label=(\roman*)"}
1. the mechanism is MAR at the realized \( (\R,\mathbf{D}_{\mathrm{obs}}) \) in the
   sense of @def-mis-mechanisms(b);

2. the parameters are **distinct**: the parameter space of
   \( (\boldsymbol{\uptheta},\boldsymbol{\uppsi}) \) is a product set
   \( \Theta\times\Psi \), so that knowing one places no restriction on the other.
:::

Then the full likelihood factors,
\[
L_{\mathrm{full}}(\boldsymbol{\uptheta},\boldsymbol{\uppsi})
=L(\boldsymbol{\uptheta})\cdot f(\R\mid\mathbf{D}_{\mathrm{obs}},\boldsymbol{\uppsi}),
\]{#eq-mis-factorization}

and consequently:

::: {.enumerate options="label=(\alph*)"}
1. the maximizer of \( L_{\mathrm{full}} \) in \( \boldsymbol{\uptheta} \) is the
   maximizer of \( L \), and likelihood ratio statistics for hypotheses about
   \( \boldsymbol{\uptheta} \) alone are the same under both;

2. if the prior is \( p(\boldsymbol{\uptheta},\boldsymbol{\uppsi})
   =p(\boldsymbol{\uptheta})p(\boldsymbol{\uppsi}) \), the posterior of
   \( \boldsymbol{\uptheta} \) given \( (\mathbf{D}_{\mathrm{obs}},\R) \) equals its
   posterior given \( \mathbf{D}_{\mathrm{obs}} \) alone, computed from
   @eq-mis-observed-likelihood;

3. if in addition the MAR condition holds at every \( (\R,\mathbf{D}_{\mathrm{obs}}) \),
   the observed-data score has mean zero and the usual maximum likelihood
   asymptotics (@thm-opt-mle) apply to \( L \) under the ordinary regularity
   conditions.
:::

A mechanism satisfying (i) and (ii) is called **ignorable**.
:::

::: {.proof}
Under (i) the factor \( f(\R\mid\mathbf{D},\boldsymbol{\uppsi}) \) in
@eq-mis-full-likelihood takes the same value for every \( \mathbf{D}_{\mathrm{mis}} \),
namely \( f(\R\mid\mathbf{D}_{\mathrm{obs}},\boldsymbol{\uppsi}) \), so it comes out
of the integral:
\[
\begin{aligned}
L_{\mathrm{full}}(\boldsymbol{\uptheta},\boldsymbol{\uppsi})
&=f(\R\mid\mathbf{D}_{\mathrm{obs}},\boldsymbol{\uppsi})
\int f(\mathbf{D}_{\mathrm{obs}},\mathbf{D}_{\mathrm{mis}}\mid\boldsymbol{\uptheta})
 \,d\mathbf{D}_{\mathrm{mis}}\\
&=f(\R\mid\mathbf{D}_{\mathrm{obs}},\boldsymbol{\uppsi})\,L(\boldsymbol{\uptheta}),
\end{aligned}
\]
which is @eq-mis-factorization.

(a) By (ii), \( \boldsymbol{\uptheta} \) ranges over \( \Theta \) whatever
\( \boldsymbol{\uppsi} \) is, so
\( \max_{\boldsymbol{\uppsi}}L_{\mathrm{full}}
=L(\boldsymbol{\uptheta})\max_{\boldsymbol{\uppsi}}f(\R\mid\mathbf{D}_{\mathrm{obs}},\boldsymbol{\uppsi}) \),
and the second factor does not involve \( \boldsymbol{\uptheta} \). Maximizing over
\( \boldsymbol{\uptheta} \) therefore maximizes \( L \). For a hypothesis
\( \boldsymbol{\uptheta}\in\Theta_0 \), both the numerator and the denominator of
the likelihood ratio carry the same second factor, which cancels.

(b) The posterior is proportional to
\( p(\boldsymbol{\uptheta})p(\boldsymbol{\uppsi})L_{\mathrm{full}} \). By
@eq-mis-factorization this is
\( \{p(\boldsymbol{\uptheta})L(\boldsymbol{\uptheta})\}\cdot
\{p(\boldsymbol{\uppsi})f(\R\mid\mathbf{D}_{\mathrm{obs}},\boldsymbol{\uppsi})\} \),
a product of a function of \( \boldsymbol{\uptheta} \) and a function of
\( \boldsymbol{\uppsi} \). Integrating out \( \boldsymbol{\uppsi} \) multiplies the
first bracket by a constant, so the marginal posterior of
\( \boldsymbol{\uptheta} \) is proportional to
\( p(\boldsymbol{\uptheta})L(\boldsymbol{\uptheta}) \).

(c) If (i) holds at every realization, then at every realization the joint density
of the observed data \( (\mathbf{D}_{\mathrm{obs}},\R) \) is exactly
\( L(\boldsymbol{\uptheta})\,f(\R\mid\mathbf{D}_{\mathrm{obs}},\boldsymbol{\uppsi}) \),
whose logarithm is \( \log L(\boldsymbol{\uptheta}) \) plus a term free of
\( \boldsymbol{\uptheta} \). So
\( \partial\log L/\partial\boldsymbol{\uptheta} \) is the \( \boldsymbol{\uptheta} \)
component of the score of a correctly specified density for the data actually seen,
and has mean zero under the ordinary regularity conditions of @thm-opt-mle.
:::

::: {.remark}
[Where each hypothesis is used]

Condition (i) alone gives @eq-mis-factorization and part (b): a Bayesian who
conditions on the data in hand needs nothing more. Condition (ii) is what makes the
factorization useful to a frequentist, who maximizes; without it the two factors
share parameters and neither can be maximized alone. Part (c) needs the everywhere
version of MAR, because a sampling distribution averages over data sets that did not
occur. Textbooks usually state "MAR plus distinctness" without separating these
roles; the separation is due to Seaman, Galati, Jackson and Carlin (2013).
:::

::: {.idea}
Under MAR with distinct parameters, the analyst may forget the mechanism and
maximize the likelihood of what was recorded. That is the entire licence for the
methods of the next two sections — and it is a licence to use a likelihood that is
usually an integral, not a product.
:::

## What the licence costs

Ignorability licenses ignoring the mechanism, not ignoring the missing values. Three
consequences follow.

**A model is now needed for variables previously conditioned on.** In a complete data
set, least squares never models the distribution of \( \X \): the likelihood
conditions on it. When a covariate is missing, @eq-mis-observed-likelihood integrates
over it, and that integral needs a density. Missing responses cost nothing extra;
missing covariates cost a model for the covariates.

**When only the response is missing, the partial records add nothing.** If every
covariate is recorded, the integral in @eq-mis-observed-likelihood over a missing
\( y_i \) is one, and the observed-data likelihood is exactly the complete-case
likelihood. This is @prp-mis-complete-case(b) again, as a statement about
likelihoods, and it is the situation of @thm-dsn-missing.

**The likelihood is no longer a product of simple factors.** Each incomplete record
contributes a marginal density, obtained by integrating the complete-data density
over the entries it lacks. For the multivariate normal those integrals are marginals
of a normal; for almost anything else they have no closed form.
[Section 41.4](04-em.html) is the standard response.

For a **monotone** pattern the integrals are free: integrating out a record's
unrecorded tail deletes exactly its conditional factors, so the likelihood becomes a
chain of conditional models on nested subsamples, each maximized on its own; that is
@exr-mis-monotone of [Section 41.1](01-mechanisms.html). The bivariate case is the
subject in miniature. With \( u \) always recorded and \( v \) recorded on the
first \( k \) units, the estimate of \( \mu_v \) is
\( \bar v_{(k)}+\hat\beta\{\bar u_{(n)}-\bar u_{(k)}\} \) by
@exr-mis-normal-monotone below: the complete-case mean of \( v \), corrected by the
shift in the mean of \( u \), using the \( n-k \) records that carry a \( u \)
and no \( v \). Those records were useless to the complete-case analysis and are
decisive here — and they repair the damage only to the extent that \( u \) predicts
\( v \), a remark that returns, with teeth, in
[Section 41.5](05-multiple-imputation.html). Anderson (1957) solved the general
monotone multivariate normal case this way.

## What changes under MNAR

Drop the MAR assumption and @eq-mis-factorization fails: the integral in
@eq-mis-full-likelihood cannot release the mechanism, and the two parameters are
entangled. The mechanism must then be modelled, and the model is not checkable, since
@prp-mis-untestable says that varying the assumed conditional distribution of the
missing values traces out a family with identical observed-data likelihoods. Whatever
identification such a model achieves comes from functional-form assumptions rather
than from the data, and the likelihood surface is often nearly flat in the directions
that matter. [Section 41.6](06-sensitivity.html) takes this up as a problem of
reporting rather than of estimation.

## Exercises

### A. Check your understanding

::: {#exr-mis-distinct-fails}
[A1]

Give a mechanism under which MAR holds but the parameters are not distinct, and
show that @eq-mis-factorization still holds while part (a) of @thm-mis-ignorable
can fail. (Hint: let the probability of recording \( y_i \) be a known function of
\( \bbeta \).)
:::

::: {.solution}
Take \( \Pr(R_i=1\mid\x_{(i)})=\Phi(\x_{(i)}\T\bbeta) \), with no free
\( \boldsymbol{\uppsi} \) at all. The mechanism depends only on recorded covariates,
so MAR holds and the factorization @eq-mis-factorization is valid, but the second
factor now depends on \( \bbeta \). Maximizing the product is not the same as
maximizing \( L(\bbeta) \), and the indicators carry information about
\( \bbeta \): the parameter spaces are not a product set.
:::

### B. Practice

::: {#exr-mis-normal-monotone}
[B1]

Let \( (u_i,v_i) \) be a bivariate normal sample with \( u \) recorded for all
\( n \) units and \( v \) for the first \( k \), the mechanism depending on
\( u \) but not on \( v \). Reparameterizing by \( (\mu_u,\sigma_u^2) \) and the
regression parameters of \( v \) on \( u \), derive the maximum likelihood
estimates and the formula for \( \hat\mu_v \) quoted above. Extend to three
variables with a monotone pattern: \( u \) recorded \( n \) times, \( v \) the
first \( k_2 \), \( z \) the first \( k_3\le k_2 \).
:::

::: {.solution}
Factor as \( f(u)f(v\mid u) \); by @exr-mis-monotone an unrecorded \( v_i \) loses
its conditional factor, so \( L=\prod_{i\le n}f(u_i)\prod_{i\le k}f(v_i\mid u_i) \).
The two parameter groups are variation independent, so each product is maximized
alone: a normal sample of size \( n \), and a normal linear model on \( k \)
observations (@thm-opt-mle). Since \( \mu_v=\alpha+\beta\mu_u \) and
\( \hat\alpha=\bar v_{(k)}-\hat\beta\bar u_{(k)} \), and since a maximum
likelihood estimate transforms with the parameterization, \( \hat\mu_v \) is as
stated. For three variables the same argument on
\( f(u)f(v\mid u)f(z\mid u,v) \) gives a univariate fit on all \( n \), the
regression of \( v \) on \( u \) over the first \( k_2 \), and the regression of
\( z \) on \( (u,v) \) over the first \( k_3 \).
:::

::: {#exr-mis-covariate-model}
[B2]

A linear model has a single covariate \( x \) that is sometimes missing, with
\( x\sim\Normal(\mu,\tau^2) \) and \( y\mid x\sim\Normal(\beta_0+\beta_1x,\sigma^2) \).
Write the observed-data likelihood contribution of a record with \( y \) recorded
and \( x \) missing, and show that it is the density of a normal with mean
\( \beta_0+\beta_1\mu \) and variance \( \beta_1^2\tau^2+\sigma^2 \). Why does this
contribution carry information about \( \beta_1 \) even though it contains no
\( x \)?
:::

### C. Going deeper

::: {#exr-mis-ignorable-design}
[C1]

Show that @thm-dsn-missing is the special case of @thm-mis-ignorable in which the
missing entries are responses, the mechanism depends only on the design, and the
complete-data model is the normal linear model. Which part of @thm-dsn-missing
corresponds to @eq-mis-observed-likelihood, and which to the warning that the
degrees of freedom must be reduced?
:::

