# Bayesian thread IV: what carries over

The Bayesian thread has three instalments behind it.
[Section 7.6](../ch07-optimality/06-bayes-conjugate.html) put a normal-inverse-gamma prior
on \( (\bbeta,\sigma^2) \) in the normal linear model and found the posterior in closed
form (@thm-opt-bayes-conjugate): precisions add, the posterior mean is a
precision-weighted average of the prior mean and the least squares
estimate (@prp-opt-shrinkage), and the improper prior
\( \pi(\bbeta,\sigma^2)\propto1/\sigma^2 \) reproduces least squares and its \( t \)
intervals exactly (@cor-opt-flat-prior).
[Section 12.6](../ch12-intervals-and-bands/06-bayes.html) summarized that posterior by
credible sets (@def-ci-credible), and
[Section 32.6](../ch32-linear-mixed-models/06-bayes.html) read a random effect as a prior,
showed that the BLUP is a posterior mean and REML a marginal posterior
mode (@prp-mix-bayes), and ran a Gibbs sampler to do what the closed form could not. This
short section says which of those survive the move to a generalized linear model, and why
the others do not.

## What survives

Three things carry over unchanged, being definitions rather than computations. **The
posterior is still prior times likelihood**: Bayes' theorem does not know what family
\( \y \) came from. **Credible sets still mean what they meant**: @def-ci-credible is a
statement about a posterior distribution, not about a normal one. **A random effect is
still a prior**: the generalized linear mixed models of
[Chapter 40](../ch40-glmm-gee/index.html) are @eq-mix-hierarchy with a nonnormal first
level, and the empirical Bayes reading (@prp-shr-empirical-bayes) survives with them.

## What does not, and why

Everything that made [Section 7.6](../ch07-optimality/06-bayes-conjugate.html) *closed
form* rests on one fact: the normal log-likelihood is a quadratic in \( \bbeta \), so
prior times likelihood is again the exponential of a quadratic and completing the square
finishes the job. Nothing weaker will do.

::: {#prp-prc-quadratic}
[Why the conjugate machinery stops at the normal family]

Let \( \y \) follow the exponential dispersion family of @def-glm-model with the canonical
link, so that \( \theta_i=\x_{(i)}\T\bbeta \) and
\( \ell(\bbeta)=\phi^{-1}\sum_iw_i\{y_i\theta_i-b(\theta_i)\}+c \). Suppose \( b \) is
twice continuously differentiable on an open interval \( \Theta \). Then the following are
equivalent.

::: {.enumerate options="label=(\alph*)"}
1. For every \( n \), every \( \X \) with \( \rank(\X)=p \) and every \( \y \), the
   function \( \bbeta\mapsto\ell(\bbeta) \) is a quadratic polynomial on
   \( \{\bbeta:\x_{(i)}\T\bbeta\in\Theta\} \).

2. \( b''  \) is constant on \( \Theta \).

3. The variance function is constant: \( V(\mu)=b''(\theta) \) does not depend on
   \( \mu \), so the family is normal (@exr-glm-variance-determines).
:::

Consequently, for every family except the normal, a normal prior on \( \bbeta \) gives a
posterior that is not normal and has no closed-form normalizing constant.
:::

::: {.proof}
The Hessian of \( \ell \) is \( -\phi^{-1}\X\T\diag\{w_ib''(\x_{(i)}\T\bbeta)\}\X \)
by @thm-glm-score(d). If (b) holds this is constant in \( \bbeta \), and a function with
a constant Hessian is a quadratic polynomial, giving (a). Conversely, if (a) holds, take
\( p=1 \), \( n=1 \), \( w_1=1 \), \( x_{11}=1 \): the Hessian is the scalar
\( -b''(\beta)/\phi \), which must then be constant over \( \Theta \), which is (b).
The equivalence of (b) and (c) is the definition \( V(\mu)=b''(\theta) \)
of @def-glm-variance with \( \mu=b'(\theta) \), which makes \( \theta\mapsto\mu \) a
bijection: \( b''\equiv\sigma^2 \) gives \( b(\theta)=\sigma^2\theta^2/2 \) up to an
affine term, the normal cumulant function.

For the last sentence, with \( \pi \) the \( \Normal_p(\mathbf{m},\bSigma) \) density
the log posterior is \( \ell(\bbeta)-\tfrac12(\bbeta-\mathbf{m})\T\bSigma^{-1}
(\bbeta-\mathbf{m}) \) up to a constant; were it the log of a normal density it would be
a quadratic in \( \bbeta \), hence so would \( \ell \), which forces the normal family.
:::

So the failure is not a want of cleverness. It is the same failure that made the
likelihood equations of @thm-glm-score nonlinear and forced iteratively reweighted least
squares on the frequentist side: the score is not linear in \( \bbeta \), so the
posterior is not normal, and the integrals defining the posterior mean, the posterior
variance and the marginal likelihood
\[
m(\y)=\int L(\bbeta;\y)\,\pi(\bbeta)\,d\bbeta
\]{#eq-prc-marginal}

have no closed form.

## Conjugate priors do exist

It is worth being precise, because "generalized linear models have no conjugate priors"
is false as stated. Diaconis and Ylvisaker (1979) characterized the conjugate priors of an
exponential family; the canonical-link likelihood is an exponential family in
\( \bbeta \), and their construction gives the density
\[
\pi(\bbeta\mid\mathbf{m}_0,\kappa_0)\ \propto\
\exp\Bigl[\kappa_0^{-1}\bigl\{\mathbf{m}_0\T\X\bbeta
-\textstyle\sum_iv_i\,b(\x_{(i)}\T\bbeta)\bigr\}\Bigr]
\]{#eq-prc-dy}

in the form developed for generalized linear models by Bedrick, Christensen and Johnson
(1996) and Chen and Ibrahim (2003), with a pseudo-response vector \( \mathbf{m}_0 \) and
a prior weight \( v_i \) on each observation. It has the shape of the likelihood, so the
posterior belongs to the same family with \( \mathbf{m}_0 \) and the \( v_i \) updated
by addition — the "precisions add" structure of @eq-opt-posterior-params, the prior read as
\( \sum_iv_i \) observations at \( \mathbf{m}_0 \). What is missing is everything one
wants afterwards: outside the normal family, and outside degenerate designs such as the
intercept-only Poisson model of @exr-prc-dy-prior, @eq-prc-dy has no closed-form
normalizing constant, so neither have the posterior mean, the variance, the marginal
likelihood or a credible interval. In practice the normal prior of the next section — not
conjugate, but interpretable — is what is used.

::: {.idea}
Everything in the Bayesian thread that was a *definition* survives; everything that was a
*closed form* does not. The next two sections replace the closed form by an approximation
accurate when \( n \) is large and a simulation accurate when it is run long enough.
:::

## Exercises

### A. Check your understanding

::: {#exr-prc-quadratic-check}
[A1]

Verify directly, from the form of the Poisson log-likelihood with a log link, that it is
not a quadratic in \( \bbeta \) for \( p=1 \), and identify the term that spoils it.
:::

::: {.solution}
With \( \mu_i=e^{x_i\beta} \) the log-likelihood is
\( \sum_i\{y_ix_i\beta-e^{x_i\beta}\} \) up to a constant. The first term is linear; the
second has second derivative \( -\sum_ix_i^2e^{x_i\beta} \), which depends on
\( \beta \). That is the cumulant term \( b(\theta)=e^{\theta} \), whose \( b'' \) is
not constant, as @prp-prc-quadratic(b) requires.
:::

::: {#exr-prc-credible-invariance}
[A2]

An equal-tailed credible interval for \( \beta_j \) transforms to an equal-tailed credible
interval for \( e^{\beta_j} \) by exponentiating its endpoints. Is the same true of a
highest-posterior-density interval? Explain.
:::

### B. Practice

::: {#exr-prc-dy-prior}
[B1]

Take the Poisson family with a log link, \( p=1 \) and \( x_{i1}=1 \) for all \( i \), so
the model has a single parameter \( \beta \) and \( \mu=e^\beta \). Write out
@eq-prc-dy explicitly, show that it is the density of \( \beta=\log\mu \) when
\( \mu \) has a gamma distribution, and identify the shape and rate in terms of
\( \mathbf{m}_0 \) and the \( v_i \). Why does the closed form appear here and not in
general?
:::

::: {.solution}
With \( b(\theta)=e^\theta \), \( \theta=\beta \) and \( \kappa_0=1 \), @eq-prc-dy is
\( \pi(\beta)\propto\exp\{a\beta-ve^\beta\} \) with \( a=\bone\T\mathbf{m}_0 \) and
\( v=\sum_iv_i \). Substituting \( \mu=e^\beta \), so \( d\beta=d\mu/\mu \), the
density of \( \mu \) is proportional to \( \mu^{a-1}e^{-v\mu} \): a gamma with shape
\( a \) and rate \( v \), normalizing constant \( v^{a}/\Gamma(a) \). That constant is
available because the single linear predictor lets the integral factor into one dimension.
With \( p>1 \) the exponent
\( \mathbf{m}_0\T\X\bbeta-\sum_iv_ie^{\x_{(i)}\T\bbeta} \) is a sum of exponentials of
different linear forms, and no change of variable makes it a product of gamma integrals.
:::

### C. Going deeper

::: {#exr-prc-quadratic-noncanonical}
[C1]

@prp-prc-quadratic assumed the canonical link. Show that with a non-canonical link the
conclusion still holds in the following sense: if the log-likelihood is a quadratic in
\( \bbeta \) for every design and every response, then the working weights of
@thm-glm-irls do not depend on \( \bbeta \) and the mean is affine in \( \bbeta \).
Deduce that the model is the normal linear model up to a reparameterization.
:::
