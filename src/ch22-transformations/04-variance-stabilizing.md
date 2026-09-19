# Variance-stabilizing transformations

For many responses the variance is tied to the mean by the nature of the measurement: counts, proportions,
amounts measured with a fixed *relative* precision. Then \( \Var(Y)=\phi\,v(\mu) \) for a known or plausible
variance function \( v \), where \( \mu=\E(Y) \) and \( \phi \) is a constant. The transformation that makes the variance
constant turns out to be determined by \( v \) up to an affine change.

## The idea

If \( g \) is smooth and \( Y \) is concentrated near \( \mu \), then \( g(Y)\approx g(\mu)+g'(\mu)(Y-\mu) \), so
\[
\Var g(Y)\approx g'(\mu)^2\,\phi\,v(\mu).
\]
This is constant in \( \mu \) exactly when \( g'(\mu) \) is proportional to \( v(\mu)^{-1/2} \), that is, when
\[
g(\mu)=\int^{\mu}\frac{dt}{\sqrt{v(t)}} .
\]{#eq-tr-stabilizing-integral}

Part (a) of the theorem below makes this precise, with \( m \) measuring the information; it is
@exr-cor-variance-stabilizing, where Fisher's \( z \) arose this way. Part (b) adds uniqueness.

::: {#thm-tr-variance-stabilizing}
[Variance-stabilizing transformations]

Let \( I \) be an open interval and \( v:I\to(0,\infty) \) continuous. Suppose that for every \( \theta\in I \) the random
variables \( T_m \) satisfy \( \sqrt m\,(T_m-\theta)\to\Normal(0,v(\theta)) \) in distribution as \( m\to\infty \). Fix
\( t_0\in I \) and let \( g(t)=\int_{t_0}^{t}v(s)^{-1/2}\,ds \) for \( t\in I \), extended in any way outside \( I \).

::: {.enumerate options="label=(\alph*)"}
1. For every \( \theta\in I \), \( \sqrt m\,\{g(T_m)-g(\theta)\}\to\Normal(0,1) \) in distribution.

2. If \( h:I\to\Real \) is continuously differentiable and \( \sqrt m\,\{h(T_m)-h(\theta)\}\to\Normal(0,1) \) for every \( \theta\in I \),
   then \( h=g+c \) or \( h=-g+c \) on \( I \) for a constant \( c \).
:::

:::

::: {.proof}
(a) By the fundamental theorem of calculus \( g'=v^{-1/2} \) on \( I \), so this is @exr-cor-variance-stabilizing with
\( s=v^{1/2} \); the values of \( g \) outside \( I \) do not matter because \( T_m\to\theta \) in probability.
(b) By the delta method
(stated in [Section 14.4](../ch14-correlation-lack-of-fit-prediction/04-correlation-coefficient.html)), \( \sqrt m\{h(T_m)-h(\theta)\}\to\Normal(0,h'(\theta)^2v(\theta)) \). Limits in distribution are
unique, so \( h'(\theta)^2v(\theta)=1 \) and \( \lvert h'(\theta)\rvert=v(\theta)^{-1/2} \) for every \( \theta\in I \). The continuous
function \( h' \) never vanishes on the interval \( I \), so by the intermediate value theorem it has one sign
throughout. Hence \( h'=g' \) or \( h'=-g' \) on \( I \), and integrating gives the result.
:::

Convergence in distribution does not by itself imply that the *variance* of \( g(T_m) \) is close to
\( 1/m \); that needs uniform integrability, which the theorem does not assert (@exr-tr-poisson-variance-limit
supplies it for Poisson counts).

::: {#cor-tr-classical}
[The classical transformations]

::: {.enumerate options="label=(\alph*)"}
1. **Counts.** If \( Y\sim\text{Poisson}(\mu) \), then \( 2(\sqrt{Y}-\sqrt{\mu})\to\Normal(0,1) \) as \( \mu\to\infty \).

2. **Proportions.** If \( S\sim\text{Binomial}(m,p) \) with \( 0<p<1 \) fixed and \( \hat{p}=S/m \), then
   \( 2\sqrt m\,\{\arcsin\sqrt{\hat{p}}-\arcsin\sqrt{p}\}\to\Normal(0,1) \) as \( m\to\infty \).

3. **Constant coefficient of variation.** If \( Y=\mu U \) for a positive random variable \( U \) whose distribution does
   not depend on \( \mu \), then \( \Var(\log Y)=\Var(\log U) \) exactly, for every \( \mu \), whenever it is finite.

4. **Power variance functions.** If \( v(\mu)=\mu^{2a} \) on \( (0,\infty) \), then @eq-tr-stabilizing-integral gives
   \( \mu^{1-a}/(1-a) \) for \( a\ne1 \) and \( \log\mu \) for \( a=1 \). Up to an affine change this is the Box–Cox
   transformation with \( \lambda=1-a \).
:::

:::

::: {.proof}
(a) The characteristic function of \( (Y-\mu)/\sqrt{\mu} \) is
\( \exp\{\mu(e^{it/\sqrt{\mu}}-1-it/\sqrt{\mu})\} \), which tends to \( e^{-t^2/2} \) because
\( e^{iu}-1-iu=-u^2/2+O(u^3) \). By the continuity theorem of Lévy (van der Vaart 1998, chapter 2),
\( (Y-\mu)/\sqrt{\mu}\to\Normal(0,1) \). Now write
\[
2(\sqrt{Y}-\sqrt{\mu})=\frac{Y-\mu}{\sqrt{\mu}}\cdot\frac{2\sqrt{\mu}}{\sqrt{Y}+\sqrt{\mu}} .
\]
Since \( \Var(Y/\mu)=1/\mu \), Chebyshev's inequality gives \( Y/\mu\to1 \) in probability, so the second factor tends
to \( 1 \) in probability, and Slutsky's lemma gives the limit.
(b) By the central limit theorem \( \sqrt m(\hat{p}-p)\to\Normal(0,p(1-p)) \). With \( v(p)=p(1-p) \) on \( (0,1) \),
the derivative of \( \arcsin\sqrt p \) is \( \tfrac12\{p(1-p)\}^{-1/2} \), so @eq-tr-stabilizing-integral gives
\( g(p)=2\arcsin\sqrt p \), and @thm-tr-variance-stabilizing(a) applies.
(c) \( \log Y=\log\mu+\log U \), and adding a constant does not change a variance.
(d) Integrate \( t^{-a} \).
:::

Part (c) is exact: for gamma data with shape \( k \), \( \Var(\log Y) \) is the trigamma function at \( k \) whatever the
mean, \( 0.2213 \) for \( k=5 \) against the first-order value \( 1/k=0.2 \). Part (d) is the rule quoted in
[Section 22.1](01-transforming-the-response.html). A known variance function determines the power; otherwise
the Box–Cox likelihood estimates one, chosen partly to stabilize the variance.

## How well they work

The limits say nothing about moderate means, but the variances can be computed exactly by summing
probability mass functions. [Figure 22.4.1](#fig-tr-stabilizing) and the table show them, scaled so that the target is \( 1 \).

| Poisson mean \( \mu \) | \( 1 \) | \( 2 \) | \( 5 \) | \( 10 \) | \( 20 \) |
|---|---|---|---|---|---|
| \( 4\Var\sqrt{Y} \) | \( 1.609 \) | \( 1.560 \) | \( 1.144 \) | \( 1.045 \) | \( 1.020 \) |
| \( 4\Var\sqrt{Y+3/8} \) | \( 0.717 \) | \( 0.924 \) | \( 1.002 \) | \( 1.001 \) | \( 1.000 \) |

Without transformation the variance of a Poisson count grows in proportion to the mean, a factor of \( 20 \)
between \( \mu=1 \) and \( \mu=20 \). After the square root the scaled variance stays between \( 1 \) and \( 1.7 \) over
the same range, and from \( \mu=5 \) on it is within \( 15 \) percent of the target. Anscombe's (1948) modification
\( \sqrt{Y+3/8} \) does better still for every mean above about \( 0.67 \). The reason is a second-order effect. A formal expansion
gives \( \Var\sqrt{Y+c}=\tfrac14+(3-8c)/(32\mu)+O(\mu^{-2}) \) (@exr-tr-anscombe), and \( c=3/8 \) removes the
\( 1/\mu \) term. For binomial proportions from \( m=20 \) trials, \( 4m\Var(\arcsin\sqrt{\hat{p}}) \) is \( 1.056 \) at
\( p=0.5 \), \( 1.118 \) at \( p=0.25 \), \( 1.546 \) at \( p=0.1 \) and \( 1.632 \) at \( p=0.05 \). Anscombe's version
\( \arcsin\sqrt{(S+3/8)/(m+3/4)} \) gives \( 0.977 \), \( 0.979 \), \( 0.919 \) and \( 0.720 \). Both
transformations fail near the boundary, where most of the probability sits on a few values and no smooth
function can equalize anything.

::: {when-format="html"}
![**Figure 22.4.1.** Exact variances of transformed counts and proportions, scaled so that the target is \( 1 \)
(dashed). (a) Poisson counts: \( \sqrt{y} \) and Anscombe's \( \sqrt{y+3/8} \), against the mean on a log scale.
(b) Binomial proportions with \( m=20 \) trials: the arcsine square root, plain and with Anscombe's constants.](variance_stabilizing.svg){#fig-tr-stabilizing width=100%}
:::

::: {when-format="pdf"}
![Exact variances of transformed counts and proportions, scaled so that the target is \( 1 \)
(dashed). (a) Poisson counts: \( \sqrt{y} \) and Anscombe's \( \sqrt{y+3/8} \), against the mean on a log scale.
(b) Binomial proportions with \( m=20 \) trials: the arcsine square root, plain and with Anscombe's constants.](variance_stabilizing.pdf){width=100%}
:::

```{.python .run #cell-variance-stabilizing-exact}
import numpy as np
from scipy import stats


def poisson_var(h, mu):
    """Var h(Y) for Y ~ Poisson(mu), by summing the probability mass function."""
    k = np.arange(0, int(mu + 40 * np.sqrt(mu) + 60))
    p = stats.poisson.pmf(k, mu)
    m1 = p @ h(k)
    return p @ (h(k) - m1) ** 2


def binomial_var(h, m, prob):
    """Var h(S) for S ~ Binomial(m, prob)."""
    s = np.arange(m + 1)
    p = stats.binom.pmf(s, m, prob)
    m1 = p @ h(s)
    return p @ (h(s) - m1) ** 2


for mu in (1, 2, 5, 10, 20):
    print(f"Poisson mean {mu:2d}:  4 Var sqrt(Y) = {4 * poisson_var(np.sqrt, mu):.3f}"
          f"   4 Var sqrt(Y + 3/8) = {4 * poisson_var(lambda k: np.sqrt(k + 3 / 8), mu):.3f}")
m = 20
for prob in (0.05, 0.1, 0.25, 0.5):
    plain = 4 * m * binomial_var(lambda s: np.arcsin(np.sqrt(s / m)), m, prob)
    anscombe = 4 * m * binomial_var(lambda s: np.arcsin(np.sqrt((s + 3 / 8) / (m + 3 / 4))), m, prob)
    print(f"binomial m = {m}, p = {prob:.2f}:  4m Var = {plain:.3f}   with 3/8: {anscombe:.3f}")
```

A stabilized scale also gives a check: for Poisson counts the residual variance on the square-root scale
should be near \( \tfrac14 \). A much larger value says the counts are *overdispersed*, and then the square root no
longer stabilizes (@exr-tr-negative-binomial).

## Estimating the variance function

With replicates at several settings, plot the log standard deviation of each group against the log of its
mean. If \( s\propto\mu^{a} \) the points lie near a line of slope \( a \), and the power \( 1-a \) stabilizes the variance (@exr-tr-spread-level), an idea that goes back to Bartlett (1947). Without replicates, residuals can be grouped
by fitted value, but misfit of the mean then inflates some spreads, so the plot should come *after* the mean
function is roughly right. Whether the variance depends on the regressors at all is tested formally by @thm-het-breusch-pagan.

## Transform, weight, or model the variance

Transformation is one of three remedies for a variance that depends on the mean.

1. **Transform.** Fit the linear model to \( g(y) \). The variance becomes nearly constant, but the *mean* is now
   modelled on the new scale, and statements about
   \( y \) need the retransformation of [Section 22.5](05-interpreting.html).

2. **Weight.** Keep the scale and fit by weighted least squares with weights \( 1/v(\hat{\mu}_i) \) (@thm-het-wls and @def-proj-gls),
   or keep ordinary least squares with heteroscedasticity-consistent standard errors (@thm-het-sandwich).

3. **Model.** Specify the mean through a link function and the variance through \( v \): a generalized linear model,
   such as Poisson or logistic regression, whose coefficients refer to the mean of \( y \) itself. Chapter 34 develops
   this modern default.

The first two can be combined. The *logit* of a binomial proportion is not variance-stabilizing, but its
variance is about \( 1/\{m\,p(1-p)\} \), so weighted least squares on empirical logits with weights
\( m\hat{p}(1-\hat{p}) \) gives a linear model for the log odds (@exr-tr-empirical-logit), the approach of Grizzle,
Starmer and Koch (1969) for categorical data.

## Exercises

### A. Check your understanding

::: {#exr-tr-cubic-variance}
[A1]

Find the variance-stabilizing transformation for \( v(\mu)=\mu^3 \), the variance function of the inverse Gaussian
distribution. Which Box–Cox power is it?
:::

::: {#exr-tr-exponential-mean}
[A2]

Let \( T_m \) be the mean of \( m \) independent exponential waiting times with mean \( \theta>0 \). Find the
variance-stabilizing transformation. Use @thm-tr-variance-stabilizing(b) to show that the square root, although
it stabilizes Poisson counts, does not stabilize \( T_m \), and find the limiting variance of
\( \sqrt m\,(\sqrt{T_m}-\sqrt\theta) \).
:::

### B. Practice

::: {#exr-tr-negative-binomial}
[B1]

Negative binomial counts have \( v(\mu)=\mu+\mu^2/\kappa \) for a shape parameter \( \kappa>0 \). Show that
@eq-tr-stabilizing-integral gives \( g(\mu)=2\sqrt{\kappa}\,\sinh^{-1}\sqrt{\mu/\kappa} \). Show that \( g \) behaves like \( 2\sqrt\mu \) for
\( \mu\ll\kappa \) and like \( \sqrt\kappa\,\log\mu \) plus a constant for \( \mu\gg\kappa \), and interpret.
:::

::: {.solution}
Substitute \( t=\kappa u^2 \): \( dt=2\kappa u\,du \) and \( t+t^2/\kappa=\kappa u^2(1+u^2) \), so the integrand becomes
\( 2\kappa u\,du/\{\sqrt\kappa\,u\sqrt{1+u^2}\}=2\sqrt\kappa\,du/\sqrt{1+u^2} \), whose integral is
\( 2\sqrt\kappa\sinh^{-1}u \) with \( u=\sqrt{\mu/\kappa} \). For small \( u \), \( \sinh^{-1}u\approx u \), giving \( 2\sqrt\mu \). For large
\( u \), \( \sinh^{-1}u=\log(2u)+o(1) \), giving \( \sqrt\kappa\log\mu \) plus a constant. Small counts behave like Poisson
counts and need a square root; large counts have a nearly constant coefficient of variation and need a logarithm.
:::

::: {#exr-tr-empirical-logit}
[B2]

Let \( S\sim\text{Binomial}(m,p) \) with \( 0<p<1 \). Show by the delta method that
\( \sqrt m\,\{\text{logit}(S/m)-\text{logit}\,p\}\to\Normal(0,1/\{p(1-p)\}) \), where \( \text{logit}\,p=\log\{p/(1-p)\} \).
For independent groups \( i=1,\dots,k \) with \( S_i\sim\text{Binomial}(m_i,p_i) \) and \( \text{logit}\,p_i=\x_{(i)}\T\bbeta \), write down the
weighted least squares estimate of \( \bbeta \) from the empirical logits, and its approximate covariance matrix.
:::

::: {.solution}
The derivative of \( \text{logit}\,p \) is \( 1/\{p(1-p)\} \), and \( \sqrt m(S/m-p)\to\Normal(0,p(1-p)) \), so the limiting variance is
\( p(1-p)/\{p(1-p)\}^2 \). With \( \ell_i=\text{logit}(S_i/m_i) \) and \( \mathbf{W}=\diag\{m_i\hat{p}_i(1-\hat{p}_i)\} \), the estimate is
\( \hbeta=(\X\T\mathbf{W}\X)^{-1}\X\T\mathbf{W}\boldsymbol{\ell} \) (@def-proj-gls), with approximate covariance
\( (\X\T\mathbf{W}\X)^{-1} \), since the weights are the approximate inverse variances and no \( \sigma^2 \) needs estimating.
The empirical logit is infinite when \( S_i=0 \) or \( m_i \), which is why \( \tfrac12 \) is often added to \( S_i \) and to
\( m_i-S_i \).
:::

::: {#exr-tr-spread-level}
[B3]

Replicate observations at four settings have means \( 10 \), \( 40 \), \( 160 \), \( 640 \) and standard deviations
\( 2 \), \( 5.7 \), \( 16 \), \( 45 \). Fit the spread–level line \( \log s=c+a\log\mu \) described above by eye or by least squares, and use
@cor-tr-classical(d) to choose the power.
:::

### C. Going deeper

::: {#exr-tr-anscombe}
[C1]

Let \( Y\sim\text{Poisson}(\mu) \) and \( c\ge0 \). Expand \( \sqrt{Y+c}=\sqrt{\mu}\,(1+u)^{1/2} \) with
\( u=(Y-\mu+c)/\mu \) to third order in \( u \), use the Poisson central moments \( \E(Y-\mu)^2=\mu \), \( \E(Y-\mu)^3=\mu \),
\( \E(Y-\mu)^4=3\mu^2+\mu \), and show formally that \( \Var\sqrt{Y+c}=\tfrac14+(3-8c)/(32\mu)+O(\mu^{-2}) \). (The
expansion is formal: making the remainder rigorous needs a bound on the tail of \( Y \).)
:::

::: {#exr-tr-poisson-variance-limit}
[C2]

@cor-tr-classical(a) is a statement about distributions. Show that the variance converges too:
\( 4\Var\sqrt{Y}\to1 \) as \( \mu\to\infty \) for \( Y\sim\text{Poisson}(\mu) \). Use the bound
\( (\sqrt{Y}-\sqrt{\mu})^2\le(Y-\mu)^2/\mu \), and the fact that a family of nonnegative random variables that converges in
distribution, and whose means converge to the mean of the limit, is uniformly integrable.
:::

::: {.solution}
The bound holds because \( \sqrt{Y}-\sqrt{\mu}=(Y-\mu)/(\sqrt{Y}+\sqrt{\mu}) \) and \( \sqrt{Y}+\sqrt{\mu}\ge\sqrt{\mu} \). The variables
\( V_\mu=(Y-\mu)^2/\mu \) converge in distribution to \( \chi^2(1) \) (by the proof of @cor-tr-classical(a)) and have
mean exactly \( 1 \), the mean of the limit, so they are uniformly integrable, and so is the dominated family
\( W_\mu=4(\sqrt{Y}-\sqrt{\mu})^2 \). Since \( W_\mu\to\chi^2(1) \) in distribution, uniform integrability gives
\( \E W_\mu\to1 \). The same argument applied to \( 2(\sqrt{Y}-\sqrt{\mu}) \), which is dominated by \( 1+W_\mu \), gives
\( \E\{2(\sqrt{Y}-\sqrt{\mu})\}\to0 \). Hence \( 4\Var\sqrt{Y}=\E W_\mu-[\E\{2(\sqrt{Y}-\sqrt{\mu})\}]^2\to1 \).
:::
