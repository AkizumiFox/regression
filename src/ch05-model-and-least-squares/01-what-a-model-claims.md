# What a regression model claims

A data set for regression consists of \( n \) cases. Each case records a response \( y_i \)
and a few other measurements, the regressors, on the same unit, such as a patient or a month. The question is how the response depends on the
regressors. Two cases with identical regressors almost always have different
responses, so a regression model cannot be an equation linking
\( y \) to the regressors. It is a statement about the *distribution* of the
response given the regressors, and usually only about one feature of that distribution,
its mean.

This section separates what is automatically true of any response and regressors
from what a model actually assumes. The linear model's assumptions come in layers, and
each layer buys a different part of the theory.

## The regression function

Let \( Y \) be a response with finite variance and \( \bm X \) a random vector of
regressors, defined on the same probability space. We need only elementary
properties of conditional expectation: the tower property \( \E[\E(Y\mid\bm X)]=\E(Y) \), and
\( \E[g(\bm X)Y\mid\bm X]=g(\bm X)\E(Y\mid\bm X) \) for functions \( g \) of the regressors.

::: {#def-lm-regression-function}
[Regression function]

The **regression function** of \( Y \) on \( \bm X \) is the conditional mean
\[
m(\bm x)=\E(Y\mid\bm X=\bm x).
\]
The **error** is \( \varepsilon=Y-m(\bm X) \), so that \( Y=m(\bm X)+\varepsilon \).
:::

The decomposition \( Y=m(\bm X)+\varepsilon \) looks like a model, but it is not one.
It holds for every pair \( (\bm X,Y) \) with \( \E Y^2<\infty \), and the error it defines
has two properties that are often listed as "assumptions" even though they are
automatic.

::: {#prp-lm-error-decomposition}
[The error of the regression function]

Let \( \E Y^2<\infty \), and let \( m \) and \( \varepsilon \) be as in @def-lm-regression-function.

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\varepsilon\mid\bm X)=0 \), and hence \( \E(\varepsilon)=0 \).

2. \( \E[\varepsilon\,g(\bm X)]=0 \) for every function \( g \) with \( \E g(\bm X)^2<\infty \). In
   particular \( \varepsilon \) is uncorrelated with every such function of the regressors.

3. For every such \( g \),
   \[
   \E\bigl[(Y-g(\bm X))^2\bigr]=\E(\varepsilon^2)+\E\bigl[(m(\bm X)-g(\bm X))^2\bigr].
   \]
   So \( m(\bm X) \) is the best predictor of \( Y \) from \( \bm X \) in mean squared
   error, and it is the only one up to events of probability zero.
:::

:::

::: {.proof}
*(a)* \( \E(\varepsilon\mid\bm X)=\E(Y\mid\bm X)-m(\bm X)=0 \), and the tower property gives
\( \E(\varepsilon)=\E[\E(\varepsilon\mid\bm X)]=0 \). *(b)* The product \( \varepsilon g(\bm X) \) is integrable
by the Cauchy–Schwarz inequality, because \( \E\varepsilon^2\le 2\E Y^2+2\E m(\bm X)^2 \) and
\( \E m(\bm X)^2\le\E Y^2 \) by Jensen's inequality for conditional expectations. Then
\( \E[\varepsilon g(\bm X)]=\E\bigl[g(\bm X)\E(\varepsilon\mid\bm X)\bigr]=0 \). *(c)* Write
\( Y-g(\bm X)=\varepsilon+(m(\bm X)-g(\bm X)) \) and expand the square. The cross term is
\( 2\E[\varepsilon\,(m(\bm X)-g(\bm X))] \), which is zero by (b) applied to the function
\( m-g \). The second term on the right is zero only if \( g(\bm X)=m(\bm X) \) with
probability one.
:::

Part (c) is the reason the conditional mean is the natural target. It is what an ideal
forecaster would report under squared-error loss. It is also the reason regression
is sometimes misread. The function \( m \) summarizes the conditional distribution by
its centre. Two populations can have the same regression function and very different
spreads, shapes and tails around it.

::: {.warning}
Part (b) says the error is uncorrelated with the regressors *by construction*.
The lack of correlation between errors and regressors therefore cannot be tested from the data or
defended as a scientific assumption about the regression function. The real
assumption hides elsewhere: in the claim that \( m \) has a particular form, or
that the quantity of interest is \( m \) at all. When someone says "the error is
correlated with the regressor", the error they have in mind is a different
object, usually the deviation from a *causal* or *structural* relation rather
than from the conditional mean. Chapter 24 and Chapter 25 take this up.
:::

## What the linear model adds

A linear regression model makes claims about \( m \) and about the errors that
\( \E Y^2<\infty \) does not guarantee. Write \( \bm x_{(i)} \) for the regressors of case
\( i \), already transformed into the columns the model uses (a constant \( 1 \), the
regressors themselves, their squares, indicators of categories, and so on;
[Section 5.3](03-multiple-regression-model.html) gives many examples). The
assumptions come in four layers.

::: {.enumerate options="label=(L\arabic*)"}
1. **Linearity.** \( \E(Y_i)=\bm x_{(i)}\T\bbeta \) for an unknown \( \bbeta\in\Real^p \). The
   regression function belongs to a fixed \( p \)-dimensional family of functions.

2. **Constant variance.** \( \Var(Y_i)=\sigma^2 \) for every \( i \), with \( \sigma^2 \) unknown.

3. **Uncorrelated errors.** \( \Cov(Y_i,Y_j)=0 \) for \( i\ne j \).

4. **Normality.** The errors are jointly normal.
:::

The first three layers are the **second-moment assumptions**. They mention only
means, variances and covariances, and they are all that [Section 5.5](05-properties.html)
and [Section 5.6](06-error-variance.html) use. Adding (L4) gives the normal linear
model, and with it exact sampling distributions for estimates and test statistics,
which [Chapter 7](../ch07-optimality/index.html) and [Chapter 11](../ch11-general-linear-hypothesis/index.html) develop. Each layer is used by a specific part of the
theory:

| Assumptions | What they deliver | Where |
|---|---|---|
| (L1) | the least squares estimator is unbiased | [Section 5.5](05-properties.html) |
| (L1)–(L3) | its covariance matrix; an unbiased estimate of \( \sigma^2 \); Gauss–Markov optimality | [Section 5.5](05-properties.html), [Section 5.6](06-error-variance.html), [Chapter 7](../ch07-optimality/index.html) |
| (L1)–(L4) | exact \( t \) and \( F \) distributions; maximum likelihood; exact intervals | [Chapter 7](../ch07-optimality/index.html), [Chapter 11](../ch11-general-linear-hypothesis/index.html), [Chapter 12](../ch12-intervals-and-bands/index.html) |

Least squares itself needs none of them. The estimator of
[Section 5.4](04-least-squares.html) is defined for any data at all. The assumptions
only tell us how good it is.

::: {#exm-lm-two-processes}
[Same regression function, different models]

[Figure 5.1.1](#fig-lm-conditional) shows two simulated data sets with the same
regression function \( m(x)=1+0.5x \). In panel (a) the errors are normal with constant
standard deviation \( 0.8 \). In panel (b) they are skewed to the right, and their
standard deviation grows linearly in \( x \), from \( 0.15 \) at \( x=0 \) to \( 2.65 \) at
\( x=10 \). The curves show the conditional density of \( Y \) at three values of
\( x \). Both processes satisfy (L1), so least squares estimates the line without
bias in both. Only process (a) satisfies (L2), so only there do the variance formulas of
[Section 5.5](05-properties.html) apply as stated. In process (b) most points lie below
the line, because the median of a right-skewed distribution is below its mean. The line
is still the conditional mean.
:::

::: {when-format="html"}
![**Figure 5.1.1.** Two data-generating processes with the same regression function
\( 1+0.5x \) (brown line). The curves are conditional densities of \( Y \) given
\( x \), drawn sideways at \( x=2,5,8 \). (a) Normal errors with constant variance: all four
layers of assumptions hold. (b) Right-skewed errors whose spread grows with
\( x \): only linearity holds.](conditional_mean.svg){#fig-lm-conditional width=100%}
:::

::: {when-format="pdf"}
![Two data-generating processes with the same regression function
\( 1+0.5x \) (brown line). The curves are conditional densities of \( Y \) given
\( x \), drawn sideways at \( x=2,5,8 \). (a) Normal errors with constant variance: all four
layers of assumptions hold. (b) Right-skewed errors whose spread grows with
\( x \): only linearity holds.](conditional_mean.pdf){width=100%}
:::

The assumptions concern the errors \( \varepsilon_i \), which are never observed. What
we see are residuals, the differences between the responses and a fitted function.
Whether the assumptions are credible is judged from residuals, from the design of the
study, and from knowledge of the subject. Part V of the book is devoted to the first of
these and to what goes wrong when the assumptions fail.

## Fixed and random regressors

The layers (L1)–(L4) were written with \( \E(Y_i) \) and \( \Var(Y_i) \), as if the
regressors were constants. In a designed experiment they are. The experimenter
chooses the temperatures, doses or fertilizer levels, and only the responses are random.
In an observational study the regressors are measured on units drawn from a
population, so the pairs \( (\bm x_{(i)},Y_i) \) are random together.

This book's convention, like that of most linear-model theory, is to treat the
regressors as fixed, and to read (L1)–(L4) in the random case as statements
*conditional on the regressors*: \( \E(Y_i\mid\X)=\bm x_{(i)}\T\bbeta \),
\( \Var(Y_i\mid\X)=\sigma^2 \), and so on, where \( \X \) collects all the regressors. There
are two reasons for this choice.

First, conditional statements imply unconditional ones where it matters. If the
estimator is unbiased given every possible \( \X \), it is unbiased when \( \X \) is random,
by the tower property. [Section 5.5](05-properties.html) makes this precise
(@cor-lm-random-regressors). Second, conditioning describes the precision of the data
actually in hand. A sample whose regressors happen to be spread out estimates a slope
better than one whose regressors happen to be bunched together. The conditional
variance says so, and the unconditional variance averages the two cases.

::: {#exm-lm-fixed-random}
[One slope, two sampling schemes]

Take \( Y=1+0.5x+\varepsilon \) with standard normal errors and \( n=15 \) cases. In the
*fixed design*, the regressor values are drawn once from the uniform distribution on
\( [0,10] \) and then kept. They have \( S_{xx}=\sum_i(x_i-\bar{x})^2=103.73 \). In the
*random design*, a new set of \( 15 \) regressor values is drawn for every sample.
Over 50,000 simulated samples, the least squares slope
(formula in [Section 5.2](02-simple-regression.html)) averages \( 0.5007 \)
under the fixed design and \( 0.5006 \) under the random one. Both are unbiased.
Their variances differ in meaning. Under the fixed design the simulated variance is
\( 0.00960 \), matching \( \sigma^2/S_{xx}=0.00964 \) for this particular
design. Under the random design it is \( 0.00933 \), matching the *average*
\( \E(\sigma^2/S_{xx})=0.00923 \) over designs, estimated from the same
simulation. Replacing \( S_{xx} \) by its expectation \( 116.67 \) would give
\( 0.00857 \), which is too small. By Jensen's inequality \( \E(1/S_{xx})\ge1/\E(S_{xx}) \),
so random bunching of the regressors costs precision on average.
:::

```{.python .run #cell-conditional-mean-fixed-random}
import numpy as np

rng = np.random.default_rng(2005)
n, beta0, beta1, sigma = 15, 1.0, 0.5, 1.0
reps = 50_000

def slope(x, y):
    """Least squares slope of y on x (with an intercept)."""
    xc = x - x.mean()
    return xc @ (y - y.mean()) / (xc @ xc)

x_fixed = rng.uniform(0, 10, size=n)            # drawn once, then held fixed
Sxx_fixed = np.sum((x_fixed - x_fixed.mean()) ** 2)
b_fixed = np.empty(reps)
b_random = np.empty(reps)
Sxx_random = np.empty(reps)
for r in range(reps):
    e = sigma * rng.standard_normal(n)
    b_fixed[r] = slope(x_fixed, beta0 + beta1 * x_fixed + e)
    x = rng.uniform(0, 10, size=n)              # a new design for every sample
    Sxx_random[r] = np.sum((x - x.mean()) ** 2)
    b_random[r] = slope(x, beta0 + beta1 * x + e)

print(f"fixed design:  mean {b_fixed.mean():.4f}  var {b_fixed.var():.5f}"
      f"  sigma^2/Sxx = {sigma**2 / Sxx_fixed:.5f}")
print(f"random design: mean {b_random.mean():.4f}  var {b_random.var():.5f}"
      f"  E(sigma^2/Sxx) = {np.mean(sigma**2 / Sxx_random):.5f}")
```

Conditioning on the regressors is not always harmless. It fails when the regressors
are measured with error, because then the regression on the *measured* values is not the
one of scientific interest (Chapter 24). It fails when a regressor is a lagged value of
the response, as in autoregressive time series, because then \( \X \) is not
determined before the errors are drawn. It also changes the question when the target is a
population quantity, such as the best linear predictor of one variable from others in a
population where the regression function is not linear
([Section 6.11](../ch06-projections/11-population.html) and [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html)).

## What regression is used for

The same fitted equation can serve different purposes, and the purpose decides which
assumptions matter. For **prediction** of new cases only the accuracy of the fitted mean
matters, and a regressor can be useful without any causal role. For **description** the
fitted mean summarizes how the response varies across the population sampled. For
**estimation of a parameter**, such as an elasticity or a dose-response slope, one
coefficient and its uncertainty are the point. For **control**, regressor values are
chosen to achieve a desired response, which requires the fitted relation to hold when the
regressors are *set* rather than observed: a causal claim.

::: {.warning}
The regression function describes how the mean response differs between cases that
*are observed* to have different regressor values. It does not by itself describe what
would happen if a regressor were *changed*. Units with high values of one regressor
usually differ in other ways too. Randomized assignment of the regressor, as in a
designed experiment, is the standard way to make the two questions coincide. Without
it, the step from association to effect needs assumptions that the data cannot check.
[Section 5.7](07-interpreting-coefficients.html) returns to this point, and
Chapter 25 treats it in full.
:::

## Exercises

### A. Check your understanding

::: {#exr-lm-hetero-tautology}
[A1]

Let \( X \) take the values \( -1,0,1 \) with probability \( 1/3 \) each, and let
\( Y=X^2+XZ \), where \( Z\sim\Normal(0,1) \) is independent of \( X \). Find the regression
function \( m \) and the error \( \varepsilon \). Show that \( \E(\varepsilon\mid X)=0 \) and
\( \Cov(\varepsilon,X)=0 \), but that \( \varepsilon \) and \( X \) are not independent. Which of the
layers (L1)–(L4) fail for a sample from this distribution, with the columns
\( 1,x,x^2 \)?
:::

::: {.solution}
Since \( \E(XZ\mid X)=X\,\E(Z)=0 \), \( m(x)=x^2 \) and \( \varepsilon=XZ \). Then
\( \E(\varepsilon\mid X)=0 \), and \( \Cov(\varepsilon,X)=\E(X^2Z)=\E(X^2)\E(Z)=0 \). But
\( \Var(\varepsilon\mid X)=X^2 \), which is \( 0 \) when \( X=0 \) and \( 1 \) otherwise. The
conditional distribution of \( \varepsilon \) depends on \( X \), so they are not
independent. With columns \( 1,x,x^2 \), (L1) holds with \( \bbeta=(0,0,1)\T \) and (L3)
holds for independent cases. (L2) fails, and so does (L4) as stated, because the
error of a case with \( x=0 \) is identically zero and the joint distribution is
not normal with covariance \( \sigma^2\I \).
:::

::: {#exr-lm-bivariate-normal}
[A2]

Let \( (X,Y) \) be bivariate normal with means \( \mu_X,\mu_Y \), standard deviations
\( \sigma_X,\sigma_Y>0 \) and correlation \( \rho \). Use @thm-mvn-conditional to find
\( m(x) \) and \( \Var(Y\mid X=x) \). Conclude that a random sample from this
distribution satisfies (L1)–(L4) conditionally on the \( x_i \), with columns \( 1,x \).
What are \( \bbeta \) and \( \sigma^2 \)?
:::

### B. Practice

::: {#exr-lm-variance-explained}
[B1]

Under the conditions of @prp-lm-error-decomposition, show that
\( \Var(Y)=\Var\bigl(m(\bm X)\bigr)+\E\bigl[\Var(Y\mid\bm X)\bigr] \) and that
\( \E(\varepsilon^2)=\E[\Var(Y\mid\bm X)] \). The ratio
\( \eta^2=\Var(m(\bm X))/\Var(Y) \) is the proportion of the variance of \( Y \) that the
regression function explains. Compute it for the distribution of @exr-lm-hetero-tautology.
:::

::: {.solution}
The first identity is the scalar case of @prp-rv-total-covariance. For the second,
\( \E(\varepsilon^2\mid\bm X)=\E[(Y-m(\bm X))^2\mid\bm X]=\Var(Y\mid\bm X) \); take expectations.
In @exr-lm-hetero-tautology, \( m(X)=X^2 \) takes the value \( 1 \) with probability
\( 2/3 \) and \( 0 \) with probability \( 1/3 \), so \( \Var(m(X))=2/9 \). Also
\( \E[\Var(Y\mid X)]=\E(X^2)=2/3 \). Hence \( \Var(Y)=8/9 \) and \( \eta^2=1/4 \).
:::

::: {#exr-lm-jensen-design}
[B2]

In @exm-lm-fixed-random the regressors are independent and uniform on \( [0,10] \).
Show that \( \E(S_{xx})=(n-1)\cdot100/12 \). Use Jensen's inequality to show that
\( \E(\sigma^2/S_{xx})\ge\sigma^2/\E(S_{xx}) \), with strict inequality here. Why is the
right side of no use for describing the precision of the slope from a given sample?
:::

### C. Going deeper

::: {#exr-lm-blp-not-regression}
[C1]

Let \( X\sim\Normal(0,1) \) and \( Y=X^2 \). Find the best linear predictor
\( \beta_0+\beta_1X \) of \( Y \) (@thm-rv-blp) and its error \( u \). Show that \( \E(u)=0 \)
and \( \Cov(u,X)=0 \), but \( \E(u\mid X)\ne0 \). Explain why fitting a straight line to a
large sample from this distribution estimates a well-defined population quantity even
though (L1) fails for the columns \( 1,x \), and why the second-moment theory of this
chapter does not describe the sampling variability of the fitted slope.
:::
