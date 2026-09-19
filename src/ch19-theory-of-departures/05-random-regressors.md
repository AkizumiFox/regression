# Random regressors and a curved mean

When the regressors are random, as in almost every observational study, the linear model is read
conditionally on \( \X \). [Section 14.1](../ch14-correlation-lack-of-fit-prediction/01-random-regressors.html)
settled the case in which the conditional model is right. If \( \Y\mid\X\sim\Normal_n(\X\bbeta,\sigma^2\I) \), every \( t \) and \( F \)
procedure keeps its exact level unconditionally, whatever the distribution of the regressors (@thm-cor-conditional), and \( \hbeta \) is unbiased with covariance \( \sigma^2\E\bigl[(\X\T\X)^{-1}\bigr] \) (@cor-lm-random-regressors). Nothing further needs to be proved about that case.

Random regressors create new problems only when the conditional model fails, and it can fail in two
essentially different ways. The regressors can be correlated with the errors, as when they are measured
with error or when a common cause of regressor and response is omitted. Then least squares does not even
estimate a coefficient of interest, a problem taken up in [Chapter 24](../ch24-errors-in-variables/index.html) and
[Chapter 25](../ch25-causal-interpretation/index.html). Or the conditional mean \( \E(Y\mid\x) \) is not linear in \( \x \). This
section treats the second case. Its message is that a curved mean, fitted by a linear model with random
regressors, behaves like heteroscedastic error, and the usual standard errors are wrong even when the true
errors are homoscedastic.

## What least squares estimates, conditionally

::: {#prp-dep-design-target}
[The design-dependent target]

Let \( (\X,\Y) \) be random with \( \rank(\X)=p<n \) almost surely, \( \E(\Y\mid\X)=\bmu(\X) \) for some function \( \bmu \)
with values in \( \Real^n \), not necessarily in \( \C(\X) \), and \( \Cov(\Y\mid\X)=\sigma^2\I_n \). Put
\( \bbeta_{\X}=(\X\T\X)^{-1}\X\T\bmu(\X) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\hbeta\mid\X)=\bbeta_{\X} \) and \( \Cov(\hbeta\mid\X)=\sigma^2(\X\T\X)^{-1} \).

2. If the second moments below are finite,
   \( \Cov(\hbeta)=\sigma^2\E\bigl[(\X\T\X)^{-1}\bigr]+\Cov(\bbeta_{\X}) \).

3. \( \E(s^2\mid\X)=\sigma^2+\norm{(\I-\M)\bmu(\X)}^2/(n-p) \).
:::

:::

::: {.proof}
(a) Given \( \X \), apply @prp-lm-misspecified with \( \boldsymbol{\updelta}=\bmu(\X)-\X\bbeta_{\X} \), which satisfies
\( \X\T\boldsymbol{\updelta}=\bzero \), and @thm-lm-moments(b). (b) is the law of total covariance (@prp-rv-total-covariance)
applied to (a). (c) is @eq-rv-rss-mean given \( \X \).
:::

Given the design, \( \hbeta \) scatters around \( \bbeta_{\X} \) exactly as the usual formula says. But \( \bbeta_{\X} \), the
least squares coefficient of the true mean over the design drawn, changes from design to design. When
\( \bmu(\X)=\X\bbeta \) it is constant (@cor-lm-random-regressors). Otherwise (b) adds a *between-design* component
that no conditional calculation can see, while (c) inflates \( s^2 \) by the misfit, which compensates only partly.

## The large-sample picture

The natural fixed target is the population projection of [Section 6.11](../ch06-projections/11-population.html),
\( \bbeta^*=\bSigma^{-1}\E(\x Y) \), which least squares estimates consistently whatever the form of the mean (@prp-proj-consistency). Its sampling distribution is the following.

::: {#thm-dep-random-x}
[Least squares for the population projection]

Let \( (\x_i,Y_i) \), \( i=1,\dots,n \), be independent copies of \( (\x,Y) \), where \( \x\in\Real^p \) (whose first component may be
the constant \( 1 \)) and \( Y \) have finite fourth moments and \( \bSigma=\E(\x\x\T) \) is positive definite. Let
\( \bbeta^*=\bSigma^{-1}\E(\x Y) \), \( U=Y-\x\T\bbeta^* \) and \( \boldsymbol{\Omega}=\E(U^2\x\x\T) \). Then \( \E(\x U)=\bzero \), and as \( n\to\infty \):

::: {.enumerate options="label=(\alph*)"}
1. \( \sqrt n\,(\hbeta-\bbeta^*)\to\Normal_p\bigl(\bzero,\bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1}\bigr) \) in distribution;

2. \( s^2\to\E(U^2) \) and \( n\,s^2(\X\T\X)^{-1}\to\E(U^2)\,\bSigma^{-1} \) in probability;

3. if \( \E(U^2\mid\x)=\E(U^2) \) almost surely, the two limiting covariances agree. Otherwise the usual interval
   for \( \mathbf{a}\T\bbeta^* \) of nominal level \( 1-\alpha \) has coverage tending to \( 2\Phi\bigl(z_{\alpha/2}\sqrt{\rho_{\mathbf{a}}}\bigr)-1 \), with
   \( \rho_{\mathbf{a}}=\E(U^2)\,\mathbf{a}\T\bSigma^{-1}\mathbf{a}\,/\,\mathbf{a}\T\bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1}\mathbf{a} \).
:::

:::

::: {.proof}
\( \E(\x U)=\E(\x Y)-\bSigma\bbeta^*=\bzero \). Write \( \mathbf{U}=\Y-\X\bbeta^* \), so that
\[
\sqrt n\,(\hbeta-\bbeta^*)=\Bigl(\frac{\X\T\X}{n}\Bigr)^{-1}\frac{\X\T\mathbf{U}}{\sqrt n} .
\]
By the weak law of large numbers \( \X\T\X/n=n^{-1}\sum_i\x_i\x_i\T\to\bSigma \) in probability, so it is invertible with
probability tending to one and its inverse tends to \( \bSigma^{-1} \), inversion being continuous at \( \bSigma \). The vectors
\( \x_iU_i \) are independent with mean \( \bzero \) and covariance \( \boldsymbol{\Omega} \), which is finite because
\( \lvert\E(U^2x_jx_k)\rvert\le\bigl(\E U^4\,\E x_j^2x_k^2\bigr)^{1/2} \) and \( U \) has a finite fourth moment. By the central limit
theorem \( \X\T\mathbf{U}/\sqrt n=\sqrt n\cdot n^{-1}\sum_i\x_iU_i\to\Normal_p(\bzero,\boldsymbol{\Omega}) \). If a random matrix \( \A_n \) tends in probability to
a constant \( \A \) and a random vector \( \mathbf{T}_n \) tends in distribution to \( \mathbf{T} \), then \( \A_n\mathbf{T}_n\to\A\mathbf{T} \): the difference
\( (\A_n-\A)\mathbf{T}_n \) tends to zero in probability because \( \mathbf{T}_n \) is bounded in probability, \( \A\mathbf{T}_n\to\A\mathbf{T} \) by
continuous mapping, and Slutsky's lemma combines the two. This proves (a).

(b) \( (n-p)s^2=\mathbf{U}\T(\I-\M)\mathbf{U}=\mathbf{U}\T\mathbf{U}-Q_n \), where
\( Q_n=(\X\T\mathbf{U}/\sqrt n)\T(\X\T\X/n)^{-1}(\X\T\mathbf{U}/\sqrt n) \) converges in distribution by the argument for (a) and
continuous mapping. Hence \( Q_n/n\to0 \) in probability, while \( \mathbf{U}\T\mathbf{U}/n\to\E(U^2) \) by the law of large numbers. Then
\( n\,s^2(\X\T\X)^{-1}=s^2(\X\T\X/n)^{-1}\to\E(U^2)\bSigma^{-1} \).

(c) If \( \E(U^2\mid\x)=c \), then \( \boldsymbol{\Omega}=\E\bigl[\E(U^2\mid\x)\x\x\T\bigr]=c\,\bSigma \), with \( c=\E(U^2) \), and
\( \bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1}=c\,\bSigma^{-1} \). In general, the \( t \) statistic for \( \mathbf{a}\T\bbeta^* \) is
\( \sqrt n\,\mathbf{a}\T(\hbeta-\bbeta^*) \) divided by \( \bigl(n\,s^2\mathbf{a}\T(\X\T\X)^{-1}\mathbf{a}\bigr)^{1/2} \). By (a), (b) and Slutsky's lemma it tends to
\( \Normal(0,1/\rho_{\mathbf{a}}) \), and the coverage follows as in @thm-dep-nonnormal(c).
:::

The condition in (c) holds when the mean is linear and the errors are homoscedastic. It typically fails when either
assumption fails (@exr-dep-hetero-harmless gives an exception), and a curved mean alone is enough. If \( \E(Y\mid\x)=m(\x) \) and \( \Var(Y\mid\x)=\sigma^2 \), then
\( U=\bigl(m(\x)-\x\T\bbeta^*\bigr)+\varepsilon \) and
\[
\E(U^2\mid\x)=\bigl(m(\x)-\x\T\bbeta^*\bigr)^2+\sigma^2,
\]
which varies with \( \x \) wherever the linear approximation is poor, exactly like a heteroscedastic error. The correct covariance \( \bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1} \) is a population version of the sandwich
@eq-lm-sandwich, and [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) shows how to estimate it.

::: {#exm-dep-curved-mean}
[A line through a parabola]

Let \( X \) be uniform on \( (0,2) \), \( \E(Y\mid X)=X^2 \) and \( \Var(Y\mid X)=0.04 \), with normal errors. A straight line is fitted.
The population line has slope \( \beta^*=\Cov(X,X^2)/\Var(X)=2 \) and intercept \( -2/3 \), and the approximation error is
\( U_0^2-1/3 \) with \( U_0=X-1 \) uniform on \( (-1,1) \). The slope's component of \( \bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1} \) is
\( \E\bigl[(X-1)^2U^2\bigr]/\Var(X)^2 \) (@exr-dep-slope-sandwich), which equals
\( 9\bigl(44/945+0.04/3\bigr)=0.5390 \). The usual formula estimates
\( \E(U^2)/\Var(X)=3\,(4/45+0.04)=0.3867 \) instead. The true standard error is
\( 1.181 \) times what the usual one estimates, and by @thm-dep-random-x(c) the nominal \( 95\% \)
interval for \( \beta^* \) has limiting coverage \( 0.903 \).

Simulation agrees. With \( 40000 \) samples the interval covered \( \beta^*=2 \) in
\( 0.896 \) of samples of size \( 50 \) and \( 0.900 \) of size \( 400 \); with \( 10000 \)
samples of size \( 1600 \), in \( 0.908 \). Over the same samples, \( n \) times the variance of the
slope was \( 0.584 \), \( 0.546 \) and \( 0.530 \), approaching the sandwich
value \( 0.5390 \), not the usual \( 0.3867 \). As a statement about the design-dependent target \( \bbeta_{\X} \), the
same intervals are *conservative*. They covered it in \( 0.999 \) of the samples of size \( 50 \),
because \( s^2 \) estimates \( \E(U^2)=0.1289 \) rather than \( \sigma^2=0.04 \). In a larger simulation of \( 200000 \) samples of size \( 50 \), the variance of the
slope, \( 0.01163 \), splits as in @prp-dep-design-target(b) into
\( 0.00249 \) from the errors and \( 0.00913 \) from the variation of
\( \beta_{\X} \) between designs.
:::

```{.python .run #cell-random-regressors-theory}
import numpy as np
from scipy import stats
from fractions import Fraction as Fr
sigma2 = Fr(1, 25)                                   # sigma = 0.2
EU2, EU4, EU6 = Fr(1, 3), Fr(1, 5), Fr(1, 7)         # moments of U = X - 1 ~ Uniform(-1, 1)
Eu2 = EU4 - Fr(2, 3) * EU2 + Fr(1, 9) + sigma2       # u = U^2 - 1/3 + error
sandwich = (EU6 - Fr(2, 3) * EU4 + Fr(1, 9) * EU2 + sigma2 * EU2) / EU2**2
usual = Eu2 / EU2
print("n Var(slope): sandwich", float(sandwich), " usual", float(usual))
z = stats.norm.ppf(0.975)
print("limiting coverage of the usual interval:",
      2 * stats.norm.cdf(z * np.sqrt(float(usual / sandwich))) - 1)
```

```{.python .run #cell-random-regressors-simulate}
def one_run(n, reps, rng):
    x = rng.uniform(0, 2, size=(reps, n))
    y = x**2 + 0.2 * rng.normal(size=(reps, n))
    xc = x - x.mean(axis=1, keepdims=True)
    Sxx = np.sum(xc**2, axis=1)
    b = np.sum(xc * y, axis=1) / Sxx                          # least squares slope
    b_X = np.sum(xc * x**2, axis=1) / Sxx                     # its conditional mean given X
    a = y.mean(axis=1) - b * x.mean(axis=1)
    s2 = np.sum((y - a[:, None] - b[:, None] * x) ** 2, axis=1) / (n - 2)
    half = stats.t.ppf(0.975, n - 2) * np.sqrt(s2 / Sxx)
    return np.mean(np.abs(b - 2.0) <= half), np.mean(np.abs(b - b_X) <= half), n * b.var()

rng = np.random.default_rng(1905)
for n in [50, 400]:
    c_star, c_X, nvar = one_run(n, 4000, rng)
    print(f"n = {n}: covers beta* {c_star:.3f}, covers beta_X {c_X:.3f}, n Var(slope) {nvar:.3f}")
```

Which target is right depends on the question: \( \bbeta^* \) summarizes the population, \( \bbeta_{\X} \) this design. With a curved mean, neither describes \( \E(Y\mid\x) \) itself.

## Controlled regressors

One situation with random regressors is more benign than it looks. In a controlled experiment the analyst sets a
nominal value \( x_i \), but the value actually delivered is \( x_i+\delta_i \), with an error \( \delta_i \) independent of \( x_i \).
If the response follows \( \beta_0+\beta_1(x_i+\delta_i)+\varepsilon_i \), then regressing on the nominal values gives
\[
Y_i=\beta_0+\beta_1x_i+(\beta_1\delta_i+\varepsilon_i),
\]
a correct linear model with a larger error variance (@exr-dep-berkson). Least squares on the nominal settings is
unbiased, as Berkson (1950) pointed out. Contrast the *classical* measurement error model, in which the true value is
fixed and the *recorded* value carries the error. There the error is correlated with the recorded regressor, and least
squares is biased towards zero (@thm-eiv-attenuation).

Rounding can fall on either side. If the grid value is what was set, as with a dial, the model is Berkson's. If a continuous true value is
*recorded* rounded to a grid of width \( w \), the rounding error is roughly uniform with variance \( w^2/12 \) and roughly
independent of the true value, so the error is classical. In simple regression, classical error of variance \( \tau^2 \)
multiplies the slope by \( \Var(x)/(\Var(x)+\tau^2) \) ([Chapter 24](../ch24-errors-in-variables/index.html)). So measurement error, rounding
included, can be ignored when its variance is a small fraction of \( \Var(x) \).

::: {.idea}
With random regressors, what matters is whether the conditional model holds. If it does, nothing changes; if
the mean is curved, \( s^2(\X\T\X)^{-1} \) is not the covariance.
:::

## Exercises

### A. Check your understanding

::: {#exr-dep-slope-sandwich}
[A1]

With \( \x=(1,X)\T \), show that the slope entries of \( \bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1} \) and \( \E(U^2)\bSigma^{-1} \) are
\( \E\bigl[(X-\mu_X)^2U^2\bigr]/\Var(X)^2 \) and \( \E(U^2)/\Var(X) \).
:::

::: {.solution}
\( \bSigma=\begin{pmatrix}1&\mu_X\\\mu_X&\E X^2\end{pmatrix} \) has determinant \( \Var(X) \), and the second row of \( \bSigma^{-1} \) is
\( (-\mu_X,1)/\Var(X) \). So \( \mathbf{e}_2\T\bSigma^{-1}\x=(X-\mu_X)/\Var(X) \), and
\( \mathbf{e}_2\T\bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1}\mathbf{e}_2=\E\bigl[U^2(X-\mu_X)^2\bigr]/\Var(X)^2 \). Similarly
\( [\bSigma^{-1}]_{22}=1/\Var(X) \).
:::

::: {#exr-dep-linear-target}
[A2]

In @prp-dep-design-target, show that \( \Cov(\bbeta_{\X})=\mathbf{0} \) when \( \bmu(\X)=\X\bbeta \), and give an example of a
nonlinear \( \bmu \) for which \( \Cov(\bbeta_{\X})=\mathbf{0} \) anyway.
:::

### B. Practice

::: {#exr-dep-berkson}
[B1]

In the controlled-regressor model, let \( \delta_i \) and \( \varepsilon_i \) be independent, with variances \( \tau^2 \) and \( \sigma^2 \),
independent of the settings. Show that least squares on the nominal \( x_i \) is unbiased for \( (\beta_0,\beta_1) \), that the
usual standard errors are valid with error variance \( \sigma^2+\beta_1^2\tau^2 \), and that \( s^2 \) estimates this variance.
:::

::: {.solution}
The composite errors \( \eta_i=\beta_1\delta_i+\varepsilon_i \) are independent of the fixed settings, with mean \( 0 \) and
variance \( \sigma^2+\beta_1^2\tau^2 \), and uncorrelated across \( i \). The model \( Y_i=\beta_0+\beta_1x_i+\eta_i \) therefore satisfies
the second-moment assumptions of [Chapter 5](../ch05-model-and-least-squares/index.html), so @thm-lm-moments and @thm-lm-sigma2
apply with \( \sigma^2 \) replaced by \( \sigma^2+\beta_1^2\tau^2 \). If \( \delta_i \) and \( \varepsilon_i \) are normal, so is \( \eta_i \), and the exact
theory of Part III applies too.
:::

::: {#exr-dep-hetero-harmless}
[B2]

In @thm-dep-random-x, show that \( \bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1}=\E(U^2)\bSigma^{-1} \) iff
\( \E\bigl[\bigl(\E(U^2\mid\x)-\E(U^2)\bigr)\x\x\T\bigr]=\mathbf{0} \). Give an example of a nonconstant \( \E(U^2\mid\x) \) for which the usual
standard errors are nevertheless consistent.
:::

### C. Going deeper

::: {#exr-dep-two-components}
[C1]

In @thm-dep-random-x, let \( \E(Y\mid\x)=m(\x) \), \( \Var(Y\mid\x)=\sigma^2 \), and put \( d(\x)=m(\x)-\x\T\bbeta^* \). Show that
\( \sqrt n\,(\bbeta_{\X}-\bbeta^*)\to\Normal_p\bigl(\bzero,\bSigma^{-1}\E(d^2\x\x\T)\bSigma^{-1}\bigr) \), and that
\( \bSigma^{-1}\boldsymbol{\Omega}\bSigma^{-1}=\sigma^2\bSigma^{-1}+\bSigma^{-1}\E(d^2\x\x\T)\bSigma^{-1} \). Interpret the two terms in the light of
@prp-dep-design-target(b).
:::

::: {.solution}
\( \bbeta_{\X}-\bbeta^*=(\X\T\X)^{-1}\X\T\mathbf{d} \) with \( \mathbf{d}=(d(\x_1),\dots,d(\x_n))\T \), and \( \E\bigl(\x\,d(\x)\bigr)=\E(\x U)=\bzero \)
because \( \E(U\mid\x)=d(\x) \). The argument of @thm-dep-random-x(a) with \( U \) replaced by \( d(\x) \) gives the first limit.
Since \( U=d(\x)+\varepsilon \) with \( \E(\varepsilon\mid\x)=0 \), \( \E(U^2\mid\x)=d(\x)^2+\sigma^2 \), and
\( \boldsymbol{\Omega}=\sigma^2\bSigma+\E(d^2\x\x\T) \). The first term, \( \sigma^2\bSigma^{-1} \), is the asymptotic version of
\( n\,\sigma^2\E\bigl[(\X\T\X)^{-1}\bigr] \), the within-design part; the second is the asymptotic version of
\( n\Cov(\bbeta_{\X}) \), the between-design part. (They are limits of those moments only under uniform
integrability, and for some designs the finite-sample moments do not even exist.)
:::
