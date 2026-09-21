# Bayesian thread II: credible intervals and predictive distributions

[Section 7.6](../ch07-optimality/06-bayes-conjugate.html) derived the posterior distribution of \( (\bbeta,\sigma^2) \) under the
conjugate normal-inverse-gamma prior, and showed that the flat prior \( \pi(\bbeta,\sigma^2)\propto1/\sigma^2 \) reproduces least
squares. Here the posterior is summarized by credible regions and used to predict. In the conjugate model these have the shapes of
the frequentist sets, and under the flat prior they coincide with them. We end by asking how often a
credible interval covers in repeated sampling.

## Credible sets

::: {#def-ci-credible}
[Credible set]

Let \( \theta \) be a function of the parameters and let \( \pi(\cdot\mid\y) \) be its posterior distribution. A set \( C(\y) \) with
\( \Pr\{\theta\in C(\y)\mid\y\}\ge1-\alpha \) is a \( 1-\alpha \) **credible set** for \( \theta \). If the posterior has a density
\( \pi(\theta\mid\y) \), a set of the form \( \{\theta:\pi(\theta\mid\y)\ge k\} \) with posterior probability \( 1-\alpha \) is a
**highest posterior density** (HPD) set.
:::

Here the statement concerns the realized set. An HPD set has the smallest volume among sets with its posterior probability. We
need three facts about the multivariate \( t \) distribution of @def-opt-mvt.

::: {#lem-ci-mvt}
[Multivariate \( t \)]

Let \( \bu\sim t_\nu(\mathbf{m},\bS) \) in \( \Real^p \), with \( \bS \) positive definite.

::: {.enumerate options="label=(\alph*)"}
1. If \( \A \) is \( q\times p \) of rank \( q \), then \( \A\bu\sim t_\nu(\A\mathbf{m},\A\bS\A\T) \).

2. \( (\bu-\mathbf{m})\T\bS^{-1}(\bu-\mathbf{m})/p\sim F(p,\nu) \).

3. \( \bu \) has density proportional to \( \bigl\{1+(\bu-\mathbf{m})\T\bS^{-1}(\bu-\mathbf{m})/\nu\bigr\}^{-(\nu+p)/2} \).
:::
:::

::: {.proof}
By @def-opt-mvt, \( \bu=\mathbf{m}+\Z/\sqrt{W/\nu} \) with \( \Z\sim\Normal_p(\bzero,\bS) \) and \( W\sim\chi^2(\nu) \) independent.

(a) \( \A\bu=\A\mathbf{m}+\A\Z/\sqrt{W/\nu} \), and \( \A\Z\sim\Normal_q(\bzero,\A\bS\A\T) \) (@thm-mvn-linear) is independent of \( W \), with
\( \A\bS\A\T \) positive definite.

(b) \( (\bu-\mathbf{m})\T\bS^{-1}(\bu-\mathbf{m})=\Z\T\bS^{-1}\Z/(W/\nu) \), and \( \Z\T\bS^{-1}\Z\sim\chi^2(p) \) independently of \( W \) (@cor-qf-mahalanobis(a)).
Divide by \( p \) and use @def-qf-noncentral-f.

(c) Given \( W=w \), \( \bu\sim\Normal_p(\mathbf{m},(\nu/w)\bS) \). Write \( D=(\bu-\mathbf{m})\T\bS^{-1}(\bu-\mathbf{m}) \). Integrating the normal density (@thm-mvn-density)
against the \( \chi^2(\nu) \) density of \( W \), and dropping factors free of \( \bu \),
\[
\begin{aligned}
\int_0^\infty w^{p/2}e^{-wD/(2\nu)}\,w^{\nu/2-1}e^{-w/2}\,dw
&=\int_0^\infty w^{(\nu+p)/2-1}e^{-w(1+D/\nu)/2}\,dw\\
&\propto(1+D/\nu)^{-(\nu+p)/2},
\end{aligned}
\]
by the substitution \( w'=w(1+D/\nu) \).
:::

## Credible sets and predictive distributions in the conjugate model

::: {#thm-ci-bayes-credible}
[Credible sets and posterior predictive distribution]

Let the prior be \( \text{NIG}(\mathbf{m}_0,\V_0,a_0,b_0) \), so that by @thm-opt-bayes-conjugate the posterior is \( \text{NIG}(\mathbf{m}_n,\V_n,a_n,b_n) \).
Put \( \nu=2a_n \) and \( \tau^2=b_n/a_n \).

::: {.enumerate options="label=(\alph*)"}
1. For \( \blambda\ne\bzero \), the interval
   \[
\blambda\T\mathbf{m}_n\ \pm\ t_{\nu,\alpha/2}\,\tau\sqrt{\blambda\T\V_n\blambda}
\]
   is the \( 1-\alpha \) HPD credible interval for \( \blambda\T\bbeta \); it is also equal-tailed.

2. For a \( p\times q \) matrix \( \bLambda \) of rank \( q \), the set
   \[
\begin{aligned}
\bigl\{\boldsymbol{\upphi}_0:\ &(\boldsymbol{\upphi}_0-\bLambda\T\mathbf{m}_n)\T(\bLambda\T\V_n\bLambda)^{-1}(\boldsymbol{\upphi}_0-\bLambda\T\mathbf{m}_n)\\
&\le q\,\tau^2F_\alpha(q,\nu)\bigr\}
\end{aligned}
\]
   is the \( 1-\alpha \) HPD credible region for \( \bLambda\T\bbeta \).

3. Let \( \mathbf{Y}_0=\X_0\bbeta+\be_0 \) be \( k \) future observations, with \( \be_0\sim\Normal_k(\bzero,\sigma^2\I_k) \) independent of \( \Y \) and of
   \( \bbeta \) given \( \sigma^2 \). The posterior predictive distribution of \( \mathbf{Y}_0 \) is
   \[
\mathbf{Y}_0\mid\y\ \sim\ t_\nu\bigl(\X_0\mathbf{m}_n,\ \tau^2(\I_k+\X_0\V_n\X_0\T)\bigr).
\]
   In particular, for one future observation at \( \x_0 \), the interval \( \x_0\T\mathbf{m}_n\pm t_{\nu,\alpha/2}\,\tau\sqrt{1+\x_0\T\V_n\x_0} \) has
   posterior predictive probability \( 1-\alpha \).

4. \( \bigl[2b_n/\chi^2_{\alpha/2}(\nu),\ 2b_n/\chi^2_{1-\alpha/2}(\nu)\bigr] \) is an equal-tailed \( 1-\alpha \) credible interval for \( \sigma^2 \).
:::
:::

::: {.proof}
By @thm-opt-bayes-conjugate and @lem-opt-nig-marginal, \( \bbeta\mid\y\sim t_\nu(\mathbf{m}_n,\tau^2\V_n) \) and \( 2b_n/\sigma^2\mid\y\sim\chi^2(\nu) \).

(a) By @lem-ci-mvt(a) with \( \A=\blambda\T \), \( \blambda\T\bbeta=\blambda\T\mathbf{m}_n+\tau\sqrt{\blambda\T\V_n\blambda}\,T \) with \( T\sim t(\nu) \). The
\( t(\nu) \) density is symmetric and decreasing in \( \lvert x\rvert \), so the symmetric interval is both equal-tailed and HPD.

(b) By @lem-ci-mvt(a), \( \bLambda\T\bbeta\mid\y\sim t_\nu(\bLambda\T\mathbf{m}_n,\tau^2\bLambda\T\V_n\bLambda) \). By (b) of the lemma the quadratic form divided
by \( q\tau^2 \) is \( F(q,\nu) \), so the set has posterior probability \( 1-\alpha \). By (c) the density is a decreasing function of the quadratic
form, so the set has the form \( \{\text{density}\ge k\} \).

(c) Given \( \sigma^2 \) and \( \y \), \( \bbeta\sim\Normal_p(\mathbf{m}_n,\sigma^2\V_n) \) and \( \be_0\sim\Normal_k(\bzero,\sigma^2\I_k) \) independently, so
\( \mathbf{Y}_0\sim\Normal_k\bigl(\X_0\mathbf{m}_n,\sigma^2(\X_0\V_n\X_0\T+\I_k)\bigr) \) (@thm-mvn-linear). Together with \( \sigma^2\mid\y\sim\text{IG}(a_n,b_n) \), this says
that \( (\mathbf{Y}_0,\sigma^2)\mid\y \) is \( \text{NIG}\bigl(\X_0\mathbf{m}_n,\I_k+\X_0\V_n\X_0\T,a_n,b_n\bigr) \) (@def-opt-nig), the matrix being positive definite.
@lem-opt-nig-marginal gives the \( t_\nu \) law, and (a) of the present theorem, applied to it with \( k=1 \), gives the interval.

(d) \( \Pr\{\chi^2_{1-\alpha/2}(\nu)\le2b_n/\sigma^2\le\chi^2_{\alpha/2}(\nu)\mid\y\}=1-\alpha \); invert the inequalities.
:::

The formulas mirror the frequentist ones: \( \hbeta \) becomes \( \mathbf{m}_n \), \( \G \) becomes \( \V_n=(\V_0^{-1}+\X\T\X)^{-1} \), \( s^2 \)
becomes \( \tau^2 \), and \( n-r \) becomes \( \nu=2a_0+n \), which counts every observation plus \( 2a_0 \) "prior observations" about
\( \sigma^2 \). The \( 1 \) in \( 1+\x_0\T\V_n\x_0 \) is again the new error.

::: {#cor-ci-flat-agreement}
[The flat prior reproduces the frequentist sets]

Let \( \rank(\X)=p<n \), \( \text{SSE}>0 \), and let the prior be \( \pi(\bbeta,\sigma^2)\propto1/\sigma^2 \). Then the credible interval (a), the region (b),
the predictive interval (c) and the interval for \( \sigma^2 \) (d) of @thm-ci-bayes-credible are, respectively, the \( t \) interval of
@thm-ci-estimable-interval, the confidence ellipsoid of @thm-ci-ellipsoid, the prediction interval of @thm-ci-prediction-interval and the interval
of @prp-ci-sigma.
:::

::: {.proof}
By @cor-opt-flat-prior the posterior is \( \text{NIG}\bigl(\hbeta,(\X\T\X)^{-1},(n-p)/2,\text{SSE}/2\bigr) \). The proof of
@thm-ci-bayes-credible used only that the posterior is normal-inverse-gamma, so its conclusions hold with \( \mathbf{m}_n=\hbeta \),
\( \V_n=(\X\T\X)^{-1}=\G \), \( \nu=n-p=n-r \), \( \tau^2=\text{SSE}/(n-p)=s^2 \) and \( 2b_n=\text{SSE} \). Substituting gives the four frequentist sets.
:::

The same set then carries two readings: "this procedure covers the fixed truth in \( 95\% \) of samples" and "given these data, the
parameter lies here with probability \( 0.95 \)". The agreement is special; it fails for proper priors and for most other models.

## A worked example

::: {#exm-ci-stackloss-predictive}
[Predicting a day of stack loss]

Take the stack loss regression with the engineer's prior of @exm-opt-stackloss-bayes, and consider a new day with air flow
\( 65 \), water temperature \( 23 \) and acid concentration \( 87 \). The posterior has \( a_n=13.5 \), so \( \nu=27 \) and
\( t_{27,0.025}=2.052 \). By @thm-ci-bayes-credible, the \( 95\% \) credible interval for the mean stack loss at these conditions is
\( 23.12\pm1.51 \), and the \( 95\% \) posterior predictive interval for the day's stack loss is \( 23.12\pm6.11 \).
The least squares answers are \( 23.15\pm1.77 \) for the mean and \( 23.15\pm7.07 \) for the new day, with
\( t_{17,0.025}=2.110 \). The Bayesian intervals are shorter because the prior informs the slopes, shrinking
\( \x_0\T\V_n\x_0 \), and \( \sigma^2 \), raising the degrees of freedom from \( 17 \) to \( 27 \).

The \( 2.5\% \) and \( 97.5\% \) points of \( 400000 \) posterior predictive draws are \( 17.00 \)
and \( 29.21 \), matching the \( t \) interval within Monte Carlo error. The simulated posterior probability of the credible ellipse
of (b) for the air-flow and water-temperature coefficients is \( 0.9503 \).
:::

```{.python .run #cell-bayes-intervals-posterior}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape
XtX, Xty = X.T @ X, X.T @ y


def nig_posterior(m0, P0, a0, b0, y=y):
    """Conjugate update; P0 is the prior precision V0^{-1} (per unit sigma^2)."""
    Vn = np.linalg.inv(P0 + XtX)
    mn = Vn @ (P0 @ m0 + X.T @ y)
    an = a0 + n / 2
    bn = b0 + (y @ y + m0 @ P0 @ m0 - mn @ (P0 + XtX) @ mn) / 2
    return mn, Vn, an, bn


m0 = np.array([0.0, 1.0, 1.0, 0.0])                     # the prior of Chapter 7
P0 = np.linalg.inv(np.diag([100.0**2, 0.1**2, 0.2**2, 0.1**2]))
a0, b0 = 3.0, 18.0
mn, Vn, an, bn = nig_posterior(m0, P0, a0, b0)
```

```{.python .run #cell-bayes-intervals-predictive}
x0 = np.array([1.0, 65.0, 23.0, 87.0])                  # a new operating condition
qb = stats.t.ppf(0.975, 2 * an)
centre = x0 @ mn
half_mean = qb * np.sqrt(bn / an * (x0 @ Vn @ x0))     # credible interval for x0^T beta
half_new = qb * np.sqrt(bn / an * (1 + x0 @ Vn @ x0))  # posterior predictive interval for Y0

beta_hat = np.linalg.solve(XtX, Xty)                   # frequentist intervals, for comparison
s2 = np.sum((y - X @ beta_hat) ** 2) / (n - p)
qf = stats.t.ppf(0.975, n - p)
h0 = x0 @ np.linalg.solve(XtX, x0)
print(f"Bayes: mean response {centre:.2f} +- {half_mean:.2f}; new day {centre:.2f} +- {half_new:.2f}")
print(f"least squares: mean response {x0 @ beta_hat:.2f} +- {qf * np.sqrt(s2 * h0):.2f}; "
      f"new day +- {qf * np.sqrt(s2 * (1 + h0)):.2f}")
```

## How often does a credible interval cover?

One coverage statement holds for every proper prior.

::: {#prp-ci-average-coverage}
[Coverage averaged over the prior]

Let the prior \( \pi \) be proper, and let \( C(\y) \) be a \( 1-\alpha \) credible set for \( \theta \) for every \( \y \). If \( (\bbeta,\sigma^2) \) is drawn from
\( \pi \) and then \( \Y \) from the model, \( \Pr\{\theta\in C(\Y)\}\ge1-\alpha \), with equality if every \( C(\y) \) has posterior probability exactly
\( 1-\alpha \). Equivalently, the frequentist coverage \( \Pr_{\bbeta,\sigma^2}\{\theta\in C(\Y)\} \), averaged over \( \pi \), is at least \( 1-\alpha \).
:::

::: {.proof}
Under the joint distribution of parameters and data, the conditional distribution of the parameters given \( \Y=\y \) is the posterior.
By the tower property, \( \Pr\{\theta\in C(\Y)\}=\E\bigl[\Pr\{\theta\in C(\Y)\mid\Y\}\bigr]\ge1-\alpha \). Conditioning instead on the parameters
gives the second form.
:::

The average hides variation: coverage is high where the prior is right and low where it is wrong. Consider the credible interval for the air-flow coefficient in @exm-ci-stackloss-predictive, and simulate its frequentist coverage when
the true air-flow coefficient is varied, the other coefficients are held at their least squares values, and \( \sigma^2=s^2 \).
[Figure 12.6.1](#fig-ci-coverage) shows the result. With the prior of [Chapter 7](../ch07-optimality/index.html), the coverage is \( 0.959 \) when the truth is the prior
mean \( 1 \), and \( 0.936 \) at \( 0.7 \) and at \( 1.5 \). It ranges from \( 0.877 \) to \( 0.965 \) over the plotted
range. This prior is mild, and prior–data conflict inflates \( b_n \), which widens the interval. A prior with slope standard deviations five times smaller behaves very differently: its coverage is \( 1.000 \) at the prior
mean but \( 0.007 \) at \( 0.7 \), close to the least squares estimate, and \( 0.000 \) at \( 1.5 \). Averaged over the prior of
[Chapter 7](../ch07-optimality/index.html), the simulated coverage is \( 0.9512 \) (\( 40000 \) draws), as @prp-ci-average-coverage requires. The flat-prior interval is
the \( t \) interval and covers exactly \( 95\% \) at every parameter value.

::: {when-format="html"}
![**Figure 12.6.1.** Frequentist coverage of the \( 95\% \) credible interval for the air-flow coefficient against the true coefficient, for the prior of Chapter 7 and for a prior with slope standard deviations divided by five (\( 4000 \) data sets per point). The flat-prior interval covers \( 0.95 \) everywhere.](credible_coverage.svg){#fig-ci-coverage width=66%}
:::

::: {when-format="pdf"}
![Frequentist coverage of the \( 95\% \) credible interval for the air-flow coefficient against the true coefficient, for the prior of Chapter 7 and for a prior with slope standard deviations divided by five (\( 4000 \) data sets per point). The flat-prior interval covers \( 0.95 \) everywhere.](credible_coverage.pdf){width=66%}
:::

An informative prior makes intervals shorter, and better when the prior is sound; when it is not, they are confidently wrong. The
frequentist interval makes no bet and pays in length. [Chapter 27](../ch27-shrinkage/index.html) meets the same
trade-off as bias against variance.

## Exercises

### A. Check your understanding

::: {#exr-ci-hpd-sigma}
[A1]

Is the equal-tailed interval of @thm-ci-bayes-credible(d) the HPD interval for \( \sigma^2 \)? Describe the HPD interval and explain the difference.
:::

::: {.solution}
No. The inverse gamma density is skewed to the right, so the HPD interval \( [l,u] \), with \( \pi(l\mid\y)=\pi(u\mid\y) \), is
shifted toward smaller values. It also depends on the parameterization: the HPD interval for \( \log\sigma^2 \) does not transform into the one for
\( \sigma^2 \), while equal-tailed intervals do.
:::

::: {#exr-ci-predictive-moments}
[A2]

For \( \nu>2 \), find the posterior predictive mean and variance of a single future observation at \( \x_0 \). Why is the variance larger than
\( \tau^2(1+\x_0\T\V_n\x_0) \)?
:::

::: {.solution}
From @thm-ci-bayes-credible(c), \( Y_0\mid\y=\x_0\T\mathbf{m}_n+\tau\sqrt{1+\x_0\T\V_n\x_0}\,T \) with \( T\sim t(\nu) \), so the mean is
\( \x_0\T\mathbf{m}_n \) and the variance is \( \tau^2(1+\x_0\T\V_n\x_0)\,\nu/(\nu-2) \). The extra factor comes from the uncertainty about
\( \sigma^2 \): given \( \sigma^2 \) the predictive variance is \( \sigma^2(1+\x_0\T\V_n\x_0) \), and averaging over the posterior gives
\( \E(\sigma^2\mid\y)=b_n/(a_n-1) \), which exceeds \( \tau^2=b_n/a_n \) by the factor \( a_n/(a_n-1)=\nu/(\nu-2) \).
:::

### B. Practice

::: {#exr-ci-g-prior-interval}
[B1]

Under the \( g \)-prior of @prp-opt-shrinkage(c) with \( \mathbf{m}_0=\bzero \), show that the credible interval for \( \blambda\T\bbeta \) is centred at
\( \{g/(1+g)\}\blambda\T\hbeta \) and has half-width \( t_{\nu,\alpha/2}\,\tau\sqrt{g/(1+g)}\sqrt{\blambda\T(\X\T\X)^{-1}\blambda} \). Compare it with the
\( t \) interval as \( g\to\infty \), paying attention to \( \nu \) and \( \tau^2 \).
:::

::: {.solution}
With \( \V_0=g(\X\T\X)^{-1} \), \( \V_n=\{g/(1+g)\}(\X\T\X)^{-1} \) and \( \mathbf{m}_n=\{g/(1+g)\}\hbeta \) (@prp-opt-shrinkage(c)). Substitute in
@thm-ci-bayes-credible(a). As \( g\to\infty \) the centre and the matrix tend to those of least squares, but \( \nu=2a_0+n \) stays larger than \( n-p \) and
\( \tau^2=b_n/a_n \) tends to \( (2b_0+\text{SSE})/(2a_0+n) \), not to \( s^2 \). The interval does not become the \( t \) interval: that needs the flat prior on
\( \sigma^2 \) as well, with \( a_0=-p/2 \) and \( b_0=0 \) (@cor-opt-flat-prior).
:::

::: {#exr-ci-joint-predictive}
[B2]

Two future observations are taken at the same \( \x_0 \). Show that, under the posterior predictive distribution, they are correlated, and find the
correlation. Why are they correlated when their errors are independent? Are they independent under the frequentist model?
:::

::: {.solution}
By @thm-ci-bayes-credible(c) with \( \X_0 \) having two equal rows \( \x_0\T \), the scale matrix is \( \tau^2\bigl(\I_2+h\bone\bone\T\bigr) \) with
\( h=\x_0\T\V_n\x_0 \), so the correlation is \( h/(1+h) \). Both observations share the unknown \( \bbeta \) (and \( \sigma^2 \)), and uncertainty about the
shared mean correlates them. Under the frequentist model with fixed parameters the two observations are independent, but the two prediction errors
\( Y_{0j}-\hat{Y}_0 \) share \( \hat{Y}_0 \) and have the same correlation \( h_0/(1+h_0) \) under the flat prior.
:::

### C. Going deeper

::: {#exr-ci-flat-rank-deficient}
[C1]

Let \( \X \) have rank \( r<p \) and take the flat prior \( \pi(\bbeta,\sigma^2)\propto1/\sigma^2 \). Show that the posterior of \( \bbeta \) is improper, but
that for estimable \( \blambda\T\bbeta \) the posterior of \( (\blambda\T\bbeta,\sigma^2) \) can be defined by reparameterizing with \( \bgamma=\bU\T\bbeta \) and
\( \boldsymbol{\updelta}=\mathbf{N}\T\bbeta \), where the columns of \( \bU \) and \( \mathbf{N} \) are orthonormal bases of \( \C(\X\T) \) and \( \Null(\X) \), and integrating over nothing
but \( \bgamma \). Show that the resulting credible interval is the \( t \) interval of @thm-ci-estimable-interval with \( n-r \) degrees of freedom, provided
the flat prior is taken as \( 1/\sigma^2 \) in \( (\bgamma,\sigma^2) \) alone.
:::

::: {.solution}
The likelihood depends on \( \bbeta \) only through \( \X\bbeta=\X\bU\bgamma \), since \( \X\mathbf{N}=\bzero \). So the posterior is flat in \( \boldsymbol{\updelta} \), which is not
integrable: the posterior of \( \bbeta \) is improper. In the model \( \Y\sim\Normal_n(\X\bU\bgamma,\sigma^2\I) \) the matrix \( \X\bU \) has full column rank \( r \),
and @cor-opt-flat-prior with \( p \) replaced by \( r \) gives a proper posterior with \( t(n-r) \) marginals. An estimable \( \blambda \) lies in
\( \C(\X\T)=\C(\bU) \), so \( \blambda=\bU\bU\T\blambda \) and \( \blambda\T\bbeta=(\bU\T\blambda)\T\bgamma \) involves \( \bgamma \) only. Its credible interval is centred at
\( (\bU\T\blambda)\T\hat{\bgamma}=\blambda\T\hbeta \) with scale \( s^2(\bU\T\blambda)\T(\bU\T\X\T\X\bU)^{-1}(\bU\T\blambda) \), and \( \bU(\bU\T\X\T\X\bU)^{-1}\bU\T \) is a
generalized inverse of \( \X\T\X \) (check \( \X\T\X\,\bU(\bU\T\X\T\X\bU)^{-1}\bU\T\X\T\X=\X\T\X \) using \( \X=\X\bU\bU\T \)), so the scale is \( s^2\blambda\T\G\blambda \).
:::

::: {#exr-ci-bayes-calibration}
[C2]

In the calibration model of [Section 12.5](05-calibration.html), give \( x_0 \) a normal prior independent of the regression parameters, and use the
flat prior \( 1/\sigma^2 \) for \( (\beta_0,\beta_1,\sigma^2) \). Write down the posterior density of \( x_0 \) up to a constant, and show that it is proper even when
the slope is not significant. Compare the behaviour of the resulting credible interval with the Fieller set in the weak assay of @exm-ci-assay.
:::
