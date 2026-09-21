# Variance functions

Fitting a generalized linear model needs a family, a link and a linear predictor.
[Section 34.3](../ch34-exponential-families-glm/03-likelihood-equations.html) noticed that
the family enters the likelihood equations
\[
\sum_{i=1}^n x_{ij}\,\frac{w_i\,h'(\eta_i)}{V(\mu_i)}\,(y_i-\mu_i)=0,
\qquad j=1,\dots,p,
\]{#eq-ql-equations}

only through the variance function \( V \): two densities with the same \( V \) produce
the same estimate. The apparatus of @def-glm-edf — a natural parameter, a cumulant
function, a density that integrates to one — *derives* @eq-ql-equations and then
disappears from it. This section puts the variance function first, and the model becomes a
mean and a variance and nothing else.

## Specifying the second moment only

::: {#def-ql-variance-function}
[A variance specification]

Let \( Y_1,\dots,Y_n \) be independent, let \( w_1,\dots,w_n>0 \) be known weights, let
\( g \) be a strictly monotone differentiable link with inverse \( h \), and let \( V \) be a
positive, continuously differentiable function on an open interval \( \mathcal M\subseteq\Real \).
A **variance specification** is the pair of assumptions
\[
\E(Y_i)=\mu_i=h\bigl(\x_{(i)}\T\bbeta\bigr),
\qquad
\Var(Y_i)=\frac{\phi\,V(\mu_i)}{w_i},
\]{#eq-ql-specification}

with \( \bbeta\in\Real^p \) and \( \phi>0 \) unknown. The function \( V \) is the **variance
function** and \( \phi \) the **dispersion**. No distribution is assumed for \( Y_i \) beyond
these two moments.
:::

The definition is weaker than @def-glm-model in saying nothing about higher moments, the
support, or whether any distribution with these moments exists, and stronger in leaving
\( \phi \) free even for a \( V \) that came from a family in which it was fixed. The
specifications in constant use are

| name | \( V(\mu) \) | \( w_i \) | \( \phi \) | \( \mathcal M \) |
|---|---|---|---|---|
| quasi-normal | \( 1 \) | \( 1 \) | free | \( \Real \) |
| quasi-Poisson | \( \mu \) | \( 1 \) | free | \( (0,\infty) \) |
| quasi-binomial | \( \mu(1-\mu) \) | \( m_i \) | free | \( (0,1) \) |
| quasi-gamma | \( \mu^2 \) | \( 1 \) | free | \( (0,\infty) \) |
| power | \( \mu^{\xi} \) | \( 1 \) | free | \( (0,\infty) \) |
| NB2 | \( \mu+\mu^2/\kappa \) | \( 1 \) | \( 1 \) | \( (0,\infty) \) |

where \( y_i \) is a group proportion in the quasi-binomial row and \( m_i \) the number
of trials behind it. A value \( \phi>1 \) is **overdispersion** relative to the family that
supplied \( V \), and \( \phi<1 \) **underdispersion**. The first four rows are variance
functions of @def-glm-variance with the dispersion set free; the last is that
of @def-cnt-negbin(a), for which a genuine likelihood exists.

::: {.remark}
[When is a variance specification a distribution?]

For \( \phi \) free and \( V(\mu)=\mu^{\xi} \), Jørgensen's classification of the Tweedie
families (@exr-glm-tweedie) settles whether some exponential dispersion family has exactly
these moments: one exists for \( \xi\le0 \) and for \( \xi\ge1 \), and for no \( \xi \)
strictly between. Existence is not usefulness. At \( \xi=1 \) the Tweedie family with
dispersion \( \phi \) is \( \phi\,W \) with \( W \) Poisson of mean \( \mu/\phi \): the
right two moments, but on the multiples of \( \phi \) rather than the integers. So no
*exponential dispersion family* on the counts has mean \( \mu \) and variance
\( \phi\mu \) for one \( \phi \) and all \( \mu \). A family of counts with exactly those
two moments does exist — the NB1 model of @def-cnt-negbin(b) — but it is not an exponential
dispersion family, and its likelihood equations are not @eq-ql-equations; that is the point
of the third subsection of this section.
:::

::: {.warning}
The quasi-binomial specification is empty when \( m_i=1 \). An ungrouped binary response
has \( Y_i^2=Y_i \), so \( \Var(Y_i)=\mu_i(1-\mu_i) \) exactly for every distribution on
\( \{0,1\} \) with mean \( \mu_i \), and there is no room for a \( \phi \). A dispersion
estimate far from \( 1 \) in an ungrouped logistic fit is a statement about the mean
model — an omitted term, a wrong link, a few gross outliers — never about the variance,
the second-moment counterpart of @thm-glm-deviance(c).
:::

## Two ways to inflate a binomial

Where a dispersion is available it can take more than one shape, and the shapes are
distinguishable. Suppose observation \( i \) averages \( m_i \) dependent binary
trials.

::: {#lem-ql-exchangeable}
[Exchangeable trials]

Let \( Z_1,\dots,Z_m \) be binary variables with common mean \( \pi \) and
\( \Cov(Z_s,Z_t)=\rho\,\pi(1-\pi) \) for all \( s\ne t \), so that \( \rho \) is
their common correlation, and let \( \bar Z=m^{-1}\sum_tZ_t \). Then
\[
\E(\bar Z)=\pi,\qquad
\Var(\bar Z)=\bigl\{1+\rho(m-1)\bigr\}\,\frac{\pi(1-\pi)}{m},
\]{#eq-ql-exchangeable}

and \( \rho\ge-1/(m-1) \), with \( \rho\ge0 \) whenever the trials are conditionally
independent and identically distributed given some unobserved variable.
:::

::: {.proof}
Each \( Z_t \) has variance \( \pi(1-\pi) \), so the stated covariance is \( \rho \)
times it and
\( \Var(\bar Z)=m^{-2}\{m\pi(1-\pi)+m(m-1)\rho\pi(1-\pi)\} \), which
is @eq-ql-exchangeable; nonnegativity gives \( 1+\rho(m-1)\ge0 \). If the \( Z_t \) are
independent given \( U \) with \( \Pr(Z_t=1\mid U)=\Pi \) the same for each \( t \), the
law of total covariance gives \( \Cov(Z_s,Z_t)=\Var(\Pi)\ge0 \).
:::

The variance specification suggested by @eq-ql-exchangeable is
\[
V_{\rho}(\mu)=\bigl\{1+\rho(m_i-1)\bigr\}\,\mu(1-\mu),\qquad w_i=m_i,\qquad \phi=1,
\]
in which the inflation grows with the group size, against the quasi-binomial
\( V(\mu)=\mu(1-\mu) \) with \( w_i=m_i \) and a constant \( \phi \). The two agree at
one group size and nowhere else
([Figure 38.1.1](01-variance-functions.html#fig-ql-binomial-inflation)): at
\( \phi=1.8 \) and \( \rho=0.05 \) they cross at \( m=17.0 \) and
differ by a factor \( 1.639 \) at \( m=40 \). Spread-out group sizes
separate them, by a plot of squared Pearson residuals against \( m_i \); equal ones do
not.

::: {when-format="html"}
![**Figure 38.1.1.** Two variance functions for a proportion from \( m \) trials, as
multiples of \( \pi(1-\pi)/m \). A constant dispersion inflates every group equally;
exchangeable trials inflate large groups more, by @eq-ql-exchangeable.](binomial_inflation.svg){#fig-ql-binomial-inflation width=58%}
:::

::: {when-format="pdf"}
![Two variance functions for a proportion from \( m \) trials, as
multiples of \( \pi(1-\pi)/m \). A constant dispersion inflates every group equally;
exchangeable trials inflate large groups more, by @eq-ql-exchangeable.](binomial_inflation.pdf){width=58%}
:::

```{.python .run #cell-variance-shapes-binomial}
import numpy as np

m_grid = np.arange(1, 41)
phi, rho = 1.8, 0.05
inflated = np.full_like(m_grid, phi, dtype=float)      # V(pi) = phi pi(1-pi)/m
exchangeable = 1 + rho * (m_grid - 1)                  # V(pi) = {1+rho(m-1)} pi(1-pi)/m
crossing = 1 + (phi - 1) / rho                         # the one group size where they agree
print(f"the two variance functions agree at m = {crossing:.1f}")
print(f"at m = 40 they differ by a factor {exchangeable[-1] / phi:.3f}")
```

The second specification is not an invention of the quasi-likelihood literature: the
**beta-binomial** distribution, a binomial mixed over a beta success probability, has
exactly this variance (@exr-ql-beta-binomial), so a mixture and a correlation leave the
same trace in the second moment.

## A variance function does not determine a fit

A variance function settles the estimating equations, which is less than settling the fit:
two models can share a variance function and disagree about \( \hbeta \), because the
second model's likelihood equations need not be @eq-ql-equations. The NB1 model
of @def-cnt-negbin(b) has variance \( \phi\mu_i \) with \( \phi=1+\alpha^{-1} \), the
quasi-Poisson variance function exactly; but its score @eq-cnt-nb1-score does not contain
the residual \( y_i-\mu_i \) at all, and the two estimates differ in
general (@exr-cnt-quasi-vs-nb), as [Section 38.4](04-overdispersion.html) shows in
numbers. A *model* is a set of distributions; a *variance specification* is two moment
assumptions plus the rule @eq-ql-equations for turning them into an estimate.

## What the missing density costs

Dropping the density is not free. There is no maximized likelihood, so
@prp-glm-three-tests loses one of its three members and
[Section 38.3](03-sandwich.html) must supply substitutes; no Akaike or Bayesian criterion,
since both are \( -2\ell(\hbeta) \) plus a penalty (@def-sel-aic-bic), though
cross-validation (@def-sel-cv) survives; and no predictive distribution, so a fitted
quasi-Poisson model cannot say what fraction of units will record zero, simulate a data
set, or be checked by a predictive comparison, where a negative binomial or beta-binomial
model does all three. What survives is the estimate, its consistency, its asymptotic
normality, and a covariance estimate that does not need the variance function to be
right — the subject of the next two sections.

## Exercises

### A. Check your understanding

::: {#exr-ql-which-phi}
[A1]

For each of the following, say whether a dispersion parameter \( \phi\ne1 \) is
meaningful, and why: (a) counts of insurance claims per policy; (b) the proportion of a
litter of \( 12 \) pups showing a defect; (c) an indicator of whether a patient survived;
(d) the proportion of \( 500 \) sampled voters supporting a candidate, the voters sampled in
clusters of \( 10 \) households.
:::

::: {.solution}
(a) Yes: nothing ties a count's variance to its mean, and heterogeneity between policies
inflates it (@lem-cnt-gamma-mixture(a)). (b) Yes, with \( w_i=12 \): correlation between
littermates, or variation between litters, inflates the variance
by @lem-ql-exchangeable. (c) No: a binary response has variance \( \mu(1-\mu) \) whatever
else is true. (d) Yes, but the group size is the cluster, not \( 500 \); using
\( m=500 \) in @eq-ql-exchangeable would overstate the inflation.
:::

::: {#exr-ql-crossing}
[A2]

Two analysts fit grouped binomial data with group sizes between \( 5 \) and \( 60 \), one
using a constant dispersion \( \phi \) and one using the exchangeable form
@eq-ql-exchangeable with correlation \( \rho \). They report \( \hat\phi=2.5 \) and
\( \hat\rho=0.04 \). At which group size do the two fitted variances agree? For which groups
does the first analyst's model claim the larger variance?
:::

::: {.solution}
Setting \( 1+\rho(m-1)=\phi \) gives \( m=1+(\phi-1)/\rho=1+1.5/0.04=38.5 \); below it
the constant-dispersion model claims the larger variance, above it the smaller.
:::

### B. Practice

::: {#exr-ql-beta-binomial}
[B1]

Let \( \Pi \) have a beta distribution with parameters \( a,b>0 \), so that
\( \E\Pi=\mu=a/(a+b) \) and \( \Var\Pi=\mu(1-\mu)/(a+b+1) \), and let \( m\bar Y \) be
binomial \( (m,\Pi) \) given \( \Pi \). Show that \( \E\bar Y=\mu \) and
\[
\Var(\bar Y)=\bigl\{1+\rho(m-1)\bigr\}\frac{\mu(1-\mu)}{m},
\qquad \rho=\frac1{a+b+1},
\]
so the beta-binomial has the variance function of @eq-ql-exchangeable. Deduce that the
beta-binomial and an exchangeable-trials model with the same \( \rho \) cannot be told apart
by first and second moments.
:::

::: {.solution}
The conditional mean is \( \Pi \) and the conditional variance
\( \Pi(1-\Pi)/m \), so \( \E\bar Y=\mu \) and, by the law of total variance and
\( \E\Pi^2=\Var\Pi+\mu^2 \),
\[
\Var(\bar Y)=\frac{\E\{\Pi(1-\Pi)\}}{m}+\Var\Pi
=\frac{\mu(1-\mu)-\Var\Pi}{m}+\Var\Pi
=\frac{\mu(1-\mu)}{m}\Bigl\{1+(m-1)\rho\Bigr\}
\]
with \( \rho=\Var\Pi/\{\mu(1-\mu)\}=1/(a+b+1) \). The two specifications therefore give
identical first and second moments for every \( m \).
:::

::: {#exr-ql-power-cv}
[B2]

Under the power specification \( \Var(Y)=\phi\mu^{\xi} \), find the coefficient of variation
\( \sqrt{\Var(Y)}/\E(Y) \) as a function of \( \mu \). For which \( \xi \) is it constant?
Which of the specifications in the table of this section have the property that a
multiplicative change of units, \( Y\mapsto cY \), leaves the specification in the same
family with a new \( \phi \)?
:::

### C. Going deeper

::: {#exr-ql-mixture-bound}
[C1]

Let \( Y\mid\Lambda \) have mean \( \Lambda \) and variance \( a\Lambda^{\xi} \), and let
\( \Lambda \) have mean \( \mu \) and finite variance \( \tau^2 \). Show that the
marginal variance is \( a\,\E(\Lambda^{\xi})+\tau^2 \), equal to \( a\mu+\tau^2 \) when
\( \xi=1 \). For \( \xi=2 \), use Jensen's inequality to show that it is at least
\( a\mu^2+\tau^2 \), so a mixture of quasi-gamma responses is overdispersed relative to
the quasi-gamma specification with the same \( a \).
:::
