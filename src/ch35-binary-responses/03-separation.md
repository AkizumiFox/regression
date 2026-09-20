# Separation

@prp-bin-irls says the logistic log-likelihood is strictly concave when the model
matrix has full column rank. That guarantees at most one maximum, not at least one,
and for binary data the failure is common: it happens whenever some linear combination
of the covariates predicts the response perfectly, as in small samples, with rare
outcomes, and whenever a categorical covariate has a level with only one outcome in
it.

## Three configurations

Write \( s_i=2y_i-1\in\{-1,1\} \) for the signed response of an ungrouped binary
observation. The logistic log-likelihood is then
\[
\ell(\bbeta)=\sum_{i=1}^n\bigl\{y_i\log\pi_i+(1-y_i)\log(1-\pi_i)\bigr\}
=-\sum_{i=1}^n\log\bigl\{1+\exp(-s_i\x_{(i)}\T\bbeta)\bigr\},
\]{#eq-bin-signed-loglik}

because \( \pi_i=1/\{1+\exp(-\x_{(i)}\T\bbeta)\} \) and \( 1-\pi_i \) is the same expression
with the sign of the linear predictor reversed. Every term is negative and rises
towards zero as \( s_i\x_{(i)}\T\bbeta \) increases, so the whole question is whether
some direction raises all of them at once.

::: {#def-bin-separation}
[Separation and overlap]

Let \( \X \) be \( n\times p \) with rows \( \x_{(i)}\T \), and let \( s_i=2y_i-1 \). The data are

::: {.enumerate options="label=(\alph*)"}
1. **completely separated** if there is a \( \bb\ne\bzero \) with \( s_i\x_{(i)}\T\bb>0 \) for
   every \( i \);

2. **quasi-completely separated** if they are not completely separated but there is a
   \( \bb\ne\bzero \) with \( s_i\x_{(i)}\T\bb\ge0 \) for every \( i \);

3. in **overlap** if \( \bb=\bzero \) is the only \( \bb \) with
   \( s_i\x_{(i)}\T\bb\ge0 \) at every observation.
:::

A vector \( \bb \) of the kind in (a) or (b) is a **separating direction**.
:::

The three cases are mutually exclusive and exhaustive, and a linear programme decides
which holds: maximize \( t \) over \( \bb \) in the unit cube subject to
\( s_i\x_{(i)}\T\bb\ge t \). A positive optimum is complete separation; an optimum of
zero attained only at \( \bb=\bzero \) is overlap, and otherwise quasi-complete
separation.

## Existence

::: {#thm-bin-separation}
[Existence and uniqueness of the maximum likelihood estimate]

Let \( \rank(\X)=p \) and let \( \ell \) be the logistic log-likelihood of @eq-bin-signed-loglik.
Then:

::: {.enumerate options="label=(\alph*)"}
1. if the data are completely or quasi-completely separated, \( \ell \) has no
   maximizer: for every \( \bbeta \) there is a \( \bbeta' \) with \( \ell(\bbeta')>\ell(\bbeta) \);

2. under complete separation \( \sup_{\bbeta}\ell(\bbeta)=0 \), and the supremum is
   approached along any separating ray \( \bbeta=t\bb \), \( t\to\infty \), with every fitted
   probability tending to \( 0 \) or \( 1 \);

3. if the data are in overlap, \( \ell \) attains its supremum at exactly one
   \( \hbeta\in\Real^p \), and \( \hbeta \) solves @eq-bin-score.
:::

:::

::: {.proof}
(a) Let \( \bb\ne\bzero \) satisfy \( s_i\x_{(i)}\T\bb\ge0 \) for all \( i \). Not every one of
these numbers can be zero: if \( \x_{(i)}\T\bb=0 \) for all \( i \) then \( \X\bb=\bzero \), and
\( \rank(\X)=p \) forces \( \bb=\bzero \). So \( s_k\x_{(k)}\T\bb>0 \) for at least one \( k \).
Each summand of @eq-bin-signed-loglik is a strictly decreasing function of
\( -s_i\x_{(i)}\T\bbeta \), hence a nondecreasing function of \( t \) along
\( \bbeta+t\bb \), and the \( k \)th is strictly increasing. Therefore
\( \ell(\bbeta+\bb)>\ell(\bbeta) \) for every \( \bbeta \), and no maximizer exists.

(b) Under complete separation put \( \delta=\min_i s_i\x_{(i)}\T\bb>0 \). Then
\( \ell(t\bb)\ge-n\log(1+e^{-t\delta})\to0 \), while \( \ell<0 \) everywhere, so the
supremum is \( 0 \) and is not attained. Since \( s_i\x_{(i)}\T(t\bb)\to\infty \), the
fitted probability \( \pi_i \) tends to \( 1 \) when \( y_i=1 \) and to \( 0 \) when \( y_i=0 \).

(c) First, \( \ell \) is strictly concave. Its Hessian is \( -\X\T\W\X \) with
\( \W=\diag\{\pi_i(1-\pi_i)\} \) positive definite (@prp-bin-irls), and
\( \bb\T\X\T\W\X\bb=0 \) forces \( \X\bb=\bzero \), hence \( \bb=\bzero \). Second, \( \ell \) is
coercive: \( \ell(\bbeta)\to-\infty \) as \( \norm{\bbeta}\to\infty \). Suppose not. Then
there are \( \bbeta_k \) with \( \norm{\bbeta_k}\to\infty \) and \( \ell(\bbeta_k)\ge-c \) for
some \( c<\infty \). Put \( \bb_k=\bbeta_k/\norm{\bbeta_k} \); since the unit sphere is
compact, a subsequence converges to some \( \bb \) with \( \norm{\bb}=1 \). If
\( s_i\x_{(i)}\T\bb<0 \) for some \( i \), then along that subsequence
\( s_i\x_{(i)}\T\bbeta_k=\norm{\bbeta_k}\,s_i\x_{(i)}\T\bb_k\to-\infty \), so the \( i \)th
term of @eq-bin-signed-loglik tends to \( -\infty \); as all the other terms are at
most \( 0 \), this gives \( \ell(\bbeta_k)\to-\infty \), a contradiction. Hence
\( s_i\x_{(i)}\T\bb\ge0 \) for every \( i \) with \( \bb\ne\bzero \), contradicting overlap.
So \( \ell \) is coercive, and the set \( \{\bbeta:\ell(\bbeta)\ge\ell(\bzero)\} \) is
closed and bounded. A continuous function attains its maximum on a compact set, and
strict concavity makes the maximizer unique. It is an interior maximum of a
differentiable function, so the gradient vanishes there, which is @eq-bin-score.
:::

Parts (a) and (c) say that **the maximum likelihood estimate exists if and only if
the data are in overlap**, a characterization due to Silvapulle (1981) and, in this
form, to Albert and Anderson (1984).

::: {.remark}
[What quasi-complete separation leaves behind]

Under quasi-complete separation let \( \bb \) be a separating direction and
\( S=\{i:\x_{(i)}\T\bb=0\} \). Along \( \bbeta+t\bb \) the terms outside \( S \) tend
to zero and those inside do not move, so the supremum of \( \ell \) is that of the
subsample \( S \) alone, approached with the coefficients along \( \bb \) running off to
infinity while the rest settle down: which is why software reports some coefficients
that look sensible alongside others in the tens with standard errors in the millions.
:::

## What the computer does

Fisher scoring does not know any of this: given separated data it climbs the ray until
it hits the iteration limit or the weights underflow, so the size of the reported
coefficients is set by the iteration count, and the standard errors are larger still,
because \( (\X\T\W\X)^{-1} \) explodes as the weights vanish.

::: {#exm-bin-separated}
[A separated subgroup]

Of the \( 944 \) respondents of @exm-bin-anes,
\( 13 \) reported that their schooling stopped at or before the
eighth grade.
Among them party identification predicts the intended vote perfectly: the ten at
\( 0 \), \( 1 \) or \( 2 \) on the party scale all expected to vote Clinton, the three
at \( 4 \) or \( 6 \) all Dole. The direction \( \bb=(-3,1)\T \), the linear predictor
\( \text{party}-3 \), satisfies \( s_i\x_{(i)}\T\bb>0 \) for all thirteen, so the data
are completely separated and no maximum likelihood estimate exists.

Asked to fit \( \operatorname{logit}\pi=\beta_0+\beta_1\,\text{party} \), software
stopped at ten iterations reports \( \hat\beta_1=8.75 \) with standard
error \( 29.3 \); at twenty-five iterations it reports
\( 23.75 \) with a standard error of order \( 10^{4} \). The deviance falls
below \( 10^{-3} \) after ten steps and below \( 10^{-9} \) after twenty-five, and the
Wald test is nowhere near significant at either stopping point.

Inference is nevertheless possible, from the likelihood ratio. The profile
log-likelihood rises to its supremum of zero, so the set of \( \beta_1 \) not rejected
at the five per cent level is the half-line \( (0.902,\infty) \): a
one-sided interval, correctly reporting that the data bound the effect from below and
not from above.
:::

```{.python .run #cell-separation-data}
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
sub = anes[anes["educ"] == 1]                      # schooling stopped at or before grade eight
X = sm.add_constant(sub[["PID"]]).to_numpy()
y = sub["vote"].to_numpy()

print(pd.crosstab(sub["PID"], sub["vote"]))
for maxiter in (5, 10, 25, 50):
    f = sm.GLM(y, X, family=sm.families.Binomial()).fit(maxiter=maxiter, tol=1e-16)
    print(f"{maxiter:3d} steps: slope {f.params[1]:10.3f}   standard error {f.bse[1]:14.3f}"
          f"   deviance {f.deviance:.2e}")
```

::: {#exm-bin-quasi}
[An empty cell]

Among the \( 127 \) respondents with a doctorate, model the vote by
party identification as a seven-level factor with the independent-Republican
category (\( \text{PID}=4 \)) as baseline, that being the only middle level with both
outcomes present. Three
levels have no variation at all, and exactly \( 3 \) of the six
contrasts run off to infinity. The others behave normally: the contrast for weak
Democrats is \( -3.140 \) with standard error
\( 0.916 \). The deviance, \( 53.808 \), is well
defined, because the diverging coefficients drive their own cells' fitted
probabilities to \( 0 \) or \( 1 \), contributing nothing to it.
:::

## Penalized likelihood

The estimate fails to exist because the likelihood has no interior maximum, so every
remedy adds something that does; three are in general use. Firth (1993) proposed adding
half the log determinant of the Fisher information to the log-likelihood, the log of
the Jeffreys prior density, as a way of removing the leading \( O(n^{-1}) \) bias of
the maximum likelihood estimate; that it also cures separation was noticed later, by
Heinze and Schemper (2002). The modified score equations take a memorable form.

::: {#prp-bin-firth}
[Firth's modified score]

For ungrouped logistic regression with \( \rank(\X)=p \), let
\( \ell^*(\bbeta)=\ell(\bbeta)+\tfrac12\log\det(\X\T\W\X) \) with
\( \W=\diag\{\pi_i(1-\pi_i)\} \). Let
\( h_i=\pi_i(1-\pi_i)\,\x_{(i)}\T(\X\T\W\X)^{-1}\x_{(i)} \) be the \( i \)th diagonal entry of
the weighted hat matrix \( \W^{1/2}\X(\X\T\W\X)^{-1}\X\T\W^{1/2} \). Then
\[
\frac{\partial\ell^*}{\partial\bbeta}
=\X\T\Bigl\{\y-\boldsymbol{\uppi}+\mathbf{h}\circ\bigl(\tfrac12\bone-\boldsymbol{\uppi}\bigr)\Bigr\},
\]{#eq-bin-firth-score}

where \( \circ \) is the entrywise product. Equivalently, the Firth estimate is a fixed
point: it is the ordinary maximum likelihood estimate for pseudo-data in which
observation \( i \) has \( y_i+h_i/2 \) successes out of \( 1+h_i \) trials, with \( h_i \)
evaluated at that same estimate. Along any ray \( t\bb \) that completely
separates the data, \( \ell^*(t\bb)\to-\infty \).
:::

::: {.proof}
Write \( F(\bbeta)=\X\T\W\X \) and \( w_i=\pi_i(1-\pi_i) \). Since
\( \partial\pi_i/\partial\beta_j=x_{ij}w_i \) and \( dw/d\pi=1-2\pi \),
\[
\frac{\partial w_i}{\partial\beta_j}=x_{ij}w_i(1-2\pi_i),
\qquad
\frac{\partial F}{\partial\beta_j}=\X\T\bD_j\X ,
\]
with \( \bD_j=\diag\{x_{ij}w_i(1-2\pi_i)\} \). By the derivative of a log determinant,
\[
\begin{aligned}
\frac{\partial\log\det F}{\partial\beta_j}
&=\tr\bigl(F^{-1}\X\T\bD_j\X\bigr)
=\sum_ix_{ij}w_i(1-2\pi_i)\,\x_{(i)}\T F^{-1}\x_{(i)}\\
&=\sum_ix_{ij}(1-2\pi_i)h_i .
\end{aligned}
\]
Halving and adding the ordinary score \( \sum_ix_{ij}(y_i-\pi_i) \) gives
@eq-bin-firth-score, since \( \tfrac12(1-2\pi_i)=\tfrac12-\pi_i \). For the pseudo-data statement, the ordinary
score for \( y_i^*=y_i+h_i/2 \) successes out of \( m_i^*=1+h_i \) trials is
\( \sum_ix_{ij}(y_i^*-m_i^*\pi_i) \), which is the same expression.

For the last claim, let \( \delta=\min_is_i\x_{(i)}\T\bb>0 \). At \( \bbeta=t\bb \) every
fitted probability is within \( e^{-t\delta} \) of \( 0 \) or \( 1 \), so
\( w_i\le e^{-t\delta} \) and \( \X\T\W\X\preceq e^{-t\delta}\X\T\X \) in the
nonnegative definite order. Hence
\( \log\det(\X\T\W\X)\le-pt\delta+\log\det(\X\T\X) \), while \( \ell\le0 \). Therefore
\( \ell^*(t\bb)\le-\tfrac12pt\delta+\text{constant}\to-\infty \).
:::

The penalty pulls the estimate towards zero exactly where the likelihood is flat:
along the offending ray it falls linearly while the likelihood has nothing left to
gain. That the Firth estimate is finite for *every* configuration with
\( \rank(\X)=p \) was proved by Kosmidis and Firth (2021), by a longer argument.

::: {#exm-bin-firth}
[The separated subgroup, penalized]

Applied to @exm-bin-separated, Firth's method gives \( \hat\beta_1=1.362 \)
with standard error \( 0.708 \) and intercept
\( -4.049 \). The Wald statistic is \( 1.92 \), short of
the conventional threshold, but the interval from inverting the *penalized* likelihood
ratio is \( (0.364,\,4.045) \), which excludes zero
comfortably. With thirteen observations and a log-likelihood as asymmetric as
[Figure 35.3.1](03-separation.html#fig-bin-separation) shows, the quadratic
approximation behind the Wald statistic is simply not available.
:::

::: {when-format="html"}
![**Figure 35.3.1.** Profile log-likelihoods for the slope in @exm-bin-separated, each
relative to its own supremum. The ordinary profile (solid) increases to zero and is
never attained, so the interval cut out by the horizontal line at
\( -\chi^2_{0.95}(1)/2 \) is a half-line. Firth's penalized profile (dashed) has an
interior maximum, marked by the vertical
line.](separation_profile.svg){#fig-bin-separation width=68%}
:::

::: {when-format="pdf"}
![Profile log-likelihoods for the slope in @exm-bin-separated, each
relative to its own supremum. The ordinary profile (solid) increases to zero and is
never attained, so the interval cut out by the horizontal line at
\( -\chi^2_{0.95}(1)/2 \) is a half-line. Firth's penalized profile (dashed) has an
interior maximum, marked by the vertical line.](separation_profile.pdf){width=68%}
:::

```{.python .run #cell-separation-firth}
def firth(X, y, steps=200):
    """Maximize the log-likelihood penalized by (1/2) log det(X'WX): Firth's estimate."""
    beta = np.zeros(X.shape[1])
    for _ in range(steps):
        mu = 1 / (1 + np.exp(-X @ beta))
        w = mu * (1 - mu)
        F = X.T @ (w[:, None] * X)                          # the expected information
        h = w * np.einsum("ij,jk,ik->i", X, np.linalg.inv(F), X)
        beta = beta + np.linalg.solve(F, X.T @ (y - mu + h * (0.5 - mu)))
    return beta


beta_firth = firth(X, y)
print("Firth estimate:", beta_firth.round(3))
```

## Ridge and Cauchy priors

A quadratic penalty works too, for a simpler reason. Maximizing
\( \ell(\bbeta)-\tfrac12\lambda\norm{\bbeta}^2 \) is the logistic analogue of ridge
regression (@def-shr-ridge), and the objective is strictly concave and coercive
whatever the data, so a unique maximizer always exists (@exr-bin-ridge-exists). It is
the posterior mode under independent \( \Normal(0,\lambda^{-1}) \) priors, which is how
[Chapter 30](../ch30-regularization-boosting/index.html) reads every penalty in @prp-reg-map.
In practice, and in the code below, the intercept is left out of the penalty, so that
the prior reading is a flat one there; the maximizer still exists, provided both
outcomes occur (@exr-bin-ridge-exists). The cost is bias towards zero that depends on
the scale of the covariates, and a \( \lambda \) that must come from somewhere, usually
cross-validation ([Chapter 29](../ch29-model-selection/index.html)).

```{.python .run #cell-separation-ridge}
def ridge_logistic(X, y, lam, steps=200):
    """Maximize the log-likelihood minus (lam/2) times the squared norm of the slopes."""
    d = np.ones(X.shape[1])
    d[0] = 0.0                                               # leave the intercept unpenalized
    beta = np.zeros(X.shape[1])
    for _ in range(steps):
        mu = 1 / (1 + np.exp(-X @ beta))
        w = mu * (1 - mu)
        H = X.T @ (w[:, None] * X) + lam * np.diag(d)
        beta = beta + np.linalg.solve(H, X.T @ (y - mu) - lam * d * beta)
    return beta


for lam in (0.1, 1.0, 10.0):
    print(f"ridge, lambda = {lam:4.1f}: {ridge_logistic(X, y, lam).round(3)}")
```

On the separated subgroup the ridge estimate of the slope is
\( 2.601 \) at \( \lambda=0.1 \), \( 1.356 \) at
\( \lambda=1 \) and \( 0.502 \) at \( \lambda=10 \): the whole range is
"the answer" for some penalty, which is the point. Gelman, Jakulin, Pittau and Su
(2008) proposed instead independent Cauchy priors of scale \( 2.5 \) on standardized
covariates, whose heavy tails leave a genuinely large coefficient nearly alone. Firth's
penalty alone is invariant to how the covariates are parameterized.

## Exact conditional inference

A third route abandons the normal approximation. The sufficient statistic for
\( \bbeta \) is \( \X\T\y \), by @lem-opt-factorization applied
to @eq-bin-signed-loglik, so the conditional distribution of the component belonging to
one coefficient, given the rest, depends on that coefficient alone. Enumerating it
gives exact \( p \)-values and a median-unbiased estimate that is finite even under
separation. The idea is Cox's (1970) and the network algorithms are those of Hirji,
Mehta and Patel (1987); the limit is computational.

::: {.idea}
Separation is not a defect of the data; it is the data saying that the effect is
larger than the sample can measure, and every remedy decides how much larger. A
one-sided likelihood ratio interval is the most honest of them, because it does not
decide.
:::

## Exercises

### A. Check your understanding

::: {#exr-bin-separation-cell}
[A1]

A model contains an indicator for one level of a categorical covariate, and every
observation at that level has \( y_i=1 \). Exhibit a separating direction, say which of
the three cases of @def-bin-separation holds, and describe the fit.
:::

::: {.solution}
Let \( \bb \) have \( 1 \) in the coordinate of that indicator and zeros elsewhere. For
observations at that level \( s_i\x_{(i)}\T\bb=+1>0 \); for the others
\( \x_{(i)}\T\bb=0 \). So \( \bb \) separates with some margins zero: quasi-complete
separation, unless the remaining data are themselves completely separated. The
coefficient of that indicator diverges, its fitted probabilities go to \( 1 \), and the
other coefficients converge to the estimates from the remaining observations.
:::

### B. Practice

::: {#exr-bin-ridge-exists}
[B1]

Show that for any data and any \( \lambda>0 \) the function
\( \bbeta\mapsto\ell(\bbeta)-\tfrac12\lambda\norm{\bbeta}^2 \) has exactly one
maximizer, without assuming \( \rank(\X)=p \). Show also that if \( \rank(\X)=p \) and
both outcomes occur the same is true when the intercept is left unpenalized.
:::

::: {.solution}
By @eq-bin-signed-loglik, \( \ell\le0 \), so the objective is at most
\( -\tfrac12\lambda\norm{\bbeta}^2 \), which tends to \( -\infty \): the objective is
coercive and continuous, so it attains its maximum on a compact sublevel set. It is
the sum of the concave \( \ell \) (@exr-bin-concavity) and a strictly concave quadratic,
hence strictly concave, so the maximizer is unique. Neither step used the rank of
\( \X \). With the intercept \( \beta_0 \) unpenalized, write \( \bgamma \) for the
slopes. If \( \norm{\bgamma}\to\infty \) the penalty alone sends the objective to
\( -\infty \), since \( \ell\le0 \); if \( \bgamma \) stays bounded and
\( |\beta_0|\to\infty \), then some observation has
\( s_i\x_{(i)}\T\bbeta\to-\infty \), because an observation with \( y_i=1 \) supplies one
as \( \beta_0\to-\infty \) and one with \( y_i=0 \) as \( \beta_0\to+\infty \). So the
objective is again coercive, and strict concavity now comes from \( \ell \) itself.
:::

::: {#exr-bin-firth-2x2}
[B2]

Take a \( 2\times2 \) table with a single binary covariate, coded as an intercept and
an indicator. Compute the \( h_i \) of @prp-bin-firth, and show that Firth's pseudo-data
amount to adding \( \tfrac12 \) to each of the four cell counts. What is the resulting
estimate of the log odds ratio?
:::

::: {#exr-bin-hauck-donner}
[B3]

For a single binary covariate with fitted probabilities \( \hat\pi_1 \) and
\( \hat\pi_0 \) in groups of sizes \( m_1 \) and \( m_0 \), show that as
\( \hat\pi_1\to1 \), with everything else fixed, the Wald statistic tends to zero while
the likelihood ratio statistic increases to a positive limit.
:::

### C. Going deeper

::: {#exr-bin-separation-high-dimensional}
[C1]

Suppose the rows \( \x_{(i)} \) are independent draws from a distribution symmetric
about the origin in \( \Real^p \) and the responses are independent coin flips. Cover's
counting argument gives the probability of complete separation as
\( 2^{1-n}\sum_{k=0}^{p-1}\binom{n-1}{k} \). Deduce that separation is almost certain
when \( n<2p \) and almost impossible when \( n\gg2p \), and relate this to
[Chapter 28](../ch28-high-dimensional/index.html).
:::

::: {.solution}
The sum is \( \Pr\{B\le p-1\} \) with \( B\sim\text{Binomial}(n-1,\tfrac12) \), so with
\( \E B=(n-1)/2 \) it is near \( 1 \) when \( p-1 \) exceeds \( (n-1)/2 \), roughly
\( n<2p \), and near \( 0 \) when \( p \) is well below \( n/2 \), the transition having
width of order \( \sqrt n \). So with pure noise maximum likelihood fails to exist below
about \( n=2p \): twice the sample size at which least squares begins to interpolate in
[Chapter 28](../ch28-high-dimensional/index.html). Both thresholds mark the point past
which an unpenalized fit reproduces the labels exactly, and both are why penalization
is the default in high dimensions.
:::
