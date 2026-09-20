# Logistic regression

A binary response records whether something happened: a patient recovered, a
voter turned out, a loan defaulted. Then \( \E(Y_i)=\Pr(Y_i=1) \), so modelling the
mean *is* modelling a probability, and a linear model for it is in trouble before it
starts: \( \x_{(i)}\T\bbeta \) ranges over the whole line while a probability may not
leave \( [0,1] \), and the variance \( \pi_i(1-\pi_i) \) moves with the mean. The
generalized linear model of
[Chapter 34](../ch34-exponential-families-glm/index.html) removes both problems by
putting a link between the mean and the linear predictor.

## The model

Let \( Y_1,\dots,Y_n \) be independent with \( Y_i\sim\text{Binomial}(m_i,\pi_i) \);
\( m_i=1 \) throughout is **ungrouped** data, and larger \( m_i \) arise when several
units share a covariate vector. Dropping the index,
\[
\binom{m}{y}\pi^y(1-\pi)^{m-y}
=\binom{m}{y}\exp\Bigl\{y\log\frac{\pi}{1-\pi}+m\log(1-\pi)\Bigr\},
\]
an exponential dispersion family (@def-glm-edf) for the mean \( \mu=m\pi \), with
natural parameter \( \theta=\log\{\pi/(1-\pi)\} \), cumulant function
\( b(\theta)=\log(1+e^{\theta}) \), dispersion \( \phi=1 \) and weight \( m \).
Differentiating \( b \) twice recovers the moments of @prp-glm-moments:
\( b'(\theta)=\pi \) and \( b''(\theta)=\pi(1-\pi) \), so the variance function is
\( V(\pi)=\pi(1-\pi) \) and \( \Var(Y)=m\pi(1-\pi) \). The natural parameter is the
**log odds**, and the canonical link is therefore the logit.

::: {#def-bin-logistic}
[Logistic regression]

The **logistic regression model** for independent \( Y_i\sim\text{Binomial}(m_i,\pi_i) \) is
\[
\operatorname{logit}\pi_i=\log\frac{\pi_i}{1-\pi_i}=\x_{(i)}\T\bbeta,
\qquad\text{equivalently}\qquad
\pi_i=\frac{\exp(\x_{(i)}\T\bbeta)}{1+\exp(\x_{(i)}\T\bbeta)} .
\]{#eq-bin-logistic}

The quantity \( \pi/(1-\pi) \) is the **odds** of the event, and the ratio of the odds at
two covariate vectors is an **odds ratio**.
:::

The logit is one choice among many ([Section 35.2](02-probit-cloglog.html)). Its claim
to first place is that it is canonical, so @thm-glm-score takes its simplest form, and
that its coefficients say something people want to hear.

## Reading the coefficients

::: {#prp-bin-interpretation}
[Three readings of a logistic coefficient]

In @eq-bin-logistic, write \( \eta_i=\x_{(i)}\T\bbeta \) and let \( x_j \) be a covariate
that appears in the linear predictor only through its own coefficient \( \beta_j \). Then:

::: {.enumerate options="label=(\alph*)"}
1. increasing \( x_j \) by \( c \) with the other covariates unchanged multiplies the odds
   by \( e^{c\beta_j} \), whatever the values of the other covariates;

2. \( \partial\pi/\partial x_j=\beta_j\,\pi(1-\pi) \), so
   \( |\partial\pi/\partial x_j|\le|\beta_j|/4 \), with equality exactly at \( \pi=1/2 \);

3. the average of that derivative over the sample, \( n^{-1}\sum_i\beta_j\hat\pi_i(1-\hat\pi_i) \),
   is the **average marginal effect**: the mean change in fitted probability per unit of
   \( x_j \).
:::

:::

::: {.proof}
(a) The log odds at \( x_j+c \) minus the log odds at \( x_j \) is \( c\beta_j \), and everything
else in \( \eta \) cancels. (b) Differentiating
\( \pi=e^{\eta}/(1+e^{\eta}) \) gives \( d\pi/d\eta=\pi(1-\pi) \), and the chain rule supplies the
factor \( \beta_j \). The function \( \pi(1-\pi) \) on \( [0,1] \) is maximized at \( \pi=1/2 \) with
value \( 1/4 \), so the derivative has modulus at most \( |\beta_j|/4 \). (c) Immediate from (b).
:::

Part (a) is why logistic regression is written as it is: one number summarizes the
effect on the odds at every setting of the other covariates. Part (b) is Gelman and
Hill's "divide by four" rule, sharp near \( \pi=1/2 \); part (c) reports the same
coefficient on the probability scale.

::: {.warning}
None of the three readings is causal: each compares two covariate vectors *within the
fitted model*, and whether that estimates an intervention is the question of
[Chapter 25](../ch25-causal-interpretation/index.html).
[Section 35.5](05-case-control.html) adds a wrinkle special to odds ratios.
:::

::: {#exm-bin-anes}
[Voting intention in the 1996 American election]

The 1996 American National Election Studies survey recorded, for
\( 944 \) respondents, the candidate they expected to vote for
(\( 393 \) said Dole), party identification from strong Democrat (0)
to strong Republican (6), age, education on a seven-point scale and household income
on a 24-point scale. Fitting @eq-bin-logistic with these four covariates and age in
decades gives

| term | \( \hat\beta_j \) | standard error |
|---|---|---|
| intercept | \( -5.471 \) | 0.619 |
| party identification | 1.2188 | 0.0709 |
| age (decades) | 0.1222 | 0.0714 |
| education | 0.0147 | 0.0771 |
| income | 0.0311 | 0.0210 |

One step along the party scale multiplies the odds of a Dole vote by
\( e^{1.2188}=3.383 \); three steps multiply them by
\( 38.72 \). Dividing by four bounds the effect on the probability by
\( 0.305 \) per step, and the average marginal effect is
\( 0.100 \): one step is worth about ten percentage points, averaged over
the respondents ([Figure 35.1.1](01-logistic-regression.html#fig-bin-anes)).
:::

::: {when-format="html"}
![**Figure 35.1.1.** Expected vote for Dole in the 1996 ANES survey. (a) The fitted
logistic curve against party identification, the other covariates held at their
means; the discs are the observed proportions in the seven party groups, with area
proportional to group size. (b) The tangent at the steepest point, of slope
\( \hat\beta_j/4 \).](anes_logistic.svg){#fig-bin-anes width=100%}
:::

::: {when-format="pdf"}
![Expected vote for Dole in the 1996 ANES survey. (a) The fitted
logistic curve against party identification, the other covariates held at their
means; the discs are the observed proportions in the seven party groups, with area
proportional to group size. (b) The tangent at the steepest point, of slope
\( \hat\beta_j/4 \).](anes_logistic.pdf){width=100%}
:::

```{.python .run #cell-anes-logistic-fit}
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
X = pd.DataFrame({"const": 1.0, "PID": anes["PID"], "age": anes["age"] / 10,
                  "educ": anes["educ"], "income": anes["income"]})
y = anes["vote"]                                   # 1 = Dole, 0 = Clinton

fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
pi_hat = fit.fittedvalues

print(fit.params.round(4).to_string())
print("odds ratio for one step on the party scale:", round(np.exp(fit.params["PID"]), 3))
print("divide by four:", round(fit.params["PID"] / 4, 3))
print("average marginal effect:", round((fit.params["PID"] * pi_hat * (1 - pi_hat)).mean(), 3))
```

## Grouped and ungrouped data

Two data sets recording the same information can look very different: one lists
\( 944 \) individual zeros and ones, the other one binomial count per
distinct covariate vector.

::: {#prp-bin-grouping}
[Grouping]

Let \( \y \) be ungrouped binary data, and suppose the rows of \( \X \) take \( N \) distinct
values \( \x_{(1)}^*,\dots,\x_{(N)}^* \), the \( k \)th occurring \( m_k \) times with \( s_k \) successes.
Then:

::: {.enumerate options="label=(\alph*)"}
1. the two log-likelihoods differ by a constant free of \( \bbeta \), so they have the same
   maximizers, the same score, the same information and the same likelihood ratio
   statistics;

2. the two deviances of @def-glm-deviance differ by \( -2\ell_{\mathrm{sat}} \), where
   \( \ell_{\mathrm{sat}}=\sum_k\{s_k\log(s_k/m_k)+(m_k-s_k)\log(1-s_k/m_k)\} \) is the
   maximized log-likelihood of the grouped saturated model;

3. if \( \bone\in\C(\X) \), then \( \sum_i m_i\hat\pi_i=\sum_i y_i \) in either form: the fit
   reproduces the total number of successes.
:::

:::

::: {.proof}
(a) With \( \eta_k=\x_{(k)}^{*\T}\bbeta \), the ungrouped log-likelihood is
\[
\sum_{k=1}^N\bigl\{s_k\eta_k-m_k\log(1+e^{\eta_k})\bigr\},
\]
because the individual contributions within a group depend on the data only through
\( s_k \). The grouped log-likelihood is the same expression plus
\( \sum_k\log\binom{m_k}{s_k} \), a constant. (b) The deviance is twice the
difference between the saturated and fitted log-likelihoods. The fitted
log-likelihood is the same in both forms by (a); the ungrouped saturated model sets
\( \hat\pi_i=y_i \) and attains log-likelihood \( 0 \), while the grouped saturated model
sets \( \hat\pi_k=s_k/m_k \) and attains \( \ell_{\mathrm{sat}}\le0 \). (c) Under the canonical
link the likelihood equations are \( \X\T(\y-\bmu)=\bzero \) (@thm-glm-score); the
column \( \bone \) gives \( \sum_i(y_i-m_i\hat\pi_i)=0 \).
:::

Part (b) is worth pausing over: the deviance of a logistic fit depends on how the
data were typed in. *Differences* of deviances between nested models do not, because
the saturated term cancels, and those are what @thm-glm-deviance licenses.

::: {#exm-bin-grouping}
[The same fit, two deviances]

Fit \( \operatorname{logit}\pi=\beta_0+\beta_1(\text{party})+\beta_2(\text{education}) \)
to the survey of @exm-bin-anes. The \( 944 \) respondents occupy
\( 47 \) distinct covariate cells. Grouped, the deviance is
\( 45.99 \) on \( 44 \) degrees of freedom; ungrouped, it is
\( 533.34 \) on \( 941 \). The coefficients and their standard errors
agree to machine precision.
:::

```{.python .run #cell-anes-logistic-grouped}
# the same model in grouped form: one binomial count per distinct covariate pattern
cols = ["PID", "educ"]
cells = anes.groupby(cols)["vote"].agg(["sum", "count"]).reset_index()
Xg = sm.add_constant(cells[cols])
counts = np.column_stack([cells["sum"], cells["count"] - cells["sum"]])
grouped = sm.GLM(counts, Xg, family=sm.families.Binomial()).fit()

Xu = sm.add_constant(anes[cols])
ungrouped = sm.GLM(anes["vote"], Xu, family=sm.families.Binomial()).fit()

print("grouped  ", grouped.params.round(4).to_numpy(), "deviance", round(grouped.deviance, 2))
print("ungrouped", ungrouped.params.round(4).to_numpy(), "deviance", round(ungrouped.deviance, 2))
```

## Estimation

There is no closed form for \( \hbeta \). Under the canonical link @thm-glm-score
reduces to
\[
\X\T(\y-\bmu)=\bzero,\qquad \mu_i=m_i\pi_i=\frac{m_i\exp(\x_{(i)}\T\bbeta)}{1+\exp(\x_{(i)}\T\bbeta)},
\]{#eq-bin-score}

a system of \( p \) smooth equations in \( p \) unknowns. Specializing
@thm-glm-irls gives the algorithm every implementation uses.

::: {#prp-bin-irls}
[Fisher scoring for logistic regression]

For @eq-bin-logistic the Fisher scoring step from \( \bb \) is the weighted least
squares fit of the working response \( \bz \) on \( \X \) with weights \( \W \), where
\[
w_i=m_i\pi_i(1-\pi_i),\qquad
z_i=\eta_i+\frac{y_i-m_i\pi_i}{m_i\pi_i(1-\pi_i)},
\]
both evaluated at \( \bb \). The observed and expected information coincide and equal
\( \X\T\W\X \); the log-likelihood is concave in \( \bbeta \), and strictly concave when
\( \rank(\X)=p \) and \( 0<\pi_i<1 \).
:::

::: {.proof}
The weights and working response are the general formulae of @thm-glm-irls with
\( V(\pi)=\pi(1-\pi) \) and \( d\eta/d\mu=1/\{m\pi(1-\pi)\} \). Differentiating
@eq-bin-score once more gives the Hessian \( -\X\T\W\X \), so the observed information
equals the expected information, as @thm-glm-score gives for every canonical link. The
matrix \( \X\T\W\X \) is nonnegative definite, and positive definite when the \( w_i \) are
positive and \( \X \) has full column rank, since \( \bb\T\X\T\W\X\bb=\sum_iw_i(\x_{(i)}\T\bb)^2 \)
vanishes only if \( \X\bb=\bzero \).
:::

Concavity means at most one local maximum, and the step is the weighted least squares
of @thm-het-wls with response and weights refreshed each iteration. Convergence from
\( \bbeta=\bzero \) is usually fast, and took \( 6 \) iterations
in @exm-bin-anes; it is not guaranteed, and [Section 35.3](03-separation.html)
describes the configuration in which there is nothing to converge to.

## Inference

With the estimate in hand, @thm-glm-asymptotics supplies its law: under the usual
regularity conditions, with \( p \) fixed and the information growing without bound,
\( \hbeta \) is treated as \( \Normal_p\bigl(\bbeta,(\X\T\hat{\W}\X)^{-1}\bigr) \).
The growth condition is real: it fails if a covariate is almost constant, and it fails
if the fitted probabilities are driven to \( 0 \) or \( 1 \), since then
\( w_i=m_i\pi_i(1-\pi_i)\to0 \). Of the three statistics of @prp-glm-three-tests, the Wald statistic is what software
prints, and it is the one the Hauck–Donner effect spoils; the likelihood ratio
statistic is invariant and survives an infinite estimate, and the score statistic
needs only the null fit.

::: {#exm-bin-anes-inference}
[Tests and intervals in the voting model]

In @exm-bin-anes, education has \( \hat\beta=0.0147 \) with standard error
\( 0.0771 \); its Wald statistic is \( 0.036 \) and so is the
likelihood ratio statistic from dropping it, while age has \( z=1.71 \)
and \( p=0.087 \). For the party coefficient the Wald interval is
\( (1.080,\,1.358) \) and the profile likelihood interval
\( (1.086,\,1.364) \): the log-likelihood is nearly quadratic here.
@exm-bin-separated shows what happens when it is not.
:::

The machinery of Part III survives in weakened form: simultaneous inference uses a
quadratic form in \( \bLambda\T\hbeta-\mathbf{d} \) with matrix
\( \{\bLambda\T(\X\T\hat{\W}\X)^{-1}\bLambda\}^{-1} \), referred to \( \chi^2(q) \) in place
of the \( F \) statistic of @thm-glh-general-f. What does *not* survive is exactness.
The right measure of how much there is to work with is the total weight
\( \sum_im_i\pi_i(1-\pi_i) \): a study of \( 500 \) subjects with twelve
events has total weight about \( 12 \), the same as a balanced study of fifty, which is
the arithmetic behind the events-per-variable rules of thumb in the applied
literature.

## Exercises

### A. Check your understanding

::: {#exr-bin-odds-arithmetic}
[A1]

A logistic model has \( \hat\beta=0.7 \) for a binary treatment indicator. Give the
estimated odds ratio, the divide-by-four bound on the change in probability, and the
actual change when the control probability is \( 0.5 \), \( 0.1 \) and \( 0.01 \).
:::

::: {.solution}
The odds ratio is \( e^{0.7}=2.014 \) and the bound is \( 0.7/4=0.175 \). Starting from
odds \( 1 \), \( 1/9 \) and \( 1/99 \), the treated odds are \( 2.014 \), \( 0.2238 \) and
\( 0.02034 \), so the probabilities become \( 0.668 \), \( 0.183 \) and \( 0.0199 \) and the
changes are \( 0.168 \), \( 0.083 \) and \( 0.0099 \). The bound is nearly attained at
\( \pi=1/2 \) and grossly overstates the change for rare events.
:::

::: {#exr-bin-canonical}
[A2]

Show from @eq-bin-score that a logistic model containing an indicator column for a
group reproduces that group's observed proportion exactly. What follows for a model
whose covariates are a full set of indicators for a categorical variable?
:::

### B. Practice

::: {#exr-bin-ld50}
[B1]

In \( \operatorname{logit}\pi=\beta_0+\beta_1x \) with \( \beta_1\ne0 \), the level at
which \( \pi=1/2 \) is \( x_{1/2}=-\beta_0/\beta_1 \) (in bioassay, the median effective
dose). Derive the delta-method standard error of \( \hat x_{1/2} \), and explain why it
becomes useless as \( \hat\beta_1 \) approaches zero.
:::

::: {#exr-bin-empirical-logit}
[B2]

For grouped data with all \( m_i \) large, @exr-tr-empirical-logit fits \( \bbeta \) by
weighted least squares on the empirical logits
\( \ell_i=\operatorname{logit}(s_i/m_i) \) with weights \( m_i\hat\pi_i(1-\hat\pi_i) \).
Show that one Fisher scoring step of @prp-bin-irls, started from the fit that
reproduces the observed proportions, is exactly that estimate.
:::

::: {.solution}
Start from a \( \bbeta \) with \( \pi_i=s_i/m_i \), so \( \eta_i=\ell_i \). The working
residual \( (y_i-m_i\pi_i)/\{m_i\pi_i(1-\pi_i)\} \) vanishes, the working response is
\( z_i=\ell_i \) and the weights are \( m_i\hat\pi_i(1-\hat\pi_i) \), so one step returns
the weighted least squares fit. The empirical-logit fit is therefore a legitimate
starting value; it is undefined when \( s_i\in\{0,m_i\} \).
:::

::: {#exr-bin-interaction}
[B3]

Two binary covariates enter with an interaction:
\( \operatorname{logit}\pi=\beta_0+\beta_1x_1+\beta_2x_2+\beta_{12}x_1x_2 \). Show that
\( \beta_{12} \) is the log of a ratio of two odds ratios, and that \( \beta_{12}=0 \) does
*not* make the effect of \( x_1 \) on the probability scale the same at both levels of
\( x_2 \).
:::

### C. Going deeper

::: {#exr-bin-concavity}
[C1]

Prove that \( \ell(\bbeta)=-\sum_i\log\{1+\exp(-s_i\x_{(i)}\T\bbeta)\} \), with
\( s_i=2y_i-1 \), is concave in \( \bbeta \), without computing the Hessian.
:::

::: {.solution}
The map \( u\mapsto\log(1+e^{u}) \) is convex, since its second derivative
\( e^{u}/(1+e^{u})^2 \) is positive. A convex function composed with an affine map is
convex, so each \( \bbeta\mapsto\log\{1+\exp(-s_i\x_{(i)}\T\bbeta)\} \) is convex, and a
sum of convex functions is convex. Hence \( \ell \) is concave.
:::
