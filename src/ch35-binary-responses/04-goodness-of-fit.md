# Goodness of fit

In the normal linear model the residual sum of squares both estimates \( \sigma^2 \)
and measures fit, and its known distribution makes a fit too poor to be chance
announce itself ([Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html)).
Binary data have no free dispersion parameter, so the deviance has only the second
role, and often cannot fill it either.

## Two statistics for grouped data

::: {#def-bin-fit-statistics}
[Pearson and deviance statistics]

For grouped binomial data \( y_i\sim\text{Binomial}(m_i,\pi_i) \), \( i=1,\dots,N \), with
fitted probabilities \( \hat\pi_i \) from a model with \( p \) parameters, the **Pearson
statistic** and the **deviance** are
\[
\begin{aligned}
X^2&=\sum_{i=1}^N\frac{(y_i-m_i\hat\pi_i)^2}{m_i\hat\pi_i(1-\hat\pi_i)},\\
D&=2\sum_{i=1}^N\Bigl\{y_i\log\frac{y_i}{m_i\hat\pi_i}
+(m_i-y_i)\log\frac{m_i-y_i}{m_i-m_i\hat\pi_i}\Bigr\},
\end{aligned}
\]{#eq-bin-fit-statistics}

with the convention \( 0\log0=0 \). Both compare the fit with the **saturated** model that
gives each group its own probability.
:::

\( X^2 \) is the sum over all \( 2N \) cells of successes and failures of
(observed \( - \) fitted)\( ^2 \)/fitted, and \( D \) is twice the sum of
observed \( \times\log \)(observed/fitted): the classical categorical-data
statistics, specialized to the binomial.

::: {#prp-bin-fit}
[What the two statistics measure]

::: {.enumerate options="label=(\alph*)"}
1. *(Grouped data.)* Fix \( N \) and the model, let \( m_i\to\infty \) for every \( i \) with
   \( \pi_i \) bounded away from \( 0 \) and \( 1 \), and suppose the model holds. Then
   \( X^2 \) and \( D \) both converge in distribution to \( \chi^2(N-p) \), and
   \( X^2-D\to0 \) in probability.

2. *(Ungrouped data, canonical link.)* If \( m_i=1 \) for every \( i \) and \( \hbeta \) solves
   the logistic likelihood equations @eq-bin-score, then
   \[
D=-2\sum_{i=1}^n\bigl\{\hat\pi_i\log\hat\pi_i+(1-\hat\pi_i)\log(1-\hat\pi_i)\bigr\}.
\]{#eq-bin-entropy-deviance}

   The deviance is therefore a function of the fitted probabilities alone: it does not
   depend on which observations were successes, it is bounded above by \( 2n\log2 \), and
   it has no \( \chi^2(n-p) \) limit.

3. *(Differences.)* For nested models \( \omega\subset\Omega \) with \( q \) parameters
   between them, \( D_\omega-D_\Omega \) does not depend on the grouping and converges to
   \( \chi^2(q) \) under the smaller model, grouped or not.
:::

:::

::: {.proof}
(a) The statement about \( D \) is @thm-glm-deviance specialized to the binomial family:
with \( N \) fixed and every \( m_i\to\infty \), the model is a smooth submodel of an
\( N \)-parameter family whose maximum likelihood estimates are asymptotically normal, so
the likelihood ratio statistic against the saturated model has \( N-p \) degrees of
freedom. For \( X^2 \), put \( b=m_i\hat\pi_i \) and
\( r_i=(y_i-b)/\{m_i\hat\pi_i(1-\hat\pi_i)\}^{1/2} \), so that \( r_i=O_p(1) \) and
\( y_i-b=O_p(m_i^{1/2}) \). The \( i \)th deviance summand is \( g(y_i) \) for
\( g(a)=2\{a\log(a/b)+(m_i-a)\log((m_i-a)/(m_i-b))\} \), and
\( g(b)=g'(b)=0 \) and \( g''(b)=2/\{m_i\hat\pi_i(1-\hat\pi_i)\} \), while the third
derivative is \( O(m_i^{-2}) \) near that point. Taylor's theorem therefore gives
\( g(y_i)=r_i^2+O_p(m_i^{-1/2}) \), using \( \hat\pi_i \) bounded away from \( 0 \) and
\( 1 \). Summing over the \( N \) fixed groups, \( D-X^2\to0 \) in probability, so \( X^2 \)
has the same limit as \( D \).

(b) Write \( \hat\eta_i=\x_{(i)}\T\hbeta \), so that
\( \hat\eta_i=\log\hat\pi_i-\log(1-\hat\pi_i) \). For ungrouped data the saturated model
sets \( \hat\pi_i=y_i \) and has log-likelihood zero, so
\[
\begin{aligned}
D&=-2\sum_i\bigl\{y_i\log\hat\pi_i+(1-y_i)\log(1-\hat\pi_i)\bigr\}\\
&=-2\sum_i\bigl\{y_i\hat\eta_i+\log(1-\hat\pi_i)\bigr\}\\
&=-2\Bigl\{\hbeta\T\X\T\y+\sum_i\log(1-\hat\pi_i)\Bigr\}.
\end{aligned}
\]
The likelihood equations say \( \X\T\y=\X\T\hat{\boldsymbol{\uppi}} \), so \( \X\T\y \) may be
replaced by \( \X\T\hat{\boldsymbol{\uppi}} \), giving
\( D=-2\sum_i\{\hat\pi_i\hat\eta_i+\log(1-\hat\pi_i)\} \), which is
@eq-bin-entropy-deviance after substituting for \( \hat\eta_i \). Each summand of
@eq-bin-entropy-deviance is at most \( 2\log2 \), attained at \( \hat\pi_i=1/2 \), which gives
the bound. For the last claim write \( H(\pi)=-\{\pi\log\pi+(1-\pi)\log(1-\pi)\} \) and
\( f_n(\bbeta)=2\sum_iH(\pi_i) \) with \( \pi_i=\pi_i(\bbeta) \), so that \( D=f_n(\hbeta) \)
identically: \( D \) is a smooth deterministic function of \( \hbeta \), with no other
dependence on the data. Now \( dH/d\eta=-\eta\,\pi(1-\pi) \), and
\( \pi(1-\pi)\le e^{-|\eta|} \), so \( |dH/d\eta|\le|\eta|e^{-|\eta|}\le e^{-1} \) and
\( \norm{\nabla f_n}\le2e^{-1}\sum_i\norm{\x_{(i)}}=O(n) \) for bounded covariates. With
\( \hbeta-\bbeta=O_p(n^{-1/2}) \) from @thm-glm-asymptotics, the mean value theorem gives
\( D-f_n(\bbeta)=O_p(n^{1/2}) \), so
\[
n^{-1}D-2n^{-1}\sum_{i=1}^nH(\pi_i)\to0\ \text{ in probability},
\]
whereas \( n^{-1}\chi^2(n-p)\to1 \). The two agree only by accident: the simulation of
@exm-bin-deviance-sim has \( n^{-1}D \) averaging \( 206.8/200=1.03 \) under one correct
model and \( 263.6/200=1.32 \) under another. So \( D \) has no \( \chi^2(n-p) \) limit.
The same bound leaves \( D \) fluctuating about \( f_n(\bbeta) \) by \( O_p(n^{1/2}) \),
which is the order of the standard deviations in @exm-bin-deviance-sim; what fails is
not the spread of \( D \) but its location.

(c) The saturated log-likelihood is common to both deviances and cancels, so the
difference is twice a log-likelihood ratio between two nested smooth models with \( q \)
free parameters between them, and @thm-glm-deviance applies with \( q \) fixed.
:::

::: {#exm-bin-grouped-fit}
[Is the party scale linear on the logit scale?]

Collapse the survey of @exm-bin-anes onto its seven party groups, of between
\( 37 \) and \( 200 \) respondents each, and fit
\( \operatorname{logit}\pi=\beta_0+\beta_1\,\text{party} \). The deviance is
\( 12.765 \) on \( 5 \) degrees of freedom
(\( p=0.0257 \)) and the Pearson statistic is
\( 11.899 \) (\( p=0.0362 \)). The
Pearson residuals locate the misfit: \( -2.067 \) at
"independent-Democrat", where the model predicts \( 0.132 \) against an
observed \( 0.065 \), and \( 2.122 \) at
"independent-Republican", where it predicts \( 0.640 \) against
\( 0.745 \). The transition across the middle of the scale is sharper
than a straight line in the log odds allows, and a quadratic term does not help: the
deviance falls only to \( 12.762 \), so the misfit is a step, not a
curve.
:::

```{.python .run #cell-goodness-of-fit-grouped}
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

anes = sm.datasets.anes96.load_pandas().data
cells = anes.groupby("PID")["vote"].agg(["sum", "count"]).reset_index()
counts = np.column_stack([cells["sum"], cells["count"] - cells["sum"]])
Xg = sm.add_constant(cells[["PID"]])

g = sm.GLM(counts, Xg, family=sm.families.Binomial()).fit()
print(cells.to_string(index=False))
print(f"deviance {g.deviance:.3f} on {g.df_resid} df, p = {stats.chi2.sf(g.deviance, g.df_resid):.4f}")
print(f"Pearson  {g.pearson_chi2:.3f} on {g.df_resid} df, "
      f"p = {stats.chi2.sf(g.pearson_chi2, g.df_resid):.4f}")
```

The conditions in part (a) are not decoration. If \( N \) grows with \( n \) instead,
as it does whenever a continuous covariate is in the model, the approximation fails,
worst of all in the extreme case \( m_i\equiv1 \) of part (b).

## The ungrouped deviance measures nothing

@eq-bin-entropy-deviance deserves to be read twice. The left-hand side looks like a
discrepancy between data and fit; the right-hand side contains no data at all. It is
twice the total Bernoulli entropy of the fitted probabilities, a measure of how
*undecided* the model is: probabilities near \( 1/2 \) give a deviance near
\( 2n\log2 \) and extreme ones give a small deviance, whether or not the model is
correct.

::: {#exm-bin-deviance-sim}
[Three deviances that mean nothing]

Simulate \( n=200 \) observations with one covariate and an intercept in three ways: a
correct model with a strong slope, a correct model with a weak slope, and a model that
leaves a genuine quadratic effect out of the fitted predictor. All three fit two
parameters, so \( n-p=198 \) and a chi-squared reference would have mean \( 198 \) and
standard deviation \( 19.9 \). Over \( 2000 \) replications:

| data-generating mechanism | mean deviance | standard deviation |
|---|---|---|
| correct, strong slope | 206.8 | 14.3 |
| correct, weak slope | 263.6 | 7.0 |
| quadratic omitted | 252.2 | 8.5 |

Two *correct* models with the same \( n \) and \( p \) differ by nearly sixty in mean
deviance, and the badly misspecified model has a *smaller* deviance than one of them:
no reference distribution can serve all three
([Figure 35.4.1](04-goodness-of-fit.html#fig-bin-fit)).

The survey model of @exm-bin-anes has ungrouped deviance \( 528.29 \)
on \( 939 \) degrees of freedom; dropping every covariate but education
gives an obviously wrong model with deviance \( 1275.16 \) on
\( 942 \). Neither number can be referred to a table.
:::

```{.python .run #cell-goodness-of-fit-entropy}
# the ungrouped deviance is a function of the fitted probabilities alone
X = pd.DataFrame({"const": 1.0, "PID": anes["PID"], "age": anes["age"] / 10,
                  "educ": anes["educ"], "income": anes["income"]})
y = anes["vote"].to_numpy()
fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
p = fit.fittedvalues.to_numpy()

entropy = -2 * np.sum(p * np.log(p) + (1 - p) * np.log(1 - p))
print("deviance reported  ", round(fit.deviance, 6))
print("entropy of the fit ", round(entropy, 6))
```

```{.python .run #cell-goodness-of-fit-devsim}
def deviance_null(n, B, seed=0):
    """Ungrouped deviances under three data-generating mechanisms: a correct model with a
    strong slope, a correct model with a weak slope, and a model that omits a real quadratic
    effect. All three fit an intercept and one slope to n observations."""
    rng = np.random.default_rng(seed)
    out = np.empty((B, 3))
    for b in range(B):
        x = rng.normal(size=n)
        Z = np.column_stack([np.ones(n), x])
        for j, eta in enumerate([0.4 + 1.5 * x, 0.4 + 0.3 * x, 0.4 + 1.5 * x - 1.8 * x ** 2]):
            yb = (rng.random(n) < 1 / (1 + np.exp(-eta))).astype(float)
            out[b, j] = sm.GLM(yb, Z, family=sm.families.Binomial()).fit().deviance
    return out


demo = deviance_null(200, 100, seed=1)
print("n = 200 and p = 2, so a chi-squared reference would have mean 198 and "
      f"standard deviation {np.sqrt(2 * 198):.1f}")
for j, label in enumerate(["correct, strong slope", "correct, weak slope  ",
                           "quadratic omitted    "]):
    print(f"{label}: mean {demo[:, j].mean():6.1f}, "
          f"standard deviation {demo[:, j].std():4.1f}")
```

::: {when-format="html"}
![**Figure 35.4.1.** (a) Simulated null distributions of the ungrouped deviance with
\( n=200 \) and two fitted parameters, under the three mechanisms of
@exm-bin-deviance-sim, against the \( \chi^2(198) \) density. (b) Cross-validated
calibration of the survey model in ten bins, with pointwise intervals. (c) ROC curves
for the four-regressor model and for party identification
alone.](fit_diagnostics.svg){#fig-bin-fit width=100%}
:::

::: {when-format="pdf"}
![(a) Simulated null distributions of the ungrouped deviance with
\( n=200 \) and two fitted parameters, under the three mechanisms of
@exm-bin-deviance-sim, against the \( \chi^2(198) \) density. (b) Cross-validated
calibration of the survey model in ten bins, with pointwise intervals. (c) ROC curves
for the four-regressor model and for party identification
alone.](fit_diagnostics.pdf){width=100%}
:::

The Pearson statistic escapes @eq-bin-entropy-deviance, since
\( X^2=\sum_i\{y_i/\hat\pi_i+(1-y_i)/(1-\hat\pi_i)\}-n \) does depend on \( \y \), but
gains nothing by it (@exr-bin-pearson-identity). What survives is part (c): for
ungrouped data, test the model against enlargements of it rather than against the
saturated alternative, and use AIC or cross-validation
([Chapter 29](../ch29-model-selection/index.html)) when the comparisons are not
nested.

## Grouping by fitted values

Hosmer and Lemeshow (1980) proposed making grouped data out of ungrouped: sort by
\( \hat\pi_i \), cut into \( G \) groups of roughly equal size, and compare
\[
\hat C=\sum_{k=1}^G\frac{(O_k-E_k)^2}{E_k(1-E_k/n_k)},\qquad
O_k=\sum_{i\in k}y_i,\quad E_k=\sum_{i\in k}\hat\pi_i ,
\]
referred to \( \chi^2(G-2) \). The reference distribution has no proof behind it: it
was established by simulation, and the two lost degrees of freedom are an empirical
finding. The statistic is a useful descriptive check and an unstable test.

::: {#exm-bin-hosmer}
[The same model, three verdicts]

For the four-regressor model of @exm-bin-anes:

| groups \( G \) | \( \hat C \) | \( p \) |
|---|---|---|
| 8 | 8.670 | 0.193 |
| 10 | 12.232 | 0.141 |
| 12 | 16.119 | 0.096 |

The \( p \)-value halves as \( G \) goes from \( 8 \) to \( 12 \), and different software makes
different choices about \( G \), about how to break ties in \( \hat\pi_i \), and about
whether to use equal-sized groups or equally spaced cut-points. None of these
choices is more principled than any other.
:::

```{.python .run #cell-goodness-of-fit-hl}
def hosmer_lemeshow(y, p, groups):
    """The Hosmer-Lemeshow statistic and its nominal p-value with the given number of groups."""
    order = np.argsort(p)
    chunks = np.array_split(order, groups)
    stat = 0.0
    for c in chunks:
        e = p[c].sum()
        stat += (y[c].sum() - e) ** 2 / (e * (1 - e / len(c)))
    return stat, stats.chi2.sf(stat, groups - 2)


for gr in (8, 10, 12):
    stat, pv = hosmer_lemeshow(y, p, gr)
    print(f"{gr:2d} groups: statistic {stat:6.3f}, p = {pv:.3f}")
```

## Calibration and discrimination

Two questions hide inside "does the model fit". **Calibration** asks whether the
fitted probabilities are right: among observations with \( \hat\pi_i\approx0.3 \), do
about thirty per cent have \( y_i=1 \)? **Discrimination** asks whether they order the
observations correctly, and a model can have either without the other. Classification
tables of sensitivity against specificity answer a third question, and depend on a
threshold and on a loss they rarely state; pseudo-\( R^2 \) measures such as
\( 1-D/D_{\text{null}} \) have no interpretation of their own, and by
@prp-bin-fit(b) the ungrouped version of that ratio is a ratio of entropies.

::: {#prp-bin-calibration}
[In-sample calibration is automatic]

Let \( \hbeta \) solve @eq-bin-score with \( \bone\in\C(\X) \), let
\( \hat\eta_i=\x_{(i)}\T\hbeta \), and suppose \( \hat\eta \) is not constant. Fit the
two-parameter logistic model \( \operatorname{logit}\pi_i=a+b\hat\eta_i \) to the same
responses. Then \( \hat a=0 \) and \( \hat b=1 \) exactly.
:::

::: {.proof}
The likelihood equations of the two-parameter model at \( (a,b)=(0,1) \) are
\( \sum_i(y_i-\hat\pi_i)=0 \) and \( \sum_i\hat\eta_i(y_i-\hat\pi_i)=0 \). The first holds
because \( \bone\in\C(\X) \) and \( \X\T(\y-\hat{\boldsymbol{\uppi}})=\bzero \); the second
because \( \hat{\boldsymbol{\upeta}}=\X\hbeta\in\C(\X) \), so
\( \hat{\boldsymbol{\upeta}}\T(\y-\hat{\boldsymbol{\uppi}})=\hbeta\T\X\T(\y-\hat{\boldsymbol{\uppi}})=0 \).
So \( (0,1) \) solves the two-parameter likelihood equations; it remains to see that a
maximizer exists. The design \( [\bone,\hat{\boldsymbol{\upeta}}] \) has rank two,
since \( \hat\eta \) is not constant. Suppose some \( (a,b)\ne(0,0) \) had
\( s_i(a+b\hat\eta_i)\ge0 \) for every \( i \). As \( \bone\in\C(\X) \) there is an
\( \mathbf{e} \) with \( \X\mathbf{e}=\bone \), and
\( \bb^*=a\mathbf{e}+b\hbeta \) satisfies
\( \X\bb^*=a\bone+b\hat{\boldsymbol{\upeta}}\ne\bzero \), so \( \bb^*\ne\bzero \)
separates \( \X \) — contradicting the overlap that the finiteness of \( \hbeta \)
guarantees (@thm-bin-separation). The two-parameter data are therefore in overlap, and
@thm-bin-separation(c) gives a unique maximizer, which must be the stationary point
\( (0,1) \).
:::

The **calibration slope** \( \hat b \) is therefore worthless on the data that produced
\( \hbeta \). It becomes informative when \( \hat\eta \) comes from a held-out fold or a
new sample, where \( \hat b<1 \) is the signature of overfitting. On the survey model of
@exm-bin-anes, ten-fold cross-validated predictions give a slope of
\( 0.983 \) and an intercept of \( 0.005 \): four
regressors and \( 944 \) observations is not a regime in which
overfitting bites.

::: {#prp-bin-auc}
[What the area under the ROC curve is, and is not]

Let \( \hat s_i \) be any scores, let \( I \) be a uniformly chosen index with \( y_I=1 \) and
\( J \) an independent uniformly chosen index with \( y_J=0 \). The area under the curve
traced by (false positive rate, true positive rate) as a threshold moves through the
scores equals
\[
\text{AUC}=\Pr(\hat s_I>\hat s_J)+\tfrac12\Pr(\hat s_I=\hat s_J),
\]
the normalized Mann–Whitney statistic. Consequently the AUC depends on the scores only
through their ranks, so it is unchanged by every strictly increasing transformation of
\( \hat s \), and in particular it says nothing whatever about calibration.
:::

::: {.proof}
Order the distinct score values downwards. Moving the threshold past one of them
advances the curve by the fraction \( \Delta\text{FPR} \) of controls and the fraction
\( \Delta\text{TPR} \) of cases at that value, adding the trapezoid
\( \Delta\text{FPR}\{\text{TPR}+\tfrac12\Delta\text{TPR}\} \). Summing, the first part
counts pairs with \( \hat s_I>\hat s_J \) and the second counts ties with weight
\( \tfrac12 \), each divided by \( n_1n_0 \). Invariance is immediate, since an
increasing transformation preserves every comparison \( \hat s_I>\hat s_J \).
:::

So a model whose probabilities are all wrong by a factor of ten has the same AUC as
the correctly calibrated model with the same ordering. The AUC is insensitive in the
other direction too: for @exm-bin-anes party identification alone gives
\( \text{AUC}=0.942 \), and adding age, education and income raises it
to \( 0.944 \). An increment of \( 0.002 \) in a rank-based
summary is the wrong scale on which to judge the three covariates.

## Residuals and influence

The diagnostics of [Chapter 20](../ch20-residuals-leverage-influence/index.html) carry
over, with the working weights as the metric.

::: {#def-bin-residuals}
[Residuals for binomial data]

With \( \hat w_i=m_i\hat\pi_i(1-\hat\pi_i) \), \( \hat{\W}=\diag(\hat w_i) \) and
\( \hat h_i \) the \( i \)th diagonal entry of
\( \hat{\W}^{1/2}\X(\X\T\hat{\W}\X)^{-1}\X\T\hat{\W}^{1/2} \), the **Pearson** and
**deviance** residuals are
\[
r_i=\frac{y_i-m_i\hat\pi_i}{\sqrt{\hat w_i}},\qquad
d_i=\operatorname{sign}(y_i-m_i\hat\pi_i)\sqrt{D_i},
\]
where \( D_i \) is the \( i \)th summand of the deviance in @eq-bin-fit-statistics. Their
**standardized** forms divide by \( (1-\hat h_i)^{1/2} \).
:::

By construction \( \sum_ir_i^2=X^2 \) and \( \sum_id_i^2=D \). The matrix in the
definition is the hat matrix of the last weighted least squares step of
@prp-bin-irls, so it is symmetric and idempotent with trace \( \rank(\X) \), and
[Section 6.8](../ch06-projections/08-leverage.html) applies to it unchanged. The
analogue of Cook's distance (@def-res-cooks) is
\( r_i^2\hat h_i/\{p(1-\hat h_i)^2\} \).

Two things change. For ungrouped data \( r_i \) takes only two values once
\( \hat\pi_i \) is given, so plots against a covariate show two curves rather than a
cloud and must be smoothed or binned. And \( \hat h_i \) carries the factor
\( \hat w_i \), so high leverage now means a covariate vector that is extreme *and* a
fitted probability near \( 1/2 \); a point whose outcome the model regards as certain
has low leverage wherever its covariates lie, the opposite of the linear-model
intuition.

## Exercises

### A. Check your understanding

::: {#exr-bin-deviance-bound}
[A1]

Use @eq-bin-entropy-deviance to find the largest and smallest possible values of the
ungrouped deviance of a logistic fit, and describe the fits that attain them.
:::

::: {.solution}
Each summand \( -2\{\hat\pi\log\hat\pi+(1-\hat\pi)\log(1-\hat\pi)\} \) lies in
\( [0,2\log2] \), so \( 0\le D\le2n\log2 \). The upper bound is attained when every
\( \hat\pi_i=1/2 \), as for the intercept-only model fitted to equal numbers of
successes and failures. The lower bound is approached only as the fitted probabilities
go to \( 0 \) and \( 1 \), which by @thm-bin-separation needs complete separation, so
\( 0 \) is an infimum and not a minimum.
:::

### B. Practice

::: {#exr-bin-pearson-identity}
[B1]

Show that for ungrouped data
\( X^2=\sum_i\{y_i/\hat\pi_i+(1-y_i)/(1-\hat\pi_i)\}-n \), and deduce that for the
intercept-only model \( X^2=n \) exactly, whatever the data. What does that say about
\( X^2 \) as a test of fit?
:::

::: {.solution}
With \( m_i=1 \), \( (y_i-\hat\pi_i)^2/\{\hat\pi_i(1-\hat\pi_i)\} \) equals
\( (1-\hat\pi_i)/\hat\pi_i \) when \( y_i=1 \) and \( \hat\pi_i/(1-\hat\pi_i) \) when
\( y_i=0 \), and both equal \( y_i/\hat\pi_i+(1-y_i)/(1-\hat\pi_i)-1 \); summing gives the
identity. For the intercept-only model \( \hat\pi_i=\bar y \) for every \( i \), so
\( X^2=n\bar y/\bar y+n(1-\bar y)/(1-\bar y)-n=n \). A statistic that equals the sample
size for every data set cannot test anything.
:::

::: {#exr-bin-calibration-out-of-sample}
[B2]

Let \( \hat\eta \) be fitted on one sample and evaluated on an independent sample from
the same population. Show that the calibration slope \( b \) has probability limit
\( 1 \) if the model is correct and \( \hbeta \) is consistent, and explain why a value
below one indicates overfitting.
:::

::: {#exr-bin-auc-limits}
[B3]

Suppose the covariate is \( \Normal(\mu_1,1) \) among cases and \( \Normal(\mu_0,1) \)
among controls. Compute the AUC of the correctly specified logistic model in terms of
\( \Delta=\mu_1-\mu_0 \), and find the \( \Delta \) needed for an AUC of \( 0.9 \).
:::

::: {.solution}
The fitted score is monotone in the covariate, so the AUC is
\( \Pr(X_1>X_0)=\Phi(\Delta/\sqrt2) \), since \( X_1-X_0\sim\Normal(\Delta,2) \). Setting
this to \( 0.9 \) gives \( \Delta=\sqrt2\,\Phi^{-1}(0.9)=1.812 \). An AUC of \( 0.65 \)
needs only \( \Delta\approx0.55 \), a difference highly significant in a large sample
yet nearly useless for classifying an individual.
:::


