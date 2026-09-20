# Offsets and rates

A count on its own is rarely interesting. Twelve accidents mean one thing on a
road carrying a thousand vehicles a day and another on one carrying a hundred
thousand. What is comparable is the **rate**: the count divided by the exposure
that produced it. The Poisson log-linear model handles rates without new
machinery, by moving a known quantity into the linear predictor and refusing to
give it a coefficient.

## Exposure

Let \( t_i>0 \) be the **exposure** of observation \( i \): person-years of
follow-up, vehicle-kilometres, population size, area surveyed. The natural
assumption is that the expected count is proportional to exposure,
\[
\mu_i=t_i\lambda_i ,
\]
where \( \lambda_i \) is the expected number of events per unit of exposure, the
quantity of interest. Modelling \( \lambda_i \) log-linearly,
\( \log\lambda_i=\x_{(i)}\T\bbeta \), gives
\[
\log\mu_i=\log t_i+\x_{(i)}\T\bbeta ,
\qquad
\mu_i=t_i\exp\bigl(\x_{(i)}\T\bbeta\bigr).
\]{#eq-cnt-offset}

The term \( \log t_i \) is an **offset**: a column of the linear predictor whose
coefficient is known to be one rather than estimated. Coefficients are now log
rate ratios, \( e^{\beta_j} \) multiplying the rate per unit of exposure rather
than the raw count. This is *not* a regression of the observed rates
\( y_i/t_i \), which would throw away the information that a rate based on a
thousand person-years is better determined than one based on five: the response
stays the count and the exposure enters as a known constant.

::: {#prp-cnt-offsets}
[Offsets]

Consider the model @eq-cnt-offset with \( \rank(\X)=p \), and write
\( \bmu=\bmu(\bbeta) \) for the mean vector it defines. Then:

::: {.enumerate options="label=(\alph*)"}
1. the log-likelihood is
   \( \ell(\bbeta)=\sum_i\{y_i(\log t_i+\x_{(i)}\T\bbeta)-t_ie^{\x_{(i)}\T\bbeta}\}-\sum_i\log y_i! \);
   its score is \( \X\T(\y-\bmu) \) and its information is \( \X\T\W\X \) with
   \( \W=\diag(\mu_1,\dots,\mu_n) \), so @thm-cnt-score holds verbatim with the
   offset in place;

2. consequently @thm-cnt-marginals holds verbatim: the fitted counts, not the
   fitted rates, reproduce the observed margins of the model matrix, and with an
   intercept \( \sum_i\hat\mu_i=\sum_iy_i \);

3. **(grouping)** suppose the rows of \( \X \) take only \( G \) distinct values
   \( \x_{(g)} \), with \( n_g \) rows equal to the \( g \)th, all with exposure
   \( t_i=1 \). Let \( s_g \) be the total count of group \( g \). Then fitting
   @eq-cnt-model to the \( n \) individual records and fitting @eq-cnt-offset to
   the \( G \) totals \( s_g \) with exposures \( t_g=n_g \) give the same
   \( \hbeta \), the same information matrix and the same deviance differences
   between nested models;

4. **(a free exposure coefficient)** the model
   \( \log\mu_i=\psi\log t_i+\x_{(i)}\T\bbeta \) contains @eq-cnt-offset as the
   submodel \( \psi=1 \). Its fit coincides with the offset fit if and only if
   \( \hat\psi=1 \), and the likelihood ratio statistic for \( H_0:\psi=1 \) is
   asymptotically \( \chi^2(1) \).
:::

:::

::: {.proof}
(a) Substituting \( \mu_i=t_ie^{\x_{(i)}\T\bbeta} \) into the Poisson
log-likelihood gives the stated expression, and
\( \partial\mu_i/\partial\bbeta=\mu_i\x_{(i)} \) exactly as before, because the
factor \( t_i \) is a constant. The derivative of \( \sum_iy_i\log t_i \) in
\( \bbeta \) is zero. Everything in the proof of @thm-cnt-score then goes through
with the new \( \mu_i \), which is (a); and (b) follows from (a) as before.

(c) Within group \( g \) all the means equal \( \mu_g=e^{\x_{(g)}\T\bbeta} \), so
\[
\ell_{\text{ind}}(\bbeta)
=\sum_{g=1}^{G}\bigl\{s_g\log\mu_g-n_g\mu_g\bigr\}-\sum_i\log y_i! .
\]
The grouped model has means \( n_g\mu_g \), so its log-likelihood is
\[
\ell_{\text{grp}}(\bbeta)
=\sum_{g=1}^{G}\bigl\{s_g(\log n_g+\log\mu_g)-n_g\mu_g\bigr\}-\sum_g\log s_g! .
\]
The two differ by a quantity free of \( \bbeta \). Log-likelihoods differing by a
constant have the same maximizer and the same derivatives of every order, hence
the same information; and differences of deviances are differences of
log-likelihoods, so the constant cancels there too.

(d) Setting \( \psi=1 \) recovers @eq-cnt-offset, so the models are nested with one
extra parameter. The larger model's log-likelihood is strictly concave in
\( (\bbeta\T,\psi)\T \) whenever \( [\X,\log\mathbf{t}] \) has full column rank, by
@thm-cnt-score(c) applied to that model matrix, so its maximizer is unique; it
therefore agrees with the offset fit exactly when \( \hat\psi=1 \). The
distribution statement is @thm-glm-deviance for one degree of freedom.
:::

Part (c) says that grouping records with identical regressors and carrying the
group size as an exposure is *exactly* a bookkeeping operation. It is a property
of the Poisson, not of generalized linear models in general, and it comes from
sums of independent Poissons being Poisson. The use is both practical — twenty
thousand records collapse to a few dozen — and theoretical: the grouped cells may
have large counts, so the deviance of a single fitted model is calibrated by
@thm-glm-deviance even though it is not for the individual records.

::: {#exm-cnt-grouping}
[Twenty thousand records in thirty-two rows]

Take the physician-visit data of @exm-cnt-visits and keep only categorical
regressors: the individual-deductible indicator, the three self-rated health
indicators, and the chronic-disease score banded into four groups. There are
\( 32 \) distinct covariate patterns, holding between
\( 9 \) and \( 3832 \) person-years. Collapse the
\( 20190 \) records to these rows, recording the total visits and the
number of person-years in each, and fit the rate model with offset the log of the
person-years. The coefficients agree with the person-level fit to fourteen decimal
places, as @prp-cnt-offsets(c) promises.

Read as rates, a person-year in excellent health, in the lowest disease band and
on a plan without an individual deductible, is expected to produce
\( 1.9758 \) visits; an individual deductible multiplies that by
\( 0.8247 \), poor self-rated health by
\( 1.6875 \), and the highest disease band by
\( 2.2621 \).

The grouped deviance is \( 665.60 \) on \( 24 \)
degrees of freedom, and now the \( \chi^2 \) calibration applies, since the
smallest fitted cell mean is \( 24.7 \). The model is rejected.
Where does the misfit sit? The fitted model has main effects only, so its
\( 24 \) residual degrees of freedom are the interactions among the
three factors; adding every two-factor interaction among them,
\( 23 \) parameters in all, drops the deviance to
\( 106.87 \) on \( 9 \) degrees of freedom, a fall
of \( 558.7 \) on \( 15 \). Five sixths of the
misfit is therefore mean model, not variance function. What survives the repair —
\( 106.87 \) on nine degrees of freedom, still far beyond any
plausible \( \chi^2 \) — is the Poisson variance function, which
[Section 37.4](04-negative-binomial.html) takes up.
:::

```{.python .run #cell-offsets-grouping}
import numpy as np
import pandas as pd
import statsmodels.api as sm


def scoring(X, y, offset=0.0, tol=1e-11, maxit=60):
    """Fisher scoring for a Poisson log-linear model with a fixed offset."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 1e-6))
    for _ in range(maxit):
        mu = np.exp(X @ beta + offset)
        z = X @ beta + (y - mu) / mu                 # working response
        XW = X * mu[:, None]                         # weights w_i = mu_i
        step = np.linalg.solve(X.T @ XW, XW.T @ z)
        if np.max(np.abs(step - beta)) < tol:
            return step
        beta = step
    return beta


data = sm.datasets.randhie.load_pandas().data.copy()
data["band"] = pd.cut(data["disea"], [-0.1, 6, 10, 15, 60], labels=[0, 1, 2, 3])
cells = ["idp", "hlthg", "hlthf", "hlthp", "band"]   # all categorical
y = data["mdvis"].to_numpy(float)


def design(frame):
    """Intercept, the plan and health indicators, and three disease-band indicators."""
    columns = [np.ones(len(frame))] + [frame[v].to_numpy(float)
                                       for v in ["idp", "hlthg", "hlthf", "hlthp"]]
    band = frame["band"].to_numpy(int)
    return np.column_stack(columns + [(band == k).astype(float) for k in (1, 2, 3)])


grouped = data.groupby(cells, as_index=False, observed=True).agg(
    visits=("mdvis", "sum"), person_years=("mdvis", "size"))
s = grouped["visits"].to_numpy(float)                # total visits in the group
t = grouped["person_years"].to_numpy(float)          # exposure: person-years
X, Xg = design(data), design(grouped)

beta_person = scoring(X, y)                          # one row per person-year
beta_group = scoring(Xg, s, offset=np.log(t))        # one row per pattern, offset log t

print("groups:", len(s), " largest coefficient difference:",
      np.abs(beta_person - beta_group).max())
print("fitted rate, free-plan person-years in excellent health:",
      np.exp(beta_group[0]))
```

## What a free exposure coefficient buys and costs

Why fix \( \psi=1 \) instead of estimating it, as @prp-cnt-offsets(d) allows? It
depends on where the exposure came from.

When exposure is a bookkeeping device — the grouping of part (c), where the
individual records are the real data — proportionality is a theorem, not an
assumption, and a \( \hat\psi \) far from one is a symptom of misspecification
elsewhere, because \( \log n_g \) is free to proxy for anything correlated with
group size. In @exm-cnt-grouping it is \( \hat\psi=1.2266 \),
standard error \( 0.0146 \), interval
\( (1.1980,\,1.2553) \), with likelihood ratio statistic
\( 249.019 \) against \( \psi=1 \). A person-year cannot
contribute more than proportionally to its own visits; the statistic is detecting
the same lack of fit as the deviance of \( 665.60 \).

When exposure genuinely varies across units — different follow-up times, different
populations — proportionality is a substantive assumption, and there is a real
trade-off, which a simulation makes concrete.

::: {#exm-cnt-exposure}
[Three ways to handle an exposure]

Counts are generated for \( m=300 \) units from
\( \log\mu_i=0.4+0.5x_i+\psi\log t_i \), with \( x_i \) standard normal and the
exposure \( t_i \) correlated with it. Three analyses are run on each of
\( 2000 \) data sets — the correct offset, a free coefficient on
\( \log t \), and no exposure term — and [Figure 37.2.1](#fig-cnt-offset) shows
the sampling distributions of the estimated slope of \( x \), whose true value is
\( 0.5 \).

When the truth is \( \psi=1 \), both the offset and the free-coefficient model are
unbiased, with means \( 0.500 \) and
\( 0.500 \), but the free coefficient has standard deviation
\( 0.029 \) against \( 0.017 \), a loss of
precision by a factor of \( 1.772 \): estimating a parameter you
already know costs information. Dropping the exposure is a disaster, the slope
averaging \( 1.096 \), because \( x \) is left to explain the
variation in exposure as well as in the rate.

When the truth is \( \psi=0.7 \), the free-coefficient model is still centred at
\( 0.501 \), while the offset model, forcing a wrong
proportionality, is centred at \( 0.322 \) — a bias of more than
six of its own standard deviations. Fixing a coefficient you have got wrong costs
more than estimating one you knew.
:::

::: {when-format="html"}
![**Figure 37.2.1.** Sampling distributions of the estimated slope of \( x \)
over \( 2000 \) simulated data sets, under three treatments of the
exposure, when the true exposure exponent is (a) one and (b) \( 0.7 \). The
dashed line marks the true slope.](offset_simulation.svg){#fig-cnt-offset width=100%}
:::

::: {when-format="pdf"}
![Sampling distributions of the estimated slope of \( x \)
over \( 2000 \) simulated data sets, under three treatments of the
exposure, when the true exposure exponent is (a) one and (b) \( 0.7 \). The
dashed line marks the true slope.](offset_simulation.pdf){width=100%}
:::

::: {.warning}
The third analysis fails by ordinary omitted variable bias: \( \log t \) belongs
in the linear predictor, and leaving it out loads its effect onto whatever is
correlated with it. What is special about an offset is only that its coefficient
is known, so it is not a variable to be selected and must never be dropped by a
stepwise procedure or a lasso.
:::

## Rates elsewhere

Two further uses of the same device are worth naming.

**Piecewise-constant hazards.** Split each subject's follow-up into intervals on
which the hazard is constant, and let \( d_{ik} \) be subject \( i \)'s events in
interval \( k \) and \( t_{ik} \) the time at risk. Treating the \( d_{ik} \) as
independent Poisson counts with means
\( t_{ik}\exp(\alpha_k+\x_{(i)}\T\bbeta) \) gives a proportional-hazards model with
a step-function baseline, fitted as an ordinary Poisson regression with offset
\( \log t_{ik} \) and one parameter per interval (@exr-cnt-piecewise); the
connection is due to Holford and to Laird and Olivier.

**Standardization.** A standardized morbidity ratio divides the observed count by
an expected count \( E_i \) from external reference rates; taking \( \log E_i \)
as the offset makes \( \exp(\x_{(i)}\T\bbeta) \) the ratio itself.

Fixing the exposure coefficient is right when theory or bookkeeping says the count
is proportional to the exposure, and costs nothing; estimating it is right when
proportionality is in doubt, and costs precision. What is never right is leaving
the exposure out.

## Exercises

### A. Check your understanding

::: {#exr-cnt-offset-units}
[A1]

A rate model uses person-years as the exposure. A colleague re-expresses exposure
in person-months, multiplying every \( t_i \) by twelve. Which coefficients change
and by how much? Which fitted values change?
:::

::: {#exr-cnt-offset-vs-rate}
[A2]

Explain why regressing \( \log(y_i/t_i) \) on \( \x_{(i)} \) by least squares is not
equivalent to fitting @eq-cnt-offset, giving two distinct reasons.
:::

### B. Practice

::: {#exr-cnt-piecewise}
[B1]

Suppose subject \( i \) is followed until an event or censoring, with a hazard that
is constant at \( h_{ik}=\exp(\alpha_k+\x_{(i)}\T\bbeta) \) on interval \( k \). Write
\( d_{ik}\in\{0,1\} \) for the event indicator and \( t_{ik} \) for the time at risk
in that interval. Show that the likelihood contribution of subject \( i \) is
\( \prod_k h_{ik}^{d_{ik}}e^{-h_{ik}t_{ik}} \), and hence that maximizing it is the
same as fitting a Poisson log-linear model to the \( d_{ik} \) with offset
\( \log t_{ik} \).
:::

::: {.solution}
With a constant hazard \( h \) on an interval, the probability of surviving a time
\( t \) is \( e^{-ht} \) and the density of an event at the end is \( he^{-ht} \), so
the contribution of interval \( k \) is \( h_{ik}^{d_{ik}}e^{-h_{ik}t_{ik}} \), and
the contributions multiply. Its logarithm is
\( \sum_{i,k}\{d_{ik}\log h_{ik}-h_{ik}t_{ik}\} \). The Poisson log-likelihood for
counts \( d_{ik} \) with means \( t_{ik}h_{ik} \) differs from it by
\( \sum_{i,k}d_{ik}\log t_{ik} \), a constant in the parameters, so the maximizers
agree. The \( d_{ik} \) need not actually be Poisson; only the likelihoods
coincide.
:::

::: {#exr-cnt-offset-power}
[B2]

In the setting of @exm-cnt-exposure, suppose the analyst tests \( H_0:\psi=1 \) and
uses the offset model if the test does not reject. Explain why the resulting
estimator of the slope of \( x \) has a sampling distribution that is neither of the
two shown in [Figure 37.2.1](#fig-cnt-offset), and relate the phenomenon to
the pre-test problem @exm-sel-pretest-coverage.
:::

### C. Going deeper

::: {#exr-cnt-grouping-general}
[C1]

Extend @prp-cnt-offsets(c) to the case where the individual records already carry
exposures \( t_i \) that differ within a group. Show that the grouped fit with
exposure \( \sum_{i\in g}t_i \) reproduces the individual fit, and identify exactly
what has to be constant within a group.
:::

::: {.solution}
What must be constant within a group is the row \( \x_{(i)} \), not the exposure.
If \( \x_{(i)}=\x_{(g)} \) for all \( i\in g \) then \( \mu_i=t_i\mu_g^{*} \) with
\( \mu_g^{*}=e^{\x_{(g)}\T\bbeta} \), so
\( \sum_{i\in g}\mu_i=\bigl(\sum_{i\in g}t_i\bigr)\mu_g^{*} \) and the
log-likelihood is
\( \sum_g\{s_g\log\mu_g^{*}-(\sum_{i\in g}t_i)\mu_g^{*}\}+\sum_iy_i\log t_i \),
whose last term is free of \( \bbeta \). Sums of independent Poissons are Poisson,
so the grouped data have exactly the claimed distribution.
:::

::: {#exr-cnt-offset-misspecified}
[C2]

Let the truth be \( \log\mu_i=\psi\log t_i+\beta_0+\beta_1x_i \) with
\( \psi\ne1 \), and suppose the offset model is fitted. Show that the limiting
value of \( \hat\beta_1 \) solves
\( \E\{\x\,t\,(t^{\psi-1}e^{\beta_0+\beta_1x}-e^{\beta_0^{*}+\beta_1^{*}x})\}=\bzero \),
and argue that \( \beta_1^{*}=\beta_1 \) when \( t^{\psi-1} \) is independent of
\( x \) under the weighting by \( t \).
:::
