# Transforming the response

The normal linear model claims that the mean of the response is linear in the regressors, that its
variance is constant, and that its distribution is normal. [Chapter 19](../ch19-theory-of-departures/index.html)
worked out what each failure costs, and [Chapter 20](../ch20-residuals-leverage-influence/index.html) and
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) showed how to see them. Replacing \( y \)
by \( g(y) \), for a strictly increasing \( g \), acts on all three claims at once. If the failures have a common
cause, one \( g \) removes them together; if not, the \( g \) that straightens the mean may worsen the variance.

## Engel's data on two scales

::: {#exm-tr-engel-scales}
[Food expenditure on the raw and the log scale]

Engel's survey of \( 235 \) Belgian households (@exm-proj-engel) records annual income, from
\( 377 \) to \( 4958 \) in the survey's currency, and expenditure on food. The least squares line on
the raw scale has slope \( 0.4852 \), and panel (b) of [Figure 22.1.1](#fig-tr-engel-scales) shows a fan: the
residual standard deviation among the richest third of households is \( 3.55 \) times that among the
poorest third. The residuals have excess kurtosis \( 9.64 \) (zero for a normal distribution), largely a
by-product of the fan, since a mixture of normals with different spreads has heavy tails.

On the log–log scale both defects nearly disappear. The spread ratio falls to \( 1.32 \), the excess kurtosis to \( -0.02 \), and the
residual plot, panel (d), shows a horizontal band. The slope is \( 0.8559 \) with standard error
\( 0.0203 \) and \( 95\% \) interval \( [0.816,\ 0.896] \). It is an elasticity
([Section 5.7](../ch05-model-and-least-squares/07-interpreting-coefficients.html)): households with
\( 1 \) percent more income spend about \( 0.86 \) percent more on food. That the interval lies below one is
Engel's law, that the share of the budget spent on food falls as income rises.
:::

::: {when-format="html"}
![**Figure 22.1.1.** Engel's data. (a) Food expenditure against income with the least squares line.
(b) Its residuals against fitted values: the spread grows with the level. (c) The same data on the
log–log scale. (d) The log–log residuals form a band of nearly constant width.](engel_scales.svg){#fig-tr-engel-scales width=100%}
:::

::: {when-format="pdf"}
![Engel's data. (a) Food expenditure against income with the least squares line.
(b) Its residuals against fitted values: the spread grows with the level. (c) The same data on the
log–log scale. (d) The log–log residuals form a band of nearly constant width.](engel_scales.pdf){width=100%}
:::

```{.python .run #cell-engel-scale-fits}
import numpy as np
import statsmodels.api as sm
from scipy import stats

df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)

raw = sm.OLS(food, sm.add_constant(income)).fit()                   # food on income
loglog = sm.OLS(np.log(food), sm.add_constant(np.log(income))).fit() # log food on log income

def spread_ratio(fit):
    """Residual SD among the richest third over that among the poorest third."""
    order = np.argsort(income)
    low, high = order[: n // 3], order[-(n // 3):]
    return fit.resid[high].std(ddof=1) / fit.resid[low].std(ddof=1)

for name, fit in [("raw", raw), ("log-log", loglog)]:
    print(f"{name:8s} slope {fit.params[1]:.4f}  spread ratio {spread_ratio(fit):.2f}"
          f"  excess kurtosis of residuals {stats.kurtosis(fit.resid):.2f}")
```

## Multiplicative models

The logarithm worked because the variation in food expenditure behaves like a *percentage* of its level.
A model in which the regressors and the error act multiplicatively makes this precise.

::: {#prp-tr-multiplicative}
[Multiplicative errors]

Let the model matrix \( \X \) have first column \( \bone \), with \( \beta_0 \) the corresponding coefficient,
and let
\[
Y_i=\exp(\x_{(i)}\T\bbeta)\,U_i,\qquad i=1,\dots,n,
\]{#eq-tr-multiplicative}

where \( U_1,\dots,U_n \) are independent and identically distributed positive random variables with
\( \nu=\E(\log U_i) \) and \( \tau^2=\Var(\log U_i) \) finite.

::: {.enumerate options="label=(\alph*)"}
1. The logarithms follow the linear model with constant variance:
   \( \log Y_i=(\beta_0+\nu)+\sum_{j\ge1}x_{ij}\beta_j+\varepsilon_i \), where
   \( \varepsilon_i=\log U_i-\nu \) are independent with mean \( 0 \) and variance \( \tau^2 \). The slopes are
   \( \beta_1,\dots,\beta_{p-1} \), unchanged; only the intercept absorbs \( \nu \).

2. If \( \E(U_i^2)<\infty \), then \( \E(Y_i)=\exp(\x_{(i)}\T\bbeta)\E(U_1) \) and the standard deviation
   of \( Y_i \) is \( \exp(\x_{(i)}\T\bbeta) \) times that of \( U_1 \). The standard deviation is proportional
   to the mean, so the coefficient of variation is the same for every case.

3. For two regressor profiles \( \x \) and \( \x' \), every quantile of \( Y \) at \( \x' \) is \( \exp\{(\x'-\x)\T\bbeta\} \)
   times the same quantile at \( \x \). The same ratio holds for the geometric means
   \( \exp\E(\log Y) \), and for the means when they are finite.
:::

:::

::: {.proof}
(a) Take logarithms in @eq-tr-multiplicative and add and subtract \( \nu \). The \( \varepsilon_i \) are functions of
the independent \( U_i \), hence independent, and have the stated moments.
(b) \( \exp(\x_{(i)}\T\bbeta) \) is a constant, so it factors out of the mean and the standard deviation.
(c) If \( u_q \) is a \( q \)-quantile of \( U_1 \), then \( \exp(\x\T\bbeta)u_q \) is a \( q \)-quantile of
\( Y \) at \( \x \), because multiplication by a positive constant preserves order. The geometric mean at \( \x \) is
\( \exp(\x\T\bbeta+\nu) \). In each case the ratio of the value at \( \x' \) to the value at \( \x \) is
\( \exp\{(\x'-\x)\T\bbeta\} \).
:::

Part (c) makes log-scale coefficients easy to report: \( (\x'-\x)\T\bbeta \) is the log of the *same* ratio
for medians, upper quartiles, geometric means and means. Part (b) is the diagnostic signature, a raw-scale
spread proportional to the fitted value, as in panel (b) of [Figure 22.1.1](#fig-tr-engel-scales).

::: {.warning}
The proposition needs the error to *multiply* the mean. If instead
\( Y_i=\exp(\x_{(i)}\T\bbeta)+\varepsilon_i \) with errors of constant variance, then
\( \log Y_i \) is not a linear function of the regressors plus an error of mean zero. For large means
\( \log Y_i\approx\x_{(i)}\T\bbeta+\varepsilon_i\exp(-\x_{(i)}\T\bbeta) \), so the logarithm *creates* a variance that
falls with the level (@exr-tr-additive-log). A mean that is exponential in the regressors combined with errors
of constant variance calls for fitting the mean directly, by nonlinear least squares or by a
generalized linear model with a logarithmic link (Chapter 34), not for transforming the data.
:::

## The ladder of powers

The logarithm belongs to Tukey's (1977) *ladder* of powers \( y^{p} \), \( p=\dots,-1,-\tfrac12,0,\tfrac12,1,2,\dots \),
with \( \log y \) in the place of \( y^0 \) and \( -y^{p} \) for \( p<0 \) so that every rung is increasing. Rungs below
\( p=1 \) are concave: they pull in a long right tail, reduce a spread that grows with the level, and straighten
a relationship that curves upwards. Rungs above \( p=1 \) do the opposite. [Section 22.2](02-box-cox.html) makes
the ordering precise (@prp-tr-family) and estimates the rung, and
[Section 22.4](04-variance-stabilizing.html) proves the rule linking rung and variance: a standard deviation
proportional to \( \mu^{a} \) calls for the power \( 1-a \), the logarithm when \( a=1 \) as in @prp-tr-multiplicative(b).

## Three goals and one function

The three claims are not equally important. With non-normal errors, least squares and the usual tests remain
asymptotically valid under mild conditions (@thm-dep-nonnormal). A non-constant variance leaves the
coefficients unbiased but the standard errors wrong (@thm-dep-covariance). A wrong mean function biases
everything. So the mean comes first, then the variance, then the shape of the errors.

Counts show the goals pulling apart. Suppose \( Y_i \) is Poisson with mean
\( \mu_i=\exp(\x_{(i)}\T\bbeta) \). The logarithm linearizes the mean exactly, but the variance of \( \log Y_i \)
is, to first order in the delta method, about \( 1/\mu_i \) (@exr-tr-poisson-scales), which changes with the regressors. The square root
stabilizes the variance at about \( \tfrac14 \) ([Section 22.4](04-variance-stabilizing.html)), but its mean is
about \( \exp(\x_{(i)}\T\bbeta/2) \), which is not linear in the regressors. (Strictly, zeros have probability \( e^{-\mu_i}>0 \), so the statement concerns \( \log Y_i \) given
\( Y_i>0 \), or \( \log(Y_i+c) \).) No power does both. The remedy is a model that states the mean and the variance separately, the Poisson
generalized linear model of Chapter 34. A transformation can also remove an interaction, under the order
condition of @prp-tw-removable, as the logarithm did for the running times of
[Section 16.1](../ch16-multiway-layouts/01-additive.html).

## Zeros and negative values

The logarithm and the negative powers need \( y>0 \). The common fix \( \log(y+c) \) for a “start” \( c>0 \) makes the
coefficients depend on \( c \) when zeros are frequent. The start can be estimated with the power (@exr-tr-shifted), but a spike at zero usually means two processes, whether a case is positive and how large
it is, which no transformation separates; two-part models and generalized linear models (Chapter 34) are
more natural. For responses of both signs, Yeo and Johnson (2000) give a family defined on the whole line.

## Exercises

### A. Check your understanding

::: {#exr-tr-elasticity-reading}
[A1]

In the log–log fit of @exm-tr-engel-scales the slope is \( 0.8559 \). By what factor does the median
food expenditure differ between two households whose incomes differ by \( 10 \) percent? By what factor
between households whose incomes differ by a factor of two? Why is the answer the same for the median,
the geometric mean and, under @eq-tr-multiplicative, the mean?
:::

::: {.solution}
The log of the ratio of incomes is \( \log 1.1 \), so the ratio of medians is
\( 1.1^{0.8559}=\exp(0.8559\log1.1) \), about \( 1.085 \): a difference of about \( 8.5 \) percent. For a factor
of two the ratio is \( 2^{0.8559} \), about \( 1.81 \). By @prp-tr-multiplicative(c) every quantile, the geometric
mean and the mean scale by the same factor \( \exp\{(\x'-\x)\T\bbeta\} \), because the error multiplies
the whole distribution.
:::

### B. Practice

::: {#exr-tr-additive-log}
[B1]

Let \( Y=\mu+\varepsilon \) with \( \mu>0 \), \( \E\varepsilon=0 \), \( \Var\varepsilon=\sigma^2 \) and \( \sigma/\mu \) small.
Use \( \log(1+u)=u-u^2/2+O(u^3) \) to show that, to the leading order,
\( \E\log Y\approx\log\mu-\sigma^2/(2\mu^2) \) and \( \Var\log Y\approx\sigma^2/\mu^2 \). Conclude that when
errors are additive with constant variance, the log scale is heteroscedastic, with the variance
falling as the mean rises, and its mean is not exactly \( \log\mu \).
:::

::: {.solution}
Write \( \log Y=\log\mu+\log(1+u) \) with \( u=\varepsilon/\mu \). Then
\( \log Y=\log\mu+u-u^2/2+O(u^3) \). Taking expectations, \( \E u=0 \) and \( \E u^2=\sigma^2/\mu^2 \), which gives the
mean. For the variance, the leading term of \( \log Y-\E\log Y \) is \( u \), with variance \( \sigma^2/\mu^2 \). Both
depend on \( \mu \), hence on the regressors, unless \( \sigma \) is proportional to \( \mu \), which is the
multiplicative case.
:::

::: {#exr-tr-poisson-scales}
[B2]

Let \( Y \) be Poisson with a large mean \( \mu \). Use the first-order expansion of the delta method to show
that \( \Var(\log Y)\approx1/\mu \) and \( \Var(\sqrt{Y})\approx\tfrac14 \), and that
\( \E\sqrt{Y}\approx\sqrt{\mu} \). (Read \( \log Y \) given \( Y>0 \), since \( \Pr(Y=0)=e^{-\mu} \).) If \( \mu=\exp(\x\T\bbeta) \), which scale has a linear mean and which a constant
variance?
:::

::: {#exr-tr-r2-scales}
[B3]

A colleague fits the same regressors to \( y \) and to \( \log y \) and prefers the log model because its
\( R^2 \) is larger. Explain why the two values of \( R^2 \) are not comparable. Then take \( x=(0,1,2) \) and
\( y=(1,8,7) \), fit a straight line to \( y \) and to \( \log y \), and compare the two values of \( R^2 \) with the
residual sums of squares on the original scale, using \( \exp \) of the fitted values for the log model.
:::

::: {.solution}
\( R^2 \) is the fraction of the total sum of squares of *the response on its own scale* that the model
explains (@thm-proj-r2-cosine). The two fits explain different totals, so their fractions answer
different questions. A likelihood on a common scale, such as the Box–Cox profile of
[Section 22.2](02-box-cox.html), which includes the Jacobian, does compare them. In the example the
raw line has \( R^2=0.628 \) and the log line \( R^2=0.699 \). On the original scale the raw line leaves
a residual sum of squares of \( 32/3 \), about \( 10.67 \), while the exponentiated log fit leaves
\( 27.37 \), more than twice as much. The larger \( R^2 \) belongs to the fit that predicts \( y \) worse,
because it is measured on a scale where the large values count for less.
:::

### C. Going deeper

::: {#exr-tr-removable}
[C1]

In a \( 2\times2 \) layout let \( \mu_{ij}=(a_i+b_j)^2 \) with positive \( a_i,b_j \). Show that the interaction
contrast \( \mu_{11}-\mu_{12}-\mu_{21}+\mu_{22} \) is \( 2(a_1-a_2)(b_1-b_2) \), and that it is zero on the square-root
scale. Now let \( Y_{ij}=\mu_{ij}+\varepsilon_{ij} \) with independent errors of variance \( \sigma^2 \). Is the table of
\( \E\sqrt{Y_{ij}} \) exactly additive? Relate your answer to the caution in
[Section 16.2](../ch16-multiway-layouts/02-interaction.html).
:::
