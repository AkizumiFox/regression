# Power and noncentrality

A test is judged by two probabilities. Its size, the probability of rejecting a true hypothesis, is
fixed by the critical value. Its **power**, the probability of rejecting a false one, depends on how
false the hypothesis is. For the \( F \) test the whole dependence runs through the single number
\( \gamma \), and that makes power calculations unusually transparent. This section states the power
function, interprets the noncentrality, uses it to plan the size of a study, and asks how sensitive the
region test of [Section 11.1](01-reduced-models.html) actually was. It ends with a popular calculation
that should not be made.

## The power function

::: {#thm-glh-power}
[Power of the \( F \) test]

Assume @eq-opt-normal-model with \( r<n \), and let \( H:\bLambda\T\bbeta=\bm d \) be testable with \( q=\rank(\bLambda) \)
(or let \( H \) be a reduced model with \( q=r-r_0 \)). The probability that the level-\( \alpha \) \( F \) test rejects is
\[
\pi(\gamma)=\Pr\bigl\{F(q,n-r,\gamma)>F_\alpha(q,n-r)\bigr\},
\]{#eq-glh-power}

where \( \gamma \) is the noncentrality of @eq-glh-general-noncentrality, or of @eq-glh-noncentrality for a
reduced model. Moreover:

::: {.enumerate options="label=(\alph*)"}
1. \( \pi \) depends on \( (\bbeta,\sigma^2) \) only through \( \gamma \). It is continuous and strictly increasing, with
   \( \pi(0)=\alpha \) and \( \pi(\gamma)\to1 \) as \( \gamma\to\infty \). In particular the test is *unbiased*: its power
   exceeds \( \alpha \) at every alternative.

2. \( \sigma^2\gamma \) is the squared distance from the mean vector \( \X\bbeta \) to the set of mean vectors
   allowed by \( H \):
   \[
\sigma^2\gamma=\min_{\bLambda\T\bb=\bm d}\norm{\X\bbeta-\X\bb}^2
=\min_{\bLambda\T\bb=\bm d}\sum_{i=1}^n\bigl(\x_{(i)}\T\bbeta-\x_{(i)}\T\bb\bigr)^2 .
\]

3. If further observations are added, with the same \( \bbeta \), \( \sigma^2 \) and hypothesis, the noncentrality does
   not decrease. If the whole design is replicated \( k \) times, \( \gamma \) is multiplied by \( k \).
:::

:::

::: {.proof}
The distribution of \( F_H \) is @thm-glh-general-f(c), which gives @eq-glh-power. (a) is @thm-qf-f-power
with \( c=F_\alpha(q,n-r) \). For (b), @thm-glh-general-f(b) applied to the data vector \( \X\bbeta \), for which
\( \bbeta \) is a least squares estimate and the residual sum of squares is zero, shows that
\( \min_{\bLambda\T\bb=\bm d}\norm{\X\bbeta-\X\bb}^2=(\bLambda\T\bbeta-\bm d)\T\W\ginv(\bLambda\T\bbeta-\bm d)=\sigma^2\gamma \).
For (c), adding rows to \( \X \) enlarges \( \C(\X\T) \), so \( H \) remains testable (the constraint
set \( \{\bb:\bLambda\T\bb=\bm d\} \) is unchanged), and the sum in (b) acquires extra nonnegative terms for every
\( \bb \). The minimum of a larger function over the same set is at least as large. Replicating the design
\( k \) times multiplies every sum of squares in (b) by \( k \), and hence the minimum.
:::

Part (b) is the reading of \( \gamma \) to remember. The noncentrality is the smallest total squared
discrepancy, over all observations, between the true mean and a mean the hypothesis allows, in units of
the error variance. It is large when the hypothesis is badly wrong, when \( \sigma^2 \) is small, and when there
are many observations at which the discrepancy shows. Part (c) makes the last point precise. Note that
adding observations also increases \( n-r \), which helps as well, and that the sum in (b) does *not*
depend on the parameterization.

Two consequences follow at once. First, for a fixed alternative, \( \gamma \) grows in proportion to the
number of replicates, so the power tends to one: *every fixed departure from the hypothesis is detected
with certainty in large enough samples.* A significant result from a very large sample may therefore
reflect a departure too small to matter, and effect sizes, not p-values, answer the question "how
much?". Second, the ratio \( \gamma/n \) is a natural *effect size*, the average squared discrepancy per
observation relative to \( \sigma^2 \). It is closely related to the effect size \( f^2 \) of Cohen's
power tables, which are widely used in the behavioural sciences.

## Planning a study

Before the data are collected, @eq-glh-power turns the question "how many observations?" into
arithmetic. One specifies the smallest departure worth detecting, a guess at \( \sigma \), the level and the
desired power, and solves for \( n \). The difficulty is the first ingredient. A departure is a whole
vector \( \bLambda\T\bbeta-\bm d \), and different vectors of the same apparent size give different \( \gamma \). For
the comparison of group means there is a clean way out.

::: {#prp-glh-least-favourable}
[Least and most favourable configurations]

In the balanced one-way layout with \( g\ge2 \) groups of \( m \) observations, group means
\( \mu_1,\dots,\mu_g \) and common variance \( \sigma^2 \), the noncentrality of the \( F \) test of equal means is
\[
\gamma=\frac{m}{\sigma^2}\sum_{i=1}^g(\mu_i-\bar{\mu})^2,\qquad \bar{\mu}=\frac1g\sum_i\mu_i .
\]
Among configurations whose range \( \max_i\mu_i-\min_i\mu_i \) equals \( \Delta \),
\[
\frac{m\Delta^2}{2\sigma^2}\ \le\ \gamma\ \le\ \frac{mg\Delta^2}{4\sigma^2}.
\]
The lower bound is attained exactly when one mean is at each extreme and all others are at the midpoint.
For even \( g \), the upper bound is attained exactly when half the means are at each extreme.
:::

::: {.proof}
The full model has one indicator per group. The reduced model is \( \C(\bone) \), so by
@thm-glh-f-test(c), \( \sigma^2\gamma=\norm{(\I-n^{-1}\bone\bone\T)\bmu}^2 \), where \( \bmu \) repeats \( \mu_i \) on the
\( m \) observations of group \( i \). This is \( m\sum_i(\mu_i-\bar{\mu})^2 \). Let \( \mu_{\max} \) and \( \mu_{\min} \) be two
groups at the extremes and put \( a=\mu_{\max}-\bar{\mu}\ge0 \), \( b=\bar{\mu}-\mu_{\min}\ge0 \), so \( a+b=\Delta \).
Dropping the other groups from the sum,
\[
\sum_i(\mu_i-\bar{\mu})^2\ \ge\ a^2+b^2\ \ge\ \tfrac12(a+b)^2=\tfrac12\Delta^2,
\]
with equality in the first step iff every other mean equals \( \bar{\mu} \), and in the second iff \( a=b \). If all
other means equal \( \bar{\mu} \), then \( \bar{\mu} \) is also the average of the two extremes, so \( a=b \) holds
automatically, and the equality case is the stated configuration. For the upper bound, let \( c \) be the midpoint
of the range. Since \( \bar{\mu} \) minimizes \( t\mapsto\sum_i(\mu_i-t)^2 \),
\( \sum_i(\mu_i-\bar{\mu})^2\le\sum_i(\mu_i-c)^2\le g\Delta^2/4 \), because every \( \lvert\mu_i-c\rvert\le\Delta/2 \). Equality
requires every mean at an extreme and \( \bar{\mu}=c \), which for even \( g \) means half at each end.
:::

A design that has power \( 0.9 \) at the least favourable configuration has at least that power at every
configuration with range \( \Delta \) or more. That is the guarantee a planner usually wants.

::: {#exm-glh-group-size}
[How many observations per group?]

Four treatments are to be compared at level \( 0.05 \), and a spread of one standard deviation between the
best and the worst treatment (\( \Delta/\sigma=1 \)) should be detected with probability at least \( 0.9 \). The
noncentrality per observation per group is \( 0.500 \) for the least favourable
configuration, \( 0.556 \) for equally spaced means, and \( 1.000 \) for the most
favourable one (two means at each end). Searching over \( m \) with @eq-glh-power gives

| Configuration | \( \gamma/m \) | smallest \( m \) with power \( \ge0.9 \) | power at that \( m \) |
|:---|:---:|:---:|:---:|
| least favourable | \( 0.500 \) | \( 30 \) | \( 0.907 \) |
| equally spaced | \( 0.556 \) | \( 27 \) | \( 0.906 \) |
| most favourable | \( 1.000 \) | \( 16 \) | \( 0.917 \) |

To be safe against every configuration with range at least \( \sigma \), the study needs \( 30 \)
observations per group. A simulation of \( 20000 \) data sets at the least favourable configuration
with \( m=30 \) rejects in a proportion \( 0.909 \) of them, against the computed
\( 0.907 \). [Figure 11.5.1](#fig-glh-power-m) shows the whole power curves.
:::

::: {when-format="html"}
![**Figure 11.5.1.** Power of the level-\( 0.05 \) \( F \) test of equal means for four groups of \( m \)
observations, when the group means span one error standard deviation. The three curves are the least
favourable configuration of @prp-glh-least-favourable, equally spaced means, and the most favourable
configuration. Dots mark the smallest \( m \) giving power \( 0.9 \).](power_group_size.svg){#fig-glh-power-m width=70%}
:::

::: {when-format="pdf"}
![Power of the level-\( 0.05 \) \( F \) test of equal means for four groups of \( m \)
observations, when the group means span one error standard deviation. The three curves are the least
favourable configuration of @prp-glh-least-favourable, equally spaced means, and the most favourable
configuration. Dots mark the smallest \( m \) giving power \( 0.9 \).](power_group_size.pdf){width=70%}
:::

```{.python .run #cell-power-design-oneway}
import numpy as np
from scipy import stats

def power(gamma, q, nu, alpha=0.05):
    """P{F(q, nu, gamma) > F_alpha(q, nu)}."""
    crit = stats.f.ppf(1 - alpha, q, nu)
    return stats.ncf.sf(crit, q, nu, gamma) if gamma > 0 else alpha

g, delta_over_sigma = 4, 1.0
patterns = {                                   # means in units of the range Delta
    "least favourable": np.array([0, 0.5, 0.5, 1]),
    "equally spaced": np.array([0, 1 / 3, 2 / 3, 1]),
    "most favourable": np.array([0, 0, 1, 1]),
}

def gamma_oneway(m, means):
    """Noncentrality m * sum (mu_i - mu_bar)^2 / sigma^2 for m observations per group."""
    dev = (means - means.mean()) * delta_over_sigma
    return m * np.sum(dev ** 2)

for name, means in patterns.items():
    m = 2
    while power(gamma_oneway(m, means), g - 1, g * (m - 1)) < 0.9:
        m += 1
    print(f"{name:17s} m = {m} per group, power {power(gamma_oneway(m, means), g - 1, g * (m - 1)):.3f}")
```

The same logic applies to any linear hypothesis. For a regression coefficient \( \beta_j \), with the other
regressors in the model, \( \sigma^2\gamma=\beta_j^2\norm{(\I-\M_{(j)})\x_j}^2 \), where \( \M_{(j)} \) projects onto the
other columns. The design enters through the length of the part of \( \x_j \) that the other regressors
cannot explain (@thm-proj-fwl). Spreading the values of \( x_j \) increases power, and collinearity with the
other regressors destroys it. Adding a covariate that is unrelated to \( \x_j \) but explains part of the
response reduces \( \sigma^2 \), and so increases \( \gamma \), with little loss in the length of the residualized
\( \x_j \).

## How sensitive was the region test?

Power calculations are also useful *after* a nonsignificant result, provided the alternatives are fixed
by the scientific question and not by the data. For the region test of @exm-glh-region, consider the
alternative in which the Northeast differs from the other three regions by \( \delta \) murders per 100,000,
with everything else as in the reduced model. The mean vector is then
\( \X_0\bm c+\delta\,\bm d_N \), where \( \bm d_N \) is the Northeast indicator, and by @thm-glh-f-test(c),
\[
\gamma=\frac{\delta^2\norm{(\I-\Mo)\bm d_N}^2}{\sigma^2}.
\]
There are \( 9 \) Northeast states, but \( \norm{(\I-\Mo)\bm d_N}^2=5.76 \) only: the
covariates already account for part of the Northeast's distinctiveness, and only the rest counts as
evidence about regions. Setting \( \sigma \) at \( s=1.423 \) as a planning value, the power of the
\( 5\% \) test is

| \( \delta \) | \( 0.5 \) | \( 1.0 \) | \( 1.5 \) | \( 2.0 \) |
|:---|:---:|:---:|:---:|:---:|
| power | \( 0.091 \) | \( 0.241 \) | \( 0.506 \) | \( 0.778 \) |

and power \( 0.8 \) needs \( \delta=2.05 \). A Northeast shortfall of \( 1.5 \) murders per 100,000, about the size estimated in
@exm-glh-regions-rank-deficient (\( 1.379 \) below the other three regions), would have been missed about half the time. The
nonsignificant verdict of @exm-glh-region therefore says little about whether regional differences of
that size exist.

```{.python .run #cell-power-design-sensitivity}
import statsmodels.api as sm

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
X0 = np.column_stack([np.ones(n), data["poverty"], data["single"], data["urban"]])
D = np.column_stack([(region == c).astype(float) for c in "NSW"])
X = np.column_stack([X0, D])

b0, *_ = np.linalg.lstsq(X0, D[:, 0], rcond=None)
k2 = np.sum((D[:, 0] - X0 @ b0) ** 2)      # ||(I - M0) d_N||^2: gamma = delta^2 k2 / sigma^2
bf, *_ = np.linalg.lstsq(X, y, rcond=None)
sigma = np.sqrt(np.sum((y - X @ bf) ** 2) / (n - 7))   # planning value for sigma
for delta in (0.5, 1.0, 1.5, 2.0):
    print(f"Northeast lower by {delta}: gamma = {delta**2 * k2 / sigma**2:.2f}, "
          f"power = {power(delta**2 * k2 / sigma**2, 3, n - 7):.3f}")
```

The calculation uses the estimate \( s \) in place of \( \sigma \), which is itself uncertain. A more careful
version would repeat it for a range of plausible values of \( \sigma \), or replace the whole exercise by a
confidence region for the regional effects, which shows directly which differences the data rule out
([Chapter 12](../ch12-intervals-and-bands/index.html)).

## Observed power

A common calculation takes the observed data, estimates the noncentrality from them, and reports the
power "the test had" at that estimate. It is meant to show whether a nonsignificant result was due to low
power. It cannot show that.

::: {#prp-glh-observed-power}
[Observed power is a function of the p-value]

Fix \( q \), \( n-r \) and \( \alpha \), and define the **observed power** as \( \pi(\hat{\gamma}) \), with \( \pi \) as in
@eq-glh-power and \( \hat{\gamma}=qF_{\text{obs}} \). Then the observed power is a strictly decreasing function of the
p-value of the \( F \) test. The same holds for any estimate \( \hat{\gamma}=h(F_{\text{obs}}) \) with \( h \) strictly
increasing and positive.
:::

::: {.proof}
The p-value is \( 1-\Phi_F(F_{\text{obs}}) \), with \( \Phi_F \) the continuous, strictly increasing distribution function
of \( F(q,n-r) \), so \( F_{\text{obs}} \) is a strictly decreasing function of the p-value. The map
\( F_{\text{obs}}\mapsto h(F_{\text{obs}}) \) is strictly increasing, and \( \pi \) is strictly increasing by
@thm-glh-power(a). The composition is strictly decreasing in the p-value.
:::

So observed power adds nothing to the p-value. In the setting of the region test (\( q=3 \),
\( n-r=43 \)), a p-value of \( 0.2 \) always corresponds to an observed power of \( 0.394 \), a
p-value of exactly \( 0.05 \) to \( 0.637 \), and \( 0.01 \) to \( 0.830 \). "The test was not
significant and its observed power was low" says the same thing twice. The informative calculations are
the ones above: power at alternatives chosen for their scientific meaning, or a confidence region.

::: {.warning}
[Reading a nonsignificant \( F \)]

A nonsignificant \( F \) statistic is weak evidence *for* the reduced model unless the test had high power
against every alternative that matters. Whether it had is a question about the design and about which
alternatives matter, answered by @eq-glh-power at alternatives fixed in advance. It is not answered by the
observed power, which merely restates the p-value, and it is not answered by the p-value being large. A
p-value of \( 0.9 \) is no stronger evidence for \( H_0 \) than one of \( 0.3 \) when the test has low power.
:::

## Exercises

### A. Check your understanding

::: {#exr-glh-slope-power}
[A1]

For the test of \( \beta_1=0 \) in simple regression, show that \( \gamma=\beta_1^2S_{xx}/\sigma^2 \), with
\( S_{xx}=\sum_i(x_i-\bar{x})^2 \). With \( n=20 \) design points in \( [0,1] \), compare \( S_{xx} \) for equal spacing
and for ten points at each end. What does the second design give up?
:::

::: {#exr-glh-unbalanced-gamma}
[A2]

For the one-way layout with group sizes \( n_1,\dots,n_g \), show that the noncentrality of the test of equal means is
\( \gamma=\sum_in_i(\mu_i-\bar{\mu}_w)^2/\sigma^2 \), with \( \bar{\mu}_w=\sum_in_i\mu_i/n \).
:::

### B. Practice

::: {#exr-glh-unbalanced-lf}
[B1]

Extend @prp-glh-least-favourable to unequal group sizes: show that among configurations with range \( \Delta \),
\[
\gamma\ \ge\ \frac{\Delta^2}{\sigma^2}\,\frac{1}{1/n_{(1)}+1/n_{(2)}},
\]
where \( n_{(1)}\le n_{(2)} \) are the two smallest group sizes, and describe the configuration attaining it.
:::

::: {.solution}
Let groups \( a \) and \( b \) carry the largest and smallest means. Dropping the other terms of @exr-glh-unbalanced-gamma,
\( \sigma^2\gamma\ge n_a(\mu_a-\bar{\mu}_w)^2+n_b(\mu_b-\bar{\mu}_w)^2\ge\min_t\{n_a(\mu_a-t)^2+n_b(\mu_b-t)^2\} \). The
minimum is at \( t=(n_a\mu_a+n_b\mu_b)/(n_a+n_b) \) and equals \( \Delta^2n_an_b/(n_a+n_b)=\Delta^2/(1/n_a+1/n_b) \). This is
smallest when \( a \) and \( b \) are the two smallest groups. It is attained when those two groups carry the
extremes and all other means equal \( (n_a\mu_a+n_b\mu_b)/(n_a+n_b) \), which is then also \( \bar{\mu}_w \). A departure is
hardest to detect when it is confined to the groups with the fewest observations.
:::

::: {#exr-glh-replication}
[B2]

Take the straight-line model with design points \( 0,1,2,3 \), one observation each, and the hypothesis
\( \beta_1=0 \). Compute \( \gamma \) as a multiple of \( \beta_1^2/\sigma^2 \) and the power at \( \beta_1/\sigma=1 \), level \( 0.05 \). Then
compare (i) replicating the design twice, (ii) adding four observations at \( x=3 \), and (iii) adding four
observations at \( x=1.5 \). Verify @thm-glh-power(c) in each case and explain the ranking.
:::

::: {#exr-glh-observed-half}
[B3]

For \( q=1 \) and large \( n-r \), the \( F \) test is approximately the test that rejects when \( \lvert Z\rvert>z_{\alpha/2} \),
with \( Z\sim\Normal(\delta,1) \). Show that if the observed \( \lvert z\rvert \) equals \( z_{\alpha/2} \) exactly, the observed
power \( \Pr\{\lvert Z\rvert>z_{\alpha/2}\} \) at \( \delta=z_{\alpha/2} \) is just above \( 1/2 \). Conclude that a result on the edge
of significance always has observed power about one half, whatever the design.
:::

::: {.solution}
At \( \delta=z_{\alpha/2} \), \( \Pr\{\lvert Z\rvert>z_{\alpha/2}\}=\Pr\{Z>z_{\alpha/2}\}+\Pr\{Z<-z_{\alpha/2}\}
=\tfrac12+\Pr\{Z-\delta<-2z_{\alpha/2}\} \), and the last probability is \( \Phi(-2z_{\alpha/2}) \), which is
\( \Phi(-3.92)\approx4\times10^{-5} \) for \( \alpha=0.05 \). The value depends only on \( \alpha \), not on \( n \), \( \sigma \) or the design,
as @prp-glh-observed-power predicts.
:::

### C. Going deeper

::: {#exr-glh-k2}
[C1]

In the sensitivity calculation, let \( n_N \) be the number of Northeast states and \( R^2_N \) the coefficient of
determination from regressing the Northeast indicator on the reduced model's columns (which include \( \bone \)).
Show that \( \norm{(\I-\Mo)\bm d_N}^2=n_N(1-n_N/n)(1-R^2_N) \). Interpret each factor, and relate the formula to @eq-proj-vif-preview.
:::

