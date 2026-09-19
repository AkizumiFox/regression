# Balanced and unbalanced one-way layouts

Every result so far holds for any group sizes, yet experimenters aim for balance. This section says what it buys: plain
orthogonality, the largest guaranteed power for a fixed number of observations, and a test of equal means that tolerates unequal
variances. The third matters most in practice, because unbalanced data with unequal variances are common and the classical \( F \)
test can then go badly wrong.

## Orthogonality and simplicity

In a balanced layout, orthogonality of contrasts is plain orthogonality \( \sum_kc_kd_k=0 \) (@thm-aov-orthogonal-contrasts(d)), so
tabulated sets such as the integer polynomial contrasts can be used as they stand. The side conditions \( \sum_k\alpha_k=0 \) and
\( \sum_kn_k\alpha_k=0 \) define the same parameters, the pooled variance is the plain average of the group variances, and every
pairwise difference has the standard error \( s\sqrt{2/m} \). Tukey's intervals have exact simultaneous coverage (@thm-mc-tukey),
whereas with unequal sizes the Tukey–Kramer intervals (@eq-mc-tukey-kramer) are conservative by an unknown amount.

## Balance maximizes guaranteed power

The noncentrality @eq-aov-noncentrality depends on the group sizes as well as the means. No allocation of a fixed total is best for
every configuration of means, and the planner does not know the means, so the relevant criterion is the power *guaranteed* against every
configuration whose means span at least a range \( \Delta \), as in @prp-glh-least-favourable. For unequal group sizes, @exr-glh-unbalanced-lf
already identifies the least favourable configuration; what is new here is the comparison of allocations.

::: {#prp-aov-balance-power}
[Balance maximizes the guaranteed noncentrality]

In a one-way layout with group sizes \( n_1,\dots,n_g \), \( g\ge2 \), let
\[
h=\min_{k\ne l}\frac{n_kn_l}{n_k+n_l}=\min_{k\ne l}\Bigl(\frac1{n_k}+\frac1{n_l}\Bigr)^{-1}.
\]

::: {.enumerate options="label=(\alph*)"}
1. For all means with \( \max_k\mu_k-\min_k\mu_k\ge\Delta>0 \), the noncentrality of the \( F \) test satisfies
   \( \gamma\ge h\Delta^2/\sigma^2 \), and there is a configuration with range exactly \( \Delta \) at which equality holds.

2. For fixed \( n=\sum_kn_k \), \( h\le n/(2g) \), with equality iff the layout is balanced.
:::

:::

::: {.proof}
Part (a) is @exr-glh-unbalanced-lf, rewritten with \( h \): its bound is \( \Delta^2/\bigl(\sigma^2(1/n_{(1)}+1/n_{(2)})\bigr)=h\Delta^2/\sigma^2 \),
since \( 1/n_k+1/n_l \) is largest for the two smallest groups. Equality holds when those two groups carry means \( \Delta \) apart and
every other mean equals their weighted average.

For (b): as just noted, \( 1/n_k+1/n_l \) is largest for the two smallest groups. Call their sizes \( a\le b \). The other \( g-2 \) groups have at least
\( b \) observations each, so \( n\ge a+(g-1)b \), and
\[
\Bigl(\frac1a+\frac1b\Bigr)n\ \ge\ \Bigl(\frac1a+\frac1b\Bigr)\bigl(a+(g-1)b\bigr)=g+(g-1)t+\frac1t,
\]
where \( t=b/a\ge1 \). Now \( (g-1)t+1/t-g=\bigl((g-1)t-1\bigr)(t-1)/t\ge0 \), because both factors are nonnegative when \( t\ge1 \) and \( g\ge2 \). Hence
\( 1/a+1/b\ge2g/n \), that is, \( h\le n/(2g) \). Equality requires \( t=1 \) and \( n=a+(g-1)b \), so all groups have the same size.
Conversely, a balanced layout has \( h=m/2=n/(2g) \).
:::

A departure is hardest to detect when it sits in the least informative pair of groups. The balanced value \( h=m/2 \) is the bound of @prp-glh-least-favourable.

::: {#exm-aov-allocation}
[Forty observations in four groups]

With \( n=40 \) and \( g=4 \), the guaranteed noncentrality per unit of \( (\Delta/\sigma)^2 \) is \( 5.000 \) for the balanced allocation
\( 10,10,10,10 \), \( 3.938 \) for the mildly unbalanced \( 7,9,11,13 \), and \( 2.400 \) for \( 4,6,10,20 \). At \( \Delta=\sigma \) the
guaranteed powers of the level-\( 0.05 \) test are \( 0.399 \), \( 0.320 \) and \( 0.205 \). To guarantee power \( 0.8 \), the range of the
means must be at least \( 1.56\sigma \), \( 1.76\sigma \) and \( 2.25\sigma \) respectively ([Figure 15.5.1](#fig-aov-allocation)). The mild
imbalance gives the same guarantee as a balanced design with \( 8\times3.938\approx31.5 \) observations, so it wastes a fifth of the
experiment.
:::

::: {when-format="html"}
![**Figure 15.5.1.** Guaranteed power of the level-\( 0.05 \) \( F \) test with \( 40 \) observations in four groups, against the
range of the group means in units of \( \sigma \), for three allocations. The curves are the power at the least favourable configuration of @prp-aov-balance-power; dots mark where each reaches \( 0.8 \).](allocation_power.svg){#fig-aov-allocation
width=66%}
:::

::: {when-format="pdf"}
![Guaranteed power of the level-\( 0.05 \) \( F \) test with \( 40 \) observations in four groups, against the
range of the group means in units of \( \sigma \), for three allocations. The curves are the power at the least favourable configuration of @prp-aov-balance-power; dots mark where each reaches \( 0.8 \).](allocation_power.pdf){width=66%}

:::

```{.python .run #cell-allocation-power-guarantee}
import itertools

import numpy as np
from scipy import stats

def guaranteed_gamma(sizes, delta_over_sigma):
    """Smallest noncentrality over all mean configurations with range at least Delta."""
    h = min(a * b / (a + b) for a, b in itertools.combinations(sizes, 2))
    return delta_over_sigma ** 2 * h

def power_f(gamma, q, nu, alpha=0.05):
    return stats.ncf.sf(stats.f.ppf(1 - alpha, q, nu), q, nu, gamma)

allocations = {"balanced": [10, 10, 10, 10], "mild": [7, 9, 11, 13], "strong": [4, 6, 10, 20]}
for name, sizes in allocations.items():
    gam = guaranteed_gamma(sizes, 1.0)
    print(f"{name:9s} {sizes}: guaranteed gamma = {gam:.3f},"
          f" guaranteed power at Delta = sigma: {power_f(gam, 3, 36):.3f}")
```

Balance is optimal for the overall test and for the family of all pairwise comparisons. It is not optimal for every purpose.
A single planned contrast is estimated most precisely with \( n_k\propto\lvert c_k\rvert \) (@exr-aov-contrast-allocation), and
comparisons with a control favour a larger control group (@exr-aov-sqrt-rule). The allocation should follow the questions.

## Unequal variances

The one-way model assumes a common variance. Suppose instead that the observations in group \( k \) have variance \( \sigma_k^2 \), with
everything else as before. What do the two mean squares then estimate?

::: {#thm-aov-heteroscedastic}
[Expected mean squares under unequal variances]

Let the \( Y_{kj} \) be uncorrelated with \( \E(Y_{kj})=\mu_k \) and \( \Var(Y_{kj})=\sigma_k^2 \), and write \( \bar{n}=n/g \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\text{SSE})=\sum_k(n_k-1)\sigma_k^2 \), so \( \E(\text{MSE}) \) is the average of the \( \sigma_k^2 \) weighted by \( n_k-1 \);

2. \( \E(\text{SSB})=\sum_k(1-n_k/n)\sigma_k^2+\sum_kn_k(\mu_k-\bar{\mu})^2 \);

3. when the means are equal,
   \[
\E(\text{MSB})-\E(\text{MSE})=\frac{g(n-1)}{n(g-1)(n-g)}\sum_{k=1}^g(\bar{n}-n_k)\,\sigma_k^2 .
\]{#eq-aov-hetero-bias}

   In a balanced layout this is zero, and both mean squares estimate \( g^{-1}\sum_k\sigma_k^2 \). In an unbalanced layout it is positive
   when the larger variances belong to the smaller groups and negative when they belong to the larger groups.
:::

:::

::: {.proof}
The argument of @thm-ss-expected-mean-squares applies with a diagonal covariance. Now \( \Cov(\Y)=\bSigma \), diagonal with \( \sigma_k^2 \) on the observations of group \( k \), and \( \E(\Y)=\boldsymbol{\uptheta} \). By @thm-rv-quadform-mean,
\( \E(\Y\T\A\Y)=\tr(\A\bSigma)+\boldsymbol{\uptheta}\T\A\boldsymbol{\uptheta} \). (a) With \( \A=\I-\M \), whose diagonal entries in group \( k \) are
\( 1-1/n_k \), \( \tr(\A\bSigma)=\sum_kn_k(1-1/n_k)\sigma_k^2 \), and \( \A\boldsymbol{\uptheta}=\bzero \). (b) With \( \A=\M-\bP_0 \), the diagonal entries in
group \( k \) are \( 1/n_k-1/n \), so \( \tr(\A\bSigma)=\sum_k(1-n_k/n)\sigma_k^2 \), and \( \boldsymbol{\uptheta}\T\A\boldsymbol{\uptheta}=\sum_kn_k(\mu_k-\bar{\mu})^2 \) as in @thm-aov-oneway-f.
(c) With equal means, multiply the difference of the two expectations by \( (g-1)(n-g) \):
\[
\begin{aligned}
&(n-g)\sum_k(1-n_k/n)\,\sigma_k^2-(g-1)\sum_k(n_k-1)\,\sigma_k^2\\
&\qquad=\sum_k\bigl(n-1-g(n-1)n_k/n\bigr)\sigma_k^2
=\frac{g(n-1)}n\sum_k(\bar{n}-n_k)\,\sigma_k^2 ,
\end{aligned}
\]
where the middle step collects the coefficient of \( \sigma_k^2 \): \( (n-g)-(n-g)n_k/n-(g-1)n_k+(g-1)=n-1-gn_k+gn_k/n \). Since
\( \sum_k(\bar{n}-n_k)=0 \), the sum equals \( -\sum_k(n_k-\bar{n})(\sigma_k^2-\tilde{\sigma}^2) \) for any constant \( \tilde{\sigma}^2 \), which is
negative when \( n_k \) and \( \sigma_k^2 \) move together and positive when they move in opposite directions. In a balanced layout every
\( \bar{n}-n_k \) is zero, and then \( \E(\text{MSE})=\sum_k(m-1)\sigma_k^2/(g(m-1))=g^{-1}\sum_k\sigma_k^2 \).
:::

The theorem predicts the behaviour of the \( F \) test. In a balanced layout both mean squares estimate the same average variance under
the hypothesis, and although the ratio is no longer exactly \( F \)-distributed, Box (1954) showed that the effect on the level is modest.
In an unbalanced layout \( \text{MSE} \) leans towards the variances of the large groups and \( \text{MSB} \) towards those of the small
groups. When the small groups are the variable ones the test rejects too often; when the large groups are, it rejects too rarely and
loses power.

::: {#exm-aov-size-simulation}
[How far the level moves]

Four normal groups have equal means and standard deviations increasing geometrically from \( 1 \) to \( R \). Three designs of \( 40 \)
observations are compared: balanced (\( 10 \) per group), *positive pairing* (sizes \( 4,6,10,20 \), the largest group most variable) and
*negative pairing* (sizes \( 20,10,6,4 \), the smallest group most variable). For each design and each \( R \), \( 40000 \) data sets were
simulated and three level-\( 0.05 \) tests applied: the classical \( F \) test and the tests of Welch and of Brown and Forsythe described
below ([Figure 15.5.2](#fig-aov-size)).

At \( R=3 \), @eq-aov-hetero-bias gives expected mean square ratios \( \E(\text{MSB})/\E(\text{MSE}) \) of \( 1.000 \), \( 0.559 \) and
\( 1.923 \) for the three designs. The simulated levels of the classical test follow: \( 0.071 \) when balanced, \( 0.011 \) with positive
pairing, and \( 0.241 \) with negative pairing, almost five times the nominal level. At \( R=4 \) the negative pairing gives \( 0.288 \). Welch's
test stays between \( 0.048 \) and \( 0.068 \) throughout, and the Brown–Forsythe test between \( 0.048 \) and \( 0.065 \).

The simulation supports the common rule of thumb that a ratio of about two between the largest and smallest standard deviations does
little harm, but only with balance: at \( R=2 \) the balanced level is \( 0.059 \), while the two unbalanced designs give \( 0.015 \) and
\( 0.160 \).
:::

::: {when-format="html"}
![**Figure 15.5.2.** Rejection rates under equal means of three tests at nominal level \( 0.05 \) (dotted), for four normal groups whose standard
deviations increase from \( 1 \) to the ratio on the horizontal axis. Each point is based on \( 40000 \) simulated data sets.](heteroscedastic_size.svg){#fig-aov-size width=100%}
:::

::: {when-format="pdf"}
![Rejection rates under equal means of three tests at nominal level \( 0.05 \) (dotted), for four normal groups whose standard
deviations increase from \( 1 \) to the ratio on the horizontal axis. Each point is based on \( 40000 \) simulated data sets.](heteroscedastic_size.pdf){width=100%}
:::

## Tests that allow unequal variances

If the variances \( \sigma_k^2 \) were known, weighted least squares ([Section 6.9](../ch06-projections/09-inner-products.html))
would compare the group means with their precision-weighted average \( \tilde{\mu}=\sum_kw_k\bar{Y}_k/\sum_kw_k \), \( w_k=n_k/\sigma_k^2 \),
and \( \sum_kw_k(\bar{Y}_k-\tilde{\mu})^2 \) would have an exact \( \chi^2(g-1) \) law under equal means. In practice the \( \sigma_k^2 \) are
replaced by the \( s_k^2 \) and the reference distribution is approximated. We describe the two approximations in common use without
deriving them; Chapter 21 treats heteroscedasticity in general.

**Welch's test.** Welch (1951) replaced \( \sigma_k^2 \) by \( s_k^2 \) in the weights, \( w_k=n_k/s_k^2 \), and divided the weighted mean square
by a correction factor:
\[
F_W=\frac{\sum_kw_k(\bar{y}_k-\tilde{\mu})^2/(g-1)}{1+\dfrac{2(g-2)}{g^2-1}\sum_kh_k},\qquad
h_k=\frac{(1-w_k/\sum_lw_l)^2}{n_k-1},
\]
referred to \( F(g-1,\nu_W) \) with \( \nu_W=(g^2-1)/(3\sum_kh_k) \). The correction and the degrees of freedom come from an approximation to
the statistic's null distribution that is accurate to first order in the quantities \( 1/(n_k-1) \). For two groups \( F_W \) is the square of Welch's
two-sample \( t \) statistic with Satterthwaite's degrees of freedom (@exr-aov-welch-two).

**The Brown–Forsythe test.** Brown and Forsythe (1974a) kept the unweighted numerator and changed the denominator. By @thm-aov-heteroscedastic(b),
under equal means \( \text{SSB} \) has expectation \( \sum_k(1-n_k/n)\sigma_k^2 \), whatever the variances.
Replacing \( \sigma_k^2 \) by \( s_k^2 \) gives an unbiased estimate of that expectation, and the statistic
\[
F_{BF}=\frac{\text{SSB}}{\sum_k(1-n_k/n)s_k^2}
\]
compares \( \text{SSB} \) with an estimate of what it should be under the hypothesis. It is referred to \( F(g-1,f) \), where
\( 1/f=\sum_kc_k^2/(n_k-1) \) with \( c_k=(1-n_k/n)s_k^2/\sum_l(1-n_l/n)s_l^2 \), Satterthwaite's approximation to the degrees of freedom of a
linear combination of independent variance estimates. In a balanced layout the denominator is \( (g-1)s^2 \), so \( F_{BF}=F \) and only the
degrees of freedom change.

Individual contrasts need the same care, since \( s\sqrt{v_{\mathbf{c}}} \) estimates the wrong variance when the \( \sigma_k^2 \) differ, even in a
balanced layout (@exr-aov-contrast-hetero). The **Welch–Satterthwaite** interval uses
\[
\widehat{\Var}(\hat{\psi})=\sum_k\frac{c_k^2s_k^2}{n_k},\qquad
\nu=\frac{\bigl(\sum_kc_k^2s_k^2/n_k\bigr)^2}{\sum_k(c_k^2s_k^2/n_k)^2/(n_k-1)},
\]
and the multiplier \( t_{\nu,\alpha/2} \) (@exr-aov-welch-contrast), adjusted as in [Chapter 13](../ch13-multiplicity/index.html) for
simultaneous inference.

::: {#exm-aov-mortar}
[An unbalanced trial with unequal variances]

A repair mortar was made with four admixture formulations, and its setting time was measured in minutes. The standard
formulation A was tested most often, and the newer formulations least (simulated data, rounded to \( 0.1 \) minute):

| Formulation | \( n_k \) | Setting times (min) | Mean | \( s_k \) |
|:--------:|:--:|:------------------------------------------------|:-----:|:---:|
| A | 12 | 18.9, 20.2, 20.6, 21.2, 20.3, 19.0, 20.2, 19.0, 20.5, 20.6, 20.1, 20.3 | \( 20.075 \) | \( 0.73 \) |
| B | 9 | 22.6, 21.3, 22.2, 22.8, 19.3, 20.7, 21.5, 20.8, 18.8 | \( 21.111 \) | \( 1.38 \) |
| C | 6 | 20.7, 18.2, 18.6, 23.6, 21.3, 22.6 | \( 20.833 \) | \( 2.14 \) |
| D | 4 | 20.1, 22.6, 28.3, 25.3 | \( 24.075 \) | \( 3.53 \) |

The standard deviations grow as the group sizes shrink, the negative pairing of @exm-aov-size-simulation. The classical analysis gives
\( F=5.34 \) on \( (3,27) \) degrees of freedom and \( p=0.0051 \), apparently strong evidence, driven by the high mean of D. But \( \text{MSE}=3.015 \) gives the variance of D weight
\( 3/27 \), while the null expectation of \( \text{MSB} \) in @thm-aov-heteroscedastic(b) gives it weight \( (1-4/31)/3\approx0.29 \); estimated
with the \( s_k^2 \), that expectation is \( 5.407 \), far above \( \text{MSE} \). Welch's test gives
\( F_W=2.64 \) on \( (3,\ 8.3) \) degrees of freedom and \( p=0.119 \), and the Brown–Forsythe test gives \( F_{BF}=2.98 \) on \( (3,\ 6.2) \)
degrees of freedom and \( p=0.115 \). The difference between D and A is \( 4.00 \) minutes. With the pooled standard error \( 1.003 \) its
\( 95\% \) interval is \( [1.94,\ 6.06] \); with the Welch–Satterthwaite standard error \( 1.776 \) on \( 3.1 \) degrees of freedom it is
\( [-1.56,\ 9.56] \). Four observations with a standard deviation of \( 3.5 \) minutes cannot establish a difference of four minutes.

The weighted and unweighted averages of the group means are \( 21.039 \) and \( 21.524 \), so the side conditions \( \sum_kn_k\alpha_k=0 \) and
\( \sum_k\alpha_k=0 \) of [Section 15.1](01-oneway-model.html) (@thm-est-side-conditions) define overall means that differ by half a minute here.
:::

```{.python .run #cell-mortar-trial-summary}
import numpy as np
from scipy import stats
from statsmodels.stats.oneway import anova_oneway

labels = ["A", "B", "C", "D"]
n_k = np.array([12, 9, 6, 4])
rng = np.random.default_rng(1520)
samples = [np.round(rng.normal(mu, sd, k), 1)
           for mu, sd, k in zip([20.0, 20.6, 21.0, 22.2], [0.8, 1.2, 2.0, 3.0], n_k)]

g, n = len(samples), n_k.sum()
means = np.array([x.mean() for x in samples])
s2_k = np.array([x.var(ddof=1) for x in samples])        # group variances
s2 = np.sum((n_k - 1) * s2_k) / (n - g)                  # pooled within-group variance
grand = means @ n_k / n
ss_b = np.sum(n_k * (means - grand) ** 2)
F = ss_b / (g - 1) / s2
print("means:", means.round(3), " sds:", np.sqrt(s2_k).round(2))
print(f"classical F = {F:.2f}, p = {stats.f.sf(F, g - 1, n - g):.4f}")
# null expectation of MS_between if the variances really are s2_k (the theorem of this section)
print(f"MS_within = {s2:.3f};  estimated null E(MS_between) = "
      f"{np.sum((1 - n_k / n) * s2_k) / (g - 1):.3f}")

welch = anova_oneway(samples, use_var="unequal")
bf = anova_oneway(samples, use_var="bf")
print(f"Welch: F = {welch.statistic:.2f} on ({g - 1}, {welch.df_denom:.1f}) df, p = {welch.pvalue:.3f}")
print(f"Brown-Forsythe: F = {bf.statistic:.2f} on ({g - 1}, {bf.df2[1]:.1f}) df, p = {bf.pvalue2:.3f}")
lev = stats.levene(*samples, center="median")           # dispersion test on |y - median|
print(f"variance test on absolute deviations from medians: p = {lev.pvalue:.4f}")
```

(By default `statsmodels` refines the numerator degrees of freedom of the Brown–Forsythe test; `pvalue2` is the original version.)

The simulation of @exm-aov-size-simulation used the same formulas, applied to group means and variances drawn directly from their
sampling distributions:

```{.python .run #cell-heteroscedastic-size-tests}
import numpy as np
from scipy import stats

def three_tests(ybar, s2, n_k):
    """p-values of the classical F, Welch and Brown-Forsythe tests; arrays of shape (reps, g)."""
    g, n = n_k.size, n_k.sum()
    grand = ybar @ n_k / n
    ss_b = ((ybar - grand[:, None]) ** 2) @ n_k
    ss_w = s2 @ (n_k - 1)
    p_f = stats.f.sf(ss_b / (g - 1) / (ss_w / (n - g)), g - 1, n - g)
    # Welch (1951): weights n_k / s_k^2 and a Satterthwaite-type denominator df
    w = n_k / s2
    mu_w = np.sum(w * ybar, axis=1) / w.sum(axis=1)
    h = (1 - w / w.sum(axis=1)[:, None]) ** 2 / (n_k - 1)
    a = np.sum(w * (ybar - mu_w[:, None]) ** 2, axis=1) / (g - 1)
    b = 1 + 2 * (g - 2) / (g ** 2 - 1) * h.sum(axis=1)
    p_w = stats.f.sf(a / b, g - 1, (g ** 2 - 1) / (3 * h.sum(axis=1)))
    # Brown-Forsythe (1974): SS_between over its null expectation under unequal variances
    d = (1 - n_k / n) * s2
    c = d / d.sum(axis=1)[:, None]
    p_bf = stats.f.sf(ss_b / d.sum(axis=1), g - 1, 1 / np.sum(c ** 2 / (n_k - 1), axis=1))
    return p_f, p_w, p_bf

def size(n_k, sds, reps=40000, seed=0, alpha=0.05):
    """Rejection rates under equal means, from group means and variances drawn exactly."""
    rng = np.random.default_rng(seed)
    n_k, sds = np.asarray(n_k, float), np.asarray(sds, float)
    ybar = rng.normal(0, sds / np.sqrt(n_k), (reps, n_k.size))
    s2 = sds ** 2 * rng.chisquare(n_k - 1, (reps, n_k.size)) / (n_k - 1)
    return [np.mean(p < alpha) for p in three_tests(ybar, s2, n_k)]

designs = {"balanced": [10, 10, 10, 10],
           "positive pairing": [4, 6, 10, 20],      # big groups get the big variances
           "negative pairing": [20, 10, 6, 4]}      # small groups get the big variances
ratios = [1, 1.5, 2, 3, 4]
res = {}
for name, sizes in designs.items():                 # sds grow geometrically from 1 to R
    res[name] = np.array([size(sizes, R ** (np.arange(4) / 3), seed=int(10 * R)) for R in ratios])
    for R, row in zip(ratios, res[name]):
        print(f"{name:17s} R = {R:3.1f}: F, Welch, BF sizes = {np.round(row, 3)}")
```

Robustness has a small price. With equal variances and a group of only four observations, Welch's test rejects slightly too often,
\( 0.061 \) and \( 0.062 \) in the two unbalanced designs. The practical advice
is simple: with unequal group sizes, use a heteroscedasticity-robust test unless there is good reason to believe the variances are
equal; with equal sizes the classical test is acceptable, although individual contrasts still need care.

## Checking the assumptions

Independence is secured by the design, above all by randomization, and cannot be checked well from the data. Correlation within
groups, which arises when the observations in a group share a batch, an operator or a day, inflates the \( F \) statistic badly (@exr-aov-intraclass).
A common variance and normality are checked with the residuals \( e_{kj}=y_{kj}-\bar{y}_k \).

The leverage of an observation in group \( k \) is \( 1/n_k \) (@prp-proj-leverage), so \( \Var(e_{kj})=\sigma^2(1-1/n_k) \). Residuals from small
groups are shrunk towards zero, so they are compared on the common scale \( r_{kj}=e_{kj}/\bigl(s\sqrt{1-1/n_k}\bigr) \). A plot of residuals
against the group means shows whether the spread differs between groups, and a spread growing with the mean suggests a transformation
(Chapter 22). A normal quantile plot of the studentized residuals shows skewness and heavy tails. Chapter 20 develops residual
diagnostics in general.

::: {when-format="html"}
![**Figure 15.5.3.** Residual plots for the mortar trial. (a) Residuals against fitted values: the spread grows from formulation A to
formulation D. (b) Normal quantile plot of the studentized residuals: the curvature at both ends reflects the mixture of groups with
different variances, not a heavy-tailed error distribution.](mortar_residuals.svg){#fig-aov-mortar-residuals width=100%}
:::

::: {when-format="pdf"}
![Residual plots for the mortar trial. (a) Residuals against fitted values: the spread grows from formulation A to
formulation D. (b) Normal quantile plot of the studentized residuals: the curvature at both ends reflects the mixture of groups with
different variances, not a heavy-tailed error distribution.](mortar_residuals.pdf){width=100%}
:::

[Figure 15.5.3](#fig-aov-mortar-residuals) shows the two plots for the mortar trial. The fan shape in panel (a) is unmistakable. The
classical tests of equal variances, which compare the logarithms of the \( s_k^2 \), are so sensitive to nonnormality that a rejection says
as much about the tails as about the variances. Levene (1960) proposed instead to apply the one-way \( F \) test to the absolute deviations
\( \lvert y_{kj}-\bar{y}_k\rvert \), which are larger on average in a more variable group; Brown and Forsythe (1974b) used deviations from the
group *medians*, which behave better for skewed data. For the mortar trial that version gives \( p=0.0019 \). Using such a test as a gate
for the classical \( F \) test is not recommended, since it has little power in small groups, exactly where unequal variances do the most
harm. The analysis should be chosen in advance, from the design and from what is known about the response.

Nonnormality is usually the least serious problem for tests about means: the level of the \( F \) test is fairly robust to it, more so
in balanced designs (Chapter 19), though heavy tails cost power by inflating \( \text{MSE} \). Transformations (Chapter 22) and resampling
(Chapter 23), for example the Kruskal–Wallis permutation test on ranks, are the usual remedies.

## Exercises

### A. Check your understanding

::: {#exr-aov-hetero-balanced}
[A1]

For the balanced design with \( m=10 \) and variances \( 1,2,4,8 \), compute \( \E(\text{MSB}) \) under equal means and \( \E(\text{MSE}) \) from @thm-aov-heteroscedastic.
Repeat for the sizes \( 4,6,10,20 \) and for \( 20,10,6,4 \).
:::

::: {.solution}
Balanced: both equal \( (1+2+4+8)/4=3.75 \). Sizes \( 4,6,10,20 \) (\( n=40 \)): \( \E(\text{SSB})=0.9+0.85\cdot2+0.75\cdot4+0.5\cdot8=9.6 \), so
\( \E(\text{MSB})=3.2 \); \( \E(\text{SSE})=3+10+36+152=201 \), so \( \E(\text{MSE})=201/36\approx5.58 \). The test is conservative. Sizes
\( 20,10,6,4 \): \( \E(\text{SSB})=0.5+0.75\cdot2+0.85\cdot4+0.9\cdot8=12.6 \), so \( \E(\text{MSB})=4.2 \); \( \E(\text{SSE})=19+18+20+24=81 \), so
\( \E(\text{MSE})=2.25 \). The test is liberal.
:::

::: {#exr-aov-leverage}
[A2]

Show that the hat matrix of the one-way model has diagonal entries \( 1/n_k \) and that \( \Var(e_{kj})=\sigma^2(1-1/n_k) \). What is the residual of
an observation that is alone in its group, and why is it useless for checking the model?
:::

::: {.solution}
The hat matrix is \( \M \) of @eq-aov-M, whose diagonal entry for an observation of group \( k \) is \( 1/n_k \). By @prp-lm-fit-moments,
\( \Cov(\he)=\sigma^2(\I-\M) \), so \( \Var(e_{kj})=\sigma^2(1-1/n_k) \). An observation alone in its group has leverage one, its fitted value is
itself, and its residual is identically zero: it carries no information about the variance or about the model.
:::

### B. Practice

::: {#exr-aov-welch-two}
[B1]

Show that for \( g=2 \) the Welch statistic \( F_W \) equals \( t_W^2 \), where \( t_W=(\bar{y}_1-\bar{y}_2)/\sqrt{s_1^2/n_1+s_2^2/n_2} \), and that
\( \nu_W \) is Satterthwaite's \( \bigl(s_1^2/n_1+s_2^2/n_2\bigr)^2/\bigl[(s_1^2/n_1)^2/(n_1-1)+(s_2^2/n_2)^2/(n_2-1)\bigr] \).
:::

::: {.solution}
With \( g=2 \) the correction factor is \( 1 \). Write \( v_k=s_k^2/n_k=1/w_k \) and \( W=w_1+w_2 \). The weighted sum of squares of two numbers about
their weighted mean is \( w_1w_2(\bar{y}_1-\bar{y}_2)^2/W=(\bar{y}_1-\bar{y}_2)^2/(v_1+v_2) \), which is \( t_W^2 \). Next
\( 1-w_1/W=w_2/W=v_1/(v_1+v_2) \), so \( h_1=v_1^2/\{(v_1+v_2)^2(n_1-1)\} \) and similarly for \( h_2 \). Then
\( \nu_W=3/(3\sum_kh_k)=(v_1+v_2)^2/\{v_1^2/(n_1-1)+v_2^2/(n_2-1)\} \).
:::

::: {#exr-aov-bf-two}
[B2]

Show that for \( g=2 \) the Brown–Forsythe statistic also equals \( t_W^2 \), and that its degrees of freedom \( f \) equal \( \nu_W \). So the two
robust tests coincide for two groups and differ only for three or more.
:::

::: {.solution}
With two groups, \( \text{SSB}=n_1n_2(\bar{y}_1-\bar{y}_2)^2/n \), and the denominator is \( (n_2s_1^2+n_1s_2^2)/n \). Their ratio is
\( (\bar{y}_1-\bar{y}_2)^2/(s_1^2/n_1+s_2^2/n_2)=t_W^2 \). The weights are \( c_1=(n_2s_1^2/n)/\{(n_2s_1^2+n_1s_2^2)/n\}=v_1/(v_1+v_2) \) and
\( c_2=v_2/(v_1+v_2) \), so \( 1/f=\sum_kc_k^2/(n_k-1)=1/\nu_W \) by the previous exercise.
:::

::: {#exr-aov-welch-contrast}
[B3]

Let \( \hat{V}=\sum_kc_k^2s_k^2/n_k \) with independent \( s_k^2\sim\sigma_k^2\chi^2(n_k-1)/(n_k-1) \) under normality. Show that
\( \E\hat{V}=V=\Var(\hat{\psi}) \) and \( \Var\hat{V}=2\sum_k(c_k^2\sigma_k^2/n_k)^2/(n_k-1) \). Choose \( \nu \) so that \( \nu\hat{V}/V \) has the mean and
variance of a \( \chi^2(\nu) \) variable, and show that replacing \( \sigma_k^2 \) by \( s_k^2 \) in the answer gives the degrees of freedom of the
Welch–Satterthwaite interval.
:::

::: {.solution}
Since \( \E s_k^2=\sigma_k^2 \), \( \E\hat{V}=\sum_kc_k^2\sigma_k^2/n_k=V \), which is \( \Var(\hat{\psi}) \) because the group means are independent with
variances \( \sigma_k^2/n_k \). \( \Var(s_k^2)=2\sigma_k^4/(n_k-1) \), and the terms are independent, which gives \( \Var\hat{V} \). A \( \chi^2(\nu) \) variable has
mean \( \nu \) and variance \( 2\nu \). The mean of \( \nu\hat{V}/V \) is \( \nu \) for every \( \nu \), and its variance is \( \nu^2\Var(\hat{V})/V^2 \). Setting this
equal to \( 2\nu \) gives \( \nu=2V^2/\Var\hat{V}=V^2/\sum_k(c_k^2\sigma_k^2/n_k)^2/(n_k-1) \). With \( s_k^2 \) for \( \sigma_k^2 \) this is the formula in the
text.
:::

::: {#exr-aov-guarantee-thirty}
[B4]

For \( n=30 \) observations in five groups, compare the quantity \( h \) of @prp-aov-balance-power for the balanced allocation and for the
allocation \( 10,5,5,5,5 \) suggested by @exr-aov-sqrt-rule for comparisons with a control. What does each allocation optimize?
:::

::: {.solution}
Balanced: \( h=6/2=3 \). For \( 10,5,5,5,5 \), the least informative pair is two treatments of size \( 5 \), with \( h=25/10=2.5 \). The guaranteed
noncentrality of the overall test is \( 3(\Delta/\sigma)^2 \) against \( 2.5(\Delta/\sigma)^2 \). But each comparison with the control has variance
\( \sigma^2(1/10+1/5)=0.3\sigma^2 \) instead of \( \sigma^2/3 \). Balance protects the overall test and the comparisons among treatments; the unbalanced
allocation sharpens the comparisons with the control.
:::

### C. Going deeper

::: {#exr-aov-intraclass}
[C1]

In a balanced layout, suppose that under equal means the errors within a group are equicorrelated, \( \Cov(\varepsilon_{kj},\varepsilon_{kj'})=\rho\sigma^2 \)
for \( j\ne j' \) with \( \rho\ge0 \), and independent between groups. Show that \( \E(\text{MSE})=\sigma^2(1-\rho) \) and
\( \E(\text{MSB})=\sigma^2(1-\rho+m\rho) \). What happens to the level of the \( F \) test with \( m=10 \) and \( \rho=0.1 \)?
:::

::: {.solution}
The covariance matrix is \( \sigma^2\{(1-\rho)\I+\rho\Z\Z\T\} \), the random-effects form with \( \sigma_a^2=\rho\sigma^2 \) and error variance
\( (1-\rho)\sigma^2 \). By @exr-aov-random-n0 with \( n_0=m \), \( \E(\text{MSE})=(1-\rho)\sigma^2 \) and \( \E(\text{MSB})=(1-\rho)\sigma^2+m\rho\sigma^2 \).
Under normality \( \text{MSB}/\E(\text{MSB}) \) and \( \text{MSE}/\E(\text{MSE}) \) are independent scaled chi-squared variables, so \( F \) is
\( 1+m\rho/(1-\rho) \) times an \( F(g-1,n-g) \) variable. With \( m=10 \) and \( \rho=0.1 \) the factor is \( 1+1/0.9\approx2.1 \), and a nominal
\( 5\% \) test rejects far more often: with \( g=4 \) groups its level is \( 0.271 \). A small intraclass correlation does great damage when
groups are large.
:::

::: {#exr-aov-contrast-hetero}
[C2]

In a balanced layout with variances \( \sigma_k^2 \), show that the pooled variance estimate of the pairwise difference \( \bar{Y}_k-\bar{Y}_l \),
\( 2s^2/m \), has expectation \( 2\bar{\sigma}^2/m \) with \( \bar{\sigma}^2=g^{-1}\sum_i\sigma_i^2 \), while the true variance is
\( (\sigma_k^2+\sigma_l^2)/m \). Conclude that balance protects the overall test but not individual comparisons, and find the ratio of true to
estimated variance for the two most variable groups when the variances are \( 1,2,4,8 \).
:::

::: {.solution}
By @thm-aov-heteroscedastic(c), in a balanced layout \( \E(s^2)=\bar{\sigma}^2 \), so \( \E(2s^2/m)=2\bar{\sigma}^2/m \), and the true variance is
\( \sigma_k^2/m+\sigma_l^2/m \) because the means are independent. For the groups with variances \( 4 \) and \( 8 \), the true variance is \( 12/m \) and the
pooled estimate targets \( 2\cdot3.75/m=7.5/m \), a ratio of \( 1.6 \). Intervals for that difference built from \( s^2 \) are too short by a factor of about
\( \sqrt{1.6}\approx1.26 \), and intervals for the difference of the two least variable groups are too long.
:::
