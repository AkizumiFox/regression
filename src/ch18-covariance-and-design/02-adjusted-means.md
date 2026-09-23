# Adjusted treatment means

A covariance analysis is usually reported as a table of treatment means "adjusted" for the covariate. The raw mean
\( \bar{y}_k \) estimates the mean response of treatment \( k \) at whatever covariate values its units happened to have; the
adjusted mean estimates it at a covariate value common to all treatments. This section defines adjusted means, gives
their standard errors, and shows how to compare them.

## Definition

Take the one-way covariance model with \( g \) groups of sizes \( n_1,\dots,n_g \) and \( q \) covariates,
\[
\E(y_i)=\mu_k+\bz_{(i)}\T\bgamma\qquad\text{for unit } i \text{ in group } k ,
\]{#eq-dsn-oneway-ancova}

where \( \bz_{(i)}\in\Real^q \) is the covariate vector of unit \( i \), and write \( \bar{\bz}_k \) for the mean covariate
vector of group \( k \) and \( \bar{\bz} \) for the overall mean. We assume @eq-dsn-ancova-rank, and write \( \mathbf{E}_{zz} \)
and \( \hat{\bgamma} \) as in @thm-dsn-ancova.

::: {#def-dsn-adjusted-means}
[Adjusted treatment means]

For a fixed covariate vector \( \bz_0\in\Real^q \), the **adjusted mean** of treatment \( k \) at \( \bz_0 \) is
\( \mu_k+\bz_0\T\bgamma \), the mean response of a unit with covariate \( \bz_0 \) under treatment \( k \). Its estimate is
\[
\hat{\mu}_k(\bz_0)=\bar{y}_k-(\bar{\bz}_k-\bz_0)\T\hat{\bgamma}.
\]{#eq-dsn-adjusted-mean}

When no value is named, \( \bz_0=\bar{\bz} \), the overall mean of the covariates, and \( \hat{\mu}_k(\bar{\bz}) \) is called
simply the adjusted mean of treatment \( k \).
:::

The fitted line (or plane) of treatment \( k \) passes through \( (\bar{\bz}_k,\bar{y}_k) \) with slope \( \hat{\bgamma} \);
the adjusted mean is its height at \( \bz_0 \). In [Figure 18.1.1](01-ancova.html#fig-dsn-ancova-lines) the filled
symbols on the dashed line are the adjusted means.

The adjusted mean is estimable: \( \mu_k+\bz_0\T\bgamma \) is the sum of \( \mu_k+\bar{\bz}_k\T\bgamma \), the mean of the
fitted values of group \( k \), and \( (\bz_0-\bar{\bz}_k)\T\bgamma \), a multiple of the estimable \( \bgamma \) (@thm-dsn-ancova(b)). By @thm-dsn-ancova(c) with \( \boldsymbol{\uprho} \) the scaled indicator of group \( k \), the least squares
estimate of \( \mu_k \) is \( \bar{y}_k-\bar{\bz}_k\T\hat{\bgamma} \), which gives @eq-dsn-adjusted-mean.

## Standard errors

::: {#prp-dsn-adjusted-se}
[Moments of adjusted means]

In the model @eq-dsn-oneway-ancova with \( \Cov(\Y)=\sigma^2\I \), and for a fixed \( \bz_0 \):

::: {.enumerate options="label=(\alph*)"}
1. \( \bar{y}_1,\dots,\bar{y}_g \) are uncorrelated with \( \hat{\bgamma} \), and
   \[
   \Var\bigl(\hat{\mu}_k(\bz_0)\bigr)=\sigma^2\Bigl(\frac1{n_k}+(\bar{\bz}_k-\bz_0)\T\mathbf{E}_{zz}^{-1}(\bar{\bz}_k-\bz_0)\Bigr);
   \]{#eq-dsn-adjusted-var}

2. for \( k\ne l \),
   \( \Cov\bigl(\hat{\mu}_k(\bz_0),\hat{\mu}_l(\bz_0)\bigr)=\sigma^2(\bar{\bz}_k-\bz_0)\T\mathbf{E}_{zz}^{-1}(\bar{\bz}_l-\bz_0) \);

3. for a contrast \( \psi=\sum_kc_k\mu_k \), the estimate \( \hat{\psi}=\sum_kc_k\hat{\mu}_k(\bz_0)=\sum_kc_k\bar{y}_k-\bar{\bz}_c\T\hat{\bgamma} \),
   with \( \bar{\bz}_c=\sum_kc_k\bar{\bz}_k \), does not depend on \( \bz_0 \), and
   \[
   \Var(\hat{\psi})=\sigma^2\Bigl(\sum_k\frac{c_k^2}{n_k}+\bar{\bz}_c\T\mathbf{E}_{zz}^{-1}\bar{\bz}_c\Bigr);
   \]{#eq-dsn-contrast-var}

4. with \( \bz_0=\bar{\bz} \), the weighted average of the adjusted means is the grand mean:
   \( \sum_kn_k\hat{\mu}_k(\bar{\bz})=\sum_i y_i \).
:::

:::

::: {.proof}
Let \( \M \) project onto the span of the group indicators. The group mean \( \bar{y}_k \) is \( \mathbf{u}_k\T\Y \) with
\( \mathbf{u}_k \) the indicator of group \( k \) divided by \( n_k \), a vector in \( \C(\M) \). The estimate
\( \hat{\bgamma}=\mathbf{E}_{zz}^{-1}\Z\T(\I-\M)\Y \) is a linear function of \( (\I-\M)\Y \). Hence
\( \Cov(\bar{y}_k,\hat{\bgamma})=\sigma^2\mathbf{u}_k\T(\I-\M)\Z\mathbf{E}_{zz}^{-1}=\bzero \), since \( (\I-\M)\mathbf{u}_k=\bzero \). Then (a) and (b) follow from \( \Var(\bar{y}_k)=\sigma^2/n_k \), \( \Cov(\bar{y}_k,\bar{y}_l)=0 \) and
\( \Cov(\hat{\bgamma})=\sigma^2\mathbf{E}_{zz}^{-1} \) (@thm-dsn-ancova(e)). For (c), \( \sum_kc_k=0 \) removes \( \bz_0 \), and the
variance follows in the same way. For (d),
\( \sum_kn_k\hat{\mu}_k(\bar{\bz})=\sum_kn_k\bar{y}_k-\bigl(\sum_kn_k\bar{\bz}_k-n\bar{\bz}\bigr)\T\hat{\bgamma} \), and the bracket is zero.
:::

Adjusted means are correlated, because their errors share \( \hat{\bgamma} \). Their variances grow with the distance of
\( \bar{\bz}_k \) from \( \bz_0 \): an adjusted mean is a prediction from the treatment's line, and for a group whose covariate
values barely overlap the others' it rests on extrapolation, as in
[Section 12.3](../ch12-intervals-and-bands/03-prediction.html). Differences do not depend on \( \bz_0 \); only the common level
and the individual standard errors do. The overall mean \( \bar{\bz} \) keeps the adjusted means on the scale of the data (d) and
minimizes \( \sum_kn_k\Var\bigl(\hat{\mu}_k(\bz_0)\bigr) \) (@exr-dsn-best-z0), but any stated value, such as a clinically
meaningful baseline, is legitimate.

## Comparing adjusted means

Pairwise differences are contrasts, so @eq-dsn-contrast-var gives
\[
\text{se}\bigl(\hat{\mu}_k-\hat{\mu}_l\bigr)=s\sqrt{\frac1{n_k}+\frac1{n_l}+(\bar{\bz}_k-\bar{\bz}_l)\T\mathbf{E}_{zz}^{-1}(\bar{\bz}_k-\bar{\bz}_l)},
\qquad s^2=\frac{\text{SSE}}{n-g-q},
\]{#eq-dsn-diff-se}

and a single comparison uses the \( t \) interval of @thm-ci-estimable-interval with \( n-g-q \) degrees of freedom.
For simultaneous comparisons the correlations matter.

::: {.enumerate options="label=(\roman*)"}
1. *Scheffé.* The adjusted treatment \( F \) test of @thm-dsn-ancova(f) is a test of \( \mu_1=\dots=\mu_g \), a
   \( (g-1) \)-dimensional space of estimable contrasts. By @thm-mc-scheffe, the intervals
   \( \hat{\psi}\pm\sqrt{(g-1)F_\alpha(g-1,n-g-q)}\,\text{se}(\hat{\psi}) \) have simultaneous coverage exactly \( 1-\alpha \) over
   *all* contrasts \( \psi \). Some interval excludes zero iff the adjusted \( F \) test rejects (@cor-mc-scheffe-f).

2. *Bonferroni.* For a planned list of \( m \) comparisons, \( t \) intervals at level \( \alpha/m \) have coverage at least
   \( 1-\alpha \) whatever the correlations (@thm-mc-bonferroni).

3. *Tukey-type intervals.* The half-width \( \bigl(q_\alpha(g,n-g-q)/\sqrt2\bigr)\,\text{se}(\hat{\mu}_k-\hat{\mu}_l) \) copies the
   Tukey–Kramer rule @eq-mc-tukey-kramer. It is exact when the covariate means are equal and the groups balanced (@exr-dsn-tukey-exact). Otherwise its coverage is only approximately \( 1-\alpha \): Hayter's theorem covers
   unequal group sizes in the one-way layout, but not the correlations that the covariate introduces. When the
   covariate means are close, as they usually are after randomization, the departure is small (@exr-dsn-tukey-coverage asks you to measure it). When they are far apart,
   use Scheffé's intervals, or compute the exact multiplier for the given correlation matrix by simulation, as
   [Section 13.3](../ch13-multiplicity/03-tukey.html) did for Tukey–Kramer.
:::

::: {#exm-dsn-tutoring-adjusted}
[Adjusted means for the four formats]

In @exm-dsn-tutoring, the overall pretest mean is \( \bar{x}=59.92 \). Format 2 had pretest mean \( 57.2 \), so its
adjusted mean is
\[
75.15-0.8670\,(57.2-59.92)=75.15+2.36=77.51 .
\]
The four adjusted means are \( 70.97 \), \( 77.51 \), \( 69.11 \) and \( 77.59 \), averaging to the grand
mean \( 73.80 \), with standard errors \( 2.02 \), \( 2.03 \), \( 2.01 \) and \( 2.01 \). These
barely differ, because the pretest means are close to \( \bar{x} \) compared with the within-group standard deviation of the
pretest, \( 10.06 \).

The six pairwise comparisons, with \( s=6.34 \) on 35 degrees of freedom, are as follows. The Tukey-type
multiplier is \( q_{0.05}(4,35)/\sqrt2=2.697 \), and Scheffé's is \( \sqrt{3F_{0.05}(3,35)}=2.936 \).

| formats | raw difference | adjusted difference | se | Tukey-type half-width | Scheffé half-width |
|---|---|---|---|---|---|
| \( 1-2 \) | \( -2.29 \) | \( -6.54 \) | \( 2.88 \) | \( 7.77 \) | \( 8.46 \) |
| \( 1-3 \) | \( 3.77 \) | \( 1.86 \) | \( 2.85 \) | \( 7.67 \) | \( 8.36 \) |
| \( 1-4 \) | \( -5.23 \) | \( -6.62 \) | \( 2.84 \) | \( 7.66 \) | \( 8.34 \) |
| \( 2-3 \) | \( 6.06 \) | \( 8.40 \) | \( 2.85 \) | \( 7.69 \) | \( 8.37 \) |
| \( 2-4 \) | \( -2.94 \) | \( -0.08 \) | \( 2.86 \) | \( 7.71 \) | \( 8.39 \) |
| \( 3-4 \) | \( -9.00 \) | \( -8.48 \) | \( 2.84 \) | \( 7.65 \) | \( 8.33 \) |

Both methods declare formats 2 and 4 better than format 3 and resolve nothing else, although the \( 2-3 \)
comparison clears Scheffé's bound by only \( 0.03 \). Without the covariate the standard error of a difference
would be \( 10.73\sqrt{2/10}=4.80 \). The raw difference between formats 3 and 4, \( -9.00 \), would then have had a
Tukey half-width near \( 13 \), and nothing would have been resolved.

Adjustment moved comparisons in both directions: format 2's advantage over format 1 grew from \( 2.29 \) to
\( 6.54 \), while format 4's lead over format 3 shrank from \( 9.00 \) to \( 8.48 \).
:::

```{.python .run #cell-tutoring-adjusted}
import itertools
import numpy as np
from scipy import stats

rng = np.random.default_rng(84)
g, m = 4, 10                                          # four formats, ten students each
n = g * m
group = rng.permutation(np.repeat(np.arange(g), m))   # the randomization
x = np.round(rng.normal(60, 10, n))                   # pretest, measured before assignment
mu_true = np.array([70.0, 74.0, 71.0, 77.0])
y = np.round(mu_true[group] + 0.8 * (x - 60) + rng.normal(0, 6, n), 1)
Z = np.eye(g)[group]                                  # indicator matrix of the formats
n_k = Z.sum(axis=0)
x_bar, y_bar = Z.T @ x / n_k, Z.T @ y / n_k           # group means

def within(v):
    """(I - M) v: deviations from the group means (M projects onto C(Z))."""
    return v - (Z.T @ v / n_k)[group]

E_xx, E_xy, E_yy = within(x) @ within(x), within(x) @ within(y), within(y) @ within(y)
beta = E_xy / E_xx                                    # slope from within-group variation only
sse = E_yy - E_xy ** 2 / E_xx                         # ANCOVA residual sum of squares
nu = n - g - 1

s2 = sse / nu
x0 = x.mean()                                         # the covariate value for adjustment
adj = y_bar - beta * (x_bar - x0)                     # adjusted means
se_adj = np.sqrt(s2 * (1 / n_k + (x_bar - x0) ** 2 / E_xx))
pairs = list(itertools.combinations(range(g), 2))
q = stats.studentized_range.ppf(0.95, g, nu) / np.sqrt(2)   # Tukey-Kramer multiplier
c_s = np.sqrt((g - 1) * stats.f.ppf(0.95, g - 1, nu))       # Scheffe multiplier
for k, l in pairs:
    d = adj[k] - adj[l]
    se_d = np.sqrt(s2 * (1 / n_k[k] + 1 / n_k[l] + (x_bar[k] - x_bar[l]) ** 2 / E_xx))
    print(f"{k + 1}-{l + 1}: raw {y_bar[k] - y_bar[l]:6.2f}  adjusted {d:6.2f}  se {se_d:.2f}"
          f"  Tukey-type +-{q * se_d:.2f}  Scheffe +-{c_s * se_d:.2f}")
print("adjusted means", adj.round(2), " standard errors", se_adj.round(2))
```

## Adjusted means in other designs

For any design part, the adjusted version at \( \bz_0 \) of an estimable \( \blambda\T\bbeta \), such as a treatment mean
averaged over blocks, is \( \blambda\T\bbeta+\bz_0\T\bgamma \), estimated by \( \boldsymbol{\uprho}\T\M(\y-\Z\hat{\bgamma})+\bz_0\T\hat{\bgamma} \) with
variance from @thm-dsn-ancova(e). Software calls these **least squares means**. In unbalanced multi-way data the choice of
\( \blambda \), equal or count-proportional weights over cells, matters as much as \( \bz_0 \)
([Chapter 17](../ch17-unbalanced-data/index.html), @thm-ub-types); it is safest to write the intended function of the cell
means (@def-ub-cell-means) explicitly.

::: {.remark}
[The adjusted mean is not a prediction interval]

The standard error @eq-dsn-adjusted-var describes the uncertainty about the *mean* response of treatment \( k \)
at \( \bz_0 \). A new unit with covariate \( \bz_0 \) given treatment \( k \) will vary around that mean with variance
\( \sigma^2 \) as well. Its prediction interval uses \( s\sqrt{1+1/n_k+(\bar{\bz}_k-\bz_0)\T\mathbf{E}_{zz}^{-1}(\bar{\bz}_k-\bz_0)} \), as in @thm-ci-prediction-interval.
:::

## Exercises

### A. Check your understanding

::: {#exr-dsn-adjusted-hand}
[A1]

From the numbers in @exm-dsn-tutoring, compute the adjusted means of formats 1 and 3 by hand and check them
against the values in @exm-dsn-tutoring-adjusted.
:::

::: {#exr-dsn-adjusted-at-zero}
[A2]

Software that fits the model \( \E(y)=\mu_k+\beta x \) reports \( \hat{\mu}_k \), the adjusted mean at \( x_0=0 \). For a pretest
scored from 0 to 100, what does this number describe? Show that its standard error can be far larger than that of the
adjusted mean at \( \bar{x} \), although both give the same differences.
:::

### B. Practice

::: {#exr-dsn-best-z0}
[B1]

Show that \( \sum_kn_k\Var\bigl(\hat{\mu}_k(\bz_0)\bigr) \) is minimized over \( \bz_0 \) by \( \bz_0=\bar{\bz} \).
:::

::: {.solution}
By @eq-dsn-adjusted-var, the sum is \( \sigma^2\bigl(g+\sum_kn_k(\bar{\bz}_k-\bz_0)\T\mathbf{E}_{zz}^{-1}(\bar{\bz}_k-\bz_0)\bigr) \).
Writing \( \bar{\bz}_k-\bz_0=(\bar{\bz}_k-\bar{\bz})+(\bar{\bz}-\bz_0) \) and using \( \sum_kn_k(\bar{\bz}_k-\bar{\bz})=\bzero \), the cross terms vanish and
the sum becomes
\[
\sigma^2\Bigl(g+\sum_kn_k(\bar{\bz}_k-\bar{\bz})\T\mathbf{E}_{zz}^{-1}(\bar{\bz}_k-\bar{\bz})+n(\bar{\bz}-\bz_0)\T\mathbf{E}_{zz}^{-1}(\bar{\bz}-\bz_0)\Bigr).
\]
The last term is nonnegative, and zero iff \( \bz_0=\bar{\bz} \), because \( \mathbf{E}_{zz} \) is positive definite.
:::

::: {#exr-dsn-baseline-change}
[B2]

A common alternative to covariance adjustment analyses the *change* \( y-x \). Show that this is the covariance model
with \( \beta \) fixed at \( 1 \). Using @prp-dsn-precision's setting with \( \tau^2 \) the variance of \( x \), find the variance of the
difference of mean changes between two groups, and show that it exceeds the variance of the adjusted difference
(ignoring the \( 1/(N-g-2) \) term) unless \( \beta=1 \).
:::

::: {.solution}
The mean change in group \( k \) is \( \bar{y}_k-\bar{x}_k=\mu_k+(\beta-1)\bar{x}_k+\bar{\varepsilon}_k \). The difference of two
groups has variance \( (2/m)\{(\beta-1)^2\tau^2+\sigma^2\} \). The adjusted difference, with the slope known, has variance
\( 2\sigma^2/m \). The change-score analysis is therefore worse by \( (2/m)(\beta-1)^2\tau^2 \), which is zero only when
\( \beta=1 \). It beats the unadjusted analysis, which has \( (\beta-1)^2 \) replaced by \( \beta^2 \), only if
\( \beta>1/2 \). Covariance adjustment estimates the slope rather than assuming one.
:::

### C. Going deeper

::: {#exr-dsn-tukey-exact}
[C1]

Suppose all groups have size \( m \) and the covariate means \( \bar{\bz}_1,\dots,\bar{\bz}_g \) are equal. Show that the adjusted
differences \( \hat{\mu}_k-\hat{\mu}_l \) are the raw differences \( \bar{y}_k-\bar{y}_l \), that these are independent of \( s^2 \) from the
covariance model, and that the Tukey intervals built from \( s \) and \( q_\alpha(g,n-g-q) \) have simultaneous coverage exactly
\( 1-\alpha \).
:::

::: {.solution}
With equal covariate means, \( \bar{\bz}_c=\sum_kc_k\bar{\bz}_k=\bzero \) for every contrast, so the adjusted contrast equals the
raw one and has variance \( \sigma^2\sum_kc_k^2/m \) by @eq-dsn-contrast-var. The group means are functions of \( \M\Y \). The
residual vector \( (\I-\M_W)\Y \) is orthogonal to \( \C(\M_W)\supseteq\C(\M) \), so under normality \( s^2 \) is independent of
\( \bar{y}_1,\dots,\bar{y}_g \) (@thm-qf-orthogonal-projections), and \( (n-g-q)s^2/\sigma^2\sim\chi^2(n-g-q) \). The means are independent
normal with common variance \( \sigma^2/m \). The maximum of \( |(\bar{y}_k-\bar{y}_l)-(\mu_k-\mu_l)|/(s/\sqrt m) \) therefore has the
studentized range distribution with parameters \( g \) and \( n-g-q \), and the proof of @thm-mc-tukey applies unchanged.
:::

::: {#exr-dsn-tukey-coverage}
[C2]

For four groups of ten and one covariate, simulate the simultaneous coverage of the Tukey-type intervals of @exm-dsn-tutoring-adjusted when the covariate means are held at \( (62,57,60,61) \) and when they are held at
\( (70,50,60,60) \), with the within-group covariate variation as in the example. How far from \( 0.95 \) is the coverage in
each case?
:::
