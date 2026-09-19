# The analysis of variance table from projections

Whether the treatments differ at all is answered by the \( F \) test, with the analysis of variance table as its bookkeeping. Both
come from one orthogonal decomposition of the observation space into the constants, the treatment space and the error space, and
every entry of the table is a property of one piece.

## Three subspaces

Write \( \bP_0=n^{-1}\bone\bone\T \) for the projection onto the constants and \( \M \) for the projection onto \( \C(\Z) \) (@eq-aov-M).
Since \( \bone\in\C(\Z) \), @thm-proj-nested applies to the pair \( \spn(\bone)\subseteq\C(\Z) \) and gives the
orthogonal decomposition
\[
\Real^n=\underbrace{\spn(\bone)}_{\dim 1}\dirsum\underbrace{\mathcal T}_{\dim g-1}\dirsum\underbrace{\C(\Z)\perpc}_{\dim n-g},
\qquad
\I=\bP_0+(\M-\bP_0)+(\I-\M),
\]{#eq-aov-three-spaces}

with the treatment space \( \mathcal T=\C(\Z)\cap\bone\perpc \) of [Section 15.2](02-contrasts.html) and
\( \bP_{\mathcal T}=\M-\bP_0 \). The three projections act on the data in the familiar way:
\[
(\bP_0\y)_{kj}=\bar{y},\qquad
(\bP_{\mathcal T}\y)_{kj}=\bar{y}_k-\bar{y},\qquad
((\I-\M)\y)_{kj}=y_{kj}-\bar{y}_k .
\]
So every observation is split as \( y_{kj}=\bar{y}+(\bar{y}_k-\bar{y})+(y_{kj}-\bar{y}_k) \), grand mean plus treatment
deviation plus residual, and the three vectors so formed are mutually orthogonal. Pythagoras gives
\[
\sum_{k,j}y_{kj}^2=n\bar{y}^2+\sum_kn_k(\bar{y}_k-\bar{y})^2+\sum_{k,j}(y_{kj}-\bar{y}_k)^2 .
\]{#eq-aov-pythagoras}

The middle term is the **between-groups sum of squares** \( \text{SSB}=\norm{\bP_{\mathcal T}\y}^2 \), and the last is the
**within-groups** or residual sum of squares \( \text{SSE}=\norm{(\I-\M)\y}^2 \) of @eq-aov-sse. Their sum
\( \sum_{k,j}(y_{kj}-\bar{y})^2 \) is the corrected total sum of squares. The same decomposition, with the distributions
attached, was the first application of Cochran's theorem (@exm-qf-oneway-cochran).

## The \( F \) test

::: {#thm-aov-oneway-f}
[The one-way \( F \) test]

In the normal one-way model of @def-aov-oneway, with mean vector \( \boldsymbol{\uptheta}=\Z\bmu \), put
\[
\text{MSB}=\frac{\text{SSB}}{g-1},\qquad s^2=\text{MSE}=\frac{\text{SSE}}{n-g},\qquad F=\frac{\text{MSB}}{\text{MSE}} .
\]

::: {.enumerate options="label=(\alph*)"}
1. \( \bar{Y} \), \( \text{SSB}=\sum_kn_k(\bar{Y}_k-\bar{Y})^2 \) and \( \text{SSE}=\sum_{k,j}(Y_{kj}-\bar{Y}_k)^2 \) are mutually independent.

2. \( \text{SSE}/\sigma^2\sim\chi^2(n-g) \), and \( \text{SSB}/\sigma^2\sim\chi^2(g-1,\gamma) \) with
   \[
\gamma=\frac{\norm{\bP_{\mathcal T}\boldsymbol{\uptheta}}^2}{\sigma^2}=\frac1{\sigma^2}\sum_{k=1}^gn_k(\mu_k-\bar{\mu})^2,
\qquad \bar{\mu}=\frac1n\sum_{k=1}^gn_k\mu_k .
\]{#eq-aov-noncentrality}

3. \( F\sim F(g-1,n-g,\gamma) \), and \( \gamma=0 \) iff \( \mu_1=\dots=\mu_g \). In particular, when the means are equal,
   \( F\sim F(g-1,n-g) \) whatever their common value and \( \sigma^2 \).

4. \( F \) is the statistic of the \( F \) test of the reduced model \( \E(\Y)\in\spn(\bone) \) against \( \C(\Z) \), and of the
   testable hypothesis \( \alpha_1=\dots=\alpha_g \) in the effects form. The test that rejects when
   \( F>F_\alpha(g-1,n-g) \) is the likelihood ratio test of equal means. Its power is a continuous, strictly increasing
   function of \( \gamma \), equal to \( \alpha \) at \( \gamma=0 \).
:::

:::

::: {.proof}
By @eq-aov-three-spaces the matrices \( \bP_0 \), \( \bP_{\mathcal T} \) and \( \I-\M \) are symmetric, idempotent and mutually
orthogonal, with ranks \( 1 \), \( g-1 \) and \( n-g \). @thm-qf-orthogonal-projections gives (a), since \( \bar{Y} \) is a function of
\( \bP_0\Y \), and it gives the chi-squared laws of (b) with noncentralities \( \norm{\bP_{\mathcal T}\boldsymbol{\uptheta}}^2/\sigma^2 \) and
\( \norm{(\I-\M)\boldsymbol{\uptheta}}^2/\sigma^2 \). The second is zero because \( \boldsymbol{\uptheta}\in\C(\Z) \). The vector \( \boldsymbol{\uptheta} \) has entry \( \mu_k \)
in group \( k \), so \( \M\boldsymbol{\uptheta}=\boldsymbol{\uptheta} \) and \( \bP_0\boldsymbol{\uptheta}=\bar{\mu}\bone \), and \( \bP_{\mathcal T}\boldsymbol{\uptheta} \) has entry
\( \mu_k-\bar{\mu} \) in group \( k \). Its squared length is \( \sum_kn_k(\mu_k-\bar{\mu})^2 \), which proves @eq-aov-noncentrality.

For (c): by (a) and (b), \( F \) is a ratio of independent chi-squared variables divided by their degrees of freedom, the numerator
noncentral, which is the noncentral \( F \) law (@def-qf-noncentral-f). Since every \( n_k\ge1 \), \( \gamma=0 \) iff \( \mu_k=\bar{\mu} \)
for every \( k \).

For (d), take \( \X_0=\bone \) inside \( \X=\Z \). Then \( \text{SSE}_0-\text{SSE}=\norm{(\M-\bP_0)\Y}^2=\text{SSB} \) by @eq-proj-extra-ss,
with \( r-r_0=g-1 \), so \( F \) is the statistic of @thm-glh-f-test. In the effects form the hypothesis is
\( \alpha_k-\alpha_g=0 \), \( k<g \): \( g-1 \) estimable, linearly independent contrasts, so it is testable (@def-glh-testable).
The mean vectors it allows are the constant vectors, so by @thm-glh-general-f(c) its \( F \) statistic is
the same one. The likelihood ratio property is @thm-glh-lrt, and the power statement is @thm-glh-power(a).
:::

The noncentrality @eq-aov-noncentrality is a weighted sum of squared deviations of the group means from their *weighted*
mean. The weights are the group sizes, because a group of \( n_k \) observations carries \( n_k \) copies of its mean in
\( \boldsymbol{\uptheta} \). In a balanced layout \( \gamma=m\sum_k(\mu_k-\bar{\mu})^2/\sigma^2 \), the case studied for planning in @prp-glh-least-favourable
and @exm-glh-group-size. The dependence on the group sizes in unbalanced layouts is taken up in
[Section 15.5](05-balance.html).

## Expected mean squares

The distributions in @thm-aov-oneway-f need normality. The expectations of the mean squares do not, and they explain why the
ratio of two mean squares is the right statistic.

::: {#thm-aov-oneway-ems}
[Expected mean squares in the one-way layout]

In the one-way model of @def-aov-oneway (uncorrelated errors with common variance \( \sigma^2 \), not necessarily normal):

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\text{MSE})=\sigma^2 \);

2. \( \displaystyle\E(\text{MSB})=\sigma^2+\frac{1}{g-1}\sum_{k=1}^gn_k(\mu_k-\bar{\mu})^2 \);

3. in the effects form, \( \sum_kn_k(\mu_k-\bar{\mu})^2=\sum_kn_k(\alpha_k-\bar{\alpha})^2 \) with \( \bar{\alpha}=n^{-1}\sum_kn_k\alpha_k \),
   which is \( \sum_kn_k\alpha_k^2 \) under the side condition \( \sum_kn_k\alpha_k=0 \). In a balanced layout it is
   \( m\sum_k(\mu_k-\tilde{\mu})^2 \) with \( \tilde{\mu}=g^{-1}\sum_k\mu_k \), which is \( m\sum_k\alpha_k^2 \) under \( \sum_k\alpha_k=0 \);

4. \( \E(\text{MSB})\ge\E(\text{MSE}) \), with equality iff \( \mu_1=\dots=\mu_g \).
:::

:::

::: {.proof}
The spaces of @eq-aov-three-spaces form an analysis of variance decomposition of \( \C(\Z) \) in the sense of @thm-ss-expected-mean-squares,
which gives \( \E(\text{SSE})=\sigma^2(n-g)+\norm{(\I-\M)\boldsymbol{\uptheta}}^2=\sigma^2(n-g) \) and
\( \E(\text{SSB})=\sigma^2(g-1)+\norm{\bP_{\mathcal T}\boldsymbol{\uptheta}}^2 \). The squared length was computed in the proof of @thm-aov-oneway-f.
Dividing by the degrees of freedom gives (a) and (b). For (c), \( \mu_k=\mu+\alpha_k \) and
\( \bar{\mu}=\mu+\bar{\alpha} \), so \( \mu_k-\bar{\mu}=\alpha_k-\bar{\alpha} \). In a balanced layout \( \bar{\mu}=\tilde{\mu} \) and
\( \bar{\alpha}=g^{-1}\sum_k\alpha_k \). (d) follows from (a) and (b). Part (b) is also @exr-ss-oneway-ems, which computes it directly from the group means.
:::

Both mean squares estimate \( \sigma^2 \) when the treatments do not differ; otherwise \( \text{MSB} \) estimates \( \sigma^2 \) plus the
treatment variation per degree of freedom. Every \( F \) ratio in an analysis of variance follows this logic: the denominator's
expectation equals the numerator's exactly when the hypothesis holds ([Section 9.6](../ch09-sums-of-squares/06-expected-mean-squares.html)).

The same projections serve the random-effects model of [Section 15.1](01-oneway-model.html), but the expectations change. With
\( \Cov(\Y)=\sigma^2\I+\sigma_a^2\Z\Z\T \) and a constant mean, \( \E(\text{MSE})=\sigma^2 \) still, while
\( \E(\text{MSB})=\sigma^2+n_0\sigma_a^2 \) with \( n_0=(n-\sum_kn_k^2/n)/(g-1) \), which is \( m \) in a balanced layout (@exr-aov-random-n0).
The ratio \( F \) then tests \( \sigma_a^2=0 \), and the difference of the two mean squares divided by \( n_0 \)
estimates \( \sigma_a^2 \). Chapter 32 develops this.

## The table, derived

Collecting the pieces gives the one-way analysis of variance table. Each row is a subspace, its degrees of freedom are the
dimension, its sum of squares is the squared length of the projection of \( \y \), and its expected mean square is \( \sigma^2 \) plus
the squared length of the projection of \( \boldsymbol{\uptheta} \) per dimension.

| Source | Subspace | df | Sum of squares | Mean square | \( F \) | Expected mean square |
|:---|:---|:---:|:---|:---|:---|:---|
| treatments | \( \mathcal T \) | \( g-1 \) | \( \sum_kn_k(\bar{y}_k-\bar{y})^2 \) | MSB | MSB/MSE | \( \sigma^2+\frac{\sum_kn_k(\mu_k-\bar{\mu})^2}{g-1} \) |
| error | \( \C(\Z)\perpc \) | \( n-g \) | \( \sum_{k,j}(y_{kj}-\bar{y}_k)^2 \) | MSE |  | \( \sigma^2 \) |
| corrected total | \( \bone\perpc \) | \( n-1 \) | \( \sum_{k,j}(y_{kj}-\bar{y})^2 \) |  |  |  |

Some tables add a row for the mean, with \( 1 \) degree of freedom and sum of squares \( n\bar{y}^2=\norm{\bP_0\y}^2 \), and an
uncorrected total \( \sum y_{kj}^2 \); its \( F \) ratio would test \( \bar{\mu}=0 \), rarely of interest.

Hand computation used the group totals \( T_k \) and grand total \( T \): \( \norm{\M\y}^2=\sum_kT_k^2/n_k \) and
\( \norm{\bP_0\y}^2=T^2/n \), the "correction factor", so \( \text{SSB}=\sum_kT_k^2/n_k-T^2/n \). Such differences of nearly equal numbers
lose accuracy in floating point ([Chapter 10](../ch10-computation/index.html)); computing deviations first, as the projections do,
avoids this.

::: {#exm-aov-nitrogen-anova}
[The analysis of variance for the nitrogen trial]

For the thirty yields of [Section 15.1](01-oneway-model.html), with grand mean \( \bar{y}=5.6930 \) t/ha:

| Source | df | Sum of squares | Mean square | \( F \) |
|:---|:---:|:---:|:---:|:---:|
| nitrogen rates | \( 4 \) | \( 28.7074 \) | \( 7.1768 \) | \( 26.43 \) |
| plots within rates | \( 25 \) | \( 6.7887 \) | \( 0.2715 \) |  |
| corrected total | \( 29 \) | \( 35.4960 \) |  |  |

The mean row would add \( n\bar{y}^2=972.307 \) on one degree of freedom. The critical value is \( F_{0.05}(4,25)=2.759 \), and
the \( p \)-value is \( 1.2\times10^{-8} \). The table says only that the treatment component of the mean is not zero; the contrasts of [Section 15.2](02-contrasts.html) and
[Section 15.4](04-orthogonal-contrasts.html) say what it is. Since the data were simulated, the true noncentrality is known:
\( \gamma=103.4 \), far more than needed for power near one.
:::

```{.python .run #cell-nitrogen-trial-anova}
import numpy as np
from scipy import stats

rates = np.array([0, 40, 80, 120, 160])            # kg N per hectare
g, m = len(rates), 6                                # 5 rates, 6 plots per rate
true_means = 6.8 - 2.7 * np.exp(-0.018 * rates)     # unknown to the analyst
rng = np.random.default_rng(1515)
group = np.repeat(np.arange(g), m)                  # plot i received rate group[i]
y = np.round(true_means[group] + rng.normal(0, 0.5, g * m), 2)   # yield, t/ha
n = len(y)
Z = (group[:, None] == np.arange(g)).astype(float)  # cell-means model matrix
n_k = Z.sum(axis=0)

M = Z @ np.diag(1 / n_k) @ Z.T
P0 = np.full((n, n), 1 / n)                         # projection onto span(1)
ss_mean = y @ P0 @ y                                # n * ybar^2
ss_between = y @ (M - P0) @ y                       # sum_k n_k (ybar_k - ybar)^2
ss_within = y @ (np.eye(n) - M) @ y                 # sum_kj (y_kj - ybar_k)^2
df_between, df_within = g - 1, n - g
F = (ss_between / df_between) / (ss_within / df_within)
p_value = stats.f.sf(F, df_between, df_within)
print(f"{'source':10s} {'df':>3s} {'SS':>9s} {'MS':>8s}")
print(f"{'rates':10s} {df_between:3d} {ss_between:9.4f} {ss_between / df_between:8.4f}")
print(f"{'plots':10s} {df_within:3d} {ss_within:9.4f} {ss_within / df_within:8.4f}")
print(f"{'corrected':10s} {n - 1:3d} {ss_between + ss_within:9.4f}")
print(f"F = {F:.2f} on ({df_between}, {df_within}) df, p = {p_value:.2g}")
```

## Three ways to read \( F \)

**As a proportion of variation.** Let \( \eta^2=\text{SSB}/(\text{SSB}+\text{SSE}) \), the fraction of the corrected total sum of
squares that lies in the treatment space. It is the \( R^2 \) of the regression of \( \y \) on the group indicators, and the sample
version of the correlation ratio of @eq-cor-correlation-ratio. By @thm-ss-r2-null,
\[
F=\frac{\eta^2/(g-1)}{(1-\eta^2)/(n-g)},
\]
and in the normal model, under equal means, \( \eta^2 \) has the Beta law of that theorem, with \( \E(\eta^2)=(g-1)/(n-1) \), not small with many groups and few observations. In the
nitrogen trial \( \eta^2=0.809 \): four fifths of the variation among plots is variation among rates.

**As the largest contrast.** By @prp-aov-contrast-ss(d), every contrast has \( \text{SS}(\mathbf{c})\le\text{SSB} \), with equality for
\( c_k\propto n_k(\bar{y}_k-\bar{y}) \). Dividing by \( s^2 \), the largest squared \( t \) statistic over all contrasts is \( (g-1)F \).
Hence the \( F \) test rejects iff some contrast, perhaps not an interesting one, is declared nonzero by Scheffé's method (@cor-mc-scheffe-f).


**As a two-sample test, for \( g=2 \).** Then \( \mathcal T \) is a line, the only contrast is \( \mu_1-\mu_2 \) up to scale, and \( F=T^2 \)
for the pooled two-sample \( t \) statistic of @exr-aov-two-groups, whose sensitivity to unequal variances the one-way test inherits
([Section 15.5](05-balance.html)).

## Exercises

### A. Check your understanding

::: {#exr-aov-fill-table}
[A1]

A one-way experiment with four treatments and six observations per treatment gives a treatment mean square of \( 12.0 \) and a
corrected total sum of squares of \( 96.0 \). Complete the analysis of variance table, compute \( F \) and \( \eta^2 \), and state the
reference distribution of \( F \) under equal means.
:::

::: {.solution}
The degrees of freedom are \( 3 \), \( 20 \) and \( 23 \). The treatment sum of squares is \( 3\times12=36 \), so the error sum of squares
is \( 96-36=60 \) and \( \text{MSE}=3 \). Then \( F=12/3=4 \), \( \eta^2=36/96=0.375 \), and under equal means \( F\sim F(3,20) \).
:::

::: {#exr-aov-computing-formula}
[A2]

Show from @eq-aov-pythagoras that \( \text{SSB}=\sum_kT_k^2/n_k-T^2/n \) and
\( \text{SSE}=\sum_{k,j}y_{kj}^2-\sum_kT_k^2/n_k \).
:::

### B. Practice

::: {#exr-aov-two-group-power}
[B1]

For \( g=2 \) with \( n_1+n_2=n \) fixed, show that \( \gamma=(\mu_1-\mu_2)^2n_1n_2/(n\sigma^2) \), and that it is largest for
\( n_1=n_2 \). By how much does the noncentrality fall for the split \( n_1=n/4 \), \( n_2=3n/4 \)?
:::

::: {.solution}
With \( \bar{\mu}=(n_1\mu_1+n_2\mu_2)/n \), \( \mu_1-\bar{\mu}=n_2(\mu_1-\mu_2)/n \) and \( \mu_2-\bar{\mu}=-n_1(\mu_1-\mu_2)/n \). So
\( \sigma^2\gamma=(\mu_1-\mu_2)^2(n_1n_2^2+n_2n_1^2)/n^2=(\mu_1-\mu_2)^2n_1n_2/n \). The product \( n_1n_2 \) with fixed sum is largest at
\( n_1=n_2=n/2 \), where it is \( n^2/4 \). For the split \( (n/4,3n/4) \) it is \( 3n^2/16 \), so the noncentrality falls by a quarter.
:::

::: {#exr-aov-eta2-null}
[B2]

In the normal model with equal means \( \E(\eta^2)=(g-1)/(n-1) \), which is not zero, so \( \eta^2 \) overstates the treatment variation in small
experiments. Show that the numerator of
\[
\hat{\omega}^2=\frac{\text{SSB}-(g-1)\,\text{MSE}}{\text{SSB}+\text{SSE}}
\]
has expectation \( \sum_kn_k(\mu_k-\bar{\mu})^2 \), whatever the means, so that it removes the bias from the numerator. When is
\( \hat{\omega}^2 \) negative?
:::

::: {.solution}
By @thm-aov-oneway-ems, \( \E(\text{SSB})=(g-1)\sigma^2+\sum_kn_k(\mu_k-\bar{\mu})^2 \) and \( \E\{(g-1)\text{MSE}\}=(g-1)\sigma^2 \), so
the difference has the stated expectation. The numerator is negative iff \( \text{MSB}<\text{MSE} \), that is, iff \( F<1 \), which
often has null probability above one half: \( 0.57 \) for the nitrogen design's \( F(4,25) \) in the normal model. The ratio itself is not unbiased, since its
denominator is random.
:::

::: {#exr-aov-contrast-within-f}
[B3]

In the nitrogen trial, compute the Scheffé-maximal contrast \( c_k=n_k(\bar{y}_k-\bar{y}) \) from the rate means, and verify that its
squared \( t \) statistic is \( 4F \). Describe in words what it compares.
:::

### C. Going deeper

::: {#exr-aov-random-n0}
[C1]

In the random-effects model \( \Cov(\Y)=\sigma^2\I+\sigma_a^2\Z\Z\T \), \( \E(\Y)=\mu\bone \), use @thm-rv-quadform-mean to show that
\( \E(\text{SSB})=(g-1)\sigma^2+(n-\sum_kn_k^2/n)\sigma_a^2 \) and \( \E(\text{SSE})=(n-g)\sigma^2 \). Show that
\( n_0=(n-\sum_kn_k^2/n)/(g-1)\le n/g \), with equality iff the layout is balanced.
:::

::: {.solution}
By @thm-rv-quadform-mean, \( \E(\Y\T\A\Y)=\tr(\A\bSigma)+\boldsymbol{\uptheta}\T\A\boldsymbol{\uptheta} \), and \( \boldsymbol{\uptheta}=\mu\bone \) is annihilated by both
\( \bP_{\mathcal T} \) and \( \I-\M \). Now \( \tr(\bP_{\mathcal T})=g-1 \) and
\( \tr(\bP_{\mathcal T}\Z\Z\T)=\tr(\Z\T\M\Z)-\tr(\Z\T\bP_0\Z) \). Since \( \M\Z=\Z \), the first term is \( \tr(\Z\T\Z)=n \). The second is
\( n^{-1}\norm{\Z\T\bone}^2=\sum_kn_k^2/n \). For the error, \( (\I-\M)\Z=\bzero \), so only \( \sigma^2\tr(\I-\M)=(n-g)\sigma^2 \) remains. Finally
\( \sum_kn_k^2\ge n^2/g \) by the Cauchy–Schwarz inequality, with equality iff all \( n_k=n/g \), so
\( (g-1)n_0\le n-n/g=(g-1)n/g \).
:::

::: {#exr-aov-invariance}
[C2]

Show that \( F \) is unchanged when \( \y \) is replaced by \( a\y+b\bone \) with \( a\ne0 \), and when the observations are permuted within
groups. Show that it changes when a constant is added to the observations of one group only, and explain why it should.
:::
