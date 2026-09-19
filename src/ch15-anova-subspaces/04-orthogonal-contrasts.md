# Orthogonal contrasts

Once the \( F \) test of [Section 15.3](03-anova-table.html) says that the treatment component of the mean is not zero, the next step is
to split that component into interpretable pieces. A complete set of orthogonal contrasts chooses an orthogonal basis of the treatment
space, one direction per question, and splits the between-groups sum of squares along it. [Chapter 8](../ch08-estimability/index.html)
proved the arithmetic (@prp-est-contrast-ss) and [Chapter 9](../ch09-sums-of-squares/index.html) placed it in the theory of
single-degree-of-freedom decompositions (@thm-ss-single-df). This section adds the distribution theory, builds orthogonal sets, and
analyses the shape of the nitrogen response.

## The decomposition

Contrasts \( \mathbf{c} \) and \( \mathbf{d} \) are orthogonal for the design if \( \sum_kc_kd_k/n_k=0 \) (@def-est-contrast), which says exactly that
\( \bu_{\mathbf{c}}\perp\bu_{\mathbf{d}} \) (@prp-est-contrast-ss(c)). Since \( \dim\mathcal T=g-1 \), a set of \( g-1 \) pairwise orthogonal contrasts is as
large as possible; we call it **complete**.

::: {#thm-aov-orthogonal-contrasts}
[Orthogonal contrasts decompose the treatment sum of squares]

In the one-way model, let \( \mathbf{c}_1,\dots,\mathbf{c}_{g-1} \) be contrasts that are pairwise orthogonal for the design,
\[
\sum_{k=1}^g\frac{c_{ik}c_{lk}}{n_k}=0\qquad(i\ne l),
\]{#eq-aov-weighted-orthogonality}

and write \( \psi_i=\sum_kc_{ik}\mu_k \) and \( F_i=\text{SS}(\mathbf{c}_i)/s^2 \).

::: {.enumerate options="label=(\alph*)"}
1. The vectors \( \bu_{\mathbf{c}_1},\dots,\bu_{\mathbf{c}_{g-1}} \) form an orthogonal basis of the treatment space \( \mathcal T \), and
   \[
\text{SSB}=\sum_{i=1}^{g-1}\text{SS}(\mathbf{c}_i).
\]

2. In the normal one-way model, \( \text{SS}(\mathbf{c}_1),\dots,\text{SS}(\mathbf{c}_{g-1}) \) and \( \text{SSE} \) are mutually independent, and
   \( \text{SS}(\mathbf{c}_i)/\sigma^2\sim\chi^2(1,\gamma_i) \) with \( \gamma_i=\psi_i^2/(\sigma^2v_{\mathbf{c}_i}) \). The noncentralities add up to
   that of the \( F \) test: \( \sum_i\gamma_i=\gamma \) of @eq-aov-noncentrality.

3. The \( F \) statistic is the average of the single-degree-of-freedom statistics: \( F=(g-1)^{-1}\sum_iF_i \).

4. In a balanced layout, @eq-aov-weighted-orthogonality reduces to \( \sum_kc_{ik}c_{lk}=0 \).

5. Conversely, every orthogonal basis \( \bu_1,\dots,\bu_{g-1} \) of \( \mathcal T \) arises in this way: if \( \bu_i \) takes the value
   \( v_{ik} \) on group \( k \), then \( c_{ik}=n_kv_{ik} \) defines a complete set of orthogonal contrasts with \( \bu_{\mathbf{c}_i}=\bu_i \).
:::

:::

::: {.proof}
Part (a) is @prp-est-contrast-ss(d), whose proof shows that the \( \bu_{\mathbf{c}_i} \) are nonzero, pairwise orthogonal vectors of
\( \mathcal T \), hence a basis, and that the projections \( \bP_i \) onto the lines they span sum to \( \bP_{\mathcal T} \). (b) The
\( \bP_i \) are mutually orthogonal projections of rank one, and each is orthogonal to \( \I-\M \) because its range lies in
\( \C(\Z) \). By @thm-qf-orthogonal-projections the sums of squares \( \norm{\bP_i\Y}^2=\text{SS}(\mathbf{c}_i) \) (@prp-aov-contrast-ss(c))
and \( \text{SSE} \) are independent, with noncentralities \( \norm{\bP_i\boldsymbol{\uptheta}}^2/\sigma^2 \). By the computation in the proof of @prp-aov-contrast-angle,
\( \norm{\bP_i\boldsymbol{\uptheta}}^2=(\bu_{\mathbf{c}_i}\T\boldsymbol{\uptheta})^2/\norm{\bu_{\mathbf{c}_i}}^2=\psi_i^2/v_{\mathbf{c}_i} \). Since
\( \sum_i\bP_i=\bP_{\mathcal T} \) and the \( \bP_i\boldsymbol{\uptheta} \) are orthogonal, \( \sum_i\norm{\bP_i\boldsymbol{\uptheta}}^2=\norm{\bP_{\mathcal T}\boldsymbol{\uptheta}}^2 \),
which is \( \sigma^2\gamma \) by @thm-aov-oneway-f. (c) Divide (a) by \( (g-1)s^2 \). (d) Put \( n_k=m \). (e) In
[Section 15.2](02-contrasts.html) we saw that each \( \bu\in\mathcal T \) equals \( \bu_{\mathbf{c}} \) for the contrast \( c_k=n_kv_k \), and by @prp-est-contrast-ss(c),
\( \sum_kc_{ik}c_{lk}/n_k=\bu_{\mathbf{c}_i}\T\bu_{\mathbf{c}_l} \), which is zero for \( i\ne l \).
:::

By part (e) there are infinitely many complete orthogonal sets, one for each orthogonal basis of \( \mathcal T \), and by part (a) all of
them split the same total. The individual pieces depend on the basis (@exr-aov-all-bases): choosing the basis is choosing the questions,
and the design of the study, not the mathematics, should make that choice.

::: {.warning}
Orthogonality is a convenience, not a requirement. Planned contrasts that are not orthogonal are valid, each with the test and
interval of @prp-aov-contrast-ss. What is lost is additivity: their sums of squares need not total \( \text{SSB} \), and their estimates are
correlated. No contrast should be replaced by a less relevant one merely to achieve orthogonality.
:::

## Building orthogonal sets

Three constructions cover most designed experiments.

**Splitting the treatments in stages.** Many questions form a hierarchy: control against treatments, then old treatments against new,
then one new formulation against another. Each splits a set of levels into two parts, and the comparisons are orthogonal for any group
sizes provided each part is summarized by the mean of all its observations.

::: {#prp-aov-tree-contrasts}
[Contrasts from successive splits]

For a set \( S \) of levels split into disjoint nonempty parts \( A \) and \( B \), let \( N_A=\sum_{k\in A}n_k \) and \( N_B=\sum_{k\in B}n_k \),
and define the contrast
\[
c_k^{(S)}=\begin{cases}n_k/N_A, & k\in A,\\ -n_k/N_B, & k\in B,\\ 0, & k\notin S.\end{cases}
\]
Its estimate is \( \bar{y}_A-\bar{y}_B \), the difference between the means of all observations in \( A \) and in \( B \), and
\[
\text{SS}(\mathbf{c}^{(S)})=\frac{N_AN_B}{N_A+N_B}\,(\bar{y}_A-\bar{y}_B)^2 .
\]

::: {.enumerate options="label=(\alph*)"}
1. Every contrast \( \mathbf{d} \) that is zero outside \( A \) (or outside \( B \)) is orthogonal to \( \mathbf{c}^{(S)} \) for the design.

2. Split the set of all \( g \) levels into two parts, split each part with more than one level into two, and continue until every part
   is a single level. The \( g-1 \) contrasts of the splits form a complete set of orthogonal contrasts, whatever the group sizes.
:::

:::

::: {.proof}
The estimate is \( \sum_{k\in A}n_k\bar{y}_k/N_A-\sum_{k\in B}n_k\bar{y}_k/N_B=\bar{y}_A-\bar{y}_B \), and
\( v_{\mathbf{c}}=\sum_{k\in A}n_k/N_A^2+\sum_{k\in B}n_k/N_B^2=1/N_A+1/N_B \), which gives the sum of squares.
(a) If \( \mathbf{d} \) is zero outside \( A \), then
\( \sum_kc_k^{(S)}d_k/n_k=\sum_{k\in A}(n_k/N_A)d_k/n_k=N_A^{-1}\sum_{k\in A}d_k=0 \), since \( \mathbf{d} \) is a contrast; the case of \( B \) is the
same. (b) A sequence of splits ending in single levels is a binary tree with \( g \) leaves, so it has \( g-1 \) internal nodes, one split
each. Take two different splits, of sets \( S \) and \( S' \). By construction either \( S \) and \( S' \) are disjoint, in which case the two
contrasts have disjoint supports and are orthogonal, or one of them, say \( S' \), lies inside one part of the split of \( S \), in which
case (a) applies. So the \( g-1 \) contrasts are pairwise orthogonal, and they form a complete set.
:::

The weighting matters: comparing plain averages of the levels in \( A \) and \( B \) gives orthogonal contrasts, in general, only in a balanced layout,
as @exm-ss-education-contrasts showed numerically. The Helmert contrasts of @exr-ss-one-way-helmert, which compare each level with the
average of those before it, are the special case of the proposition that splits off one level at a time; with the weights \( n_k/N_A \)
they stay orthogonal for unequal group sizes.

**Factorial treatments.** Four treatments may be the combinations of two fertilizers, each present or absent: none, \( P \) only,
\( K \) only, and both. The contrasts
\[
\begin{gathered}
\text{effect of }P:\ (-1,1,-1,1),\qquad
\text{effect of }K:\ (-1,-1,1,1),\\
\text{interaction}:\ (1,-1,-1,1)
\end{gathered}
\]
are orthogonal in a balanced layout. The third asks whether the effect of \( P \) depends on \( K \). A one-way analysis with these contrasts
is a two-factor analysis in disguise; [Chapter 16](../ch16-multiway-layouts/index.html) treats such layouts directly and defines
interaction in general (@def-tw-interaction).

**Polynomial trends.** When the levels are values of a quantitative variable, such as a dose or a fertilizer rate, the orthogonal
polynomial contrasts describe the shape of the response one degree at a time. For equally
spaced levels and equal group sizes they are the integer vectors obtained from Gram–Schmidt applied to \( 1,s,s^2,\dots \) in the level
scores ([Section 8.7](../ch08-estimability/07-contrasts.html), [Section 9.4](../ch09-sums-of-squares/04-orthogonal-designs.html)). For
five levels they are
\[
\begin{gathered}
\text{linear }(-2,-1,0,1,2),\qquad
\text{quadratic }(2,-1,-2,-1,2),\\
\text{cubic }(-1,2,0,-2,1),\qquad
\text{quartic }(1,-4,6,-4,1).
\end{gathered}
\]
The degree-\( d \) contrast is orthogonal to every polynomial of lower degree in the score, so it vanishes whenever the means lie on a
polynomial of degree less than \( d \). With unequal spacing the Gram–Schmidt step uses the actual level values (@exr-aov-unequal-spacing),
and with unequal group sizes it uses the weighted inner product \( \sum_kn_kp_kq_k \) and the contrasts
\( c_k=n_kp_k \) (@exm-est-party-trend).

## Trend analysis of the nitrogen trial

::: {#exm-aov-nitrogen-trend}
[The shape of the nitrogen response]

The five rates are equally spaced and the design is balanced, so the integer contrasts above apply. The four single-degree-of-freedom
sums of squares, each tested against \( s^2=0.2715 \) on \( 25 \) degrees of freedom, are:

| Component | Sum of squares | \( F \) | \( p \) | Share of \( \text{SSB} \) |
|:---|:---:|:---:|:---:|:---:|
| linear | \( 25.9121 \) | \( 95.42 \) | \( 5.1\times10^{-10} \) | \( 90.3\% \) |
| quadratic | \( 2.7904 \) | \( 10.28 \) | \( 0.0037 \) | \( 9.7\% \) |
| cubic | \( 0.0049 \) | \( 0.02 \) | \( 0.89 \) | \( 0.0\% \) |
| quartic | \( 8.6\times10^{-6} \) | \( 3.2\times10^{-5} \) | \( 0.996 \) | \( 0.0\% \) |
| rates (total) | \( 28.7074 \) | \( 26.43 \) |  |  |

The pieces add to the between-rates sum of squares of @exm-aov-nitrogen-anova, and \( F \) is the average of the four single-degree-of-freedom
\( F \) values. Yield rises with the rate and the rise flattens: the quadratic component is clearly significant. The cubic and quartic
components together test whether a quadratic in the rate describes the five means; their sum of squares \( 0.0049 \) on two degrees of
freedom gives \( F=0.009 \) and \( p=0.99 \) ([Figure 15.4.1](#fig-aov-nitrogen-trend)).
:::

::: {when-format="html"}
![**Figure 15.4.1.** Trend analysis of the nitrogen trial. (a) Rate means with bars of two standard errors, and the least
squares linear and quadratic curves in the rate fitted to all thirty plots. (b) The four orthogonal polynomial sums of squares
(logarithmic scale); the dotted line is the value a single degree of freedom needs to be significant at level \( 0.05 \),
\( s^2F_{0.05}(1,25) \).](nitrogen_trend.svg){#fig-aov-nitrogen-trend width=100%}
:::

::: {when-format="pdf"}
![Trend analysis of the nitrogen trial. (a) Rate means with bars of two standard errors, and the least
squares linear and quadratic curves in the rate fitted to all thirty plots. (b) The four orthogonal polynomial sums of squares
(logarithmic scale); the dotted line is the value a single degree of freedom needs to be significant at level \( 0.05 \),
\( s^2F_{0.05}(1,25) \).](nitrogen_trend.pdf){width=100%}
:::

```{.python .run #cell-nitrogen-trend-trend}
import numpy as np
from scipy import stats

rates = np.array([0, 40, 80, 120, 160])            # kg N per hectare
g, m = len(rates), 6                                # 5 rates, 6 plots per rate
true_means = 6.8 - 2.7 * np.exp(-0.018 * rates)     # unknown to the analyst
rng = np.random.default_rng(1515)
group = np.repeat(np.arange(g), m)                  # plot i received rate group[i]
y = np.round(true_means[group] + rng.normal(0, 0.5, g * m), 2)   # yield, t/ha
n = len(y)

means = np.array([y[group == k].mean() for k in range(g)])
nu = n - g
s2 = np.sum((y - means[group]) ** 2) / nu
ss_between = m * np.sum((means - means.mean()) ** 2)

C_poly = np.array([[-2, -1, 0, 1, 2],               # linear
                   [2, -1, -2, -1, 2],              # quadratic
                   [-1, 2, 0, -2, 1],               # cubic
                   [1, -4, 6, -4, 1]], float)       # quartic

def contrast_ss(c):
    return (c @ means) ** 2 / (c @ c / m)           # balanced: sum c_k^2 / n_k = c'c / m

ss_poly = np.array([contrast_ss(c) for c in C_poly])
F_poly = ss_poly / s2
p_poly = stats.f.sf(F_poly, 1, nu)
for name, ss, F, p in zip(["linear", "quadratic", "cubic", "quartic"], ss_poly, F_poly, p_poly):
    print(f"{name:9s} SS = {ss:8.4f}  F = {F:7.2f}  p = {p:.2g}")
print(f"sum = {ss_poly.sum():.4f}   between rates = {ss_between:.4f}")
```

The test of the cubic and quartic components is not a new test. Every model in which the mean depends on the rate alone has its mean
space inside \( \C(\Z) \) (@lem-cor-row-structure), so the one-way residual sum of squares is the pure-error sum of squares of
[Section 14.6](../ch14-correlation-lack-of-fit-prediction/06-lack-of-fit.html). The quadratic regression on the rate misses exactly the
cubic and quartic directions of \( \C(\Z) \), so their sum of squares is its lack-of-fit sum of squares, and the trend test is the pure-error
test of @thm-cor-lack-of-fit (@exr-aov-trend-lof). The listing checks the identity: the quadratic regression has residual sum of squares
\( 6.7935 \), the lack of fit plus the within-rates sum of squares.

```{.python .run #cell-nitrogen-trend-lack}
ss_lof = ss_poly[2] + ss_poly[3]                    # cubic + quartic, 2 df
F_lof = (ss_lof / 2) / s2
X2 = np.column_stack([np.ones(n), rates[group], rates[group] ** 2])
coef, *_ = np.linalg.lstsq(X2, y, rcond=None)       # quadratic regression on the rate
sse_quad = np.sum((y - X2 @ coef) ** 2)
print(f"lack of fit of the quadratic: SS = {ss_lof:.4f}, F = {F_lof:.3f},"
      f" p = {stats.f.sf(F_lof, 2, nu):.2f}")
print(f"quadratic regression SSE = {sse_quad:.4f} = lack of fit + within = {ss_lof + s2 * nu:.4f}")
```

The fitted quadratic is \( 4.014+0.03466\,N-0.0001139\,N^2 \) t/ha at \( N \) kg N/ha. It has a maximum at
\( N=0.03466/(2\times0.0001139)\approx152 \), inside the range of rates, and it declines after that. That decline is a property of the
parabola, not of the crop: the data were generated from a response that rises at every rate, and five noisy means cannot distinguish a
parabola that turns over from a curve that levels off. A trend analysis describes the means at the rates used and should not be
extrapolated.

## One data set, two decompositions

A different set of questions gives a different complete set. Suppose the agronomist had planned to ask first whether fertilizer helps
at all, and then how the response behaves among the fertilized rates. The contrast of the control against the four fertilized rates is a
split of the kind in @prp-aov-tree-contrasts, with coefficients proportional to \( (-4,1,1,1,1) \). The linear, quadratic and cubic
polynomial contrasts among the four fertilized rates, \( (0,-3,-1,1,3) \), \( (0,1,-1,-1,1) \) and \( (0,-1,3,-3,1) \), are zero at the control,
so they are orthogonal to it by part (a) of that proposition and to each other by balance. Their sums of squares are \( 21.3701 \)
(\( F=78.70 \)), \( 6.6364 \) (\( F=24.44 \)), \( 0.7004 \) (\( F=2.58 \), \( p=0.12 \)) and \( 0.0005 \), and they add to the same
\( 28.7074 \).

```{.python .run #cell-nitrogen-trend-second}
C_ctrl = np.array([[-4, 1, 1, 1, 1],                # control against the four fertilized
                   [0, -3, -1, 1, 3],               # linear among 40, ..., 160
                   [0, 1, -1, -1, 1],               # quadratic among them
                   [0, -1, 3, -3, 1]], float)       # cubic among them
ss_ctrl = np.array([contrast_ss(c) for c in C_ctrl])
print("second set:", ss_ctrl.round(4), " sum =", round(ss_ctrl.sum(), 4))
```

The two decompositions tell different but compatible stories. In the first, curvature is significant; in the second, the curvature
among the fertilized rates is not, because most of the bend is the jump from no nitrogen to \( 40 \) kg, which the second set assigns to
the control contrast. Computing several decompositions and reporting the one whose pieces
look best is the data snooping of [Section 13.6](../ch13-multiplicity/06-choosing.html).

## Exercises

### A. Check your understanding

::: {#exr-aov-poly5-check}
[A1]

Check that the four polynomial contrasts for five levels are pairwise orthogonal. Show that the quadratic contrast vanishes when
\( \mu_k=a+bk \), and the cubic contrast vanishes when \( \mu_k=a+bk+ck^2 \).
:::

::: {#exr-aov-factorial-contrasts}
[A2]

For the four fertilizer combinations of "Factorial treatments", verify that the three contrasts are orthogonal when every group has
the same size, and that the \( P \) and \( K \) contrasts are not orthogonal for group sizes \( 3,3,3,5 \) (in the order none, \( P \), \( K \),
both). What does the lack of orthogonality mean for the estimates?
:::

::: {.solution}
With equal sizes, orthogonality is the plain inner product: \( (-1,1,-1,1)\cdot(-1,-1,1,1)=1-1-1+1=0 \), and similarly for the other
two pairs. With sizes \( 3,3,3,5 \), \( \sum_kc_kd_k/n_k=\tfrac13-\tfrac13-\tfrac13+\tfrac15=-\tfrac2{15}\ne0 \). By @prp-est-contrast-ss(c) the two
estimates are then correlated, with covariance \( -\tfrac2{15}\sigma^2 \), and their sums of squares do not add to the part of \( \text{SSB} \) that
the two contrasts span.
:::

### B. Practice

::: {#exr-aov-extend-planned}
[B1]

Show that the two planned contrasts of @exm-aov-nitrogen-planned, \( (-1,\tfrac14,\tfrac14,\tfrac14,\tfrac14) \) and \( (0,0,-1,0,1) \), are
orthogonal. Find two further contrasts that complete them to a complete orthogonal set, and say what they compare.
:::

::: {.solution}
The design is balanced, and \( -1\cdot0+\tfrac14(0-1+0+1)=0 \). A contrast \( \mathbf{e} \) orthogonal to the first satisfies
\( -e_1+\tfrac14(e_2+e_3+e_4+e_5)=0 \), and since \( e_2+e_3+e_4+e_5=-e_1 \) this forces \( e_1=0 \). Orthogonality to the second forces
\( e_3=e_5 \). So \( \mathbf{e}=(0,a,b,c,b) \) with \( a+2b+c=0 \). Two orthogonal choices are \( (0,1,0,-1,0) \), comparing \( 40 \) with \( 120 \), and
\( (0,1,-1,1,-1) \), comparing \( 40 \) and \( 120 \) together with \( 80 \) and \( 160 \) together. Their inner product is \( 1-1=0 \). The completion is
valid but not very natural: completing a planned set is optional.
:::

::: {#exr-aov-unequal-spacing}
[B2]

Five plots at each of four rates \( 0 \), \( 1 \), \( 2 \) and \( 4 \) (in units of \( 50 \) kg N/ha) are to be analysed for trend. Show that the linear
and quadratic orthogonal polynomial contrasts are proportional to \( (-7,-3,1,9) \) and \( (7,-4,-8,5) \). Check that the quadratic contrast
vanishes when the means are linear in the rate.
:::

::: {.solution}
The rates have mean \( 7/4 \), so the linear contrast is proportional to \( 4(s_k-7/4)=(-7,-3,1,9) \). For the quadratic, centre
\( s_k^2=(0,1,4,16) \) at its mean \( 21/4 \) and multiply by four: \( (-21,-17,-5,43) \). Its inner product with \( (-7,-3,1,9) \) is \( 580 \), and the
linear vector has squared length \( 140 \), so subtract \( \tfrac{580}{140}(-7,-3,1,9) \) to get \( (8,-\tfrac{32}{7},-\tfrac{64}{7},\tfrac{40}{7}) \),
which is proportional to \( (7,-4,-8,5) \). Its entries sum to zero and its inner product with the linear contrast is
\( -49+12-8+45=0 \). For means \( a+bs_k \) the quadratic contrast gives \( b(0-4-16+20)=0 \).
:::

::: {#exr-aov-symmetric-sizes}
[B3]

Four equally spaced levels have group sizes \( 2,4,4,2 \). Show that the classical linear contrast \( (-3,-1,1,3) \) is orthogonal for this
design to the quadratic \( (1,-1,-1,1) \) but not to the cubic \( (-1,3,-3,1) \). Explain the first fact by symmetry.
:::

::: {.solution}
For linear and quadratic, \( \sum_kc_kd_k/n_k=-\tfrac32+\tfrac14-\tfrac14+\tfrac32=0 \). For linear and cubic, it is
\( \tfrac32-\tfrac34-\tfrac34+\tfrac32=\tfrac32\ne0 \). The sizes are symmetric about the middle of the scale. The linear contrast is odd under the
reflection \( k\mapsto5-k \), the quadratic is even, and the weights \( 1/n_k \) are even, so the weighted inner product of an odd and an even
vector vanishes. Linear and cubic are both odd, and nothing forces their weighted inner product to vanish.
:::

::: {#exr-aov-tree-nitrogen}
[B4]

Use the formula of @prp-aov-tree-contrasts to recompute the sum of squares of the control against the fertilized rates in the nitrogen trial
from the rate means, and compare it with the value \( 21.3701 \) in the text.
:::

::: {.solution}
Here \( N_A=24 \), \( N_B=6 \), \( \bar{y}_A=\tfrac14(5.237+6.057+6.515+6.652)=6.11525 \) and \( \bar{y}_B=4.005 \). The sum of squares is
\( \tfrac{24\cdot6}{30}(6.11525-4.005)^2=4.8\times2.11025^2\approx21.375 \). The small difference from \( 21.3701 \) comes from using the
rounded rate means.
:::

### C. Going deeper

::: {#exr-aov-all-bases}
[C1]

Let \( \mathbf{U}=[\bu_1,\dots,\bu_{g-1}] \) and \( \mathbf{V}=[\bv_1,\dots,\bv_{g-1}] \) be two orthonormal bases of the treatment space. Show that
\( \mathbf{V}=\mathbf{U}\mathbf{Q} \) for an orthogonal \( (g-1)\times(g-1) \) matrix \( \mathbf{Q} \), that the vectors of squared coordinates of \( \y \) in the two bases
have the same sum \( \text{SSB} \), and that the individual squared coordinates generally differ. More generally, show that for any symmetric
idempotent \( \A \) of rank \( q \) and \( \Y\sim\Normal_n(\boldsymbol{\uptheta},\sigma^2\I) \), \( \Y\T\A\Y \) splits into \( q \) independent scaled noncentral
\( \chi^2(1) \) variables, one for each vector of an orthonormal basis of \( \C(\A) \).
:::

::: {.solution}
Both \( \mathbf{U}\mathbf{U}\T \) and \( \mathbf{V}\mathbf{V}\T \) equal \( \bP_{\mathcal T} \), the unique projection onto \( \mathcal T \). Put \( \mathbf{Q}=\mathbf{U}\T\mathbf{V} \).
Then \( \mathbf{U}\mathbf{Q}=\bP_{\mathcal T}\mathbf{V}=\mathbf{V} \) and \( \mathbf{Q}\T\mathbf{Q}=\mathbf{V}\T\bP_{\mathcal T}\mathbf{V}=\mathbf{V}\T\mathbf{V}=\I \). The coordinate vectors satisfy
\( \mathbf{V}\T\y=\mathbf{Q}\T\mathbf{U}\T\y \), and an orthogonal matrix preserves length, so both have squared length
\( \y\T\bP_{\mathcal T}\y=\text{SSB} \). The individual coordinates are mixed by \( \mathbf{Q} \) and generally change. For the general statement, write
\( \A=\sum_{i=1}^q\bw_i\bw_i\T \) with \( \bw_1,\dots,\bw_q \) an orthonormal basis of \( \C(\A) \) (@prp-proj-orthonormal-formula). The rank-one
projections \( \bw_i\bw_i\T \) are mutually orthogonal, so by @thm-qf-orthogonal-projections the terms \( (\bw_i\T\Y)^2 \) are independent, and
\( (\bw_i\T\Y)^2/\sigma^2\sim\chi^2(1,(\bw_i\T\boldsymbol{\uptheta})^2/\sigma^2) \). Eigenvectors of \( \A \) for the eigenvalue one are one such basis.
:::

::: {#exr-aov-trend-lof}
[C2]

For \( g \) equally spaced levels with equal group sizes, show that the sum of the sums of squares of the polynomial contrasts of degrees
\( d+1,\dots,g-1 \) is the lack-of-fit sum of squares of the polynomial regression of degree \( d<g-1 \) on the level score, and that the
sum for degrees \( 1,\dots,d \) is that regression's sum of squares after the mean.
:::

::: {.solution}
Let \( \mathcal P_d\subseteq\C(\Z) \) be the mean space of the degree-\( d \) regression, spanned by the vectors whose entry in group \( k \) is
\( s_k^i \), \( i=0,\dots,d \) (@lem-cor-row-structure). The contrast vector of the degree-\( i \) orthogonal polynomial \( p_i \) has entry
\( p_i(s_k)/m \) in group \( k \), and in a balanced layout orthogonality of these vectors is orthogonality of the \( p_i \) in the plain inner
product. By construction \( p_0,\dots,p_d \) span the polynomials of degree at most \( d \) on the scores, so the vectors of degrees \( 1,\dots,d \)
span \( \mathcal P_d\cap\bone\perpc \), and those of degrees \( d+1,\dots,g-1 \) span \( \C(\Z)\cap\mathcal P_d\perpc \). The projection onto the
first space gives the regression sum of squares after the mean (@thm-proj-nested), and the projection onto the second gives
\( \norm{(\M-\bP_{\mathcal P_d})\y}^2 \), the lack-of-fit sum of squares of @thm-cor-lack-of-fit. Each projection is the sum of the
projections onto the lines, by @thm-proj-sum.
:::
