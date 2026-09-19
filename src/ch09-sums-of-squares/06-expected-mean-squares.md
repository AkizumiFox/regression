# Expected mean squares

This section computes the expectations of sums of squares when the data are random, using only two moments, and their
distributions under normality. Each mean square estimates \( \sigma^2 \) plus a term that vanishes
exactly when its hypothesis is true, while the residual mean square estimates \( \sigma^2 \) alone.
That is the structure of the \( F \) ratio, on which Part III builds.

## Expected sums of squares for a decomposition

In this section \( \E(\Y)=\bmu \) and \( \Cov(\Y)=\sigma^2\I \). We do not assume at first that
\( \bmu\in\C(\X) \), because the effect of a wrong model on the mean squares is part of the story.

::: {#thm-ss-expected-mean-squares}
[Expected mean squares]

Let \( \mathcal V_1,\dots,\mathcal V_k \) be an analysis of variance decomposition of \( \C(\X) \), with
projections \( \bP_i \) and dimensions \( r_i \), and let \( \E(\Y)=\bmu \), \( \Cov(\Y)=\sigma^2\I \).

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\Y\T\bP_i\Y)=\sigma^2r_i+\norm{\bP_i\bmu}^2 \), so the mean square of \( \mathcal V_i \) has
   expectation
   \[
\E(\text{MS}_i)=\sigma^2+\frac{\norm{\bP_i\bmu}^2}{r_i}.
\]{#eq-ss-ems}

2. \( \E(\text{SSE})=\sigma^2(n-r)+\norm{(\I-\M)\bmu}^2 \), so \( \E(s^2)=\sigma^2 \) iff \( \bmu\in\C(\X) \).

3. If \( \bmu=\X\bbeta \) with \( \X=[\X_0,\dots,\X_k] \) partitioned into terms, the partial sum of
   squares of \( \X_j \) has
   \[
\E\,\text{SS}(\X_j\mid\X_{(-j)})=\sigma^2d_j+\bbeta_j\T\X_j\T(\I-\M_{(-j)})\X_j\bbeta_j ,
\]
   with \( d_j=\rank(\M-\M_{(-j)}) \), and the sequential sum of squares of \( \X_j \) has
   \( \E\,\text{SS}(\X_j\mid\X_0,\dots,\X_{j-1})=\sigma^2d_j'+\norm{(\M_j-\M_{j-1})\X\bbeta}^2 \), with
   \( d_j'=\rank(\M_j-\M_{j-1}) \).

4. If moreover \( \Y \) is normal, then \( \Var(\Y\T\bP_i\Y)=2\sigma^4r_i+4\sigma^2\norm{\bP_i\bmu}^2 \),
   and the sums of squares of different subspaces, including the error space, are uncorrelated.
:::

:::

::: {.proof}
(a) By @thm-rv-quadform-mean with \( \A=\bP_i \) and \( \bSigma=\sigma^2\I \),
\( \E(\Y\T\bP_i\Y)=\sigma^2\tr(\bP_i)+\bmu\T\bP_i\bmu \). The trace is \( r_i \) (@prp-proj-trace-rank), and
\( \bmu\T\bP_i\bmu=\bmu\T\bP_i\T\bP_i\bmu=\norm{\bP_i\bmu}^2 \). Divide by \( r_i \).

(b) The same with \( \A=\I-\M \), of trace \( n-r \). The second term vanishes iff \( \bmu\in\C(\X) \).

(c) For the partial sum of squares, \( (\M-\M_{(-j)})\X_l\bbeta_l=\bzero \) for \( l\neq j \), because
\( \X_l\bbeta_l \) lies in \( \C(\X_{(-j)}) \) and is fixed by both projections. And
\( (\M-\M_{(-j)})\X_j\bbeta_j=\X_j\bbeta_j-\M_{(-j)}\X_j\bbeta_j=(\I-\M_{(-j)})\X_j\bbeta_j \). Its squared
length is \( \bbeta_j\T\X_j\T(\I-\M_{(-j)})\X_j\bbeta_j \), since \( \I-\M_{(-j)} \) is symmetric and
idempotent. The sequential case is (a) for the chain decomposition (@def-ss-sequential).

(d) By @thm-qf-mean-var(b) with \( \bSigma=\sigma^2\I \),
\( \Var(\Y\T\bP_i\Y)=2\tr\{(\sigma^2\bP_i)^2\}+4\sigma^2\bmu\T\bP_i\bP_i\bmu=2\sigma^4r_i+4\sigma^2\norm{\bP_i\bmu}^2 \).
By part (d) of the same theorem the covariance of two forms is
\( 2\sigma^4\tr(\bP_i\bP_j)+4\sigma^2\bmu\T\bP_i\bP_j\bmu=0 \) when \( \bP_i\bP_j=\bzero \).
:::

Parts (a) to (c) use nothing beyond second moments. They hold for skewed or heavy-tailed errors,
provided the errors are uncorrelated with a common variance. Part (d) uses normality, and without
it both the variances and the zero correlations fail in general (@exr-ss-kurtosis-sse).

Part (c) confirms @prp-ss-tested-regression. For full-rank \( \X \) the term
\( \norm{(\I-\M_{(-j)})\X_j\bbeta_j}^2 \) is positive unless \( \bbeta_j=\bzero \), and it is large only if
\( \X_j \) has a large part orthogonal to the other terms. A nearly collinear term has a small
expected partial sum of squares even when its coefficients are substantial.

Under normality the expectations become full distributions.

::: {#cor-ss-ms-distributions}
[Distributions of the mean squares]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \) with \( \bmu\in\C(\X) \), and let \( \mathcal V_1,\dots,\mathcal V_k \)
be an analysis of variance decomposition of \( \C(\X) \). Then the sums of squares
\( \Y\T\bP_1\Y,\dots,\Y\T\bP_k\Y \) and \( \text{SSE} \) are mutually independent,
\[
\frac{\Y\T\bP_i\Y}{\sigma^2}\sim\chi^2\Bigl(r_i,\ \gamma_i\Bigr),\quad
\gamma_i=\frac{\norm{\bP_i\bmu}^2}{\sigma^2},
\qquad
\frac{\text{SSE}}{\sigma^2}\sim\chi^2(n-r),
\]
and \( F_i=\text{MS}_i/s^2\sim F(r_i,n-r,\gamma_i) \), which is central iff \( \bP_i\bmu=\bzero \).
The same holds for a single projection \( \bP \) onto any subspace of \( \C(\X) \), such as a partial
sum of squares, paired with \( \text{SSE} \).
:::

::: {.proof}
The projections \( \bP_1,\dots,\bP_k,\I-\M \) are mutually orthogonal (@thm-ss-decomposition(a)), so
@thm-qf-orthogonal-projections gives independence and the noncentral chi-squared laws, with
\( \norm{(\I-\M)\bmu}^2=0 \) for the error space. The ratio of two independent chi-squared variables,
each divided by its degrees of freedom and the numerator noncentral, is the noncentral \( F \)
(@def-qf-noncentral-f). A single projection \( \bP \) onto a subspace of \( \C(\X) \) satisfies
\( \bP(\I-\M)=\bzero \), and the same argument applies to the pair.
:::

The residual part is the familiar @thm-opt-sampling; the rest applies the same argument to the
other pieces.

## The logic of the \( F \) ratio

By @eq-ss-ems every mean square estimates \( \sigma^2 \) plus \( \norm{\bP_i\bmu}^2/r_i \), the squared size
per dimension of the mean's component in \( \mathcal V_i \), while under a correct model the residual
mean square estimates \( \sigma^2 \) alone. So the ratio
\[
F_i=\frac{\text{MS}_i}{s^2}
\]
compares two estimates of \( \sigma^2 \) that agree in expectation exactly when \( \bP_i\bmu=\bzero \).
This is the logic of the analysis of
variance: *choose a denominator whose expectation equals that of the numerator when the
hypothesis holds.* Independence gives the distribution exactly, and the noncentrality
\( \gamma_i=\norm{\bP_i\bmu}^2/\sigma^2 \) governs the power. The probability that \( F_i \) exceeds any
fixed critical value increases strictly with \( \gamma_i \) (@thm-qf-f-power). The mean of the ratio
follows from independence and \( \E(1/V)=1/(\nu-2) \) for \( V\sim\chi^2(\nu) \), \( \nu>2 \):
\[
\E(F_i)=\frac{n-r}{n-r-2}\cdot\frac{r_i+\gamma_i}{r_i}.
\]{#eq-ss-f-mean}

The principle is more general than the choice of \( s^2 \). With random effects the expected mean
squares contain variance components, and the right denominator may be another mean square
(@exr-ss-random-oneway; Chapters 32 and 33). Also, when the model is wrong, @thm-ss-expected-mean-squares(b) shows that the denominator
itself is inflated, by \( \norm{(\I-\M)\bmu}^2/(n-r) \) (compare @exm-rv-rss-bias). Omitted terms inflate
the denominator of every \( F \) ratio, which pushes the ratios towards zero, and they can also
shift the numerators, because the omitted part of
\( \bmu \) may have components in the \( \mathcal V_i \) (@exr-ss-omitted-term).

[Chapter 11](../ch11-general-linear-hypothesis/index.html) derives the \( F \) test from the likelihood ratio principle (@thm-glh-lrt) and studies its power (@thm-glh-power).

## Expected mean squares in an unbalanced layout

In the two-way layout of [Section 9.3](03-sequential-partial.html), @eq-ss-ems makes the hypotheses
of @thm-ss-two-way-hypotheses concrete. For the sequential sum of squares of \( A \) entered first,
\( (\bP_A-\bP_0)\bmu \) has entry \( \bar{\mu}^w_i-\bar{\mu}^w \) on each of the \( n_{i\cdot} \) observations of
row \( i \), where \( \bar{\mu}^w=\sum_in_{i\cdot}\bar{\mu}^w_i/n \). So
\[
\E\,\text{MS}(A\mid\mu)=\sigma^2+\frac{1}{a-1}\sum_in_{i\cdot}\bigl(\bar{\mu}^w_i-\bar{\mu}^w\bigr)^2 .
\]
For the Type III sum of squares, the noncentrality is the sum of squares @eq-ss-hypothesis evaluated
at the true means, divided by \( \sigma^2 \) (@exr-ss-hypothesis-noncentrality):
\[
\gamma=\frac{(\bLambda\T\bmu_c)\T\bigl(\bLambda\T\bm N^{-1}\bLambda\bigr)^{-1}(\bLambda\T\bmu_c)}{\sigma^2},
\]
where \( \bmu_c \) is the vector of \( ab \) cell means, \( \bm N=\diag(n_{ij}) \), and the columns of
\( \bLambda \) express the \( a-1 \) contrasts of unweighted row means. In the first, the design enters
through the weights \( n_{ij}/n_{i\cdot} \), which also define the hypothesis. In the second, the
hypothesis does not depend on the design, which enters only through the precision matrix
\( \bm N^{-1} \) of the cell means.

::: {#exm-ss-ems-simulation}
[Checking expected mean squares by simulation]

Take the design of @exm-ss-anes-types: \( n=944 \) respondents in the nine cells of party by
education, with the observed counts. Set the true cell means equal to the observed cell means, and
\( \sigma^2=1.21 \). The script generates \( 40{,}000 \) data sets with normal errors and
\( 40{,}000 \) with centred exponential errors of the same variance, which are strongly skewed.
For each it computes the sequential decomposition (party, then education, then interaction), the
residual, and the Type III sum of squares for party, and averages the mean squares:

| Source | df | \( \gamma \) | \( \E(\text{MS}) \) | normal | skewed | \( \Var(\text{MS}) \) | simulated |
|---|---|---|---|---|---|---|---|
| mean | 1 | 14594.91 | 17661.050 | 17660.423 | 17661.259 | 85476.55 | 86357.17 |
| party | 2 | 535.09 | 324.942 | 324.915 | 324.787 | 784.89 | 786.39 |
| education given party | 2 | 35.75 | 22.840 | 22.836 | 22.891 | 53.81 | 53.13 |
| interaction | 4 | 29.60 | 10.163 | 10.174 | 10.179 | 11.57 | 11.75 |
| residual | 935 | 0.00 | 1.210 | 1.210 | 1.210 | 0.0031 | 0.0031 |
| party, Type III | 2 | 537.53 | 326.414 | 326.362 | 326.259 | 788.46 | 789.47 |

The columns \( \E(\text{MS}) \) and \( \Var(\text{MS}) \) come from @eq-ss-ems and
@thm-ss-expected-mean-squares(d), with \( \Var(\text{MS}_i)=\sigma^4(2r_i+4\gamma_i)/r_i^2 \); the
last column is the simulated variance under normal errors. The averages agree with the
expectations under both error distributions, as parts (a) and (b) promise. The variances are another matter. Under skewed errors the variance of the residual
mean square is \( 0.0124 \), \( 3.97 \) times the normal-theory value, because the variance of
a quadratic form involves the fourth moment of the errors (@thm-rv-quadform-variance). And the
correlation between the interaction and residual mean squares, \( 0.008 \) under
normal errors (zero in theory), is \( 0.020 \) under skewed errors, against a theoretical
\( 0.025 \) from the third and fourth moments: the two mean squares are no longer independent.

For the interaction, the noncentrality \( \gamma=29.60 \) is large. By @eq-ss-f-mean, the \( F \) ratio
has mean \( 8.417 \), and the simulated mean is \( 8.427 \). Its \( 5\% \)
critical value is \( 2.381 \), and the power is \( 0.9967 \) in theory and
\( 0.9967 \) in the simulation. [Figure 9.6.1](#fig-ss-ems-f) shows the simulated
distribution with the noncentral \( F \) density.
:::

::: {when-format="html"}
![**Figure 9.6.1.** The \( F \) ratio for the interaction in @exm-ss-ems-simulation over \( 40{,}000 \)
simulated data sets with normal errors, with the noncentral \( F(4,935,29.6) \) density of
@cor-ss-ms-distributions and the central null density \( F(4,935) \) with its upper \( 5\% \) point.](ems_f_ratio.svg){#fig-ss-ems-f width=82%}
:::

::: {when-format="pdf"}
![The \( F \) ratio for the interaction in @exm-ss-ems-simulation over \( 40{,}000 \)
simulated data sets with normal errors, with the noncentral \( F(4,935,29.6) \) density of
@cor-ss-ms-distributions and the central null density \( F(4,935) \) with its upper \( 5\% \) point.](ems_f_ratio.pdf){width=82%}
:::

```{.python .run #cell-ems-simulation-pieces}
import numpy as np
import pandas as pd
import statsmodels.api as sm
d = sm.datasets.anes96.load_pandas().data
party = pd.cut(d["PID"], [-1, 1.5, 4.5, 6], labels=["Dem", "Ind", "Rep"])
educ = pd.cut(d["educ"], [0, 3.5, 5.5, 7], labels=["HS", "College", "Graduate"])
cell = party.cat.codes.to_numpy() * 3 + educ.cat.codes.to_numpy()
n = len(cell)

def basis(*blocks):
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    return U[:, s > 1e-10 * s[0]]

one = np.ones((n, 1))
A = np.eye(3)[party.cat.codes.to_numpy()]
B = np.eye(3)[educ.cat.codes.to_numpy()]
Cell = np.eye(9)[cell]
U0, UA, UAB, UC = basis(one), basis(one, A), basis(one, A, B), basis(Cell)

def complement(U_big, U_small):
    """Orthonormal basis of C(U_small)-perp within C(U_big)."""
    R = U_big - U_small @ (U_small.T @ U_big)
    return basis(R)

pieces = {"mean": U0, "party": complement(UA, U0),
          "educ | party": complement(UAB, UA), "interaction": complement(UC, UAB)}
Ub = basis(one, B, np.column_stack([(A[:, [i]] - A[:, [0]]) * (B[:, [j]] - B[:, [0]])
                                     for i in (1, 2) for j in (1, 2)]))
pieces["party (Type III)"] = complement(UC, Ub)       # unweighted-means hypothesis

sigma = 1.1
mu = d.groupby(cell)["selfLR"].mean().to_numpy()[cell]  # true means: observed cell means
for name, U in pieces.items():
    r = U.shape[1]
    gamma = np.sum((U.T @ mu) ** 2) / sigma ** 2        # noncentrality ||P mu||^2 / sigma^2
    print(f"{name:17s} rank {r}  gamma = {gamma:9.2f}"
          f"  E(MS) = {sigma ** 2 * (1 + gamma / r):9.3f}")
```

```{.python .run #cell-ems-simulation-simulate}
rng = np.random.default_rng(99)
reps, chunk = 40_000, 4_000
dfe = n - 9
sums = {"normal": [], "skewed": []}
for kind in sums:
    out = []
    for _ in range(reps // chunk):
        if kind == "normal":
            E = sigma * rng.normal(size=(n, chunk))
        else:                                             # centred exponential, variance sigma^2
            E = sigma * (rng.exponential(size=(n, chunk)) - 1.0)
        Y = mu[:, None] + E
        ms = {name: np.sum((U.T @ Y) ** 2, axis=0) / U.shape[1] for name, U in pieces.items()}
        ms["residual"] = (np.sum(Y ** 2, axis=0) - np.sum((UC.T @ Y) ** 2, axis=0)) / dfe
        out.append(pd.DataFrame(ms))
    sums[kind] = pd.concat(out, ignore_index=True)
print(pd.DataFrame({k: v.mean() for k, v in sums.items()}).round(3))
```

The simulation never forms an \( n\times n \) matrix: each sum of squares comes from the coordinates
of \( \Y \) in an orthonormal basis of its subspace (@thm-ss-decomposition(d)).

## Exercises

### A. Check your understanding

::: {#exr-ss-simple-ems}
[A1]

In simple regression with \( \E(\Y)=\beta_0\bone+\beta_1\x \), show that the regression mean square
has expectation \( \sigma^2+\beta_1^2S_{xx} \), and that the \( F \) ratio of the corrected table has
noncentrality \( \beta_1^2S_{xx}/\sigma^2 \) under normality.
:::

::: {#exr-ss-oneway-ems}
[A2]

In a one-way layout with \( g \) groups, group sizes \( n_i \) and means \( \mu_i \), show that the
between-groups mean square has expectation
\( \sigma^2+\sum_in_i(\mu_i-\bar{\mu})^2/(g-1) \), where \( \bar{\mu}=\sum_in_i\mu_i/n \). Compare
with @exm-qf-oneway-cochran.
:::

### B. Practice

::: {#exr-ss-random-oneway}
[B1]

In a balanced one-way layout with \( g \) groups of \( m \) observations, let
\( Y_{ij}=\mu+a_i+\varepsilon_{ij} \), where the \( a_i \) and \( \varepsilon_{ij} \) are uncorrelated
random variables with mean zero and variances \( \sigma_a^2 \) and \( \sigma^2 \). Show that
\( \Cov(\Y)=\sigma^2\I+\sigma_a^2\bm Z\bm Z\T \), with \( \bm Z \) the matrix of group indicators, and use
@thm-rv-quadform-mean to show that
\[
\E(\text{MS}_{\text{between}})=\sigma^2+m\sigma_a^2,\qquad
\E(\text{MS}_{\text{within}})=\sigma^2 .
\]
Which ratio tests \( \sigma_a^2=0 \), and what unbiased estimator of \( \sigma_a^2 \) do the two mean
squares suggest?
:::

::: {.solution}
\( \Y=\mu\bone+\bm Z\bm a+\be \), so \( \E(\Y)=\mu\bone \) and
\( \Cov(\Y)=\bm Z(\sigma_a^2\I_g)\bm Z\T+\sigma^2\I \) by @thm-rv-linear. The between-groups
projection is \( \M-\bP_1 \), where \( \M=m^{-1}\bm Z\bm Z\T \) averages within groups, and the
within-groups projection is \( \I-\M \). The mean \( \mu\bone \) is annihilated by both. By
@thm-rv-quadform-mean, \( \E(\Y\T(\M-\bP_1)\Y)=\tr\{(\M-\bP_1)(\sigma^2\I+\sigma_a^2\bm Z\bm Z\T)\} \).
Now \( \bm Z\bm Z\T=m\M \), and \( (\M-\bP_1)\M=\M-\bP_1 \), so the trace is
\( \sigma^2(g-1)+\sigma_a^2m(g-1) \). For the within-groups projection, \( (\I-\M)\bm Z=\bzero \), so the
trace is \( \sigma^2(n-g) \). Dividing by the degrees of freedom gives the two expectations. The ratio
\( \text{MS}_{\text{between}}/\text{MS}_{\text{within}} \) has expectations that agree iff
\( \sigma_a^2=0 \), and \( (\text{MS}_{\text{between}}-\text{MS}_{\text{within}})/m \) is unbiased for
\( \sigma_a^2 \) (though it can be negative).
:::

::: {#exr-ss-omitted-term}
[B2]

Suppose \( \E(\Y)=\X\bbeta+\bz\delta \) but the model \( \E(\Y)=\X\bbeta \) is fitted, with
\( \bone\in\C(\X) \). Show that
\[
\E(\text{SSE})=\sigma^2(n-r)+\delta^2\norm{(\I-\M)\bz}^2,
\qquad
\E(\text{SSR})=\sigma^2(r-1)+\norm{(\M-\bP_1)(\X\bbeta+\bz\delta)}^2 .
\]
Give an example in which \( \bbeta \) has all slopes zero but the regression mean square is inflated.
:::

::: {#exr-ss-hypothesis-noncentrality}
[B3]

In the setting of @thm-ss-hypothesis with \( \bm d=\bzero \), show that the expected sum of squares of
the hypothesis \( \bLambda\T\bbeta=\bzero \) is
\( \sigma^2q+(\bLambda\T\bbeta)\T(\bLambda\T\G\bLambda)^{-1}\bLambda\T\bbeta \), so the
noncentrality is obtained by evaluating @eq-ss-hypothesis at the true parameter.
:::

::: {.solution}
By @thm-ss-hypothesis(b), the sum of squares is \( \Y\T\bm K\Y \) with
\( \bm K=\M\bT(\bLambda\T\G\bLambda)^{-1}\bT\T\M \), a projection of rank \( q \), where \( \bT \) is the
\( n\times q \) matrix with \( \bLambda=\X\T\bT \). By @thm-ss-expected-mean-squares(a) its expectation
is \( \sigma^2q+\bbeta\T\X\T\bm K\X\bbeta \). Since \( \bT\T\M\X\bbeta=\bT\T\X\bbeta=\bLambda\T\bbeta \),
the second term is \( (\bLambda\T\bbeta)\T(\bLambda\T\G\bLambda)^{-1}\bLambda\T\bbeta \).
:::

### C. Going deeper

::: {#exr-ss-kurtosis-sse}
[C1]

Let the errors be independent with mean zero, variance \( \sigma^2 \) and fourth moment
\( \mu_4=\E(\varepsilon_i^4) \), and let \( \bmu\in\C(\X) \). Use @thm-rv-quadform-variance to show that
\[
\Var(\text{SSE})=(\mu_4-3\sigma^4)\sum_i(1-h_{ii})^2+2\sigma^4(n-r),
\]
where \( h_{ii} \) are the leverages. For the centred exponential errors of @exm-ss-ems-simulation,
\( \mu_4=9\sigma^4 \). Show that \( \Var(s^2)\approx8\sigma^4/(n-r) \) when all leverages are small, and
compare with the simulated variance of the residual mean square.
:::

::: {.solution}
With \( \A=\I-\M \) and \( \bm\theta=\bmu \), the terms of @thm-rv-quadform-variance involving
\( \bm\theta \) vanish because \( (\I-\M)\bmu=\bzero \). The diagonal of \( \I-\M \) is \( 1-h_{ii} \), and
\( \tr\{(\I-\M)^2\}=n-r \), which gives the formula. If the leverages are small,
\( \sum_i(1-h_{ii})^2\approx n-r \) (it equals \( n-2r+\sum_ih_{ii}^2 \)), so
\( \Var(\text{SSE})\approx(6+2)\sigma^4(n-r) \) and \( \Var(s^2)\approx8\sigma^4/(n-r) \), four times the
normal value \( 2\sigma^4/(n-r) \). This is the factor seen in the simulation.
:::

::: {#exr-ss-power-design}
[C2]

For the Type I sum of squares of \( A \) entered first in a two-way layout, show that its
noncentrality can be positive when \( \mu_{ij}=\beta_j \) depends on \( B \) only, and find the condition on
the counts \( n_{ij} \) under which it is zero for all such \( \bmu \). Relate your answer to @thm-ss-proportional.
:::
